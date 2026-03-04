!> MCB sprayers with evaporative cooling.
!!
!! @author Stephan de Roode
!! @author Annelot Broerze
module modspraying
  use fortran_support, only: nnml_output
  use modfields,       only: qtp, qt0, thlp, exnf, ql0, svp, sv0, rhobf
  use modglobal,       only: dx, dy, dzf, i1, j1, imax, jmax, kmax, ifnamopt, &
                             fname_options, checknamelisterror, cp, rlv, pi
  use modlogging,      only: profile_output
  use modmpi,          only: myid, myidx, myidy, comm3d, mpierr, d_mpi_bcast
  use modprecision,    only: field_r
  use modsprayingdata, only: i_glob_spray, j_glob_spray, k_glob_spray, &
                             i_spray, j_spray, k_spray, &
                             water_spray_rate, salt_spray_rate, &
                             lwater_spraying, lsalt_spraying, salinity, &
                             isv_salt,tracer, lsalt_sponge, lcoupled, &
                             my_process_sprays, target_mode, isv_salt_n
  use modtracers,      only: add_tracer, get_tracer_index

  implicit none

  private

  character(len=*), parameter :: modname = 'modspraying'

  public :: spraying_read_namelist
  public :: initspraying
  public :: spraying

contains

  !> Read the namelist for spraying.
  subroutine spraying_read_namelist(nml_filename)

    character(len=*), intent(in) :: nml_filename

    character(len=*), parameter :: routine = modname//'/spraying_read_namelist'

    integer :: ierr

    namelist /namspraying/ lwater_spraying, lsalt_spraying, &
                           i_glob_spray, j_glob_spray, k_glob_spray, &
                           water_spray_rate, salt_spray_rate, salinity, &
                           tracer, lsalt_sponge, lcoupled, target_mode

    if(myid==0) then
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read(ifnamopt, namspraying, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMSPRAYING')
      write(nnml_output, namspraying)
      close(ifnamopt)
    endif

    call D_MPI_BCAST(lwater_spraying ,    1,  0, comm3d, mpierr)
    call D_MPI_BCAST(lsalt_spraying  ,    1,  0, comm3d, mpierr)
    call D_MPI_BCAST(i_glob_spray    ,    1,  0, comm3d, mpierr)
    call D_MPI_BCAST(j_glob_spray    ,    1,  0, comm3d, mpierr)
    call D_MPI_BCAST(k_glob_spray    ,    1,  0, comm3d, mpierr)
    call D_MPI_BCAST(water_spray_rate,    1,  0, comm3d, mpierr)
    call D_MPI_BCAST(salt_spray_rate,     1,  0, comm3d, mpierr)
    call D_MPI_BCAST(tracer,             20,  0, comm3d, mpierr)
    call D_MPI_BCAST(lsalt_sponge,        1,  0, comm3d, mpierr)
    call D_MPI_BCAST(lcoupled,            1,  0, comm3d, mpierr)
    call D_MPI_BCAST(target_mode,         3, 0, comm3d, mpierr)

  end subroutine spraying_read_namelist

  !> Initialize spraying parameters and determine local spraying location.
  subroutine initspraying

    character(len=*), parameter :: routine = modname//'/initspraying'

    if (lwater_spraying) then
      lsalt_spraying  = .true.
      salt_spray_rate = water_spray_rate * salinity ! directy coupled to water spray rate
    else
      water_spray_rate = 0
    endif

    if (lsalt_spraying) then
      if (lcoupled) then
        isv_salt = get_tracer_index('ss_'//target_mode)
        isv_salt_n = get_tracer_index(target_mode//'_n')
      else
        call add_tracer(trim(tracer), long_name=trim(tracer)//" mixing ratio", &
                        unit="kg/kg", isv=isv_salt)
      end if
    endif

    ! determine local position of spraying from global position
    i_spray = i_glob_spray - myidx*imax
    j_spray = j_glob_spray - myidy*jmax
    k_spray = k_glob_spray

    ! are the local coordinates actually in the domain?
    if (i_spray >= 2 .and. i_spray <= i1 .and. &
        j_spray >= 2 .and. j_spray <= j1 .and. &
        k_spray >= 1 .and. k_spray <= kmax) then
      write(profile_output,*) 'spraying point at myidx = ',myidx, ' myidy = ', myidy
      write(profile_output,*) 'global locations ',i_glob_spray,j_glob_spray,k_glob_spray
      write(profile_output,*) 'local locations ',i_spray,j_spray,k_spray
      my_process_sprays = .true.
    else  ! if not, there is no sprayer here
      my_process_sprays = .false.
      i_spray = -999
      j_spray = -999
      k_spray = -999
    endif

    if (myid==0) then
      write(profile_output,*) 'Spraying data used: '
      write(profile_output,*) 'lwater_spraying     ',lwater_spraying
      write(profile_output,*) 'lsalt_spraying      ',lsalt_spraying
      write(profile_output,*) 'i_glob_spray        ',i_glob_spray
      write(profile_output,*) 'j_glob_spray        ',j_glob_spray
      write(profile_output,*) 'k_glob_spray        ',k_glob_spray
      write(profile_output,*) 'water_spray_rate    ',water_spray_rate
      write(profile_output,*) 'salt_spray_rate     ',salt_spray_rate
      write(profile_output,*) 'salt scalar number  ',isv_salt
      write(profile_output,*)
    endif

  end subroutine initspraying

  !> Apply spraying tendencies to the model fields.
  subroutine spraying()

    real(field_r) :: dqldt_spraying, dsvdt_spraying
    real(field_r) :: dm, dn
    real(field_r) :: cell_volume !< Air density times grid cell volume [kg]

    if (my_process_sprays) then

      cell_volume = dx * dy * dzf(k_spray)

      if (lwater_spraying) then
        dqldt_spraying = water_spray_rate / (rhobf(k_spray) * cell_volume)

        !$acc serial default(present) async
        qtp(i_spray,j_spray,k_spray) = qtp(i_spray,j_spray,k_spray) &
          + (1-qt0(i_spray,j_spray,k_spray)) * dqldt_spraying

        ! Evaporative cooling
        thlp(i_spray,j_spray,k_spray) = thlp(i_spray,j_spray,k_spray) & 
          - (rlv / (cp * exnf(k_spray))) &
          * (1 - ql0(i_spray,j_spray,k_spray)) * dqldt_spraying
        !$acc end serial
      end if

      if (lsalt_spraying) then
        if (lcoupled) then
          dm = salt_spray_rate / (rhobf(k_spray) * cell_volume)

          ! Increase in number concentration, assuming monodisperse aerosol
          dn = salt_spray_rate / (2165.0 * pi / 6 * (75e-9)**3)
          dn = dn / cell_volume ! Number concentrations are in #/m3

          !$acc serial default(present) async
          svp(i_spray,j_spray,k_spray,isv_salt) = &
            svp(i_spray,j_spray,k_spray,isv_salt) + dm
        
          svp(i_spray,j_spray,k_spray,isv_salt_n) = &
            svp(i_spray,j_spray,k_spray,isv_salt_n) + dn
          !$acc end serial
        else
          dsvdt_spraying = salt_spray_rate / (rhobf(k_spray) * cell_volume) &
            * (1 - sv0(i_spray,j_spray,k_spray,isv_salt) / salinity)

          !$acc serial default(present) async
          svp(i_spray,j_spray,k_spray,isv_salt) = &
            svp(i_spray,j_spray,k_spray,isv_salt) + dsvdt_spraying
          !$acc end serial
        end if
      endif
    end if

    !$acc wait

  end subroutine spraying

end module modspraying
