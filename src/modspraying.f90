!> \file modspraying.f99
!! Stephan de Roode and Annelot Broerze



module modspraying
   use modaerosol_common,  only: iACS, iCOS, iSS
   use modaerosol,         only: modes_f
   use modaerosol_sources, only: aerosol_point_source, MASS_SOURCE
   use modprecision, only: field_r
   use modsprayingdata, only: i_glob_spray,j_glob_spray,k_glob_spray,&
                              i_spray,j_spray,k_spray,&
                              water_spray_rate,salt_spray_rate,&
                              lwater_spraying,lsalt_spraying,salinity,&
                              isv_salt,tracer,lsalt_sponge, lcoupled, my_process_sprays, target_mode, &
                              isv_salt_n
   implicit none

  ! for lateral sponge
   public :: spraying

contains
  subroutine initspraying
  use modglobal,    only : i1,j1,imax,jmax,kmax,ifnamopt,fname_options,checknamelisterror
  use modmpi,       only : myid,myidx,myidy,comm3d, mpierr, d_mpi_bcast
  use modtracers,   only: add_tracer, get_tracer_index
  use fortran_support, only: nnml_output
  use modlogging,      only : profile_output
  !use modnudgeboundary, only : lnudgeboundary


  integer ierr

  namelist/NAMSPRAYING/ lwater_spraying,lsalt_spraying,&
                        i_glob_spray,j_glob_spray,k_glob_spray,&
                        water_spray_rate,salt_spray_rate,salinity,&
                        tracer,lsalt_sponge, lcoupled, target_mode

  if(myid==0) then    !first myid
    open(ifnamopt,file=fname_options,status='old',iostat=ierr)
    read (ifnamopt,NAMSPRAYING,iostat=ierr)
    call checknamelisterror(ierr, ifnamopt, 'NAMSPRAYING')
    write(nnml_output ,NAMSPRAYING)
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
  
  if (lwater_spraying) then
     lsalt_spraying  = .true.
     salt_spray_rate = water_spray_rate * salinity   !directy coupled to water spray rate
  else
     water_spray_rate = 0.
  endif

  if (lsalt_spraying) then
    if (.not. lcoupled) then
     call add_tracer(trim(tracer), long_name=trim(tracer)//" mixing ratio", &
          unit="kg/kg", isv=isv_salt)
    else
      isv_salt = get_tracer_index('ss_'//target_mode)
      isv_salt_n = get_tracer_index(target_mode//'_n')
    end if
  endif

  !determine local position of spraying from global position
  i_spray = i_glob_spray - myidx*imax
  j_spray = j_glob_spray - myidy*jmax
  k_spray = k_glob_spray

  ! are the local coordinates actually in the domain?
  if (i_spray >= 2 .and. i_spray <= i1 .and. &
       j_spray >= 2 .and. j_spray <= j1 .and. &
       k_spray >= 1 .and. k_spray <= kmax) then
     write(profile_output,*) 'spraying point at myid = ',myid
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

  subroutine spraying()
    use modfields, only: qtp, qt0, thlp, exnf, ql0, svp, sv0
    use modglobal, only: dx, dy, dzf, cp, rlv, pi
    use modtracers, only: get_tracer_index

    real(field_r) :: dqldt_spraying, dsvdt_spraying
    real(field_r) :: dm, dn

    if (lwater_spraying .and. my_process_sprays) then
      dqldt_spraying = water_spray_rate/(rhobf(k_spray)*dx*dy*dzf(k_spray)) 
      qtp(i_spray,j_spray,k_spray)  = qtp(i_spray,j_spray,k_spray)  &
                   + (1-qt0(i_spray,j_spray,k_spray)) * dqldt_spraying
      ! Evaporative cooling
      thlp(i_spray,j_spray,k_spray) = thlp(i_spray,j_spray,k_spray) & 
                   - (rlv/(cp*exnf(k_spray)))* (1-ql0(i_spray,j_spray,k_spray))* dqldt_spraying
    end if

    if (lsalt_spraying .and. my_process_sprays) then
      if (lcoupled) then
        dm = salt_spray_rate / (rhobf(k_spray)*dx*dy*dzf(k_spray))

        ! Increase in number concentration, assuming monodisperse aerosol
        dn = salt_spray_rate / (2165.0 * pi / 6 * (75e-9)**3)
        dn = dn / (dx*dy*dzf(k_spray))

        svp(i_spray,j_spray,k_spray,isv_salt) = &
          svp(i_spray,j_spray,k_spray,isv_salt) + dm
       
        svp(i_spray,j_spray,k_spray,isv_salt_n) = &
          svp(i_spray,j_spray,k_spray,isv_salt_n) + dn
      else
        dsvdt_spraying = salt_spray_rate/(rhobf(k_spray)*dx*dy*dzf(k_spray)) * &
                          (1-sv0(i_spray,j_spray,k_spray,isv_salt)/salinity)
        svp(i_spray,j_spray,k_spray,isv_salt)  = svp(i_spray,j_spray,k_spray,isv_salt) &
                     + dsvdt_spraying
      end if
    endif

  end subroutine spraying

end module modspraying
