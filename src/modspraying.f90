!> MCB sprayers with evaporative cooling.
!!Stephan de Roode and Annelot Broerze

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
                             spray_Dg, spray_sigma_g, particle_emission_rate, ldistribution, isv_salt,tracer, lsalt_sponge, lcoupled, &
                             my_process_sprays, target_mode, isv_salt_n, isv_ss_acs, isv_ss_acs_n, isv_ss_cos, isv_ss_cos_n
  use modtracers,      only: add_tracer, get_tracer_index
  use modlogging,      only: message, warning, finish
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
                           spray_dg, spray_sigma_g, ldistribution,particle_emission_rate, tracer, lsalt_sponge, lcoupled, target_mode

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
    call D_MPI_BCAST(spray_dg,           1,  0, comm3d, mpierr)
    call D_MPI_BCAST(spray_sigma_g,           1,  0, comm3d, mpierr)
    call D_MPI_BCAST(ldistribution,           1,  0, comm3d, mpierr)
    call D_MPI_BCAST(particle_emission_rate,           1,  0, comm3d, mpierr)
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
      if (.not. ldistribution) then 
        salt_spray_rate = water_spray_rate * salinity ! directy coupled to water spray rate
      else 
        salt_spray_rate = 0
      endif
    else
      water_spray_rate = 0
    endif
    print * , "Salt spray rate:", salt_spray_rate

    if (lsalt_spraying) then
      if (lcoupled) then
        if (.not. ldistribution) then
                isv_salt = get_tracer_index('ss_'//target_mode)
                isv_salt_n = get_tracer_index(target_mode//'_n')
        else
                isv_ss_acs   = get_tracer_index('ss_acs')
                isv_ss_acs_n = get_tracer_index('acs_n')
                isv_ss_cos   = get_tracer_index('ss_cos')
                isv_ss_cos_n = get_tracer_index('cos_n')
        endif
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
    if ((lwater_spraying .or. lsalt_spraying) .and. &
        i_spray >= 2 .and. i_spray <= i1 .and. &
        j_spray >= 2 .and. j_spray <= j1 .and. &
        k_spray >= 1 .and. k_spray <= kmax) then
      call message(routine,'spraying point at myidx = ',myidx, ' myidy = ', myidy)
      call message(routine,'global locations ',i_glob_spray,j_glob_spray,k_glob_spray)
      call message(routine,'local locations ',i_spray,j_spray,k_spray)
      my_process_sprays = .true.
    else  ! if not, there is no sprayer here
      my_process_sprays = .false.
      i_spray = -999
      j_spray = -999
      k_spray = -999
    endif

  end subroutine initspraying

  real function lognormal_cdf(d, dg, sigma_g)

    real, intent(in) :: d, dg, sigma_g
    real :: z

    z = (log(d) - log(dg)) / (sqrt(2.0)*log(sigma_g))

    lognormal_cdf = 0.5 * (1.0 + erf(z))

  end function lognormal_cdf

  real function lognormal_cdf_mass(d, dg, sigma_g)

    real, intent(in) :: d, dg, sigma_g
    real :: z

    z = ( log(d) - log(dg) - 3.0*log(sigma_g)**2 ) / &
        ( sqrt(2.0)*log(sigma_g) )

    lognormal_cdf_mass = 0.5 * (1.0 + erf(z))

  end function lognormal_cdf_mass

  !> Apply spraying tendencies to the model fields.
  subroutine spraying()

    real(field_r) :: dqldt_spraying, dsvdt_spraying
    real(field_r) :: dm, dn
    real(field_r) :: cell_volume !< Air density times grid cell volume [kg]
    real(field_r) :: fracn_acs, fracn_cos, fracm_acs, fracm_cos, ndot_acs, ndot_cos, dn_acs, dn_cos, dm_acs, dm_cos, mdot_total,&
           ndot_total, mdot_acs, mdot_cos, mean_particle_mass
    real(field_r) :: dacs_max = 500e-9

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
                if (.not. ldistribution) then       
                        dm = salt_spray_rate / (rhobf(k_spray) * cell_volume)

                        ! Increase in number concentration, assuming monodisperse aerosol
                        dn = salt_spray_rate / (2165.0 * pi / 6 * (spray_dg)**3)
                        dn = dn / cell_volume ! Number concentrations are in #/m3

                        !$acc serial default(present) async
                        svp(i_spray,j_spray,k_spray,isv_salt) = &
                        svp(i_spray,j_spray,k_spray,isv_salt) + dm
        
                        svp(i_spray,j_spray,k_spray,isv_salt_n) = &
                        svp(i_spray,j_spray,k_spray,isv_salt_n) + dn
                        !$acc end serial
                else

                        !Calculate fraction of number concentration for acs and cos between respective boundaries
                        !fracN_acs = lognormal_cdf(dacs_max,spray_Dg,spray_sigma_g) - lognormal_cdf(dacs_min, spray_Dg,
                        !spray_sigma_g)
                        !fracN_cos = lognormal_cdf(dcos_max, spray_Dg, spray_sigma_g) - lognormal_cdf(dcos_min, spray_Dg, spray_sigma_g

                        fracn_acs = lognormal_cdf(dacs_max,spray_dg,spray_sigma_g)
                        fracn_cos = 1.0 - fracn_acs

                        !Calculate number for acs and cos
                        ndot_acs = particle_emission_rate * fracn_acs
                        ndot_cos = particle_emission_rate * fracn_cos

                        !Do the same for mass
                        !fracM_acs = lognormal_cdf_mass(dacs_max,spray_Dg,spray_sigma_g) - lognormal_cdf_mass(dacs_min, spray_Dg,
                        !spray_sigma_g)
                        !fracM_cos = lognormal_cdf_mass(dcos_max, spray_Dg, spray_sigma_g) - lognormal_cdf_mass(dcos_min, spray_Dg,
                        !spray_sigma_g)
                        fracm_acs = lognormal_cdf_mass(dacs_max,spray_dg,spray_sigma_g)
                        fracm_cos = 1.0 - fracm_acs

                        print *, fracn_acs, fracn_cos, fracn_acs + fracn_cos

                        !Calculate mean particle mass
                        mean_particle_mass = 2165.0 * pi / 6.0 * spray_dg**3 * &
                        exp(4.5 * log(spray_sigma_g)**2)

                        !Calculate total emitted mass rate, and for acs, cos
                        mdot_total = particle_emission_rate * mean_particle_mass

                        mdot_acs = fracM_acs * mdot_total
                        mdot_cos = fracM_cos * mdot_total
                        
                        print *, "particle_emission_rate =", particle_emission_rate
                        print *, "mean_particle_mass =", mean_particle_mass
                        print *, "mdot_total =", mdot_total

                        print *, "Number fractions:", fracn_acs, fracn_cos
                        print *, "Mass fractions:  ", fracm_acs, fracm_cos

                        !Add to M7 tracers
                        dn_acs = ndot_acs / cell_volume
                        dn_cos = ndot_cos / cell_volume

                        dm_acs = mdot_acs / (rhobf(k_spray)*cell_volume)
                        dm_cos = mdot_cos / (rhobf(k_spray)*cell_volume)

                        print *, "Number emissions:", dn_acs, dn_cos
                        print *, "Mass emissions:  ", dm_acs, dm_cos

                        svp(i_spray,j_spray,k_spray,isv_ss_acs) = &
                        svp(i_spray,j_spray,k_spray,isv_ss_acs) + dm_acs

                        svp(i_spray,j_spray,k_spray,isv_ss_cos) = &
                        svp(i_spray,j_spray,k_spray,isv_ss_cos) + dm_cos

                        svp(i_spray,j_spray,k_spray,isv_ss_acs_n) = &
                        svp(i_spray,j_spray,k_spray,isv_ss_acs_n) + dn_acs

                        svp(i_spray,j_spray,k_spray,isv_ss_cos_n) = &
                        svp(i_spray,j_spray,k_spray,isv_ss_cos_n) + dn_cos
                endif
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
