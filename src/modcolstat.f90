!> Point-column statistics output.
module modcolstat

  use fortran_support,   only: nnml_output
  use modglobal,         only: ifnamopt, fname_options, checknamelisterror, &
                               i1, j1, kmax,k1, itot, jtot, imax, jmax, &
                               dtav_glob, timeav_glob, &
                               cp, rlv, rd, rv, x0, y0, dx, dy, dzf, dzh, dzhi, dxi, dyi, cu, cv, eps1, nsv, grav, tdn, tup, output_prefix
  use modfields,         only: um, vm, wm, thlm, qtm, u0, v0, w0, thl0, thl0h, qt0, qt0h, ql0, ql0h, thv0h, &
                               e12m, e120, exnf, exnh, presf, presh, rhof, rhobf, rhobh, tmp0, sv0, svm, svp, thlpcar, u0av, v0av, thvh
  use modgenstat,        only: umav, vmav, wmav, thlmav, thvmav, thmav, qtmav, qlmav, svmav, w2av, dtav, timeav
  use modsubgriddata,    only: ekm, ekh, sbshr, sbbuo, sbdiss, csz
  use modsurfdata,       only: thlflux, qtflux, ustar, thls, qts, svflux
  use modraddata,        only: lwd, lwu, swd, swu, lwdca, lwuca, swdca, swuca, thlprad
  use modthermodynamics, only: qsat_tab, thv0
  use modtracers,        only: tracer_prop, get_tracer_index
  use modmicrodata,      only: imicro, imicro_sice, imicro_sice2
  use modsimpleice_data, only: tuprsg, tdnrsg
  use modpois_data,      only: p
  use modmpi,            only: myid, myidx, myidy, comm3d, mpierr, d_mpi_bcast, &
                               nprocs, D_MPI_ALLREDUCE, mpi_sum
  use modlogging,        only: finish
  use modnetcdf_file_t,  only: multi_profile_file_t
  use modprecision,      only: field_r
  use modstat_nc_files,  only: add_output_file, is_sampling_timestep, is_writing_timestep

  implicit none

  private

  character(len=*), parameter :: modname = 'modcolstat'
  integer, parameter :: max_points = 999 !< Max number of column statistics points.

  public :: initcolstat
  public :: colstat
  public :: exitcolstat

  logical :: lcolstat = .false.
  integer :: nsamples = 1
  logical :: l_sbtke_beg_set = .false.
  real(field_r), allocatable :: sbtke_beg_col(:,:), sbtke_last_col(:,:)
  integer :: npoints = 0
  integer :: x_idx(max_points) = 0
  integer :: y_idx(max_points) = 0

  type(multi_profile_file_t) :: ofile
  integer :: ofile_id
  integer, allocatable :: i_local(:), j_local(:), global_index(:)
  real(field_r), allocatable :: locx(:), locy(:)

contains

  subroutine initcolstat

    integer :: ierr, ip, nlocal, n, il, rank_offset
    integer, allocatable :: nlocal_send(:), nlocal_all_arr(:)
    integer :: x_idx_global, y_idx_global
    real(field_r) :: locx_point, locy_point
    real(field_r), allocatable :: locx_all(:), locy_all(:)

    character(len=*), parameter :: routine = modname//'/initcolstat'

    namelist /NAMCOLSTAT/ lcolstat, npoints, x_idx, y_idx

    x_idx = 0
    y_idx = 0
    npoints = 0

    if (myid == 0) then
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read(ifnamopt, NAMCOLSTAT, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMCOLSTAT')
      write(nnml_output, NAMCOLSTAT)
      close(ifnamopt)
    end if

    call d_mpi_bcast(lcolstat, 1, 0, comm3d, mpierr)
    call d_mpi_bcast(npoints,  1, 0, comm3d, mpierr)
    call d_mpi_bcast(x_idx, max_points, 0, comm3d, mpierr)
    call d_mpi_bcast(y_idx, max_points, 0, comm3d, mpierr)

    if (.not. lcolstat) return

    nsamples = nint(timeav / dtav)
    if (nsamples < 1) then
      call finish(routine, 'NAMGENSTAT: timeav must be >= dtav when NAMCOLSTAT:lcolstat is true')
    end if
    if (abs(timeav / dtav - real(nsamples)) > 1e-4) then
      call finish(routine, 'NAMGENSTAT: timeav must be an integer multiple of dtav when NAMCOLSTAT:lcolstat is true')
    end if
    if (npoints < 1 .or. npoints > max_points) then
      call finish(routine, 'NAMCOLSTAT: npoints must be between 1 and ', max_points)
    end if

    if (dtav <= 0.0) then
      call finish(routine, 'NAMGENSTAT: dtav must be > 0 when NAMCOLSTAT:lcolstat is true')
    end if

    do ip = 1, npoints
      if (x_idx(ip) < 1 .or. x_idx(ip) > itot) then
        call finish(routine, 'NAMCOLSTAT: x_idx out of range at entry ', ip)
      end if
      if (y_idx(ip) < 1 .or. y_idx(ip) > jtot) then
        call finish(routine, 'NAMCOLSTAT: y_idx out of range at entry ', ip)
      end if
    end do

    nlocal = 0
    do ip = 1, npoints
      x_idx_global = x_idx(ip)
      y_idx_global = y_idx(ip)
      if (x_idx_global >= myidx*imax + 1 .and. x_idx_global <= myidx*imax + imax .and. &
          y_idx_global >= myidy*jmax + 1 .and. y_idx_global <= myidy*jmax + jmax) then
        nlocal = nlocal + 1
      end if
    end do

    allocate(i_local(nlocal), j_local(nlocal), global_index(nlocal))
    allocate(locx(nlocal), locy(nlocal))
    allocate(locx_all(npoints), locy_all(npoints))

    ! here we determine the global index of each local point, and the local x,y coordinates of each point
    nlocal = 0
    do ip = 1, npoints
      x_idx_global = x_idx(ip)
      y_idx_global = y_idx(ip)
      if (x_idx_global >= myidx*imax + 1 .and. x_idx_global <= myidx*imax + imax .and. &
          y_idx_global >= myidy*jmax + 1 .and. y_idx_global <= myidy*jmax + jmax) then
        nlocal = nlocal + 1
        i_local(nlocal) = x_idx_global - myidx*imax + 1
        j_local(nlocal) = y_idx_global - myidy*jmax + 1
        global_index(nlocal) = ip
        locx_point = real(x0 + dx * (real(x_idx_global,kind=field_r) - 0.5_field_r), kind=field_r)
        locy_point = real(y0 + dy * (real(y_idx_global,kind=field_r) - 0.5_field_r), kind=field_r)
        locx(nlocal) = locx_point
        locy(nlocal) = locy_point
      end if
    end do

    ! we want to write each measurement into a contiguous slice in a collective array in the output file.
    ! to do this, we need to know the offset of each rank's local points in the global array.    
    ! we create an array of length nprocs. then in index myid we store the number of local points.
    ! myid starts at 0, so we add 1 for proper indicing
    allocate(nlocal_send(nprocs), nlocal_all_arr(nprocs))
    nlocal_send = 0
    nlocal_send(myid + 1) = nlocal

    ! then we do an allreduce sum to get the total number of points on each rank.
    call D_MPI_ALLREDUCE(nlocal_send, nlocal_all_arr, nprocs, MPI_SUM, comm3d, mpierr)
    ! we sum the array up to myid-1 to get the offset for this rank UP UNTIL THIS ID.
    ! that way, we have counted all measurement points before our own in the global array
    rank_offset = sum(nlocal_all_arr(1:myid))
    deallocate(nlocal_send, nlocal_all_arr)

    ! we construct the contiguous slice of points in the global array
    global_index = [(rank_offset + il, il = 1, nlocal)]

    ! set all global x,y coordinates to 0
    locx_all = 0.0_field_r
    locy_all = 0.0_field_r
    ! hopefully, at this point all ranks are properly separated
    do il = 1, nlocal
      locx_all(rank_offset + il) = locx(il)
      locy_all(rank_offset + il) = locy(il)
    end do
    ! when we sum, and we're properly separated, each PE will be left with the correct coordinate array
    call D_MPI_ALLREDUCE(locx_all, npoints, MPI_SUM, comm3d, mpierr)
    call D_MPI_ALLREDUCE(locy_all, npoints, MPI_SUM, comm3d, mpierr)

    ofile = multi_profile_file_t(trim(output_prefix)//'colstat', npoints, nz=kmax, lgpu=.false., &
                                locx=locx_all, locy=locy_all, profile_ids=global_index)
    call add_output_file(ofile, dtav, ofile_id, dt_write=timeav)

    call ofile%add_var('rhof','Full level slab averaged density','kg/m^3','izt','tttt')
    call ofile%add_var('rhobf','Full level base-state density','kg/m^3','izt','tttt')
    call ofile%add_var('rhobh','Half level base-state density','kg/m^3','izmt','ttmt')
    call ofile%add_var('presh','Pressure at cell center','Pa','izt','tttt')
    call ofile%add_var('u','West-East velocity','m/s','izt','mttt')
    call ofile%add_var('v','South-North velocity','m/s','izt','tmtt')
    call ofile%add_var('w','Vertical velocity','m/s','izmt','ttmt')
    call ofile%add_var('thl','Liquid water potential temperature','K','izt','tttt')
    call ofile%add_var('thv','Virtual potential temperature','K','izt','tttt')
    call ofile%add_var('qt','Total water specific humidity','kg/kg','izt','tttt')
    call ofile%add_var('ql','Liquid water specific humidity','kg/kg','izt','tttt')
    call ofile%add_var('wthls','SFS-Theta_l flux','Km/s','izmt','ttmt')
    call ofile%add_var('wthlr','Resolved Theta_l flux','Km/s','izmt','ttmt')
    call ofile%add_var('wthlt','Total Theta_l flux','Km/s','izmt','ttmt')
    call ofile%add_var('wthvs','SFS-buoyancy flux','Km/s','izmt','ttmt')
    call ofile%add_var('wthvr','Resolved buoyancy flux','Km/s','izmt','ttmt')
    call ofile%add_var('wthvt','Total buoyancy flux','Km/s','izmt','ttmt')
    call ofile%add_var('wqts','SFS-moisture flux','kg/kg m/s','izmt','ttmt')
    call ofile%add_var('wqtr','Resolved moisture flux','kg/kg m/s','izmt','ttmt')
    call ofile%add_var('wqtt','Total moisture flux','kg/kg m/s','izmt','ttmt')
    call ofile%add_var('wqls','SFS-liquid water flux','kg/kg m/s','izmt','ttmt')
    call ofile%add_var('wqlr','Resolved liquid water flux','kg/kg m/s','izmt','ttmt')
    call ofile%add_var('wqlt','Total liquid water flux','kg/kg m/s','izmt','ttmt')
    call ofile%add_var('uws','SFS-momentum flux (uw)','m^2/s^2','izmt','ttmt')
    call ofile%add_var('uwr','Resolved momentum flux (uw)','m^2/s^2','izmt','ttmt')
    call ofile%add_var('uwt','Total momentum flux (uw)','m^2/s^2','izmt','ttmt')
    call ofile%add_var('vws','SFS-momentum flux (vw)','m^2/s^2','izmt','ttmt')
    call ofile%add_var('vwr','Resolved momentum flux (vw)','m^2/s^2','izmt','ttmt')
    call ofile%add_var('vwt','Total momentum flux (vw)','m^2/s^2','izmt','ttmt')
    call ofile%add_var('w2s','SFS-TKE','m^2/s^2','izmt','ttmt')
    call ofile%add_var('w2r','Resolved vertical velocity variance','m^2/s^2','izmt','ttmt')
    call ofile%add_var('skew','vertical velocity skewness','-','izmt','ttmt')
    call ofile%add_var('u2r','Resolved horizontal velocity variance (u)','m^2/s^2','izt','mttt')
    call ofile%add_var('v2r','Resolved horizontal velocity variance (v)','m^2/s^2','izt','tmtt')
    call ofile%add_var('thl2r','Resolved theta_l variance','K^2','izt','tttt')
    call ofile%add_var('thv2r','Resolved buoyancy variance','K^2','izt','tttt')
    call ofile%add_var('th2r','Resolved theta variance','K^2','izt','tttt')
    call ofile%add_var('qt2r','Resolved total water variance','(kg/kg)^2','izt','tttt')
    call ofile%add_var('ql2r','Resolved liquid water variance','(kg/kg)^2','izt','tttt')
    call ofile%add_var('cs','Smagorinsky constant','-','izt','tttt')
    call ofile%add_var('cfrac','Cloud fraction','-','izt','tttt')
    call ofile%add_var('hur','Relative humidity','%','izt','tttt')
    call ofile%add_var('hus','Specific humidity','kg/kg','izt','tttt')
    call ofile%add_var('ta', 'Temperature','K','izt','tttt')
    call ofile%add_var('clw', 'Specific cloud liquid water content','kg/kg','izt','tttt')
    call ofile%add_var('cli', 'Specific cloud ice content','kg/kg','izt','tttt')
    call ofile%add_var('plw', 'Specific precipitation liquid water content','kg/kg','izt','tttt')
    call ofile%add_var('pli', 'Specific precipitation ice content','kg/kg','izt','tttt')
    do n = 1, nsv
      call ofile%add_var(trim(tracer_prop(n)%tracname), trim(tracer_prop(n)%traclong), trim(tracer_prop(n)%unit),'izt','tttt')
      call ofile%add_var(trim(tracer_prop(n)%tracname)//'p', trim(tracer_prop(n)%traclong)//' tendency',trim(tracer_prop(n)%unit)//'/s)','izt','tttt')
      call ofile%add_var(trim(tracer_prop(n)%tracname)//'2r','Resolved '//trim(tracer_prop(n)%traclong)//' variance','('//trim(tracer_prop(n)%unit)//')^2','izt','tttt')
      call ofile%add_var('w'//trim(tracer_prop(n)%tracname)//'s','SFS '//trim(tracer_prop(n)%traclong)//' flux',trim(tracer_prop(n)%unit)//' m/s','izmt','ttmt')
      call ofile%add_var('w'//trim(tracer_prop(n)%tracname)//'r','Resolved '//trim(tracer_prop(n)%traclong)//' flux',trim(tracer_prop(n)%unit)//' m/s','izmt','ttmt')
      call ofile%add_var('w'//trim(tracer_prop(n)%tracname)//'t','Total '//trim(tracer_prop(n)%traclong)//' flux',trim(tracer_prop(n)%unit)//' m/s','izmt','ttmt')
    end do
    ! modbudget variables
    call ofile%add_var('tker','Resolved TKE','kg/ms^2','izt','tttt')
    call ofile%add_var('shr','Resolved Shear','kg/ms^3','izt','tttt')
    call ofile%add_var('buo','Resolved Buoyancy','kg/ms^3','izt','tttt')
    call ofile%add_var('trsp','Resolved Transport','kg/ms^3','izt','tttt')
    call ofile%add_var('ptrsp','Resolved Pressure transport (redistribution)','kg/ms^3','izt','tttt')
    call ofile%add_var('sbtke','Subgrid TKE','kg/ms^2','izt','tttt')
    call ofile%add_var('sbshr','Subgrid Shear','kg/m^2s^2','izt','tttt')
    call ofile%add_var('sbbuo','Subgrid Buoyancy','kg/m^2s^2','izt','tttt')
    call ofile%add_var('sbdiss','Subgrid Dissipation','kg/m^2s^2','izt','tttt')
    call ofile%add_var('sbstor','Subgrid Storage','kg/m^2s^2','izt','tttt')
    call ofile%add_var('sbbudg','Subgrid Budget = sum of contributions excl storage','kg/m^2s^2','izt','tttt')
    call ofile%add_var('sbresid','Subgrid Residual = budget - storage','kg/m^2s^2','izt','tttt')
    call ofile%add_var('ekm','Turbulent exchange coefficient momentum','m^2/s','izt','tttt')
    call ofile%add_var('khkm','Kh / Km, in post-processing used to determine filter-grid ratio','-','izt','tttt')
    ! modradstat variables
    call ofile%add_var('thltend','Total radiative tendency','K/s','izt','tttt')
    call ofile%add_var('thllwtend','Long wave radiative tendency','K/s','izt','tttt')
    call ofile%add_var('thlswtend','Short wave radiative tendency','K/s','izt','tttt')
    call ofile%add_var('thlradls','Prescribed large scale radiative tendency','K/s','izt','tttt')
    call ofile%add_var('lwu','Long wave upward radiative flux','W/m^2','izmt','ttmt')
    call ofile%add_var('lwd','Long wave downward radiative flux','W/m^2','izmt','ttmt')
    call ofile%add_var('swu','Short wave upward radiative flux','W/m^2','izmt','ttmt')
    call ofile%add_var('swd','Short wave downward radiative flux','W/m^2','izmt','ttmt')
    call ofile%add_var('lwuca','Long wave clear air upward radiative flux','W/m^2','izmt','ttmt')
    call ofile%add_var('lwdca','Long wave clear air downward radiative flux','W/m^2','izmt','ttmt')
    call ofile%add_var('swuca','Short wave clear air upward radiative flux','W/m^2','izmt','ttmt')
    call ofile%add_var('swdca','Short wave clear air downward radiative flux','W/m^2','izmt','ttmt')
    call ofile%add_var('thllwtendca','Long wave clear air radiative tendency','K/s','izt','tttt')
    call ofile%add_var('thlswtendca','Short wave clear air radiative tendency','K/s','izt','tttt')
    
    deallocate(locx_all, locy_all)

    allocate(sbtke_beg_col(kmax, nlocal), sbtke_last_col(kmax,nlocal))

  end subroutine initcolstat

  subroutine colstat

    integer :: col_idx, i, j, k, n, iqr
    real(field_r) :: ekhalf, euhalf, evhalf, wthls_loc, wthlr_loc, wqts_loc, wqtr_loc
    real(field_r) :: wqls_loc, wqlr_loc, wthvs_loc, wthvr_loc, uws_loc, uwr_loc, vws_loc, vwr_loc
    real(field_r) :: upcu, vpcv, egp, ph, uwrs, vwrs, sbtke_inst
    real(field_r) :: qs0h, t0h, den, cthl, cqt, a_dry, b_dry, a_moist, b_moist, c1, c2, tsurf, qsat
    real(field_r) :: weres_loc(k1), ptrsp_loc(k1), buo_loc(k1), trsp_loc(k1), shr_loc(k1)
    real(field_r) :: u2r_inst, v2r_inst, w_prime, w2r_inst, wsvs_loc, wsvr_loc
    real(field_r) :: scale_rn, ilratio, prof_tmp(kmax)

    ! genstat variables
    real(field_r), pointer :: rhof_col(:), rhobf_col(:), rhobh_col(:), presh_col(:)
    real(field_r), pointer :: u_col(:), v_col(:), w_col(:), thl_col(:), thv_col(:), qt_col(:), ql_col(:)
    real(field_r), pointer :: wthls_col(:), wthlr_col(:), wthlt_col(:), wthvs_col(:), wthvr_col(:), wthvt_col(:)
    real(field_r), pointer :: wqts_col(:), wqtr_col(:), wqtt_col(:), wqls_col(:), wqlr_col(:), wqlt_col(:)
    real(field_r), pointer :: uws_col(:), uwr_col(:), uwt_col(:), vws_col(:), vwr_col(:), vwt_col(:)
    real(field_r), pointer :: w2s_col(:), w2r_col(:), skew_col(:), u2r_col(:), v2r_col(:)
    real(field_r), pointer :: thl2r_col(:), thv2r_col(:), th2r_col(:), qt2r_col(:), ql2r_col(:)
    real(field_r), pointer :: cs_col(:), cfrac_col(:), hur_col(:), hus_col(:), ta_col(:)
    real(field_r), pointer :: clw_col(:), cli_col(:), plw_col(:), pli_col(:)

    real(field_r), pointer :: sv_col(:), svp_col(:), sv2r_col(:)
    real(field_r), pointer :: wsvs_col(:), wsvr_col(:), wsvt_col(:)
    ! modbudget variables
    real(field_r), pointer :: tker_col(:), shr_col(:), buo_col(:), trsp_col(:), ptrsp_col(:)
    real(field_r), pointer :: sbtke_col(:), sbshr_col(:), sbbuo_col(:), sbdiss_col(:)
    real(field_r), pointer :: sbstor_col(:), sbbudg_col(:), sbresid_col(:), ekm_col(:), khkm_col(:)
    ! radstat variables
    real(field_r), pointer :: thltend_col(:), thllwtend_col(:), thlswtend_col(:), thlradls_col(:)
    real(field_r), pointer :: lwu_col(:), lwd_col(:), swu_col(:), swd_col(:)
    real(field_r), pointer :: lwuca_col(:), lwdca_col(:), swuca_col(:), swdca_col(:)
    real(field_r), pointer :: thllwtendca_col(:), thlswtendca_col(:)

    if (.not. lcolstat) return
    if (.not. allocated(i_local)) return

    if (.not. is_sampling_timestep(ofile_id)) return
    ! ----------------------------------------------------------------
    ! Accumulation pass — every dtav, add instantaneous column values
    ! ----------------------------------------------------------------
    do col_idx = 1, size(i_local)

      i = i_local(col_idx)
      j = j_local(col_idx)
      weres_loc = 0.0_field_r
      ptrsp_loc = 0.0_field_r
      buo_loc = 0.0_field_r
      trsp_loc = 0.0_field_r
      shr_loc = 0.0_field_r

      call get_pointers(col_idx)


      do k = 1, kmax
        u_col(k)     = u_col(k)    + u0(i,j,k)
        v_col(k)     = v_col(k)    + v0(i,j,k)
        w_col(k)     = w_col(k)    + w0(i,j,k)
        thl_col(k)   = thl_col(k)  + thl0(i,j,k)
        qt_col(k)    = qt_col(k)   + qt0(i,j,k)
        ql_col(k)    = ql_col(k)   + ql0(i,j,k)
        ta_col(k)    = ta_col(k)   + tmp0(i,j,k)
        thv_col(k)   = thv_col(k)  + thv0(i,j,k)

        presh_col(k)  = presh_col(k)  + presh(k)
        rhof_col(k)   = rhof_col(k)   + rhof(k)
        rhobf_col(k)  = rhobf_col(k)  + rhobf(k)
        rhobh_col(k)  = rhobh_col(k)  + rhobh(k)

        ! Moments — use local instantaneous values for consistent skewness
        u2r_inst = (um(i,j,k) + cu - umav(k))**2
        v2r_inst = (vm(i,j,k) + cv - vmav(k))**2
        w_prime  = wm(i,j,k) - wmav(k)
        w2r_inst = w_prime**2
        u2r_col(k)   = u2r_col(k)   + u2r_inst
        v2r_col(k)   = v2r_col(k)   + v2r_inst
        w2r_col(k)   = w2r_col(k)   + w2r_inst
        w2s_col(k)   = w2s_col(k)   + e12m(i,j,k)**2
        thl2r_col(k) = thl2r_col(k) + (thlm(i,j,k) - thlmav(k))**2
        thv2r_col(k) = thv2r_col(k) + (thv0(i,j,k) - thvmav(k))**2
        ! Match modgenstat: th2 uses thlm moments around thmav.
        th2r_col(k)  = th2r_col(k)  + (thlm(i,j,k) - thmav(k))**2
        qt2r_col(k)  = qt2r_col(k)  + (qtm(i,j,k) - qtmav(k))**2
        ql2r_col(k)  = ql2r_col(k)  + (ql0(i,j,k) - qlmav(k))**2

        ! Accumulate w'^3 and normalize later with slab w2av, matching modgenstat.
        skew_col(k) = skew_col(k) + w_prime**3

        if (k == 1) then
          wthls_loc = thlflux(i,j)
          wthlr_loc = 0.0_field_r
          wqts_loc  = qtflux(i,j)
          wqtr_loc  = 0.0_field_r
          wqls_loc  = 0.0_field_r
          wqlr_loc  = 0.0_field_r

          ! qls is 0 as in modgenstat
          tsurf = thls * exnh(1)
          qsat  = qts

          a_dry   = 1.0_field_r + (rv/rd - 1.0_field_r) * qts
          b_dry   = rv/rd - 1.0_field_r

          a_moist = (1.0_field_r - qts + rv/rd * qsat * (1.0_field_r + rlv / (rv * tsurf))) / &
                    (1.0_field_r + rlv / (rv * tsurf) * rlv / (cp * tsurf) * qsat)
          b_moist = a_moist * rlv / (tsurf * cp) - 1.0_field_r

          c1 = a_dry
          c2 = b_dry

          wthvs_loc = c1 * thlflux(i,j) + c2 * thls * qtflux(i,j)
          wthvr_loc = 0.0_field_r

          upcu = um(i,j,1) + cu
          upcu = sign(1.0_field_r, upcu) * max(abs(upcu), eps1)
          uws_loc = - (0.5_field_r * (ustar(i,j) + ustar(i-1,j)))**2 * upcu / &
                    sqrt(upcu**2 + ((vm(i,j,1) + vm(i-1,j,1) + vm(i,j+1,1) + vm(i-1,j+1,1)) / 4.0_field_r + cv)**2)
          uwr_loc = 0.0_field_r

          vpcv = vm(i,j,1) + cv
          vpcv = sign(1.0_field_r, vpcv) * max(abs(vpcv), eps1)
          vws_loc = - (0.5_field_r * (ustar(i,j) + ustar(i,j-1)))**2 * vpcv / &
                    sqrt(vpcv**2 + ((um(i,j,1) + um(i+1,j,1) + um(i,j-1,1) + um(i+1,j-1,1)) / 4.0_field_r + cu)**2)
          vwr_loc = 0.0_field_r
        else
          !-----------------------------------------------------------
          ! Calculate prefactors for subgrid wthv and wql fluxes
          ! at half levels
          !-----------------------------------------------------------
          qs0h = qt0h(i,j,k) - ql0h(i,j,k)
          t0h  = exnh(k) * thl0h(i,j,k) + (rlv/cp) * ql0h(i,j,k)

          den  = 1.0_field_r + (rlv**2) * qs0h / (rv * cp * (t0h**2))
          cthl = (exnh(k) * cp / rlv) * ((1.0_field_r - den) / den)
          cqt  = 1.0_field_r / den

          a_dry   = 1.0_field_r + (rv/rd - 1.0_field_r) * qt0h(i,j,k)
          b_dry   = rv/rd - 1.0_field_r
          a_moist = (1.0_field_r - qt0h(i,j,k) + rv/rd * qs0h * (1.0_field_r + rd/rv * rlv / (rd * t0h))) / den
          b_moist = a_moist * rlv / (t0h * cp) - 1.0_field_r

          c1 = merge(a_moist, a_dry, ql0h(i,j,k) > 0.0_field_r)
          c2 = merge(b_moist, b_dry, ql0h(i,j,k) > 0.0_field_r)

          !-----------------------------------------------------------
          ! Calculate resolved and subgrid fluxes at half levels
          !-----------------------------------------------------------
          ekhalf = (ekh(i,j,k)*dzf(k-1)+ekh(i,j,k-1)*dzf(k))/(2*dzh(k))
          euhalf = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i-1,j,k  ) )  + &
                     dzf(k  ) * ( ekm(i,j,k-1) + ekm(i-1,j,k-1) ) ) * &
                      ( 0.25 * dzhi(k) )
          evhalf = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i,j-1,k  ) )  + &
                     dzf(k  ) * ( ekm(i,j,k-1) + ekm(i,j-1,k-1) ) ) * &
                      ( 0.25 * dzhi(k) )

          wthls_loc = -ekhalf * (thl0(i,j,k) - thl0(i,j,k-1)) * dzhi(k)
          wthlr_loc = (w0(i,j,k) - wmav(k)) * thl0h(i,j,k)

          wqts_loc  = -ekhalf * (qt0(i,j,k) - qt0(i,j,k-1)) * dzhi(k)
          wqtr_loc  = (w0(i,j,k) - wmav(k)) * qt0h(i,j,k)

          wqls_loc  = 0.0_field_r
          if (ql0h(i,j,k) > 0.0_field_r) wqls_loc = cthl * wthls_loc + cqt * wqts_loc
          wqlr_loc  = (w0(i,j,k) - wmav(k)) * ql0h(i,j,k)

          wthvs_loc = c1 * wthls_loc + c2 * thl0h(i,j,k) * wqts_loc
          wthvr_loc = (w0(i,j,k) - wmav(k)) * thv0h(i,j,k)

          uwr_loc = (w0(i,j,k)+w0(i-1,j,k)-2*wmav(k)) &
                    *((u0(i,j,k-1)+cu)*dzf(k)+(u0(i,j,k)+cu)*dzf(k-1))*(0.25*dzhi(k))
          vwr_loc = (w0(i,j,k)+w0(i,j-1,k)-2*wmav(k)) &
                    *((v0(i,j,k-1)+cv)*dzf(k)+(v0(i,j,k)+cv)*dzf(k-1))*(0.25*dzhi(k))
          uws_loc = -euhalf &
                    *((u0(i,j,k)-u0(i,j,k-1))/dzh(k)+(w0(i,j,k)-w0(i-1,j,k))*dxi)
          vws_loc = -evhalf &
                    *((v0(i,j,k)-v0(i,j,k-1))/dzh(k)+(w0(i,j,k)-w0(i,j-1,k))*dyi)
        end if

        wthls_col(k) = wthls_col(k) + wthls_loc
        wthlr_col(k) = wthlr_col(k) + wthlr_loc
        wthlt_col(k) = wthlt_col(k) + wthls_loc + wthlr_loc
        wthvs_col(k) = wthvs_col(k) + wthvs_loc
        wthvr_col(k) = wthvr_col(k) + wthvr_loc
        wthvt_col(k) = wthvt_col(k) + wthvs_loc + wthvr_loc
        wqts_col(k)  = wqts_col(k)  + wqts_loc
        wqtr_col(k)  = wqtr_col(k)  + wqtr_loc
        wqtt_col(k)  = wqtt_col(k)  + wqts_loc  + wqtr_loc
        wqls_col(k)  = wqls_col(k)  + wqls_loc
        wqlr_col(k)  = wqlr_col(k)  + wqlr_loc
        wqlt_col(k)  = wqlt_col(k)  + wqls_loc  + wqlr_loc
        uws_col(k)   = uws_col(k)   + uws_loc
        uwr_col(k)   = uwr_col(k)   + uwr_loc
        uwt_col(k)   = uwt_col(k)   + uws_loc   + uwr_loc
        vws_col(k)   = vws_col(k)   + vws_loc
        vwr_col(k)   = vwr_col(k)   + vwr_loc
        vwt_col(k)   = vwt_col(k)   + vws_loc   + vwr_loc


        cs_col(k)    = cs_col(k)    + csz(k)
        hus_col(k)   = hus_col(k)   + (qt0(i,j,k) - ql0(i,j,k))
        hur_col(k)   = hur_col(k)   + 100.0_field_r * (qt0(i,j,k) - ql0(i,j,k)) / qsat_tab(tmp0(i,j,k), presf(k))

        ilratio = max(0._field_r,min(1._field_r,(tmp0(i,j,k)-tdn) / (tup-tdn)))
        clw_col(k)   = clw_col(k)   + ql0(i,j,k) * ilratio
        cli_col(k)   = cli_col(k)   + ql0(i,j,k) * (1.0_field_r - ilratio)
        if (ql0(i,j,k) > 0.0_field_r) cfrac_col(k) = cfrac_col(k) + 1.0_field_r
        
        iqr = get_tracer_index("qr")
        if (iqr > 0) then
          if (imicro == imicro_sice .or. imicro == imicro_sice2) then
              ilratio = max(0._field_r,min(1._field_r,(tmp0(i,j,k)-tdnrsg)/(tuprsg-tdnrsg)))
              plw_col(k) = plw_col(k) + sv0(i,j,k,iqr) * ilratio
              pli_col(k) = pli_col(k) + sv0(i,j,k,iqr) * (1-ilratio)
          else
            plw_col(k) = plw_col(k) + sv0(i,j,k,iqr)
          end if
        end if

        do n = 1, nsv
          call get_sv_pointer(n)
          sv_col(k)   = sv_col(k)   + svm(i,j,k,n)
          svp_col(k)  = svp_col(k)  + svp(i,j,k,n)
          sv2r_col(k) = sv2r_col(k) + (svm(i,j,k,n) - svmav(k,n))**2
          if (k == 1) then
            wsvs_loc = svflux(i,j,n)
            wsvr_loc = 0.0_field_r
          else
            ekhalf   = (ekh(i,j,k)*dzf(k-1) + ekh(i,j,k-1)*dzf(k)) / (2.0_field_r * dzh(k))
            wsvs_loc = -ekhalf * (sv0(i,j,k,n) - sv0(i,j,k-1,n)) / dzh(k)
            wsvr_loc = (w0(i,j,k) - wmav(k)) * ((sv0(i,j,k,n)*dzf(k-1) + sv0(i,j,k-1,n)*dzf(k)) / (2.0_field_r*dzh(k)))
          end if
          wsvs_col(k) = wsvs_col(k) + wsvs_loc
          wsvr_col(k) = wsvr_col(k) + wsvr_loc
          wsvt_col(k) = wsvt_col(k) + wsvs_loc + wsvr_loc
        end do

        ! Budget terms
        if (k == 1) then
          ! no shear, buoyancy or pressure transport at the surface
          buo_loc(k)   = 0.0_field_r
          shr_loc(k)   = 0.0_field_r
          ptrsp_loc(k) = 0.0_field_r
        else
          ! buoyancy
          buo_loc(k) = rhobh(k) * grav / thvh(k) * w0(i,j,k) * (thv0h(i,j,k) - thvh(k))
          
          ! shear
          uwrs = rhobh(k) * (w0(i,j,k) + w0(i-1,j,k)) * (u0(i,j,k-1) + u0(i,j,k)) * 0.25_field_r
          vwrs = rhobh(k) * (w0(i,j,k) + w0(i,j-1,k)) * (v0(i,j,k-1) + v0(i,j,k)) * 0.25_field_r
          shr_loc(k) = -uwrs * (u0av(k) - u0av(k-1)) / dzh(k) - vwrs * (v0av(k) - v0av(k-1)) / dzh(k)
          
          ! pressure transport
          ph = (p(i,j,k)*dzf(k-1) + p(i,j,k-1)*dzf(k)) / (2.0_field_r * dzh(k))
          ptrsp_loc(k) = rhobh(k) * w0(i,j,k) * ph
        end if

        ! Accumulate budget terms (method from modbudget)
        shr_col(k)   = shr_col(k)   + shr_loc(k)
        buo_col(k)   = buo_col(k)   + buo_loc(k)
        ptrsp_col(k) = ptrsp_col(k) + ptrsp_loc(k)

        egp = 0.5*rhobf(k)*( (0.5*(u0(i,j,k)+u0(i+1,j,k))-(u0av(k)-cu))**2 &
                  +(0.5*(v0(i,j,k)+v0(i,j+1,k))-(v0av(k)-cv))**2 &
                  +(0.5*(w0(i,j,k)+w0(i,j,k+1))             )**2 )
        weres_loc(k) = egp * 0.5_field_r * (w0(i,j,k) + w0(i,j,k+1))
        tker_col(k)  = tker_col(k) + egp
        
        sbtke_inst = e120(i,j,k)**2 * rhobf(k)
        sbtke_col(k)  = sbtke_col(k) + sbtke_inst
        if (.not. l_sbtke_beg_set) sbtke_beg_col(k,col_idx) = sbtke_inst
        sbtke_last_col(k,col_idx) = sbtke_inst
        ekm_col(k)    = ekm_col(k)    + ekm(i,j,k)
        if (ekm(i,j,k) > eps1) then
          khkm_col(k) = khkm_col(k) + ekh(i,j,k) / ekm(i,j,k)
        end if



        thltend_col(k) = thltend_col(k) + thlprad(i,j,k)
        
        !absolute values to handle different sign conventions in radiation code
        thllwtend_col(k) = thllwtend_col(k) + &
            (abs(lwd(i,j,k+1)) - abs(lwu(i,j,k+1)) - abs(lwd(i,j,k)) + abs(lwu(i,j,k))) / &
            (rhof(k)*exnf(k)*cp*dzf(k))
        thlswtend_col(k) = thlswtend_col(k) + &
            (abs(swd(i,j,k+1)) - abs(swu(i,j,k+1)) - abs(swd(i,j,k)) + abs(swu(i,j,k))) / &
            (rhof(k)*exnf(k)*cp*dzf(k))

        thlradls_col(k) = thlradls_col(k) + thlpcar(k)

        lwu_col(k)  = lwu_col(k)  + abs(lwu(i,j,k))
        lwd_col(k)  = lwd_col(k)  + abs(lwd(i,j,k))
        swu_col(k)  = swu_col(k)  + abs(swu(i,j,k))
        swd_col(k)  = swd_col(k)  + abs(swd(i,j,k))

        ! assume clear-sky fluxes are already calculated in modradstat and available as lwuca/lwdca/swuca/swdca; if not, these will just accumulate zeros and not contribute to the final averages
        lwuca_col(k) = lwuca_col(k) + abs(lwuca(i,j,k))
        lwdca_col(k) = lwdca_col(k) + abs(lwdca(i,j,k))
        swuca_col(k) = swuca_col(k) + abs(swuca(i,j,k))
        swdca_col(k) = swdca_col(k) + abs(swdca(i,j,k))

        thllwtendca_col(k) = thllwtendca_col(k) + &
            (-lwdca(i,j,k+1) - lwuca(i,j,k+1) + lwdca(i,j,k) + lwuca(i,j,k)) / &
            (rhof(k)*exnf(k)*cp*dzf(k))
        thlswtendca_col(k) = thlswtendca_col(k) + &
            (-swdca(i,j,k+1) - swuca(i,j,k+1) + swdca(i,j,k) + swuca(i,j,k)) / &
            (rhof(k)*exnf(k)*cp*dzf(k))

      end do

      ! Post-loop: resolved transport (divergence-based)
      trsp_col(1) = trsp_col(1) - weres_loc(1) / (0.5_field_r * dzh(1))
      do k = 2, kmax
        trsp_col(k) = trsp_col(k) - (weres_loc(k) - weres_loc(k-1)) / dzh(k)
      end do
      
      ! Note: ptrsp_loc contains raw flux values accumulated into ptrsp_col.
      ! During write phase (after averaging), will convert ptrsp_col to divergence form.
      ! This matches modbudget's final transformation: ptrsp = -(d/dz)(ptrsp_flux)

      k = 1
      !sbshr, sbbuo, sbdiss are only defined at k=1 (surface layer)
      !These are multiplied by rhobf during accumulation per modbudget convention
      sbshr_col(k)  = sbshr_col(k)  + sbshr(i,j,k) * rhobf(k)
      sbbuo_col(k)  = sbbuo_col(k)  + sbbuo(i,j,k) * rhobf(k)
      sbdiss_col(k) = sbdiss_col(k) + sbdiss(i,j,k) * rhobf(k)
      ! sbbudg is computed from components: sum of shear, buoyancy, dissipation  
      ! Will be calculated during output normalization from accumulated sbshr+sbbuo+sbdiss

    end do  ! accumulation loop

    if (.not. l_sbtke_beg_set) l_sbtke_beg_set = .true.

    ! ----------------------------------------------------------------
    ! Write pass — every timeav: average, fill file buffers, and zero
    ! ----------------------------------------------------------------
    if (is_writing_timestep(ofile_id)) then
      scale_rn = 1.0_field_r / real(nsamples, kind=field_r)

      do col_idx = 1, size(i_local)

        call get_pointers(col_idx)

        rhof_col(:)   = rhof_col(:)   * scale_rn
        rhobf_col(:)  = rhobf_col(:)  * scale_rn
        rhobh_col(:)  = rhobh_col(:)  * scale_rn
        presh_col(:)  = presh_col(:)  * scale_rn
        u_col(:)      = u_col(:)      * scale_rn
        v_col(:)      = v_col(:)      * scale_rn
        w_col(:)      = w_col(:)      * scale_rn
        thl_col(:)    = thl_col(:)    * scale_rn
        thv_col(:)    = thv_col(:)    * scale_rn
        qt_col(:)     = qt_col(:)     * scale_rn
        ql_col(:)     = ql_col(:)     * scale_rn
        wthls_col(:)  = wthls_col(:)  * scale_rn
        wthlr_col(:)  = wthlr_col(:)  * scale_rn
        wthlt_col(:)  = wthlt_col(:)  * scale_rn
        wthvs_col(:)  = wthvs_col(:)  * scale_rn
        wthvr_col(:)  = wthvr_col(:)  * scale_rn
        wthvt_col(:)  = wthvt_col(:)  * scale_rn
        wqts_col(:)   = wqts_col(:)   * scale_rn
        wqtr_col(:)   = wqtr_col(:)   * scale_rn
        wqtt_col(:)   = wqtt_col(:)   * scale_rn
        wqls_col(:)   = wqls_col(:)   * scale_rn
        wqlr_col(:)   = wqlr_col(:)   * scale_rn
        wqlt_col(:)   = wqlt_col(:)   * scale_rn
        uws_col(:)    = uws_col(:)    * scale_rn
        uwr_col(:)    = uwr_col(:)    * scale_rn
        uwt_col(:)    = uwt_col(:)    * scale_rn
        vws_col(:)    = vws_col(:)    * scale_rn
        vwr_col(:)    = vwr_col(:)    * scale_rn
        vwt_col(:)    = vwt_col(:)    * scale_rn
        w2s_col(:)    = w2s_col(:)    * scale_rn
        w2r_col(:)    = w2r_col(:)    * scale_rn
        skew_col(:)   = skew_col(:)   * scale_rn
        do k = 1, kmax
          skew_col(k) = skew_col(k) / max(w2av(k)**1.5_field_r, epsilon(1.0_field_r))
        end do
        u2r_col(:)    = u2r_col(:)    * scale_rn
        v2r_col(:)    = v2r_col(:)    * scale_rn
        thl2r_col(:)  = thl2r_col(:)  * scale_rn
        thv2r_col(:)  = thv2r_col(:)  * scale_rn
        th2r_col(:)   = th2r_col(:)   * scale_rn
        qt2r_col(:)   = qt2r_col(:)   * scale_rn
        ql2r_col(:)   = ql2r_col(:)   * scale_rn
        cs_col(:)     = cs_col(:)     * scale_rn
        cfrac_col(:)  = cfrac_col(:)  * scale_rn
        hur_col(:)    = hur_col(:)    * scale_rn
        hus_col(:)    = hus_col(:)    * scale_rn
        ta_col(:)     = ta_col(:)     * scale_rn
        clw_col(:)    = clw_col(:)    * scale_rn
        cli_col(:)    = cli_col(:)    * scale_rn
        plw_col(:)    = plw_col(:)    * scale_rn
        pli_col(:)    = pli_col(:)    * scale_rn
        tker_col(:)   = tker_col(:)   * scale_rn
        shr_col(:)    = shr_col(:)    * scale_rn
        buo_col(:)    = buo_col(:)    * scale_rn
        trsp_col(:)   = trsp_col(:)   * scale_rn
        ptrsp_col(:)  = ptrsp_col(:)  * scale_rn

        ! Match modbudget staggering: report full-level values from adjacent half levels.
        prof_tmp = shr_col(:)
        do k = 1, kmax-1
          shr_col(k) = 0.5_field_r * (prof_tmp(k) + prof_tmp(k+1))
        end do
        shr_col(kmax) = prof_tmp(kmax)

        prof_tmp = buo_col(:)
        do k = 1, kmax-1
          buo_col(k) = 0.5_field_r * (prof_tmp(k) + prof_tmp(k+1))
        end do
        buo_col(kmax) = prof_tmp(kmax)

        prof_tmp = trsp_col(:)
        do k = 1, kmax-1
          trsp_col(k) = 0.5_field_r * (prof_tmp(k) + prof_tmp(k+1))
        end do
        trsp_col(kmax) = prof_tmp(kmax)
        
        ! Convert ptrsp from accumulated raw flux to divergence form (matching modbudget):
        ! ptrsp_final(k) = -(ptrsp_flux(k+1) - ptrsp_flux(k)) / dzf(k)
        prof_tmp = ptrsp_col(:)
        ptrsp_col(1) = 0.0_field_r  ! ptrsp(1) = 0 (no divergence at surface)
        do k = 2, kmax-1
          ptrsp_col(k) = -(prof_tmp(k+1) - prof_tmp(k)) / dzf(k)
        end do
        ptrsp_col(kmax) = 0.0_field_r
        
        sbtke_col(:)  = sbtke_col(:)  * scale_rn
        sbshr_col(:)  = sbshr_col(:)  * scale_rn
        sbbuo_col(:)  = sbbuo_col(:)  * scale_rn
        sbdiss_col(:) = sbdiss_col(:) * scale_rn
        ! sbbudg = sum of shear + buoyancy + dissipation (per modbudget convention)
        sbbudg_col(:) = sbshr_col(:) + sbbuo_col(:) + sbdiss_col(:)
        sbstor_col(:)  = (sbtke_last_col(:,col_idx) - sbtke_beg_col(:,col_idx)) / max(timeav, eps1)
        sbresid_col(:) = sbbudg_col(:) - sbstor_col(:)
        ekm_col(:)    = ekm_col(:)    * scale_rn
        khkm_col(:)   = khkm_col(:)   * scale_rn
        thltend_col(:)    = thltend_col(:)    * scale_rn
        thllwtend_col(:)  = thllwtend_col(:)  * scale_rn
        thlswtend_col(:)  = thlswtend_col(:)  * scale_rn
        thlradls_col(:)   = thlradls_col(:)   * scale_rn
        lwu_col(:)    = lwu_col(:)    * scale_rn
        lwd_col(:)    = lwd_col(:)    * scale_rn
        swu_col(:)    = swu_col(:)    * scale_rn
        swd_col(:)    = swd_col(:)    * scale_rn
        lwuca_col(:)      = lwuca_col(:)      * scale_rn
        lwdca_col(:)      = lwdca_col(:)      * scale_rn
        swuca_col(:)      = swuca_col(:)      * scale_rn
        swdca_col(:)      = swdca_col(:)      * scale_rn
        thllwtendca_col(:) = thllwtendca_col(:) * scale_rn
        thlswtendca_col(:) = thlswtendca_col(:) * scale_rn
        if (nsv > 0) then
          do n = 1, nsv
            call get_sv_pointer(n)
            ! get pointers here
            sv_col(:)   = sv_col(:)   * scale_rn
            svp_col(:)  = svp_col(:)  * scale_rn
            sv2r_col(:) = sv2r_col(:) * scale_rn
            wsvs_col(:) = wsvs_col(:) * scale_rn
            wsvr_col(:) = wsvr_col(:) * scale_rn
            wsvt_col(:) = wsvt_col(:) * scale_rn
          end do
        end if
      end do
    end if

    ! zeroing of buffer accumulators is done in modnetcdf_file_t

    contains


    subroutine get_pointers(col_idx)
      integer, intent(in) :: col_idx
      integer :: n
      character(len=64) :: vname

      call ofile%get_pointer('rhof', global_index(col_idx), rhof_col)
      call ofile%get_pointer('rhobf', global_index(col_idx), rhobf_col)
      call ofile%get_pointer('rhobh', global_index(col_idx), rhobh_col)
      call ofile%get_pointer('presh', global_index(col_idx), presh_col)
      call ofile%get_pointer('u', global_index(col_idx), u_col)
      call ofile%get_pointer('v', global_index(col_idx), v_col)
      call ofile%get_pointer('w', global_index(col_idx), w_col)
      call ofile%get_pointer('thl', global_index(col_idx), thl_col)
      call ofile%get_pointer('thv', global_index(col_idx), thv_col)
      call ofile%get_pointer('qt', global_index(col_idx), qt_col)
      call ofile%get_pointer('ql', global_index(col_idx), ql_col)
      call ofile%get_pointer('wthls', global_index(col_idx), wthls_col)
      call ofile%get_pointer('wthlr', global_index(col_idx), wthlr_col)
      call ofile%get_pointer('wthlt', global_index(col_idx), wthlt_col)
      call ofile%get_pointer('wthvs', global_index(col_idx), wthvs_col)
      call ofile%get_pointer('wthvr', global_index(col_idx), wthvr_col)
      call ofile%get_pointer('wthvt', global_index(col_idx), wthvt_col)
      call ofile%get_pointer('wqts', global_index(col_idx), wqts_col)
      call ofile%get_pointer('wqtr', global_index(col_idx), wqtr_col)
      call ofile%get_pointer('wqtt', global_index(col_idx), wqtt_col)
      call ofile%get_pointer('wqls', global_index(col_idx), wqls_col)
      call ofile%get_pointer('wqlr', global_index(col_idx), wqlr_col)
      call ofile%get_pointer('wqlt', global_index(col_idx), wqlt_col)
      call ofile%get_pointer('uws', global_index(col_idx), uws_col)
      call ofile%get_pointer('uwr', global_index(col_idx), uwr_col)
      call ofile%get_pointer('uwt', global_index(col_idx), uwt_col)
      call ofile%get_pointer('vws', global_index(col_idx), vws_col)
      call ofile%get_pointer('vwr', global_index(col_idx), vwr_col)
      call ofile%get_pointer('vwt', global_index(col_idx), vwt_col)
      call ofile%get_pointer('w2s', global_index(col_idx), w2s_col)
      call ofile%get_pointer('w2r', global_index(col_idx), w2r_col)
      call ofile%get_pointer('skew', global_index(col_idx), skew_col)
      call ofile%get_pointer('u2r', global_index(col_idx), u2r_col)
      call ofile%get_pointer('v2r', global_index(col_idx), v2r_col)
      call ofile%get_pointer('thl2r', global_index(col_idx), thl2r_col)
      call ofile%get_pointer('thv2r', global_index(col_idx), thv2r_col)
      call ofile%get_pointer('th2r', global_index(col_idx), th2r_col)
      call ofile%get_pointer('qt2r', global_index(col_idx), qt2r_col)
      call ofile%get_pointer('ql2r', global_index(col_idx), ql2r_col)
      call ofile%get_pointer('cs', global_index(col_idx), cs_col)
      call ofile%get_pointer('cfrac', global_index(col_idx), cfrac_col)
      call ofile%get_pointer('hur', global_index(col_idx), hur_col)
      call ofile%get_pointer('hus', global_index(col_idx), hus_col)
      call ofile%get_pointer('ta', global_index(col_idx), ta_col)
      call ofile%get_pointer('clw', global_index(col_idx), clw_col)
      call ofile%get_pointer('cli', global_index(col_idx), cli_col)
      call ofile%get_pointer('plw', global_index(col_idx), plw_col)
      call ofile%get_pointer('pli', global_index(col_idx), pli_col)

      call ofile%get_pointer('tker', global_index(col_idx), tker_col)
      call ofile%get_pointer('shr', global_index(col_idx), shr_col)
      call ofile%get_pointer('buo', global_index(col_idx), buo_col)
      call ofile%get_pointer('trsp', global_index(col_idx), trsp_col)
      call ofile%get_pointer('ptrsp', global_index(col_idx), ptrsp_col)
      call ofile%get_pointer('sbtke', global_index(col_idx), sbtke_col)
      call ofile%get_pointer('sbshr', global_index(col_idx), sbshr_col)
      call ofile%get_pointer('sbbuo', global_index(col_idx), sbbuo_col)
      call ofile%get_pointer('sbdiss', global_index(col_idx), sbdiss_col)
      call ofile%get_pointer('sbstor', global_index(col_idx), sbstor_col)
      call ofile%get_pointer('sbbudg', global_index(col_idx), sbbudg_col)
      call ofile%get_pointer('sbresid', global_index(col_idx), sbresid_col)
      call ofile%get_pointer('ekm', global_index(col_idx), ekm_col)
      call ofile%get_pointer('khkm', global_index(col_idx), khkm_col)
      call ofile%get_pointer('thltend', global_index(col_idx), thltend_col)
      call ofile%get_pointer('thllwtend', global_index(col_idx), thllwtend_col)
      call ofile%get_pointer('thlswtend', global_index(col_idx), thlswtend_col)
      call ofile%get_pointer('thlradls', global_index(col_idx), thlradls_col)
      call ofile%get_pointer('lwu', global_index(col_idx), lwu_col)
      call ofile%get_pointer('lwd', global_index(col_idx), lwd_col)
      call ofile%get_pointer('swu', global_index(col_idx), swu_col)
      call ofile%get_pointer('swd', global_index(col_idx), swd_col)
      call ofile%get_pointer('lwuca', global_index(col_idx), lwuca_col)
      call ofile%get_pointer('lwdca', global_index(col_idx), lwdca_col)
      call ofile%get_pointer('swuca', global_index(col_idx), swuca_col)
      call ofile%get_pointer('swdca', global_index(col_idx), swdca_col)
      call ofile%get_pointer('thllwtendca', global_index(col_idx), thllwtendca_col)
      call ofile%get_pointer('thlswtendca', global_index(col_idx), thlswtendca_col)

    end subroutine get_pointers

    subroutine get_sv_pointer(n)
      character(len=64) :: vname
      integer, intent(in) :: n

      vname = trim(tracer_prop(n)%tracname)
      call ofile%get_pointer(vname, global_index(col_idx), sv_col)
      call ofile%get_pointer(trim(vname)//'p', global_index(col_idx), svp_col)
      call ofile%get_pointer(trim(vname)//'2r', global_index(col_idx), sv2r_col)
      call ofile%get_pointer('w'//trim(vname)//'s', global_index(col_idx), wsvs_col)
      call ofile%get_pointer('w'//trim(vname)//'r', global_index(col_idx), wsvr_col)
      call ofile%get_pointer('w'//trim(vname)//'t', global_index(col_idx), wsvt_col)

    end subroutine get_sv_pointer
  end subroutine colstat


  subroutine exitcolstat

    if (allocated(i_local)) deallocate(i_local)
    if (allocated(j_local)) deallocate(j_local)
    if (allocated(global_index)) deallocate(global_index)
    if (allocated(locx)) deallocate(locx)
    if (allocated(locy)) deallocate(locy)
    if (allocated(sbtke_beg_col)) deallocate(sbtke_beg_col)
    if (allocated(sbtke_last_col)) deallocate(sbtke_last_col)

  end subroutine exitcolstat

end module modcolstat
