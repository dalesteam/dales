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
  integer :: npoints = 0
  integer :: x_idx(max_points) = 0
  integer :: y_idx(max_points) = 0

  type(multi_profile_file_t) :: ofile
  integer :: ofile_id
  integer, allocatable :: i_local(:), j_local(:), global_index(:)
  real(field_r), allocatable :: locx(:), locy(:)

  ! Per-profile storage, dimensioned as (kmax, nlocal, nvar) unless noted.
  real(field_r), allocatable :: rhof_col(:,:), rhobf_col(:,:), rhobh_col(:,:), presh_col(:,:)
  real(field_r), allocatable :: u_col(:,:), v_col(:,:), w_col(:,:), thl_col(:,:), thv_col(:,:), qt_col(:,:), ql_col(:,:)
  real(field_r), allocatable :: wthls_col(:,:), wthlr_col(:,:), wthlt_col(:,:), wthvs_col(:,:), wthvr_col(:,:), wthvt_col(:,:)
  real(field_r), allocatable :: wqts_col(:,:), wqtr_col(:,:), wqtt_col(:,:), wqls_col(:,:), wqlr_col(:,:), wqlt_col(:,:)
  real(field_r), allocatable :: uws_col(:,:), uwr_col(:,:), uwt_col(:,:), vws_col(:,:), vwr_col(:,:), vwt_col(:,:)
  real(field_r), allocatable :: w2s_col(:,:), w2r_col(:,:), skew_col(:,:), u2r_col(:,:), v2r_col(:,:)
  real(field_r), allocatable :: thl2r_col(:,:), thv2r_col(:,:), th2r_col(:,:), qt2r_col(:,:), ql2r_col(:,:)
  real(field_r), allocatable :: cs_col(:,:), cfrac_col(:,:), hur_col(:,:), hus_col(:,:), ta_col(:,:)
  real(field_r), allocatable :: clw_col(:,:), cli_col(:,:), plw_col(:,:), pli_col(:,:)

  real(field_r), allocatable :: sv_col(:,:,:), svp_col(:,:,:), sv2r_col(:,:,:)
  real(field_r), allocatable :: wsvs_col(:,:,:), wsvr_col(:,:,:), wsvt_col(:,:,:)

  real(field_r), allocatable :: tker_col(:,:), shr_col(:,:), buo_col(:,:), trsp_col(:,:), ptrsp_col(:,:)
  real(field_r), allocatable :: sbtke_col(:,:), sbshr_col(:,:), sbbuo_col(:,:), sbdiss_col(:,:)
  real(field_r), allocatable :: sbstor_col(:,:), sbbudg_col(:,:), sbresid_col(:,:), ekm_col(:,:), khkm_col(:,:)
  real(field_r), allocatable :: sbtke_beg_col(:,:), sbtke_last_col(:,:)

  real(field_r), allocatable :: thltend_col(:,:), thllwtend_col(:,:), thlswtend_col(:,:), thlradls_col(:,:)
  real(field_r), allocatable :: lwu_col(:,:), lwd_col(:,:), swu_col(:,:), swd_col(:,:)
  real(field_r), allocatable :: lwuca_col(:,:), lwdca_col(:,:), swuca_col(:,:), swdca_col(:,:)
  real(field_r), allocatable :: thllwtendca_col(:,:), thlswtendca_col(:,:)

  ! genstat variables
  ! modbudget variables
  ! radstat variables
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
    call allocate_profile_arrays(nlocal)
    call zero_profile_arrays()

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
    logical :: any_write

    if (.not. lcolstat) return
    if (.not. allocated(i_local)) return

    any_write = .false.

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


      do k = 1, kmax
        u_col(k,col_idx)     = u_col(k,col_idx)    + u0(i,j,k)
        v_col(k,col_idx)     = v_col(k,col_idx)    + v0(i,j,k)
        w_col(k,col_idx)     = w_col(k,col_idx)    + w0(i,j,k)
        thl_col(k,col_idx)   = thl_col(k,col_idx)  + thl0(i,j,k)
        qt_col(k,col_idx)    = qt_col(k,col_idx)   + qt0(i,j,k)
        ql_col(k,col_idx)    = ql_col(k,col_idx)   + ql0(i,j,k)
        ta_col(k,col_idx)    = ta_col(k,col_idx)   + tmp0(i,j,k)
        thv_col(k,col_idx)   = thv_col(k,col_idx)  + thv0(i,j,k)

        presh_col(k,col_idx)  = presh_col(k,col_idx)  + presh(k)
        rhof_col(k,col_idx)   = rhof_col(k,col_idx)   + rhof(k)
        rhobf_col(k,col_idx)  = rhobf_col(k,col_idx)  + rhobf(k)
        rhobh_col(k,col_idx)  = rhobh_col(k,col_idx)  + rhobh(k)

        ! Moments — use local instantaneous values for consistent skewness
        u2r_inst = (um(i,j,k) + cu - umav(k))**2
        v2r_inst = (vm(i,j,k) + cv - vmav(k))**2
        w_prime  = wm(i,j,k) - wmav(k)
        w2r_inst = w_prime**2
        u2r_col(k,col_idx)   = u2r_col(k,col_idx)   + u2r_inst
        v2r_col(k,col_idx)   = v2r_col(k,col_idx)   + v2r_inst
        w2r_col(k,col_idx)   = w2r_col(k,col_idx)   + w2r_inst
        w2s_col(k,col_idx)   = w2s_col(k,col_idx)   + e12m(i,j,k)**2
        thl2r_col(k,col_idx) = thl2r_col(k,col_idx) + (thlm(i,j,k) - thlmav(k))**2
        thv2r_col(k,col_idx) = thv2r_col(k,col_idx) + (thv0(i,j,k) - thvmav(k))**2
        ! Match modgenstat: th2 uses thlm moments around thmav.
        th2r_col(k,col_idx)  = th2r_col(k,col_idx)  + (thlm(i,j,k) - thmav(k))**2
        qt2r_col(k,col_idx)  = qt2r_col(k,col_idx)  + (qtm(i,j,k) - qtmav(k))**2
        ql2r_col(k,col_idx)  = ql2r_col(k,col_idx)  + (ql0(i,j,k) - qlmav(k))**2

        ! Accumulate w'^3 and normalize later with slab w2av, matching modgenstat.
        skew_col(k,col_idx) = skew_col(k,col_idx) + w_prime**3

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

        wthls_col(k,col_idx) = wthls_col(k,col_idx) + wthls_loc
        wthlr_col(k,col_idx) = wthlr_col(k,col_idx) + wthlr_loc
        wthlt_col(k,col_idx) = wthlt_col(k,col_idx) + wthls_loc + wthlr_loc
        wthvs_col(k,col_idx) = wthvs_col(k,col_idx) + wthvs_loc
        wthvr_col(k,col_idx) = wthvr_col(k,col_idx) + wthvr_loc
        wthvt_col(k,col_idx) = wthvt_col(k,col_idx) + wthvs_loc + wthvr_loc
        wqts_col(k,col_idx)  = wqts_col(k,col_idx)  + wqts_loc
        wqtr_col(k,col_idx)  = wqtr_col(k,col_idx)  + wqtr_loc
        wqtt_col(k,col_idx)  = wqtt_col(k,col_idx)  + wqts_loc  + wqtr_loc
        wqls_col(k,col_idx)  = wqls_col(k,col_idx)  + wqls_loc
        wqlr_col(k,col_idx)  = wqlr_col(k,col_idx)  + wqlr_loc
        wqlt_col(k,col_idx)  = wqlt_col(k,col_idx)  + wqls_loc  + wqlr_loc
        uws_col(k,col_idx)   = uws_col(k,col_idx)   + uws_loc
        uwr_col(k,col_idx)   = uwr_col(k,col_idx)   + uwr_loc
        uwt_col(k,col_idx)   = uwt_col(k,col_idx)   + uws_loc   + uwr_loc
        vws_col(k,col_idx)   = vws_col(k,col_idx)   + vws_loc
        vwr_col(k,col_idx)   = vwr_col(k,col_idx)   + vwr_loc
        vwt_col(k,col_idx)   = vwt_col(k,col_idx)   + vws_loc   + vwr_loc


        cs_col(k,col_idx)    = cs_col(k,col_idx)    + csz(k)
        hus_col(k,col_idx)   = hus_col(k,col_idx)   + (qt0(i,j,k) - ql0(i,j,k))
        hur_col(k,col_idx)   = hur_col(k,col_idx)   + 100.0_field_r * (qt0(i,j,k) - ql0(i,j,k)) / qsat_tab(tmp0(i,j,k), presf(k))

        ilratio = max(0._field_r,min(1._field_r,(tmp0(i,j,k)-tdn) / (tup-tdn)))
        clw_col(k,col_idx)   = clw_col(k,col_idx)   + ql0(i,j,k) * ilratio
        cli_col(k,col_idx)   = cli_col(k,col_idx)   + ql0(i,j,k) * (1.0_field_r - ilratio)
        if (ql0(i,j,k) > 0.0_field_r) cfrac_col(k,col_idx) = cfrac_col(k,col_idx) + 1.0_field_r
        
        iqr = get_tracer_index("qr")
        if (iqr > 0) then
          if (imicro == imicro_sice .or. imicro == imicro_sice2) then
              ilratio = max(0._field_r,min(1._field_r,(tmp0(i,j,k)-tdnrsg)/(tuprsg-tdnrsg)))
              plw_col(k,col_idx) = plw_col(k,col_idx) + sv0(i,j,k,iqr) * ilratio
              pli_col(k,col_idx) = pli_col(k,col_idx) + sv0(i,j,k,iqr) * (1-ilratio)
          else
            plw_col(k,col_idx) = plw_col(k,col_idx) + sv0(i,j,k,iqr)
          end if
        end if

        do n = 1, nsv
          sv_col(k,n,col_idx)   = sv_col(k,n,col_idx)   + svm(i,j,k,n)
          svp_col(k,n,col_idx)  = svp_col(k,n,col_idx)  + svp(i,j,k,n)
          sv2r_col(k,n,col_idx) = sv2r_col(k,n,col_idx) + (svm(i,j,k,n) - svmav(k,n))**2
          if (k == 1) then
            wsvs_loc = svflux(i,j,n)
            wsvr_loc = 0.0_field_r
          else
            ekhalf   = (ekh(i,j,k)*dzf(k-1) + ekh(i,j,k-1)*dzf(k)) / (2.0_field_r * dzh(k))
            wsvs_loc = -ekhalf * (sv0(i,j,k,n) - sv0(i,j,k-1,n)) / dzh(k)
            wsvr_loc = (w0(i,j,k) - wmav(k)) * ((sv0(i,j,k,n)*dzf(k-1) + sv0(i,j,k-1,n)*dzf(k)) / (2.0_field_r*dzh(k)))
          end if
          wsvs_col(k,n,col_idx) = wsvs_col(k,n,col_idx) + wsvs_loc
          wsvr_col(k,n,col_idx) = wsvr_col(k,n,col_idx) + wsvr_loc
          wsvt_col(k,n,col_idx) = wsvt_col(k,n,col_idx) + wsvs_loc + wsvr_loc
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
        shr_col(k,col_idx)   = shr_col(k,col_idx)   + shr_loc(k)
        buo_col(k,col_idx)   = buo_col(k,col_idx)   + buo_loc(k)
        ptrsp_col(k,col_idx) = ptrsp_col(k,col_idx) + ptrsp_loc(k)

        egp = 0.5*rhobf(k)*( (0.5*(u0(i,j,k)+u0(i+1,j,k))-(u0av(k)-cu))**2 &
                  +(0.5*(v0(i,j,k)+v0(i,j+1,k))-(v0av(k)-cv))**2 &
                  +(0.5*(w0(i,j,k)+w0(i,j,k+1))             )**2 )
        weres_loc(k) = egp * 0.5_field_r * (w0(i,j,k) + w0(i,j,k+1))
        tker_col(k,col_idx)  = tker_col(k,col_idx) + egp
        
        sbtke_inst = e120(i,j,k)**2 * rhobf(k)
        sbtke_col(k,col_idx)  = sbtke_col(k,col_idx) + sbtke_inst
        if (.not. l_sbtke_beg_set) sbtke_beg_col(k,col_idx) = sbtke_inst
        sbtke_last_col(k,col_idx) = sbtke_inst
        ekm_col(k,col_idx)    = ekm_col(k,col_idx)    + ekm(i,j,k)
        if (ekm(i,j,k) > eps1) then
          khkm_col(k,col_idx) = khkm_col(k,col_idx) + ekh(i,j,k) / ekm(i,j,k)
        end if



        thltend_col(k,col_idx) = thltend_col(k,col_idx) + thlprad(i,j,k)
        
        !absolute values to handle different sign conventions in radiation code
        thllwtend_col(k,col_idx) = thllwtend_col(k,col_idx) + &
            (abs(lwd(i,j,k+1)) - abs(lwu(i,j,k+1)) - abs(lwd(i,j,k)) + abs(lwu(i,j,k))) / &
            (rhof(k)*exnf(k)*cp*dzf(k))
        thlswtend_col(k,col_idx) = thlswtend_col(k,col_idx) + &
            (abs(swd(i,j,k+1)) - abs(swu(i,j,k+1)) - abs(swd(i,j,k)) + abs(swu(i,j,k))) / &
            (rhof(k)*exnf(k)*cp*dzf(k))

        thlradls_col(k,col_idx) = thlradls_col(k,col_idx) + thlpcar(k)

        lwu_col(k,col_idx)  = lwu_col(k,col_idx)  + abs(lwu(i,j,k))
        lwd_col(k,col_idx)  = lwd_col(k,col_idx)  + abs(lwd(i,j,k))
        swu_col(k,col_idx)  = swu_col(k,col_idx)  + abs(swu(i,j,k))
        swd_col(k,col_idx)  = swd_col(k,col_idx)  + abs(swd(i,j,k))

        ! assume clear-sky fluxes are already calculated in modradstat and available as lwuca/lwdca/swuca/swdca; if not, these will just accumulate zeros and not contribute to the final averages
        lwuca_col(k,col_idx) = lwuca_col(k,col_idx) + abs(lwuca(i,j,k))
        lwdca_col(k,col_idx) = lwdca_col(k,col_idx) + abs(lwdca(i,j,k))
        swuca_col(k,col_idx) = swuca_col(k,col_idx) + abs(swuca(i,j,k))
        swdca_col(k,col_idx) = swdca_col(k,col_idx) + abs(swdca(i,j,k))

        thllwtendca_col(k,col_idx) = thllwtendca_col(k,col_idx) + &
            (-lwdca(i,j,k+1) - lwuca(i,j,k+1) + lwdca(i,j,k) + lwuca(i,j,k)) / &
            (rhof(k)*exnf(k)*cp*dzf(k))
        thlswtendca_col(k,col_idx) = thlswtendca_col(k,col_idx) + &
            (-swdca(i,j,k+1) - swuca(i,j,k+1) + swdca(i,j,k) + swuca(i,j,k)) / &
            (rhof(k)*exnf(k)*cp*dzf(k))

      end do

      ! Post-loop: resolved transport (divergence-based)
      trsp_col(1,col_idx) = trsp_col(1,col_idx) - weres_loc(1) / (0.5_field_r * dzh(1))
      do k = 2, kmax
        trsp_col(k,col_idx) = trsp_col(k,col_idx) - (weres_loc(k) - weres_loc(k-1)) / dzh(k)
      end do
      
      ! Note: ptrsp_loc contains raw flux values accumulated into ptrsp_col.
      ! During write phase (after averaging), will convert ptrsp_col to divergence form.
      ! This matches modbudget's final transformation: ptrsp = -(d/dz)(ptrsp_flux)

      k = 1
      !sbshr, sbbuo, sbdiss are only defined at k=1 (surface layer)
      !These are multiplied by rhobf during accumulation per modbudget convention
      sbshr_col(k,col_idx)  = sbshr_col(k,col_idx)  + sbshr(i,j,k) * rhobf(k)
      sbbuo_col(k,col_idx)  = sbbuo_col(k,col_idx)  + sbbuo(i,j,k) * rhobf(k)
      sbdiss_col(k,col_idx) = sbdiss_col(k,col_idx) + sbdiss(i,j,k) * rhobf(k)
      ! sbbudg is computed from components: sum of shear, buoyancy, dissipation  
      ! Will be calculated during output normalization from accumulated sbshr+sbbuo+sbdiss

    end do  ! accumulation loop

    if (.not. l_sbtke_beg_set) l_sbtke_beg_set = .true.

    ! ----------------------------------------------------------------
    ! Write pass — every timeav: average, fill file buffers, and zero
    ! ----------------------------------------------------------------
    if (is_writing_timestep(ofile_id)) then
      any_write = .true.
      scale_rn = 1.0_field_r / real(nsamples, kind=field_r)

      do col_idx = 1, size(i_local)

        rhof_col(:,col_idx)   = rhof_col(:,col_idx)   * scale_rn
        rhobf_col(:,col_idx)  = rhobf_col(:,col_idx)  * scale_rn
        rhobh_col(:,col_idx)  = rhobh_col(:,col_idx)  * scale_rn
        presh_col(:,col_idx)  = presh_col(:,col_idx)  * scale_rn
        u_col(:,col_idx)      = u_col(:,col_idx)      * scale_rn
        v_col(:,col_idx)      = v_col(:,col_idx)      * scale_rn
        w_col(:,col_idx)      = w_col(:,col_idx)      * scale_rn
        thl_col(:,col_idx)    = thl_col(:,col_idx)    * scale_rn
        thv_col(:,col_idx)    = thv_col(:,col_idx)    * scale_rn
        qt_col(:,col_idx)     = qt_col(:,col_idx)     * scale_rn
        ql_col(:,col_idx)     = ql_col(:,col_idx)     * scale_rn
        wthls_col(:,col_idx)  = wthls_col(:,col_idx)  * scale_rn
        wthlr_col(:,col_idx)  = wthlr_col(:,col_idx)  * scale_rn
        wthlt_col(:,col_idx)  = wthlt_col(:,col_idx)  * scale_rn
        wthvs_col(:,col_idx)  = wthvs_col(:,col_idx)  * scale_rn
        wthvr_col(:,col_idx)  = wthvr_col(:,col_idx)  * scale_rn
        wthvt_col(:,col_idx)  = wthvt_col(:,col_idx)  * scale_rn
        wqts_col(:,col_idx)   = wqts_col(:,col_idx)   * scale_rn
        wqtr_col(:,col_idx)   = wqtr_col(:,col_idx)   * scale_rn
        wqtt_col(:,col_idx)   = wqtt_col(:,col_idx)   * scale_rn
        wqls_col(:,col_idx)   = wqls_col(:,col_idx)   * scale_rn
        wqlr_col(:,col_idx)   = wqlr_col(:,col_idx)   * scale_rn
        wqlt_col(:,col_idx)   = wqlt_col(:,col_idx)   * scale_rn
        uws_col(:,col_idx)    = uws_col(:,col_idx)    * scale_rn
        uwr_col(:,col_idx)    = uwr_col(:,col_idx)    * scale_rn
        uwt_col(:,col_idx)    = uwt_col(:,col_idx)    * scale_rn
        vws_col(:,col_idx)    = vws_col(:,col_idx)    * scale_rn
        vwr_col(:,col_idx)    = vwr_col(:,col_idx)    * scale_rn
        vwt_col(:,col_idx)    = vwt_col(:,col_idx)    * scale_rn
        w2s_col(:,col_idx)    = w2s_col(:,col_idx)    * scale_rn
        w2r_col(:,col_idx)    = w2r_col(:,col_idx)    * scale_rn
        skew_col(:,col_idx)   = skew_col(:,col_idx)   * scale_rn
        do k = 1, kmax
          skew_col(k,col_idx) = skew_col(k,col_idx) / max(w2av(k)**1.5_field_r, epsilon(1.0_field_r))
        end do
        u2r_col(:,col_idx)    = u2r_col(:,col_idx)    * scale_rn
        v2r_col(:,col_idx)    = v2r_col(:,col_idx)    * scale_rn
        thl2r_col(:,col_idx)  = thl2r_col(:,col_idx)  * scale_rn
        thv2r_col(:,col_idx)  = thv2r_col(:,col_idx)  * scale_rn
        th2r_col(:,col_idx)   = th2r_col(:,col_idx)   * scale_rn
        qt2r_col(:,col_idx)   = qt2r_col(:,col_idx)   * scale_rn
        ql2r_col(:,col_idx)   = ql2r_col(:,col_idx)   * scale_rn
        cs_col(:,col_idx)     = cs_col(:,col_idx)     * scale_rn
        cfrac_col(:,col_idx)  = cfrac_col(:,col_idx)  * scale_rn
        hur_col(:,col_idx)    = hur_col(:,col_idx)    * scale_rn
        hus_col(:,col_idx)    = hus_col(:,col_idx)    * scale_rn
        ta_col(:,col_idx)     = ta_col(:,col_idx)     * scale_rn
        clw_col(:,col_idx)    = clw_col(:,col_idx)    * scale_rn
        cli_col(:,col_idx)    = cli_col(:,col_idx)    * scale_rn
        plw_col(:,col_idx)    = plw_col(:,col_idx)    * scale_rn
        pli_col(:,col_idx)    = pli_col(:,col_idx)    * scale_rn
        tker_col(:,col_idx)   = tker_col(:,col_idx)   * scale_rn
        shr_col(:,col_idx)    = shr_col(:,col_idx)    * scale_rn
        buo_col(:,col_idx)    = buo_col(:,col_idx)    * scale_rn
        trsp_col(:,col_idx)   = trsp_col(:,col_idx)   * scale_rn
        ptrsp_col(:,col_idx)  = ptrsp_col(:,col_idx)  * scale_rn

        ! Match modbudget staggering: report full-level values from adjacent half levels.
        prof_tmp = shr_col(:,col_idx)
        do k = 1, kmax-1
          shr_col(k,col_idx) = 0.5_field_r * (prof_tmp(k) + prof_tmp(k+1))
        end do
        shr_col(kmax,col_idx) = prof_tmp(kmax)

        prof_tmp = buo_col(:,col_idx)
        do k = 1, kmax-1
          buo_col(k,col_idx) = 0.5_field_r * (prof_tmp(k) + prof_tmp(k+1))
        end do
        buo_col(kmax,col_idx) = prof_tmp(kmax)

        prof_tmp = trsp_col(:,col_idx)
        do k = 1, kmax-1
          trsp_col(k,col_idx) = 0.5_field_r * (prof_tmp(k) + prof_tmp(k+1))
        end do
        trsp_col(kmax,col_idx) = prof_tmp(kmax)
        
        ! Convert ptrsp from accumulated raw flux to divergence form (matching modbudget):
        ! ptrsp_final(k) = -(ptrsp_flux(k+1) - ptrsp_flux(k)) / dzf(k)
        prof_tmp = ptrsp_col(:,col_idx)
        ptrsp_col(1,col_idx) = 0.0_field_r  ! ptrsp(1) = 0 (no divergence at surface)
        do k = 2, kmax-1
          ptrsp_col(k,col_idx) = -(prof_tmp(k+1) - prof_tmp(k)) / dzf(k)
        end do
        ptrsp_col(kmax,col_idx) = 0.0_field_r
        
        sbtke_col(:,col_idx)  = sbtke_col(:,col_idx)  * scale_rn
        sbshr_col(:,col_idx)  = sbshr_col(:,col_idx)  * scale_rn
        sbbuo_col(:,col_idx)  = sbbuo_col(:,col_idx)  * scale_rn
        sbdiss_col(:,col_idx) = sbdiss_col(:,col_idx) * scale_rn
        ! sbbudg = sum of shear + buoyancy + dissipation (per modbudget convention)
        sbbudg_col(:,col_idx) = sbshr_col(:,col_idx) + sbbuo_col(:,col_idx) + sbdiss_col(:,col_idx)
        sbstor_col(:,col_idx)  = (sbtke_last_col(:,col_idx) - sbtke_beg_col(:,col_idx)) / max(timeav, eps1)
        sbresid_col(:,col_idx) = sbbudg_col(:,col_idx) - sbstor_col(:,col_idx)
        ekm_col(:,col_idx)    = ekm_col(:,col_idx)    * scale_rn
        khkm_col(:,col_idx)   = khkm_col(:,col_idx)   * scale_rn
        thltend_col(:,col_idx)    = thltend_col(:,col_idx)    * scale_rn
        thllwtend_col(:,col_idx)  = thllwtend_col(:,col_idx)  * scale_rn
        thlswtend_col(:,col_idx)  = thlswtend_col(:,col_idx)  * scale_rn
        thlradls_col(:,col_idx)   = thlradls_col(:,col_idx)   * scale_rn
        lwu_col(:,col_idx)    = lwu_col(:,col_idx)    * scale_rn
        lwd_col(:,col_idx)    = lwd_col(:,col_idx)    * scale_rn
        swu_col(:,col_idx)    = swu_col(:,col_idx)    * scale_rn
        swd_col(:,col_idx)    = swd_col(:,col_idx)    * scale_rn
        lwuca_col(:,col_idx)      = lwuca_col(:,col_idx)      * scale_rn
        lwdca_col(:,col_idx)      = lwdca_col(:,col_idx)      * scale_rn
        swuca_col(:,col_idx)      = swuca_col(:,col_idx)      * scale_rn
        swdca_col(:,col_idx)      = swdca_col(:,col_idx)      * scale_rn
        thllwtendca_col(:,col_idx) = thllwtendca_col(:,col_idx) * scale_rn
        thlswtendca_col(:,col_idx) = thlswtendca_col(:,col_idx) * scale_rn
        if (nsv > 0) then
          sv_col(:,:,col_idx)   = sv_col(:,:,col_idx)   * scale_rn
          svp_col(:,:,col_idx)  = svp_col(:,:,col_idx)  * scale_rn
          sv2r_col(:,:,col_idx) = sv2r_col(:,:,col_idx) * scale_rn
          wsvs_col(:,:,col_idx) = wsvs_col(:,:,col_idx) * scale_rn
          wsvr_col(:,:,col_idx) = wsvr_col(:,:,col_idx) * scale_rn
          wsvt_col(:,:,col_idx) = wsvt_col(:,:,col_idx) * scale_rn
        end if

      end do

      call write_all_profiles_to_file()
    end if

    ! Zero all accumulators after writing (all ifiles share the same timing)
    if (any_write) then
      call zero_profile_arrays()
      l_sbtke_beg_set = .false.
    end if

  end subroutine colstat


  subroutine write_all_profiles_to_file

    integer :: n
    character(len=64) :: vname

    call write_field_to_file('rhof', rhof_col)
    call write_field_to_file('rhobf', rhobf_col)
    call write_field_to_file('rhobh', rhobh_col)
    call write_field_to_file('presh', presh_col)
    call write_field_to_file('u', u_col)
    call write_field_to_file('v', v_col)
    call write_field_to_file('w', w_col)
    call write_field_to_file('thl', thl_col)
    call write_field_to_file('thv', thv_col)
    call write_field_to_file('qt', qt_col)
    call write_field_to_file('ql', ql_col)
    call write_field_to_file('wthls', wthls_col)
    call write_field_to_file('wthlr', wthlr_col)
    call write_field_to_file('wthlt', wthlt_col)
    call write_field_to_file('wthvs', wthvs_col)
    call write_field_to_file('wthvr', wthvr_col)
    call write_field_to_file('wthvt', wthvt_col)
    call write_field_to_file('wqts', wqts_col)
    call write_field_to_file('wqtr', wqtr_col)
    call write_field_to_file('wqtt', wqtt_col)
    call write_field_to_file('wqls', wqls_col)
    call write_field_to_file('wqlr', wqlr_col)
    call write_field_to_file('wqlt', wqlt_col)
    call write_field_to_file('uws', uws_col)
    call write_field_to_file('uwr', uwr_col)
    call write_field_to_file('uwt', uwt_col)
    call write_field_to_file('vws', vws_col)
    call write_field_to_file('vwr', vwr_col)
    call write_field_to_file('vwt', vwt_col)
    call write_field_to_file('w2s', w2s_col)
    call write_field_to_file('w2r', w2r_col)
    call write_field_to_file('skew', skew_col)
    call write_field_to_file('u2r', u2r_col)
    call write_field_to_file('v2r', v2r_col)
    call write_field_to_file('thl2r', thl2r_col)
    call write_field_to_file('thv2r', thv2r_col)
    call write_field_to_file('th2r', th2r_col)
    call write_field_to_file('qt2r', qt2r_col)
    call write_field_to_file('ql2r', ql2r_col)
    call write_field_to_file('cs', cs_col)
    call write_field_to_file('cfrac', cfrac_col)
    call write_field_to_file('hur', hur_col)
    call write_field_to_file('hus', hus_col)
    call write_field_to_file('ta', ta_col)
    call write_field_to_file('clw', clw_col)
    call write_field_to_file('cli', cli_col)
    call write_field_to_file('plw', plw_col)
    call write_field_to_file('pli', pli_col)

    do n = 1, nsv
      vname = trim(tracer_prop(n)%tracname)
      call write_field_to_file(vname, sv_col(:,n,:))
      call write_field_to_file(trim(vname)//'p', svp_col(:,n,:))
      call write_field_to_file(trim(vname)//'2r', sv2r_col(:,n,:))
      call write_field_to_file('w'//trim(vname)//'s', wsvs_col(:,n,:))
      call write_field_to_file('w'//trim(vname)//'r', wsvr_col(:,n,:))
      call write_field_to_file('w'//trim(vname)//'t', wsvt_col(:,n,:))
    end do

    call write_field_to_file('tker', tker_col)
    call write_field_to_file('shr', shr_col)
    call write_field_to_file('buo', buo_col)
    call write_field_to_file('trsp', trsp_col)
    call write_field_to_file('ptrsp', ptrsp_col)
    call write_field_to_file('sbtke', sbtke_col)
    call write_field_to_file('sbshr', sbshr_col)
    call write_field_to_file('sbbuo', sbbuo_col)
    call write_field_to_file('sbdiss', sbdiss_col)
    call write_field_to_file('sbstor', sbstor_col)
    call write_field_to_file('sbbudg', sbbudg_col)
    call write_field_to_file('sbresid', sbresid_col)
    call write_field_to_file('ekm', ekm_col)
    call write_field_to_file('khkm', khkm_col)
    call write_field_to_file('thltend', thltend_col)
    call write_field_to_file('thllwtend', thllwtend_col)
    call write_field_to_file('thlswtend', thlswtend_col)
    call write_field_to_file('thlradls', thlradls_col)
    call write_field_to_file('lwu', lwu_col)
    call write_field_to_file('lwd', lwd_col)
    call write_field_to_file('swu', swu_col)
    call write_field_to_file('swd', swd_col)
    call write_field_to_file('lwuca', lwuca_col)
    call write_field_to_file('lwdca', lwdca_col)
    call write_field_to_file('swuca', swuca_col)
    call write_field_to_file('swdca', swdca_col)
    call write_field_to_file('thllwtendca', thllwtendca_col)
    call write_field_to_file('thlswtendca', thlswtendca_col)

  end subroutine write_all_profiles_to_file

  subroutine write_field_to_file(name, local_field)

    character(len=*), intent(in) :: name
    real(field_r), intent(in) :: local_field(:,:)

    integer :: il
    real(field_r), pointer :: ptr(:)

    do il = 1, size(local_field, dim=2)
      call ofile%get_pointer(name, global_index(il), ptr)
      ptr = local_field(:,il)
    end do

  end subroutine write_field_to_file

  subroutine allocate_profile_arrays(nlocal)

    integer, intent(in) :: nlocal

    allocate(rhof_col(kmax,nlocal), rhobf_col(kmax,nlocal), rhobh_col(kmax,nlocal), presh_col(kmax,nlocal))
    allocate(u_col(kmax,nlocal), v_col(kmax,nlocal), w_col(kmax,nlocal), thl_col(kmax,nlocal), thv_col(kmax,nlocal), &
             qt_col(kmax,nlocal), ql_col(kmax,nlocal))
    allocate(wthls_col(kmax,nlocal), wthlr_col(kmax,nlocal), wthlt_col(kmax,nlocal), wthvs_col(kmax,nlocal), &
             wthvr_col(kmax,nlocal), wthvt_col(kmax,nlocal))
    allocate(wqts_col(kmax,nlocal), wqtr_col(kmax,nlocal), wqtt_col(kmax,nlocal), wqls_col(kmax,nlocal), &
             wqlr_col(kmax,nlocal), wqlt_col(kmax,nlocal))
    allocate(uws_col(kmax,nlocal), uwr_col(kmax,nlocal), uwt_col(kmax,nlocal), vws_col(kmax,nlocal), &
             vwr_col(kmax,nlocal), vwt_col(kmax,nlocal))
    allocate(w2s_col(kmax,nlocal), w2r_col(kmax,nlocal), skew_col(kmax,nlocal), u2r_col(kmax,nlocal), v2r_col(kmax,nlocal))
    allocate(thl2r_col(kmax,nlocal), thv2r_col(kmax,nlocal), th2r_col(kmax,nlocal), qt2r_col(kmax,nlocal), ql2r_col(kmax,nlocal))
    allocate(cs_col(kmax,nlocal), cfrac_col(kmax,nlocal), hur_col(kmax,nlocal), hus_col(kmax,nlocal), ta_col(kmax,nlocal))
    allocate(clw_col(kmax,nlocal), cli_col(kmax,nlocal), plw_col(kmax,nlocal), pli_col(kmax,nlocal))

    allocate(sv_col(kmax,nsv,nlocal), svp_col(kmax,nsv,nlocal), sv2r_col(kmax,nsv,nlocal))
    allocate(wsvs_col(kmax,nsv,nlocal), wsvr_col(kmax,nsv,nlocal), wsvt_col(kmax,nsv,nlocal))

    allocate(tker_col(kmax,nlocal), shr_col(kmax,nlocal), buo_col(kmax,nlocal), trsp_col(kmax,nlocal), ptrsp_col(kmax,nlocal))
    allocate(sbtke_col(kmax,nlocal), sbshr_col(kmax,nlocal), sbbuo_col(kmax,nlocal), sbdiss_col(kmax,nlocal))
    allocate(sbstor_col(kmax,nlocal), sbbudg_col(kmax,nlocal), sbresid_col(kmax,nlocal), ekm_col(kmax,nlocal), khkm_col(kmax,nlocal))
    allocate(sbtke_beg_col(kmax,nlocal), sbtke_last_col(kmax,nlocal))

    allocate(thltend_col(kmax,nlocal), thllwtend_col(kmax,nlocal), thlswtend_col(kmax,nlocal), thlradls_col(kmax,nlocal))
    allocate(lwu_col(kmax,nlocal), lwd_col(kmax,nlocal), swu_col(kmax,nlocal), swd_col(kmax,nlocal))
    allocate(lwuca_col(kmax,nlocal), lwdca_col(kmax,nlocal), swuca_col(kmax,nlocal), swdca_col(kmax,nlocal))
    allocate(thllwtendca_col(kmax,nlocal), thlswtendca_col(kmax,nlocal))

  end subroutine allocate_profile_arrays

  subroutine zero_profile_arrays

    rhof_col = 0.0_field_r
    rhobf_col = 0.0_field_r
    rhobh_col = 0.0_field_r
    presh_col = 0.0_field_r
    u_col = 0.0_field_r
    v_col = 0.0_field_r
    w_col = 0.0_field_r
    thl_col = 0.0_field_r
    thv_col = 0.0_field_r
    qt_col = 0.0_field_r
    ql_col = 0.0_field_r
    wthls_col = 0.0_field_r
    wthlr_col = 0.0_field_r
    wthlt_col = 0.0_field_r
    wthvs_col = 0.0_field_r
    wthvr_col = 0.0_field_r
    wthvt_col = 0.0_field_r
    wqts_col = 0.0_field_r
    wqtr_col = 0.0_field_r
    wqtt_col = 0.0_field_r
    wqls_col = 0.0_field_r
    wqlr_col = 0.0_field_r
    wqlt_col = 0.0_field_r
    uws_col = 0.0_field_r
    uwr_col = 0.0_field_r
    uwt_col = 0.0_field_r
    vws_col = 0.0_field_r
    vwr_col = 0.0_field_r
    vwt_col = 0.0_field_r
    w2s_col = 0.0_field_r
    w2r_col = 0.0_field_r
    skew_col = 0.0_field_r
    u2r_col = 0.0_field_r
    v2r_col = 0.0_field_r
    thl2r_col = 0.0_field_r
    thv2r_col = 0.0_field_r
    th2r_col = 0.0_field_r
    qt2r_col = 0.0_field_r
    ql2r_col = 0.0_field_r
    cs_col = 0.0_field_r
    cfrac_col = 0.0_field_r
    hur_col = 0.0_field_r
    hus_col = 0.0_field_r
    ta_col = 0.0_field_r
    clw_col = 0.0_field_r
    cli_col = 0.0_field_r
    plw_col = 0.0_field_r
    pli_col = 0.0_field_r
    sv_col = 0.0_field_r
    svp_col = 0.0_field_r
    sv2r_col = 0.0_field_r
    wsvs_col = 0.0_field_r
    wsvr_col = 0.0_field_r
    wsvt_col = 0.0_field_r
    tker_col = 0.0_field_r
    shr_col = 0.0_field_r
    buo_col = 0.0_field_r
    trsp_col = 0.0_field_r
    ptrsp_col = 0.0_field_r
    sbtke_col = 0.0_field_r
    sbshr_col = 0.0_field_r
    sbbuo_col = 0.0_field_r
    sbdiss_col = 0.0_field_r
    sbtke_beg_col = 0.0_field_r
    sbtke_last_col = 0.0_field_r
    sbstor_col = 0.0_field_r
    sbbudg_col = 0.0_field_r
    sbresid_col = 0.0_field_r
    ekm_col = 0.0_field_r
    khkm_col = 0.0_field_r
    thltend_col = 0.0_field_r
    thllwtend_col = 0.0_field_r
    thlswtend_col = 0.0_field_r
    thlradls_col = 0.0_field_r
    lwu_col = 0.0_field_r
    lwd_col = 0.0_field_r
    swu_col = 0.0_field_r
    swd_col = 0.0_field_r
    lwuca_col = 0.0_field_r
    lwdca_col = 0.0_field_r
    swuca_col = 0.0_field_r
    swdca_col = 0.0_field_r
    thllwtendca_col = 0.0_field_r
    thlswtendca_col = 0.0_field_r

  end subroutine zero_profile_arrays

  subroutine deallocate_profile_arrays

    if (allocated(rhof_col)) deallocate(rhof_col, rhobf_col, rhobh_col, presh_col)
    if (allocated(u_col)) deallocate(u_col, v_col, w_col, thl_col, thv_col, qt_col, ql_col)
    if (allocated(wthls_col)) deallocate(wthls_col, wthlr_col, wthlt_col, wthvs_col, wthvr_col, wthvt_col)
    if (allocated(wqts_col)) deallocate(wqts_col, wqtr_col, wqtt_col, wqls_col, wqlr_col, wqlt_col)
    if (allocated(uws_col)) deallocate(uws_col, uwr_col, uwt_col, vws_col, vwr_col, vwt_col)
    if (allocated(w2s_col)) deallocate(w2s_col, w2r_col, skew_col, u2r_col, v2r_col)
    if (allocated(thl2r_col)) deallocate(thl2r_col, thv2r_col, th2r_col, qt2r_col, ql2r_col)
    if (allocated(cs_col)) deallocate(cs_col, cfrac_col, hur_col, hus_col, ta_col)
    if (allocated(clw_col)) deallocate(clw_col, cli_col, plw_col, pli_col)
    if (allocated(sv_col)) deallocate(sv_col, svp_col, sv2r_col, wsvs_col, wsvr_col, wsvt_col)
    if (allocated(tker_col)) deallocate(tker_col, shr_col, buo_col, trsp_col, ptrsp_col)
    if (allocated(sbtke_col)) deallocate(sbtke_col, sbshr_col, sbbuo_col, sbdiss_col)
    if (allocated(sbtke_beg_col)) deallocate(sbtke_beg_col, sbtke_last_col)
    if (allocated(sbstor_col)) deallocate(sbstor_col, sbbudg_col, sbresid_col, ekm_col, khkm_col)
    if (allocated(thltend_col)) deallocate(thltend_col, thllwtend_col, thlswtend_col, thlradls_col)
    if (allocated(lwu_col)) deallocate(lwu_col, lwd_col, swu_col, swd_col)
    if (allocated(lwuca_col)) deallocate(lwuca_col, lwdca_col, swuca_col, swdca_col)
    if (allocated(thllwtendca_col)) deallocate(thllwtendca_col, thlswtendca_col)

  end subroutine deallocate_profile_arrays

  subroutine exitcolstat

    if (allocated(i_local)) deallocate(i_local)
    if (allocated(j_local)) deallocate(j_local)
    if (allocated(global_index)) deallocate(global_index)
    if (allocated(locx)) deallocate(locx)
    if (allocated(locy)) deallocate(locy)
    call deallocate_profile_arrays()

  end subroutine exitcolstat

end module modcolstat
