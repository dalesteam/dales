!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!

!> Dumps instantaneous cross sections of several fields.
module modcrosssection

  use fortran_support,   only: nnml_output, finish, &
#if FIELD_PRECISION==64
                               t_ptr_2d => t_ptr_2d_dp
#else
                               t_ptr_2d => t_ptr_2d_sp
#endif
  use modglobal,         only: longint, kmax, nsv, cu, cv, itot, jtot, imax, &
                               jmax, kmax, i1, j1, ifnamopt, dtav_glob, &
                               rk3step, checknamelisterror, dx, dy, zf
  use modtracers,        only: tracer_prop
  use modnetcdf_file_t,  only: cross_section_file_t
  use modstat_nc_files,  only: add_output_file, is_sampling_timestep
  use modprecision,      only: field_r
  use modthermodynamics, only: calc_virt_pot_temp
  use modmpi,            only: D_MPI_BCAST, commwrld, mpierr, myid, myidx, &
                               myidy
  use modfields,         only: um, vm, wm, thlm, qtm, ql0, thvf, e12m, exnf
#if defined(_OPENACC)
  use modgpu,            only: update_host
#endif

  implicit none

  private

  character(len=*), parameter :: modname = 'modcrosssection'

  public :: crosssection_read_namelist
  public :: initcrosssection
  public :: crosssection

  ! User settings
  logical :: lcross = .false. !< switch for doing the crosssection (on/off).
  real    :: dtav             !< Sampling interval [s].
  integer :: crossplane(100)  !< Locations of the xz cross sections.
  integer :: crossheight(100) !< Locations of the xy cross sections.
  integer :: crossortho(100)  !< Locations of the yz cross sections.
  logical :: lxz = .true.     !< Switch for doing xz crosssections.
  logical :: lxy = .true.     !< Switch for doing xy crosssections.
  logical :: lyz = .true.     !< Switch for doing yz crosssections.

  integer :: nxz !< Number of xz cross sections.
  integer :: nxy !< Number of xy cross sections.
  integer :: nyz !< Number of yz cross sections.

  ! NetCDF files.
  type(cross_section_file_t), allocatable :: xy_files(:)    !< List of xy cross files.
  type(cross_section_file_t), allocatable :: xz_files(:)    !< List of xz cross files.
  type(cross_section_file_t), allocatable :: yz_files(:)    !< List of yz cross files.
  integer,                    allocatable :: xy_file_ids(:) !< List of xy cross file ids.
  integer,                    allocatable :: yz_file_ids(:) !< List of yz cross file ids.
  integer,                    allocatable :: xz_file_ids(:) !< List of xz cross file ids.

contains

  !> Read the cross section namelist.
  subroutine crosssection_read_namelist(nml_filename)

    character(len=*), intent(in) :: nml_filename !< Filename of the namelist.

    integer :: ierr

    namelist /NAMCROSSSECTION/ lcross, dtav, crossheight, crossplane, &
                               crossortho, lxy, lxz, lyz

    crossheight(1) = 2
    crossheight(2:100) = -999
    crossplane(1) = 2
    crossplane(2:100) = -999
    crossortho(1) = 2
    crossortho(2:100) = -999

    dtav = dtav_glob

    if (myid==0) then
      open(ifnamopt, file=nml_filename, status='old', iostat=ierr)
      read(ifnamopt, NAMCROSSSECTION, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMCROSSSECTION')
      write(nnml_output, NAMCROSSSECTION)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(dtav, 1, 0, commwrld, mpierr)
    call D_MPI_BCAST(lcross, 1, 0, commwrld, mpierr)
    call D_MPI_BCAST(crossheight(1:100), 100, 0, commwrld, mpierr)
    call D_MPI_BCAST(crossplane(1:100), 100, 0, commwrld, mpierr)
    call D_MPI_BCAST(crossortho(1:100), 100, 0, commwrld, mpierr)
    call D_MPI_BCAST(lxy, 1, 0, commwrld, mpierr)
    call D_MPI_BCAST(lxz, 1, 0, commwrld, mpierr)
    call D_MPI_BCAST(lyz, 1, 0, commwrld, mpierr)

  end subroutine crosssection_read_namelist

  !> Initializing Crosssection. Read out the namelist, initializing the variables
  subroutine initcrosssection

    character(len=*), parameter :: routine = modname//'/initcrosssection'

    integer          :: k, ifile
    real(field_r)    :: loc
    character(len=4) :: cloc

    if (.not. lcross) return

    if (any(crossheight > kmax) .or. any(crossplane > jtot + 1) &
        .or. any(crossortho > itot + 1)) then
      call finish(routine, 'crosssection out of range')
    end if

    k=1
    do while (crossheight(k) > 0)
       nxy=nxy+1
       k=k+1
    end do

    k=1
    do while (crossplane(k) > 0)
       crossplane(k) = crossplane(k) - myidy*jmax  ! convert to local grid index
       nxz=nxz+1
       k=k+1
    end do

    k=1
    do while (crossortho(k) > 0)
       crossortho(k) = crossortho(k) - myidx*imax  ! convert to local grid index
       nyz=nyz+1
       k=k+1
    end do

    ! cross sections with
    ! 2  <= crossplane, <= j1
    ! 2  <= crossortho <= i1
    ! belong to this processor

    ! XY cross sections

    allocate(xy_files(nxy), xy_file_ids(nxy))

    ifile = 0
    do k = 1, size(crossheight)
      if (crossheight(k) > 0) then
        ifile = ifile + 1
        write(cloc, '(i4.4)') crossheight(k)
        loc = zf(crossheight(k))
        xy_files(ifile) = cross_section_file_t('crossxy.'//cloc, nx=itot, &
                                               ny=jtot, loc=loc)
        call add_output_file(xy_files(ifile), dtav, xy_file_ids(ifile))
      end if
    end do

    do ifile = 1, nxy
      call xy_files(ifile)%add_var('u',    'xy crosssection of the west-east velocity',                'm/s',     'mt0t')
      call xy_files(ifile)%add_var('v',    'xy crosssection of the south-north velocity',              'm/s',     'tm0t')
      call xy_files(ifile)%add_var('w',    'xy crosssection of the vertical velocity',                 'm/s',     'tt0t')
      call xy_files(ifile)%add_var('thl',  'xy crossection of the liquid water potential temperature', 'K',       'tt0t')
      call xy_files(ifile)%add_var('thv',  'xy crosssection of the virtual potential temperature',     'K',       'tt0t')
      call xy_files(ifile)%add_var('qt',   'xy crosssection of the total water specific humidity',     'kg/kg',   'tt0t')
      call xy_files(ifile)%add_var('ql',   'xy crosssection of the liquid water specific humidity',    'kg/kg',   'tt0t')
      call xy_files(ifile)%add_var('buoy', 'xy crosssection of the buoyancy',                          'K',       'tt0t')
      call xy_files(ifile)%add_var('e120', 'xy crosssection of the sqrt(turbulent kinetic energy',     'm^2/s^2', 'tt0t')
    end do

    ! XZ cross sections

    allocate(xz_files(nxz), xz_file_ids(nxz))

    ifile = 0
    do k = 1, size(crossplane)
      if (crossplane(k) > 0) then
        ifile = ifile + 1
        write(cloc, '(i4.4)') crossplane(k)
        loc = dy * (crossplane(k) - 2) + 0.5_field_r * dy
        xz_files(ifile) = cross_section_file_t('crossxz.'//cloc, &
                                               nx=itot, nz=kmax, loc=loc)
        call add_output_file(xz_files(ifile), dtav, xz_file_ids(ifile))
      end if
    end do

    do ifile = 1, nxz
      call xz_files(ifile)%add_var('u',    'xz crosssection of the west-east velocity',                'm/s',     'm0tt')
      call xz_files(ifile)%add_var('v',    'xz crosssection of the south-north velocity',              'm/s',     't0tt')
      call xz_files(ifile)%add_var('w',    'xz crosssection of the vertical velocity',                 'm/s',     't0mt')
      call xz_files(ifile)%add_var('thl',  'xz crossection of the liquid water potential temperature', 'K',       't0tt')
      call xz_files(ifile)%add_var('thv',  'xz crosssection of the virtual potential temperature',     'K',       't0tt')
      call xz_files(ifile)%add_var('qt',   'xz crosssection of the total water specific humidity',     'kg/kg',   't0tt')
      call xz_files(ifile)%add_var('ql',   'xz crosssection of the liquid water specific humidity',    'kg/kg',   't0tt')
      call xz_files(ifile)%add_var('buoy', 'xz crosssection of the buoyancy',                          'K',       't0tt')
      call xz_files(ifile)%add_var('e120', 'xz crosssection of the sqrt(turbulent kinetic energy',     'm^2/s^2', 't0tt')
    end do

    ! YZ cross sections

    allocate(yz_files(nyz), yz_file_ids(nyz))

    ifile = 0
    do k = 1, size(crossortho)
      if (crossortho(k) > 0) then
        ifile = ifile + 1
        write(cloc, '(i4.4)') crossortho(k)
        loc = dx * (crossortho(k) - 2) + 0.5_field_r * dx
        yz_files(ifile) = cross_section_file_t('crossyz.'//cloc, &
                                               ny=jtot, nz=kmax, loc=loc)
        call add_output_file(yz_files(ifile), dtav, yz_file_ids(ifile))
      end if
    end do

    do ifile = 1, nyz
      call yz_files(ifile)%add_var('u',    'yz crosssection of the west-east velocity',                'm/s',     '0ttt')
      call yz_files(ifile)%add_var('v',    'yz crosssection of the south-north velocity',              'm/s',     '0mtt')
      call yz_files(ifile)%add_var('w',    'yz crosssection of the vertical velocity',                 'm/s',     '0tmt')
      call yz_files(ifile)%add_var('thl',  'yz crossection of the liquid water potential temperature', 'K',       '0ttt')
      call yz_files(ifile)%add_var('thv',  'yz crosssection of the virtual potential temperature',     'K',       '0ttt')
      call yz_files(ifile)%add_var('qt',   'yz crosssection of the total water specific humidity',     'kg/kg',   '0ttt')
      call yz_files(ifile)%add_var('ql',   'yz crosssection of the liquid water specific humidity',    'kg/kg',   '0ttt')
      call yz_files(ifile)%add_var('buoy', 'yz crosssection of the buoyancy',                          'K',       '0ttt')
      call yz_files(ifile)%add_var('e120', 'yz crosssection of the sqrt(turbulent kinetic energy',     'm^2/s^2', '0ttt')
    end do

  end subroutine initcrosssection

  !> Run crosssection.
  subroutine crosssection

    if (lcross .and. rk3step == 3) then
#if defined(_OPENACC)
      call update_host
#endif
      if (lxz) call wrtvert
      if (lxy) call wrthorz
      if (lyz) call wrtorth
    end if

  end subroutine crosssection


  !> Do the xz crosssections and dump them to file.
  subroutine wrtvert

    integer :: i, j, k, n, cross

    ! Lists of pointers, one per file. This is so that each cross section can be
    ! taken asyncronously on GPU.
    type(t_ptr_2d) :: u(100)
    type(t_ptr_2d) :: v(100)
    type(t_ptr_2d) :: w(100)
    type(t_ptr_2d) :: thl(100)
    type(t_ptr_2d) :: thv(100)
    type(t_ptr_2d) :: qt(100)
    type(t_ptr_2d) :: ql(100)
    type(t_ptr_2d) :: buoy(100)
    type(t_ptr_2d) :: e12(100)

    if (is_sampling_timestep(xz_file_ids(1))) then
      ! Setup the pointers
      ! CJ: this could be done one time during initialization, have to make sure that
      ! the buffer is allocated though...
      do cross = 1, nxz
        call xz_files(cross)%get_pointer('u', u(cross)%p)
        call xz_files(cross)%get_pointer('v', v(cross)%p)
        call xz_files(cross)%get_pointer('w', w(cross)%p)
        call xz_files(cross)%get_pointer('thl', thl(cross)%p)
        call xz_files(cross)%get_pointer('thv', thv(cross)%p)
        call xz_files(cross)%get_pointer('qt', qt(cross)%p)
        call xz_files(cross)%get_pointer('ql', ql(cross)%p)
        call xz_files(cross)%get_pointer('buoy', buoy(cross)%p)
        call xz_files(cross)%get_pointer('e120', e12(cross)%p)
      end do

      do cross = 1, nxz
        j = crossplane(cross)

        !$acc kernels default(present) async
        u(cross)%p(:,:) = um(2:i1,j,1:kmax) + cu
        v(cross)%p(:,:) = vm(2:i1,j,1:kmax) + cv
        w(cross)%p(:,:) = wm(2:i1,j,1:kmax)
        e12(cross)%p(:,:) = e12m(2:i1,j,1:kmax)
        !$acc end kernels

        !$acc parallel loop collapse(2) default(present) async
        do k = 1, kmax
          do i = 2, i1
            qt(cross)%p(i,k) = qtm(i,j,k)
            ql(cross)%p(i,k) = ql0(i,j,k)
            thl(cross)%p(i,k) = thlm(i,j,k)
            thv(cross)%p(i,k) = calc_virt_pot_temp(thlm(i,j,k), qtm(i,j,k), &
                                                   ql0(i,j,k), exnf(k))
            buoy(cross)%p(i,k) = thv(cross)%p(i,k) - thvf(k)
          end do
        end do
      end do
    end if

  end subroutine wrtvert

  !> Do the xy crosssections and dump them to file.
  subroutine wrthorz

    integer :: i, j, k, n, cross

    type(t_ptr_2d) :: u(100)
    type(t_ptr_2d) :: v(100)
    type(t_ptr_2d) :: w(100)
    type(t_ptr_2d) :: thl(100)
    type(t_ptr_2d) :: thv(100)
    type(t_ptr_2d) :: qt(100)
    type(t_ptr_2d) :: ql(100)
    type(t_ptr_2d) :: buoy(100)
    type(t_ptr_2d) :: e12(100)

    if (is_sampling_timestep(xy_file_ids(1))) then
      do cross = 1, nxy
        call xy_files(cross)%get_pointer('u', u(cross)%p)
        call xy_files(cross)%get_pointer('v', v(cross)%p)
        call xy_files(cross)%get_pointer('w', w(cross)%p)
        call xy_files(cross)%get_pointer('thl', thl(cross)%p)
        call xy_files(cross)%get_pointer('thv', thv(cross)%p)
        call xy_files(cross)%get_pointer('qt', qt(cross)%p)
        call xy_files(cross)%get_pointer('ql', ql(cross)%p)
        call xy_files(cross)%get_pointer('buoy', buoy(cross)%p)
        call xy_files(cross)%get_pointer('e120', e12(cross)%p)
      end do

      do cross = 1, nxy
        k = crossheight(cross)

        !$acc kernels default(present) async
        u(cross)%p(:,:) = um(2:i1,2:j1,k) + cu
        v(cross)%p(:,:) = vm(2:i1,2:j1,k) + cv
        w(cross)%p(:,:) = wm(2:i1,2:j1,k)
        e12(cross)%p(:,:) = e12m(2:i1,2:j1,k)
        !$acc end kernels

        !$acc parallel loop collapse(2) default(present) async
        do j = 2, j1
          do i = 2, i1
            qt(cross)%p(i,j) = qtm(i,j,k)
            ql(cross)%p(i,j) = ql0(i,j,k)
            thl(cross)%p(i,j) = thlm(i,j,k)
            thv(cross)%p(i,j) = calc_virt_pot_temp(thlm(i,j,k), qtm(i,j,k), &
                                                   ql0(i,j,k), exnf(k))
            buoy(cross)%p(i,j) = thv(cross)%p(i,j) - thvf(k)
          end do
        end do
      end do
    end if

  end subroutine wrthorz

  !> Do the yz crosssections and dump them to file.
  subroutine wrtorth

    integer :: i, j, k, n, cross

    type(t_ptr_2d) :: u(100)
    type(t_ptr_2d) :: v(100)
    type(t_ptr_2d) :: w(100)
    type(t_ptr_2d) :: thl(100)
    type(t_ptr_2d) :: thv(100)
    type(t_ptr_2d) :: qt(100)
    type(t_ptr_2d) :: ql(100)
    type(t_ptr_2d) :: buoy(100)
    type(t_ptr_2d) :: e12(100)

    if (is_sampling_timestep(yz_file_ids(1))) then
      do cross = 1, nyz
        call yz_files(cross)%get_pointer('u', u(cross)%p)
        call yz_files(cross)%get_pointer('v', v(cross)%p)
        call yz_files(cross)%get_pointer('w', w(cross)%p)
        call yz_files(cross)%get_pointer('thl', thl(cross)%p)
        call yz_files(cross)%get_pointer('thv', thv(cross)%p)
        call yz_files(cross)%get_pointer('qt', qt(cross)%p)
        call yz_files(cross)%get_pointer('ql', ql(cross)%p)
        call yz_files(cross)%get_pointer('buoy', buoy(cross)%p)
        call yz_files(cross)%get_pointer('e120', e12(cross)%p)
      end do

      do cross = 1, nyz
        i = crossortho(cross)

        !$acc kernels default(present) async
        u(cross)%p(:,:) = um(i,2:j1,1:kmax) + cu
        v(cross)%p(:,:) = vm(i,2:j1,1:kmax) + cv
        w(cross)%p(:,:) = wm(i,2:j1,1:kmax)
        e12(cross)%p(:,:) = e12m(i,2:j1,1:kmax)
        !$acc end kernels

        !$acc parallel loop collapse(2) default(present) async
        do k = 1, kmax
          do j = 2, j1
            qt(cross)%p(j,k) = qtm(i,j,k)
            ql(cross)%p(j,k) = ql0(i,j,k)
            thl(cross)%p(j,k) = thlm(i,j,k)
            thv(cross)%p(j,k) = calc_virt_pot_temp(thlm(i,j,k), qtm(i,j,k), &
                                                   ql0(i,j,k), exnf(k))
            buoy(cross)%p(j,k) = thv(cross)%p(i,j) - thvf(k)
          end do
        end do
      end do
    end if

  end subroutine wrtorth

end module modcrosssection
