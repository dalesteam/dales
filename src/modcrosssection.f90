!> \file modcrosssection.f90
!!   Dumps an instantenous crosssection of the field

!>
!! Dumps an instantenous crosssection of the field.
!>
!! Crosssections in the yz-plane and in the xy-plane            |
    !        of u,v,w,thl,thv,qt,ql. Written to movv_*.expnr and movh_*.expnr
!! If netcdf is true, this module leads the cross.myid.expnr.nc output

!!  \par Revision list
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
module modcrosssection

  use fortran_support, only: &
#if FIELD_PRECISION==64
                               t_ptr_2d => t_ptr_2d_dp
#else
                               t_ptr_2d => t_ptr_2d_sp
#endif
  use modglobal, only  : longint, kmax, nsv
  use modtracers, only : tracer_prop
  use modlogging, only: finish
  use modnetcdf_file_t, only: cross_section_file_t
  use modstat_nc_files, only: add_output_file, is_sampling_timestep
  use modprecision, only: field_r
  use modthermodynamics, only: calc_virt_pot_temp

implicit none
character(len=*), parameter :: modname = 'modcrosssection'
private
PUBLIC :: initcrosssection, crosssection
save
  integer :: crossheight(100)
  integer :: nxy = 0, nxz = 0, nyz = 0
  integer :: cross
  character(4) :: cheight

  real    :: dtav
  integer(kind=longint) :: idtav,tnext
  logical :: lcross = .false. !< switch for doing the crosssection (on/off)
  logical :: lbinary = .false. !< switch for doing the crosssection (on/off)
  integer :: crossplane(100) !< Location of the xz crosssection
  integer :: crossortho(100) !< Location of the yz crosssection

  logical :: lxy = .true.   !< switch for doing xy crosssections
  logical :: lxz = .true.   !< switch for doing xz crosssections
  logical :: lyz = .true.   !< switch for doing yz crosssections

  type(cross_section_file_t), allocatable :: xy_files(:)    !< List of xy cross files.
  type(cross_section_file_t), allocatable :: xz_files(:)    !< List of xz cross files.
  type(cross_section_file_t), allocatable :: yz_files(:)    !< List of yz cross files.
  integer,                    allocatable :: xy_file_ids(:) !< List of xy cross file ids.
  integer,                    allocatable :: yz_file_ids(:) !< List of yz cross file ids.
  integer,                    allocatable :: xz_file_ids(:) !< List of xz cross file ids.

contains
!> Initializing Crosssection. Read out the namelist, initializing the variables
  subroutine initcrosssection
    use modmpi,   only :myid,mpierr,comm3d,cmyid,cmyidx,cmyidy,myidx,myidy,D_MPI_BCAST
    use modglobal,only :imax,jmax,itot,jtot,ifnamopt,fname_options,dtmax,dtav_glob,ladaptive,j1,kmax,i1,dt_lim,cexpnr,&
                        tres,btime,checknamelisterror,output_prefix
    use modstat_nc,only : lnetcdf,open_nc, define_nc,ncinfo,nctiminfo,writestat_dims_nc
    use fortran_support, only: nnml_output

   implicit none
    character(len=*), parameter :: routine = modname//'/initcrosssection'
    integer :: ierr,k,n
    integer :: ifile

    namelist/NAMCROSSSECTION/ &
    lcross, lbinary, dtav, crossheight, crossplane, crossortho, lxy, lxz, lyz


    crossheight(1)=2
    crossheight(2:100)=-999
    crossplane(1)=2
    crossplane(2:100)=-999
    crossortho(1)=2
    crossortho(2:100)=-999

    dtav = dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMCROSSSECTION,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMCROSSSECTION')
      write(nnml_output ,NAMCROSSSECTION)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(dtav       ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lcross     ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lbinary    ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(crossheight(1:100),100,0,comm3d,mpierr)
    call D_MPI_BCAST(crossplane(1:100),100,0,comm3d,mpierr)
    call D_MPI_BCAST(crossortho(1:100),100,0,comm3d,mpierr)
    call D_MPI_BCAST(lxy        ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lxz        ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lyz        ,1,0,comm3d,mpierr)

    if(any((crossheight(1:100).gt.kmax)) .or. any(crossplane > jtot+1) .or. any(crossortho > itot+1) ) then
      call finish(routine, 'CROSSSECTION: crosssection out of range')
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
        write(cheight, '(i4.4)') crossheight(k)
        xy_files(ifile) = cross_section_file_t('crossxy.'//cheight, nx=itot, &
                                               ny=jtot)
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
        write(cheight, '(i4.4)') crossplane(k)
        xz_files(ifile) = cross_section_file_t('crossxz.'//cheight, &
                                               nx=itot, nz=kmax)
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
        write(cheight, '(i4.4)') crossortho(k)
        yz_files(ifile) = cross_section_file_t('crossyz.'//cheight, &
                                               ny=jtot, nz=kmax)
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
!>Run crosssection. Mainly timekeeping
  subroutine crosssection
    use modglobal, only : rk3step,timee,dt_lim
    use modstat_nc, only : writestat_nc
#if defined(_OPENACC)
    use modgpu, only: update_host
#endif
    implicit none


    if (.not. lcross) return
    if (rk3step/=3) return
#if defined(_OPENACC)
    call update_host
#endif
    if (lxz) call wrtvert
    if (lxy) call wrthorz
    if (lyz) call wrtorth

  end subroutine crosssection


!> Do the xz crosssections and dump them to file
  subroutine wrtvert
  use modglobal, only : imax,i1,j1,kmax,nsv,rlv,cp,rv,rd,cu,cv,cexpnr,ifoutput,rtimee
  use modfields, only : um,vm,wm,thlm,qtm,svm,thl0,qt0,ql0,e120,exnf,thvf
  use modmpi,    only : myidy
  use modstat_nc, only : lnetcdf, writestat_nc
  implicit none

  integer i,k,n,isv, j, cross

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
      e12(cross)%p(:,:) = e120(2:i1,j,1:kmax)
      !$acc end kernels

      !$acc parallel loop collapse(2) default(present) async
      do k = 1, kmax
        do i = 2, i1
          qt(cross)%p(i,k) = qt0(i,j,k)
          ql(cross)%p(i,k) = ql0(i,j,k)
          thl(cross)%p(i,k) = thl0(i,j,k)
          thv(cross)%p(i,k) = calc_virt_pot_temp(thl0(i,j,k), qt0(i,j,k), &
                                                 ql0(i,j,k), exnf(k))
          buoy(cross)%p(i,k) = thv(cross)%p(i,k) - thvf(k)
        end do
      end do
    end do
  
  end if

  end subroutine wrtvert

!> Do the xy crosssections and dump them to file
  subroutine wrthorz
    use modglobal, only : imax,jmax,i1,j1,nsv,rlv,cp,rv,rd,cu,cv,cexpnr,ifoutput,rtimee
    use modfields, only : um,vm,wm,thlm,qtm,svm,thl0,qt0,ql0,e120,exnf,thvf
    use modmpi,    only : cmyid
    use modstat_nc, only : lnetcdf, writestat_nc
    implicit none


    ! LOCAL
    integer i,j,n,isv, cross,k

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

    ! Setup the pointers
    ! CJ: this could be done one time during initialization, have to make sure that
    ! the buffer is allocated though...
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
      e12(cross)%p(:,:) = e120(2:i1,2:j1,k)
      !$acc end kernels

      !$acc parallel loop collapse(2) default(present) async
      do j = 2, j1
        do i = 2, i1
          qt(cross)%p(i,j) = qt0(i,j,k)
          ql(cross)%p(i,j) = ql0(i,j,k)
          thl(cross)%p(i,j) = thl0(i,j,k)
          thv(cross)%p(i,j) = calc_virt_pot_temp(thl0(i,j,k), qt0(i,j,k), &
                                                 ql0(i,j,k), exnf(k))
          buoy(cross)%p(i,j) = thv(cross)%p(i,j) - thvf(k)
        end do
      end do
    end do

  end if

  end subroutine wrthorz

  ! yz cross section
  subroutine wrtorth
    use modglobal, only : jmax,kmax,i1,j1,nsv,rlv,cp,rv,rd,cu,cv,cexpnr,ifoutput,rtimee
    use modfields, only : um,vm,wm,thlm,qtm,svm,thl0,qt0,ql0,e120,exnf,thvf
    use modmpi,    only : cmyid, myidx
    use modstat_nc, only : lnetcdf, writestat_nc
    implicit none


    ! LOCAL
    integer j,k,n,isv,i
    character(21) :: name

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

    ! Setup the pointers
    ! CJ: this could be done one time during initialization, have to make sure that
    ! the buffer is allocated though...
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
      e12(cross)%p(:,:) = e120(i,2:j1,1:kmax)
      !$acc end kernels

      !$acc parallel loop collapse(2) default(present) async
      do k = 1, kmax
        do j = 2, j1
          qt(cross)%p(j,k) = qt0(i,j,k)
          ql(cross)%p(j,k) = ql0(i,j,k)
          thl(cross)%p(j,k) = thl0(i,j,k)
          thv(cross)%p(j,k) = calc_virt_pot_temp(thl0(i,j,k), qt0(i,j,k), &
                                                 ql0(i,j,k), exnf(k))
          buoy(cross)%p(j,k) = thv(cross)%p(i,j) - thvf(k)
        end do
      end do
    end do

  end if

  end subroutine wrtorth

end module modcrosssection
