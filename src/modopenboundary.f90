!> \file modopenboundary.f90
!!  Sets ghost cells at domain boundary and implements radiation
!!  conditions for the boundary normal velocity components
!>
!!  Sets ghost cells at domain boundary and implements radiation
!!  conditions for the boundary normal velocity components
!>
!!  \author Frans Liqui Lung
!  This file is part of DALES.
!  To do:
!  - Allow for different zint and adjust rhointi accordingly
!  - Correct non divergence free input more elegently
!  - Change definition uphase for division by 0
!  - When to use nextval and currentval for nudging and check rtimee
!  - How to handle vertical derivative in top boundary condition (full levels) now obtained from profile (horizontal average)
!  - Use correct velocity level to determine in or outflow in full levels, half levels are used now
!  - Check rtimee and half level nudgin and correction term
!  - Use um or u0 in half levels
!  - Add possibility for higher order integration schemes
!  - Extent turbulent pertubation generation options
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
module modopenboundary
use modglobal, only : boundary_type,boundary,lopenbc,linithetero,lboundary,lperiodic,lsynturb,dzint,dxint,dyint,ntboundary,tboundary,dxturb,dyturb
use modsynturb, only : synturb,initsynturb,exitsynturb
use netcdf
use modprecision, only : field_r
use modmpi, only : openboundary_excjs

implicit none
integer :: nxpatch, nypatch, nzpatch, nxturb, nyturb, nzturb
real, dimension(:,:), allocatable :: uturbtemp,vturbtemp,wturbtemp
real, dimension(:), allocatable :: rhointi
real, dimension(:,:,:), allocatable :: thls_hetero,ps_hetero
real :: ibuoy = 0.
logical :: lbuoytop = .false.
integer :: nx1max,nx2max

contains
  subroutine initopenboundary
    ! Initialisation routine for openboundaries
    use modmpi, only : myidx, myidy, nprocx, nprocy,myid
    use modglobal, only : imax,jmax,kmax,i1,j1,k1,dx,dy,dzf,itot,jtot,zf,zh,solver_id, &
      & iadv_mom,iadv_thl,iadv_qt,iadv_tke,iadv_sv,nsv,cu,cv
    use modsurfdata, only : isurf
    implicit none
    integer :: i,j

    if(.not.lopenbc) return
    ! Check for conflicting options
    if(solver_id == 0) stop 'Openboundaries only possible with HYPRE or FFTW pressure solver, change solver_id'
    if(iadv_mom /=2) stop 'Only second order advection scheme supported with openboundaries, change iadv_mom to 2'
    if(iadv_thl /=2) stop 'Only second order advection scheme supported with openboundaries, change iadv_thl to 2'
    if(iadv_qt  /=2) stop 'Only second order advection scheme supported with openboundaries, change iadv_qt to 2'
    if(iadv_tke /=2) stop 'Only second order advection scheme supported with openboundaries, change iadv_tke to 2'
    if(any(iadv_sv(1:nsv)/=2)) stop 'Only second order advection scheme supported with openboundaries, change iadv_sv to 2'
    if(cu/=0.) stop 'Translation velocity not allowed in combination with open boundaries, set cu to 0'
    if(cv/=0.) stop 'Translation velocity not allowed in combination with open boundaries, set cv to 0'
    ! Set buoyancy term at top boundary on or off (off default)
    if(lbuoytop) ibuoy = 1.
    ! Check if boundary is present on process
    if(myidx==0)        lboundary(1) = .true.
    if(myidx==nprocx-1) lboundary(2) = .true.
    if(myidy==0)        lboundary(3) = .true.
    if(myidy==nprocy-1) lboundary(4) = .true.
    lboundary(5) = .true.
    ! Set boundary names
    boundary(1)%name = "west"
    boundary(2)%name = "east"
    boundary(3)%name = "south"
    boundary(4)%name = "north"
    boundary(5)%name = "top"
    ! Set dimensions for each boundary
    boundary(1)%nx1  = jmax; boundary(1)%nx2  = kmax
    boundary(2)%nx1  = jmax; boundary(2)%nx2  = kmax
    boundary(3)%nx1  = imax; boundary(3)%nx2  = kmax
    boundary(4)%nx1  = imax; boundary(4)%nx2  = kmax
    boundary(5)%nx1  = imax; boundary(5)%nx2  = jmax
    boundary(1)%nx1u = jmax; boundary(1)%nx2u = kmax
    boundary(2)%nx1u = jmax; boundary(2)%nx2u = kmax
    boundary(3)%nx1u = i1;   boundary(3)%nx2u = kmax
    boundary(4)%nx1u = i1;   boundary(4)%nx2u = kmax
    boundary(5)%nx1u = i1;   boundary(5)%nx2u = jmax
    boundary(1)%nx1v = j1;   boundary(1)%nx2v = kmax
    boundary(2)%nx1v = j1;   boundary(2)%nx2v = kmax
    boundary(3)%nx1v = imax; boundary(3)%nx2v = kmax
    boundary(4)%nx1v = imax; boundary(4)%nx2v = kmax
    boundary(5)%nx1v = imax; boundary(5)%nx2v = j1
    boundary(1)%nx1w = jmax; boundary(1)%nx2w = k1
    boundary(2)%nx1w = jmax; boundary(2)%nx2w = k1
    boundary(3)%nx1w = imax; boundary(3)%nx2w = k1
    boundary(4)%nx1w = imax; boundary(4)%nx2w = k1
    boundary(5)%nx1w = imax; boundary(5)%nx2w = jmax
    ! Set number of patches for correction factor for radiation boundary conditions
    if(dxint == -1.) dxint = real(itot)*dx ! Set dxint to entire width as default
    if(dyint == -1.) dyint = real(jtot)*dy ! Set dyint to entire width as default
    if(dxturb == -1.) dxturb = real(itot)*dx ! Set dxint to entire width as default
    if(dyturb == -1.) dyturb = real(jtot)*dy ! Set dyint to entire width as default
    nxpatch = int(dx/dxint*real(itot));
    nypatch = int(dy/dyint*real(jtot));
    nzpatch = kmax ! For now vertical integration scale is set equal to dz
    nxturb = int(dx/dxturb*real(itot));
    nyturb = int(dy/dyturb*real(jtot));
    nzturb = kmax ! For now vertical resolution turbulence input must equal dz
    if(mod(dxint,dx)/=0 .or. mod(dyint,dy)/=0) stop 'dxint and dyint should be multiples of dx and dy respectively.'
    if(mod(dxturb,dx)/=0 .or. mod(dyturb,dy)/=0) stop 'dxturb and dyturb should be multiples of dx and dy respectively.'
    boundary(1)%nx1patch = nypatch; boundary(1)%nx2patch = nzpatch
    boundary(2)%nx1patch = nypatch; boundary(2)%nx2patch = nzpatch
    boundary(3)%nx1patch = nxpatch; boundary(3)%nx2patch = nzpatch
    boundary(4)%nx1patch = nxpatch; boundary(4)%nx2patch = nzpatch
    boundary(5)%nx1patch = nxpatch; boundary(5)%nx2patch = nypatch
    boundary(1)%nx1turb = nyturb; boundary(1)%nx2turb = nzturb
    boundary(2)%nx1turb = nyturb; boundary(2)%nx2turb = nzturb
    boundary(3)%nx1turb = nxturb; boundary(3)%nx2turb = nzturb
    boundary(4)%nx1turb = nxturb; boundary(4)%nx2turb = nzturb
    boundary(5)%nx1turb = nxturb; boundary(5)%nx2turb = nyturb
    if(myid==0) print *,"dxint/dx,dyint/dy,nxpatch,nypatch",int(dxint/dx),int(dyint/dy),nxpatch,nypatch
    ! Allocate phase velocity, correction term radiation boundaries and pertubation fields
    do i = 1,5
      if(.not.lboundary(i) .or. lperiodic(i)) cycle ! Open boundary not present
      allocate(boundary(i)%radcorr(boundary(i)%nx1patch,boundary(i)%nx2patch), &
        boundary(i)%radcorrsingle(boundary(i)%nx1patch,boundary(i)%nx2patch), &
        boundary(i)%uphase(boundary(i)%nx1patch,boundary(i)%nx2patch), &
        boundary(i)%uphasesingle(boundary(i)%nx1patch,boundary(i)%nx2patch), &
        boundary(i)%uturb(boundary(i)%nx1u,boundary(i)%nx2u), &
        boundary(i)%vturb(boundary(i)%nx1v,boundary(i)%nx2v), &
        boundary(i)%wturb(boundary(i)%nx1w,boundary(i)%nx2w), &
        boundary(i)%thlturb(boundary(i)%nx1,boundary(i)%nx2), &
        boundary(i)%qtturb(boundary(i)%nx1,boundary(i)%nx2), &
        boundary(i)%e12turb(boundary(i)%nx1,boundary(i)%nx2))
        boundary(i)%uturb = 0.; boundary(i)%vturb = 0.; boundary(i)%wturb = 0.
        boundary(i)%thlturb = 0.; boundary(i)%qtturb = 0.; boundary(i)%e12turb = 0.
        if(nsv>0) then
          allocate(boundary(i)%svturb(boundary(i)%nx1,boundary(i)%nx2,nsv))
          boundary(i)%svturb = 0.
        endif
    end do
    call initsynturb
  end subroutine initopenboundary

  subroutine exitopenboundary
    ! Exit routine for openboundaries
    use modglobal, only : nsv
    implicit none
    integer :: i
    if(.not.lopenbc) return
    deallocate(tboundary)
    do i = 1,5
      if(.not.lboundary(i) .or. lperiodic(i)) cycle
      deallocate(boundary(i)%thl,boundary(i)%qt,boundary(i)%e12, &
        boundary(i)%u,boundary(i)%v,boundary(i)%w,boundary(i)%uphasesingle,boundary(i)%uphase, &
        boundary(i)%radcorr,boundary(i)%radcorrsingle,boundary(i)%uturb,boundary(i)%vturb, &
          boundary(i)%wturb,boundary(i)%thlturb,boundary(i)%qtturb,boundary(i)%name,boundary(i)%e12turb)
      if(nsv>0) deallocate(boundary(i)%sv,boundary(i)%svturb)
    end do
    deallocate(rhointi)
    call exitsynturb
  end subroutine exitopenboundary

  subroutine openboundary_initfields
    ! Routine that reads the fields for a heterogneous initialisation
    ! Variables not present in the input file are initialised to their profiles
    use modglobal, only : imax,jmax,kmax,i1,j1,cexpnr,nsv,i2,j2,k1
    use modfields, only : u0,um,v0,vm,w0,wm,thl0,thlm,qt0,qtm,e120,e12m, sv0, svm, &
      uprof,vprof,thlprof,qtprof,e12prof,svprof
    use modmpi, only : myidx,myidy,myid
    implicit none
    integer :: VARID,STATUS,NCID,mpierr,timeID,n
    integer, dimension(3) :: istart
    if(.not.lopenbc) return
    if(.not.linithetero) return
    u0 = 0.; um = 0.; v0 = 0.; vm = 0.; w0 = 0.; wm = 0.;
    thl0 = 0.; thlm = 0.; qt0 = 0; qtm = 0; e120 = 0.; e12m = 0.
    !--- open initfields.input.xxx.nc ---
    STATUS = NF90_OPEN('initfields.inp.'//cexpnr//'.nc', nf90_nowrite, NCID)
    if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
    istart = (/myidx*imax+1,myidy*jmax+1,1/)
    STATUS = NF90_INQ_VARID(NCID,'u0', VARID)
    if(STATUS == NF90_ENOTVAR) then ! If not available initialise with profile
      call take_prof(u0,um,uprof)
      if(myid==0) print *, "u initialized with prof.inp"
    else
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID,u0(2:i2,2:j1,1:kmax),start=istart,count=(/i1,jmax,kmax/))
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      um = u0
    endif
    STATUS = NF90_INQ_VARID(NCID,'v0', VARID)
    if(STATUS == NF90_ENOTVAR) then ! If not available initialise with profile
      call take_prof(v0,vm,vprof)
      if(myid==0) print *, "v initialized with prof.inp"
    else
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID,v0(2:i1,2:j2,1:kmax),start=istart,count=(/imax,j1,kmax/))
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      vm = v0
    endif
    STATUS = NF90_INQ_VARID(NCID,'w0', VARID)
    if(STATUS == NF90_ENOTVAR) then ! If not available initialise with 0.
      wm = 0.
      w0 = 0.
      if(myid==0) print *, "w initilised with 0."
    else
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID,w0(2:i1,2:j1,1:k1),start=istart,count=(/imax,jmax,k1/))
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      wm = w0
    endif
    STATUS = NF90_INQ_VARID(NCID,'thl0', VARID)
    if(STATUS == NF90_ENOTVAR) then ! If not available initialise with profile
      call take_prof(thl0,thlm,thlprof)
      if(myid==0) print *, "thl initialized with prof.inp"
    else
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID,thl0(2:i1,2:j1,1:kmax),start=istart,count=(/imax,jmax,kmax/))
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      thlm = thl0
    endif
    STATUS = NF90_INQ_VARID(NCID,'qt0', VARID)
    if(STATUS == NF90_ENOTVAR) then ! If not available initialise with profile
      call take_prof(qt0,qtm,qtprof)
      if(myid==0) print *, "qt initialized with prof.inp"
    else
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID,qt0(2:i1,2:j1,1:kmax),start=istart,count=(/imax,jmax,kmax/))
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      qtm = qt0
    endif
    STATUS = NF90_INQ_VARID(NCID,'e120', VARID)
    if(STATUS == NF90_ENOTVAR) then ! If not available initialise with profile
      call take_prof(e120,e12m,e12prof)
      if(myid==0) print *, "e12 initialized with prof.inp"
    else
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID,e120(2:i1,2:j1,1:kmax),start=istart,count=(/imax,jmax,kmax/))
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      e12m = e120
    endif
    if(nsv>0) then
      STATUS = NF90_INQ_VARID(NCID,'sv0', VARID)
      if(STATUS == NF90_ENOTVAR) then ! If not available initialise with profile
        do n = 1,nsv
          call take_prof(sv0(:,:,:,n),svm(:,:,:,n),svprof(:,n))
        end do
        if(myid==0) print *, "sv initialized with scalar.inp"
      else ! Initial field taken from initfields.inp
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID,sv0(2:i1,2:j1,1:kmax,1),start=istart,count=(/imax,jmax,kmax,nsv/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        svm = sv0
      endif
    endif
    status = nf90_close(ncid)
    if (status /= nf90_noerr) call handle_err(status)
  end subroutine openboundary_initfields

  subroutine openboundary_readboundary
    ! Routine reads the boundary input for all time steps
    use modglobal, only : kmax,cexpnr,imax,jmax,itot,jtot,k1,ntboundary, &
      tboundary,i1,j1,i2,j2,kmax,nsv,iturb
    use modmpi, only : myidx,myidy
    implicit none
    integer :: it,i,j,k,ib
    character(len = nf90_max_name) :: RecordDimName
    integer :: VARID,STATUS,NCID,mpierr,timeID
    integer, dimension(3) :: istart

    if(.not.lopenbc) return
    !--- open openboundaries.inp.xxx.nc ---
    STATUS = NF90_OPEN('openboundaries.inp.'//cexpnr//'.nc', nf90_nowrite, NCID)
    if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
    !--- get time dimensions
    status = nf90_inq_dimid(ncid, "time", timeID)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inquire_dimension(NCID, timeID, len=ntboundary, name=RecordDimName)
    if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
    !--- read time
    allocate(tboundary(ntboundary))
    STATUS = NF90_INQ_VARID(NCID, 'time', VARID)
    if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
    STATUS = NF90_GET_VAR (NCID, VARID, tboundary, start=(/1/), count=(/ntboundary/) )
    if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
    do ib = 1,5 ! loop over boundaries
      ! Allocate input fields
      if(.not.lboundary(ib) .or. lperiodic(ib)) cycle ! Open boundary not present
      allocate(boundary(ib)%thl(boundary(ib)%nx1,boundary(ib)%nx2,ntboundary), &
        boundary(ib)%qt(boundary(ib)%nx1,boundary(ib)%nx2,ntboundary),  &
        boundary(ib)%e12(boundary(ib)%nx1,boundary(ib)%nx2,ntboundary), &
        boundary(ib)%u(boundary(ib)%nx1u,boundary(ib)%nx2u,ntboundary), &
        boundary(ib)%v(boundary(ib)%nx1v,boundary(ib)%nx2v,ntboundary), &
        boundary(ib)%w(boundary(ib)%nx1w,boundary(ib)%nx2w,ntboundary))
      if(nsv>0) then
        allocate(boundary(ib)%sv(boundary(ib)%nx1,boundary(ib)%nx2,ntboundary,nsv))
      endif
      if(lsynturb .and. iturb<10) then ! Allocate turbulent input fields
        allocate(boundary(ib)%u2(boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary), &
        & boundary(ib)%v2(boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary), &
        & boundary(ib)%w2(boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary), &
        & boundary(ib)%uv(boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary), &
        & boundary(ib)%uw(boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary), &
        & boundary(ib)%vw(boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary), &
        & boundary(ib)%thl2(boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary),&
        & boundary(ib)%qt2(boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary), &
        & boundary(ib)%wthl(boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary), &
        & boundary(ib)%wqt(boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary))
      endif
      ! Read fields
      select case(ib)
      case(1,2)
        istart = (/myidy*jmax+1,1,1/)
      case(3,4)
        istart = (/myidx*imax+1,1,1/)
      case(5)
        istart = (/myidx*imax+1,myidy*jmax+1,1/)
      end select
      ! Read u
      STATUS = NF90_INQ_VARID(NCID,'u'//boundary(ib)%name, VARID)
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%u, start=istart, &
        & count=(/boundary(ib)%nx1u,boundary(ib)%nx2u,ntboundary/))
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      ! Read v
      STATUS = NF90_INQ_VARID(NCID, 'v'//boundary(ib)%name, VARID)
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%v, start=istart, &
        & count=(/boundary(ib)%nx1v,boundary(ib)%nx2v,ntboundary/))
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      ! Read w
      STATUS = NF90_INQ_VARID(NCID, 'w'//boundary(ib)%name, VARID)
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%w, start=istart, &
        & count=(/boundary(ib)%nx1w,boundary(ib)%nx2w,ntboundary/))
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      ! Read thl
      STATUS = NF90_INQ_VARID(NCID, 'thl'//boundary(ib)%name, VARID)
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%thl, start=istart, &
        & count=(/boundary(ib)%nx1,boundary(ib)%nx2,ntboundary/))
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      ! Read qt
      STATUS = NF90_INQ_VARID(NCID, 'qt'//boundary(ib)%name, VARID)
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%qt, start=istart, &
        & count=(/boundary(ib)%nx1,boundary(ib)%nx2,ntboundary/))
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      ! Read e12
      STATUS = NF90_INQ_VARID(NCID, 'e12'//boundary(ib)%name, VARID)
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%e12, start=istart, &
        & count=(/boundary(ib)%nx1,boundary(ib)%nx2,ntboundary/))
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      ! Read sv
      if(nsv>0) then
        STATUS = NF90_INQ_VARID(NCID, 'sv'//boundary(ib)%name, VARID)
        if(STATUS == NF90_ENOTVAR) then
          boundary(ib)%sv = 0.
          if(myidx==0 .and. myidy==0) print *, "No boundary information for sv at boundary",ib,"Values set to 0"
        else
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%sv, start=(/istart,1/), &
            & count=(/boundary(ib)%nx1,boundary(ib)%nx2,ntboundary,nsv/))
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        endif
      endif
      ! Read input for turbulent pertubations
      if(lsynturb .and. iturb <10) then
        ! Read u2
        STATUS = NF90_INQ_VARID(NCID, 'u2'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%u2, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! Read v2
        STATUS = NF90_INQ_VARID(NCID, 'v2'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%v2, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! Read w2
        STATUS = NF90_INQ_VARID(NCID, 'w2'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%w2, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! Read uv
        STATUS = NF90_INQ_VARID(NCID, 'uv'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%uv, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! Read uw
        STATUS = NF90_INQ_VARID(NCID, 'uw'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%uw, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! Read vw
        STATUS = NF90_INQ_VARID(NCID, 'vw'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%vw, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! Read thl2
        STATUS = NF90_INQ_VARID(NCID, 'thl2'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%thl2, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! Read qt2
        STATUS = NF90_INQ_VARID(NCID, 'qt2'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%qt2, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! Read wthl
        STATUS = NF90_INQ_VARID(NCID, 'wthl'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%wthl, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! Read wqt
        STATUS = NF90_INQ_VARID(NCID, 'wqt'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%wqt, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1turb,boundary(ib)%nx2turb,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      endif
    end do
    status = nf90_close(ncid)
    if (status /= nf90_noerr) call handle_err(status)
  end subroutine openboundary_readboundary

  subroutine openboundary_divcorr()
    ! Correct for any integrated divergence present in the boundary input
    use modmpi, only : myid,comm3d,mpierr,d_mpi_allreduce,MPI_SUM
    use modglobal, only : imax,jmax,kmax,dzf,dy,dx,xsize,ysize,zh,k1,lwarmstart, &
      i1,i2,j1,j2,k1,dzh
    use modfields, only : u0,um,v0,vm,w0,wm,rhobf,rhobh
    implicit none
    real(field_r) :: sumdiv,div,divpart,divnew,divold
    integer :: i,j,k,it,icalc
    ! Create 1/int(rho)
    allocate(rhointi(k1))
    rhointi = 1./(rhobf*dzf)
    ! Divergence correction
    if(myid==0) print *, "Start divergence correction boundaries"
    do it = 1,ntboundary
      do icalc=1,2
        ! Calculate divergence
        div = 0.
        if(lboundary(1)) then
          do j = 1,jmax
            do k = 1,kmax
              div = div - rhobf(k)*boundary(1)%u(j,k,it)*dzf(k)*dy
            end do
          end do
        endif
        if(lboundary(2)) then
          do j = 1,jmax
            do k = 1,kmax
              div = div + rhobf(k)*boundary(2)%u(j,k,it)*dzf(k)*dy
            end do
          end do
        endif
        if(lboundary(3)) then
          do i = 1,imax
            do k = 1,kmax
              div = div - rhobf(k)*boundary(3)%v(i,k,it)*dzf(k)*dx
            end do
          end do
        endif
        if(lboundary(4)) then
          do i = 1,imax
            do k = 1,kmax
              div = div + rhobf(k)*boundary(4)%v(i,k,it)*dzf(k)*dx
            end do
          end do
        endif
        if(lboundary(5)) then
          do i = 1,imax
            do j = 1,jmax
              div = div + rhobh(k1)*boundary(5)%w(i,j,it)*dx*dy
            end do
          end do
        endif
        call D_MPI_ALLREDUCE(div,sumdiv,1,MPI_SUM,comm3d,mpierr)
        if(icalc==1) then
           divold=sumdiv
         else
           divnew=sumdiv
           if(myid==0) print *, 't,input,corrected;',tboundary(it),divold,divnew
           exit
         endif
        ! Apply correction, spread divergence over lateral boundaries
        if(lboundary(1)) then
          do k = 1,kmax
            divpart = sumdiv*ysize*dzf(k)/(2*xsize*zh(k1)+2*ysize*zh(k1))
            boundary(1)%u(:,k,it)=boundary(1)%u(:,k,it)+divpart/(rhobf(k)*ysize*dzf(k))
          end do
        endif
        if(lboundary(2)) then
          do k = 1,kmax
            divpart = sumdiv*ysize*dzf(k)/(2*xsize*zh(k1)+2*ysize*zh(k1))
            boundary(2)%u(:,k,it)=boundary(2)%u(:,k,it)-divpart/(rhobf(k)*ysize*dzf(k))
          end do
        endif
        if(lboundary(3)) then
          do k = 1,kmax
            divpart = sumdiv*xsize*dzf(k)/(2*xsize*zh(k1)+2*ysize*zh(k1))
            boundary(3)%v(:,k,it)=boundary(3)%v(:,k,it)+divpart/(rhobf(k)*dzf(k)*xsize)
          end do
        endif
        if(lboundary(4)) then
          do k = 1,kmax
            divpart = sumdiv*xsize*dzf(k)/(2*xsize*zh(k1)+2*ysize*zh(k1))
            boundary(4)%v(:,k,it)=boundary(4)%v(:,k,it)-divpart/(rhobf(k)*xsize*dzf(k))
          enddo
        endif
      end do
    end do
    if(myid==0) print *, "Finished divergence correction boundaries"
    ! Copy data to boundary information
    if(.not.lwarmstart) then
      if(lboundary(1).and..not.lperiodic(1)) then
        do j = 2,j1
          do k = 1,kmax
            u0(2,j,k) = boundary(1)%u(j-1,k,1)
            um(2,j,k) = boundary(1)%u(j-1,k,1)
          end do
        end do
      endif
      if(lboundary(2).and..not.lperiodic(2)) then
        do j = 2,j1
          do k = 1,kmax
            u0(i2,j,k) = boundary(2)%u(j-1,k,1)
            um(i2,j,k) = boundary(2)%u(j-1,k,1)
          end do
        end do
      endif
      if(lboundary(3).and..not.lperiodic(3)) then
        do i = 2,i1
          do k = 1,kmax
            v0(i,2,k) = boundary(3)%v(i-1,k,1)
            vm(i,2,k) = boundary(3)%v(i-1,k,1)
          end do
        end do
      endif
      if(lboundary(4).and..not.lperiodic(4)) then
        do i = 2,i1
          do k = 1,kmax
            v0(i,j2,k) = boundary(4)%v(i-1,k,1)
            vm(i,j2,k) = boundary(4)%v(i-1,k,1)
          end do
        end do
      endif
      if(lboundary(5).and..not.lperiodic(5)) then
        do i = 2,i1
          do j = 2,j1
            w0(i,j,k1) = boundary(5)%w(i-1,j-1,1)
            wm(i,j,k1) = boundary(5)%w(i-1,j-1,1)
          end do
        end do
      endif
    endif
  end subroutine openboundary_divcorr

  subroutine openboundary_ghost
    ! Subroutine that fills the ghost cells for the cell centred variables at the boundary
    use modglobal, only : i1,j1,k1,ih,jh,nsv
    use modfields, only : um,u0,vm,v0,wm,w0,e12m,e120,thlm,thl0,qtm,qt0,svm,sv0,thl0av,qt0av,u0av,v0av
    implicit none
    integer :: i,n
    if(.not.lopenbc) return
    ! Apply non domain boundaries and periodic boundaries
    call openboundary_excjs(um   , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(u0   , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(vm   , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(v0   , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(wm   , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(w0   , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(e12m , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(e120 , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(thlm , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(thl0 , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(qtm  , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(qt0  , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    do n = 1,nsv
      call openboundary_excjs(svm(:,:,:,n), 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
      call openboundary_excjs(sv0(:,:,:,n), 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    end do
    ! Apply open boundaries for domain boundaries for full levels (ghost cells)
    do i = 1,5 ! Loop over boundaries
      if(.not.lboundary(i).or.lperiodic(i)) cycle
      call applyboundaryf(thlm ,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%thl,boundary(i)%nx1,boundary(i)%nx2,0,boundary(i)%thlturb,profile=thl0av)
      call applyboundaryf(thl0 ,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%thl,boundary(i)%nx1,boundary(i)%nx2,0,boundary(i)%thlturb,profile=thl0av)
      call applyboundaryf(qtm  ,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%qt,boundary(i)%nx1,boundary(i)%nx2,1,boundary(i)%qtturb,profile=qt0av)
      call applyboundaryf(qt0  ,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%qt,boundary(i)%nx1,boundary(i)%nx2,1,boundary(i)%qtturb,profile=qt0av)
      call applyboundaryf(e12m ,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%e12,boundary(i)%nx1,boundary(i)%nx2,1,boundary(i)%e12turb)
      call applyboundaryf(e120 ,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%e12,boundary(i)%nx1,boundary(i)%nx2,1,boundary(i)%e12turb)
      do n = 1,nsv
        call applyboundaryf(svm(:,:,:,n) ,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%sv(:,:,:,n),boundary(i)%nx1,boundary(i)%nx2,1,boundary(i)%svturb(:,:,n))
        call applyboundaryf(sv0(:,:,:,n) ,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%sv(:,:,:,n),boundary(i)%nx1,boundary(i)%nx2,1,boundary(i)%svturb(:,:,n))
      end do
      if(i/=1.and.i/=2) call applyboundaryf(um,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%u,boundary(i)%nx1u,boundary(i)%nx2u,0,boundary(i)%uturb,profile=u0av)
      if(i/=1.and.i/=2) call applyboundaryf(u0,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%u,boundary(i)%nx1u,boundary(i)%nx2u,0,boundary(i)%uturb,profile=u0av)
      if(i/=3.and.i/=4) call applyboundaryf(vm,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%v,boundary(i)%nx1v,boundary(i)%nx2v,0,boundary(i)%vturb,profile=v0av)
      if(i/=3.and.i/=4) call applyboundaryf(v0,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%v,boundary(i)%nx1v,boundary(i)%nx2v,0,boundary(i)%vturb,profile=v0av)
      if(i/=5) call applyboundaryf(wm,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%w,boundary(i)%nx1w,boundary(i)%nx2w,0,boundary(i)%wturb)
      if(i/=5) call applyboundaryf(w0,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%w,boundary(i)%nx1w,boundary(i)%nx2w,0,boundary(i)%wturb)
    end do
  end subroutine openboundary_ghost

  subroutine openboundary_tend
    ! Subroutine that handles the tendencies of the boundary normal velocity components
    ! Radiation boundary conditions are used for outflow boundaries and nudging
    ! boundary conditions for the inflow boundaries.
    ! Outflow:
    ! du/dt = -uphase*du/dx
    ! Inflow:
    ! du/dt = (u-uboundary)/tau
    implicit none
    integer :: ib

    if(.not.lopenbc) return
    if(lboundary(1).and..not.lperiodic(1)) call applyboundaryh(1,boundary(1)%nx1,boundary(1)%nx2,boundary(1)%uturb)
    if(lboundary(2).and..not.lperiodic(2)) call applyboundaryh(2,boundary(2)%nx1,boundary(2)%nx2,boundary(2)%uturb)
    if(lboundary(3).and..not.lperiodic(3)) call applyboundaryh(3,boundary(3)%nx1,boundary(3)%nx2,boundary(3)%vturb)
    if(lboundary(4).and..not.lperiodic(4)) call applyboundaryh(4,boundary(4)%nx1,boundary(4)%nx2,boundary(4)%vturb)
    if(lboundary(5).and..not.lperiodic(5)) call applyboundaryh(5,boundary(5)%nx1,boundary(5)%nx2,boundary(5)%wturb)
    ! Calculate and add correction term to guarantee conservation of mass
    do ib = 1,5
      if(.not. lboundary(ib).or.lperiodic(ib)) cycle
      call radcorrection(ib)
    end do
  end subroutine openboundary_tend

  subroutine openboundary_phasevelocity
    ! Subroutine that calculates the phase velocity that is required for the
    ! radiation outflow boundary. The phase velocity is calculated from the phase
    ! velocity of one gridcell to the interior at the prior time and is averaged
    ! over the integration length scales. Lower limit is given by input boundary
    ! velocity and upper limit by cfd criterium.
    use modmpi, only : comm3d,commrow,commcol,myidx,myidy,mpierr,D_MPI_ALLREDUCE,MPI_SUM
    use modglobal, only : imax,jmax,kmax,i1,j1,dx,dy,dzf
    use modfields, only : u0,up,v0,vp,w0,wp,rhobh,rhobf
    implicit none
    integer :: ib,i,j,k,ipatch,jpatch,kpatch
    real :: ipos,jpos

    if(.not.lopenbc) return
    do ib = 1,5 ! Loop over boundaries
      if(.not.lboundary(ib).or.lperiodic(ib)) cycle
      select case(ib) ! Select boundary
      case(1) ! West
        boundary(1)%uphasesingle=0.
        do j = 1,jmax
          jpos = j + (myidy * jmax) - 1
          jpatch = int((jpos-0.5)*dy/dyint)+1
          do k = 1,kmax
            kpatch = k
            boundary(1)%uphasesingle(jpatch,kpatch) = boundary(1)%uphasesingle(jpatch,kpatch) + &
              (-up(3,j+1,k)*dx/sign(max(abs(u0(4,j+1,k)-u0(3,j+1,k)),real(1e-10, field_r)),u0(4,j+1,k)-u0(3,j+1,k))) &
              *dy*rhobf(k)*dzf(k)/dyint*rhointi(kpatch)
          end do
        end do
        ! Integrate over processes
        call D_MPI_ALLREDUCE(boundary(1)%uphasesingle,boundary(1)%uphase,nypatch*nzpatch, &
                           MPI_SUM, commcol,mpierr)
      case(2) ! East
        boundary(2)%uphasesingle=0.
        do j = 1,jmax
          jpos = j + (myidy * jmax) - 1
          jpatch = int((jpos-0.5)*dy/dyint)+1
          do k = 1,kmax
            kpatch = k
            boundary(2)%uphasesingle(jpatch,kpatch) = boundary(2)%uphasesingle(jpatch,kpatch) + &
              (-up(i1,j+1,k)*dx/sign(max(abs(u0(i1,j+1,k)-u0(i1-1,j+1,k)),real(1e-10, field_r)),u0(i1,j+1,k)-u0(i1-1,j+1,k))) &
              *dy*rhobf(k)*dzf(k)/dyint*rhointi(kpatch)
          end do
        end do
        ! Integrate over processes
        call D_MPI_ALLREDUCE(boundary(2)%uphasesingle,boundary(2)%uphase,nypatch*nzpatch, &
                           MPI_SUM, commcol,mpierr)
      case(3) ! South
        boundary(3)%uphasesingle=0.
        do i = 1,imax
          ipos = i + (myidx * imax) - 1
          ipatch = int((ipos-0.5)*dx/dxint)+1
          do k = 1,kmax
            kpatch = k
            boundary(3)%uphasesingle(ipatch,kpatch) = boundary(3)%uphasesingle(ipatch,kpatch) + &
              (-vp(i+1,3,k)*dy/sign(max(abs(v0(i+1,4,k)-v0(i+1,3,k)),real(1e-10, field_r)),v0(i+1,4,k)-v0(i+1,3,k))) &
              *dx*rhobf(k)*dzf(k)/dxint*rhointi(kpatch)
          end do
        end do
        ! Integrate over processes
        call D_MPI_ALLREDUCE(boundary(3)%uphasesingle,boundary(3)%uphase,nxpatch*nzpatch, &
                           MPI_SUM, commrow,mpierr)
      case(4) ! North
        boundary(4)%uphasesingle=0.
        do i = 1,imax
          ipos = i + (myidx * imax) - 1
          ipatch = int((ipos-0.5)*dx/dxint)+1
          do k = 1,kmax
            kpatch = k
            boundary(4)%uphasesingle(ipatch,kpatch) = boundary(4)%uphasesingle(ipatch,kpatch) + &
              (-vp(i+1,j1,k)*dy/sign(max(abs(v0(i+1,j1,k)-v0(i+1,j1-1,k)),real(1e-10, field_r)),v0(i+1,j1,k)-v0(i+1,j1-1,k))) &
              *dx*rhobf(k)/dxint*rhointi(kpatch)
          end do
        end do
        ! Integrate over processes
        call D_MPI_ALLREDUCE(boundary(4)%uphasesingle,boundary(4)%uphase,nxpatch*nzpatch, &
                           MPI_SUM, commrow,mpierr)
      case(5) ! Top
        boundary(5)%uphasesingle=0.
        do i = 1,imax
          ipos = i + (myidx * imax) - 1
          ipatch = int((ipos-0.5)*dx/dxint)+1
          do j = 1,jmax
            jpos = j + (myidy * jmax) - 1
            jpatch = int((jpos-0.5)*dy/dyint)+1
            boundary(5)%uphasesingle(ipatch,jpatch) = boundary(5)%uphasesingle(ipatch,jpatch) + &
              (-rhobh(kmax)*wp(i+1,j+1,kmax)*dzf(kmax-1)/sign(max(abs(rhobh(kmax)*w0(i+1,j+1,kmax)-rhobh(kmax-1)*w0(i+1,j+1,kmax-1)),real(1e-10, field_r)), &
              & rhobh(kmax)*w0(i+1,j+1,kmax)-rhobh(kmax-1)*w0(i+1,j+1,kmax-1)))*dx*dy/(dxint*dyint)
          end do
        end do
        ! Integrate over processes
        call D_MPI_ALLREDUCE(boundary(5)%uphasesingle,boundary(5)%uphase,nxpatch*nypatch, &
                           MPI_SUM, comm3d,mpierr)
      end select
    end do
  end subroutine openboundary_phasevelocity

  subroutine openboundary_turb
    ! Subroutine that calls the synthetic turbulence routine for the generation
    ! of synthetic turbulence at the dirichlet inflow boundaries.
    use modglobal, only : rk3step
    implicit none
    if(rk3step == 1) call synturb()
  end subroutine openboundary_turb

  subroutine applyboundaryf(a,sx,ex,sy,ey,sz,ez,ih,jh,ib,val,nx1,nx2,lmax0,turb,profile)
    ! Routine fills ghost cells based on robin (inflow) or
    ! homogeneous neumann (outflow) boundary conditions. Adds turbulent
    ! pertubations to inflow condition if lsynturb=true.
    use modglobal, only : dzh,dx,dy,imax,jmax,kmax,rtimee,rdt,i2,j2,k1,i1,j1,tauh,pbc
    use modfields, only : u0,v0,w0,e120
    use modmpi, only : myid
    implicit none
    integer, intent(in) :: sx,ex,sy,ey,sz,ez,ih,jh,ib,nx1,nx2,lmax0
    real(field_r), intent(in), dimension(nx1,nx2,ntboundary) :: val
    real(field_r), intent(in), dimension(nx1,nx2) :: turb
    real(field_r), intent(in), dimension(k1), optional :: profile ! optional for top boundary to take gradient into account
    real(field_r), intent(inout), dimension(sx-ih:ex+ih,sy-jh:ey+jh,sz:ez) :: a
    integer :: i,j,k,itp,itm,kav=5,itpn,itmn
    real :: coefdir,coefneu,tp,tm,fp,fm,fpn,fmn,ddz,valtarget,un,e

    ! Get interpolation coefficients for boundary input
    itm=1
    if(ntboundary>1) then
      do while(rtimee-rdt>tboundary(itm))
        itm=itm+1
      end do
      if (rtimee-rdt>tboundary(1)) then
        itm=itm-1
      end if
      itp = itm+1
      tm = tboundary(itm)
      tp = tboundary(itp)
      fm = (tp-rtimee+rdt)/(tp-tm)
      fp = (rtimee-rdt-tm)/(tp-tm)
    else
      itp = 1
      fp  = 0.
      fm  = 1.
    endif
    select case(ib) ! Select domain boundary
    case(1) ! West
      do k = 1,nx2
        do j = 1,nx1
          un = u0(sx,min(j+1,j1),min(k,kmax))
          if(un<=0) then ! Homogeneous Neumann outflow
            a(sx-1,j+1,k)=a(sx,j+1,k)
          else ! Robin inflow conditions
            e = e120(sx,min(j+1,j1),min(k,kmax))
            coefdir = abs(un)**pbc
            coefneu = -tauh*un*(abs(un)**pbc+e**pbc)
            valtarget = (fp*val(j,k,itp)+fm*val(j,k,itm)+turb(j,k))*coefdir
            a(sx-1,j+1,k) = ( 2.*dx*valtarget - &
              a(sx,j+1,k)*(coefdir*dx+2.*coefneu) ) / (coefdir*dx-2.*coefneu)
            if(lmax0==1) a(sx-1,j+1,k) = max(0.,a(sx-1,j+1,k))
          endif
        end do
      end do
    case(2) ! East
      do k = 1,nx2
        do j = 1,nx1
          un = u0(ex+1,min(j+1,j1),min(k,kmax))
          if(un>=0) then ! Homogeneous Neumann outflow
            a(ex+1,j+1,k)=a(ex,j+1,k)
          else ! Robin inflow conditions
            e = e120(ex,min(j+1,j1),min(k,kmax))
            coefdir = abs(un)**pbc
            coefneu = -tauh*un*(abs(un)**pbc+e**pbc)
            valtarget = (fp*val(j,k,itp)+fm*val(j,k,itm)+turb(j,k))*coefdir
            a(ex+1,j+1,k) = ( 2.*dx*valtarget - &
              a(ex,j+1,k)*(coefdir*dx-2.*coefneu) ) / (coefdir*dx+2.*coefneu)
            if(lmax0==1) a(ex+1,j+1,k) = max(a(ex+1,j+1,k),0.)
          endif
        end do
      end do
    case(3) ! South
      do k = 1,nx2
        do i = 1,nx1
          un = v0(min(i+1,i1),sy,min(k,kmax))
          if(un<=0) then ! Homogeneous Neumann outflow
            a(i+1,sy-1,k)=a(i+1,sy,k)
          else ! Robin inflow conditions
            e = e120(min(i+1,i1),sy,min(k,kmax))
            coefdir = abs(un)**pbc
            coefneu = -tauh*un*(abs(un)**pbc+e**pbc)
            valtarget = (fp*val(i,k,itp)+fm*val(i,k,itm)+turb(i,k))*coefdir
            a(i+1,sy-1,k) = ( 2.*dy*valtarget - &
              a(i+1,sy,k)*(coefdir*dy+2.*coefneu) ) / (coefdir*dy-2.*coefneu)
            if(lmax0==1) a(i+1,sy-1,k) = max(a(i+1,sy-1,k),0.)
          endif
        end do
      end do
    case(4) ! North
      do k = 1,nx2
        do i = 1,nx1
          un = v0(min(i+1,i1),ey+1,min(k,kmax))
          if(un>=0) then ! Homogeneous Neumann outflow
            a(i+1,ey+1,k)=a(i+1,ey,k)
          else ! Robin inflow conditions
            e = e120(min(i+1,i1),ey,min(k,kmax))
            coefdir = abs(un)**pbc
            coefneu = -tauh*un*(abs(un)**pbc+e**pbc)
            valtarget = (fp*val(i,k,itp)+fm*val(i,k,itm)+turb(i,k))*coefdir
            a(i+1,ey+1,k) = ( 2.*dy*valtarget - &
              a(i+1,ey,k)*(coefdir*dy-2.*coefneu) ) / (coefdir*dy+2.*coefneu)
            if(lmax0==1) a(i+1,ey+1,k) = max(a(i+1,ey+1,k),0.)
          endif
        end do
      end do
    case(5) ! Top
      ! Obtain verticle gradient if slab averaged profile is given
      if(present(profile)) then
        ddz = sum((profile(kmax-kav+1:kmax)-profile(kmax-kav:kmax-1))/ &
                   dzh(kmax-kav+1:kmax))/kav
      else
        ddz = 0.
      endif
      do i = 1,nx1
        do j = 1,nx2
          un = w0(min(i+1,i1),min(j+1,j1),ez)
          if(un>=0) then ! Neumann outflow
            a(i+1,j+1,ez)=ddz*dzh(ez)+a(i+1,j+1,ez-1)
          else ! Robin inflow conditions
            e = e120(min(i+1,i1),min(j+1,j1),ez-1)
            coefdir = abs(un)**pbc
            coefneu = -tauh*un*(abs(un)**pbc+e**pbc)
            valtarget = (fp*val(i,j,itp)+fm*val(i,j,itm)+turb(i,j))*coefdir+ddz*coefneu
            a(i+1,j+1,ez) = ( 2.*dzh(ez)*valtarget - &
              a(i+1,j+1,ez-1)*(coefdir*dzh(ez)-2.*coefneu) ) / (coefdir*dzh(ez)+2.*coefneu)
            if(lmax0==1) a(i+1,j+1,ez) = max(a(i+1,j+1,ez),0.)
          endif
        end do
      end do
    end select
  end subroutine applyboundaryf

  subroutine applyboundaryh(ib,nx1,nx2,turb)
    ! Subroutine that applies the radiation and dirichlet boundary conditions
    ! for the boundary-normal velocity components. Adds turbulence to
    ! the inflow dirichlet boundaries if lsynturb=.true.
    use modmpi, only : myidx,myidy,myid
    use modglobal, only : dx,dy,dzf,dxi,dyi,rdt,i2,j2,k1,i1,j1,kmax,rtimee,rdt,itot,jtot,imax,jmax,grav,taum
    use modfields, only : um,u0,up,vm,v0,vp,wm,w0,wp,rhobf,rhobh,thvh,thv0h
    implicit none
    integer, intent(in) :: nx1,nx2,ib
    real(field_r), intent(in), dimension(nx1,nx2) :: turb
    integer :: i,j,k,itmc,itmn,itpc,itpn,ipatch,jpatch,kpatch
    real :: tm,tp,fpc,fmc,fpn,fmn,unext,uwallcurrent,ipos,jpos,tau
    itmc=1
    itmn=1
    if(ntboundary>1) then
      do while(rtimee-rdt>tboundary(itmc))
        itmc=itmc+1
      end do
      if (rtimee-rdt>tboundary(1)) then
        itmc=itmc-1
      end if
      do while(tboundary(itmn)<rtimee)
        itmn=itmn+1
      end do
      if (rtimee>tboundary(1)) then
        itmn=itmn-1
      end if
      itpc = itmc+1
      itpn = itmn+1
      tm = tboundary(itmc)
      tp = tboundary(itpc)
      fmc = (tp-rtimee+rdt)/(tp-tm)
      fpc = (rtimee-rdt-tm)/(tp-tm)
      tm = tboundary(itmn)
      tp = tboundary(itpn)
      fmn = (tp-rtimee)/(tp-tm)
      fpn = (rtimee-tm)/(tp-tm)
    else
      itpc = 1
      itpn = 1
      fpc  = 0.
      fmc  = 1.
      fpn  = 0.
      fmn  = 1.
    endif
    ! Apply domain boundaries
    select case(ib) ! Select boundary
    case(1) ! West
      tau = max(taum,rdt)
      do j = 1,nx1
        jpos = j + (myidy * jmax) - 1
        jpatch = int((jpos-0.5)*dy/dyint)+1
        do k = 1,nx2
          kpatch = k
          uwallcurrent = fpc*boundary(1)%u(j,k,itpc)+fmc*boundary(1)%u(j,k,itmc)
          if(uwallcurrent<=0.) then ! Outflow (Radiation)
            up(2,j+1,k) = -max(min(boundary(1)%uphase(jpatch,kpatch),uwallcurrent),-dx/rdt) * &
              (u0(3,j+1,k)-u0(2,j+1,k))*dxi
          else ! Inflow nudging
            unext = fpn*boundary(1)%u(j,k,itpn)+fmn*boundary(1)%u(j,k,itmn)
            up(2,j+1,k) = ((unext+turb(j,k)) - u0(2,j+1,k))/tau
          endif
        end do
      end do
    case(2) ! East
      tau = max(taum,rdt)
      do j = 1,nx1
        jpos = j + (myidy * jmax) - 1
        jpatch = int((jpos-0.5)*dy/dyint)+1
        do k = 1,nx2
          kpatch = k
          uwallcurrent = fpc*boundary(2)%u(j,k,itpc)+fmc*boundary(2)%u(j,k,itmc)
          if(uwallcurrent>=0.) then ! Outflow (Radiation)
            up(i2,j+1,k) = -min(max(boundary(2)%uphase(jpatch,kpatch),uwallcurrent),dx/rdt) * &
              (u0(i2,j+1,k)-u0(i1,j+1,k))*dxi
          else ! Inflow (Dirichlet)
            unext = fpn*boundary(2)%u(j,k,itpn)+fmn*boundary(2)%u(j,k,itmn)
            up(i2,j+1,k) = ((unext+turb(j,k)) - u0(i2,j+1,k))/tau
          endif
        end do
      end do
    case(3) ! South
      tau = max(taum,rdt)
      do i = 1,nx1
        ipos = i + (myidx * imax) - 1
        ipatch = int((ipos-0.5)*dx/dxint)+1
        do k = 1,nx2
          kpatch = k
          uwallcurrent = fpc*boundary(3)%v(i,k,itpc)+fmc*boundary(3)%v(i,k,itmc)
          if(uwallcurrent<=0.) then ! Outflow (Radiation)
            vp(i+1,2,k) = -max(min(boundary(3)%uphase(ipatch,kpatch),uwallcurrent),-dy/rdt) * &
              (v0(i+1,3,k)-v0(i+1,2,k))*dyi
          else ! Inflow (Dirichlet)
            unext = fpn*boundary(3)%v(i,k,itpn)+fmn*boundary(3)%v(i,k,itmn)
            vp(i+1,2,k) = ((unext+turb(i,k)) - v0(i+1,2,k))/tau
          endif
        end do
      end do
    case(4) ! North
      tau = max(taum,rdt)
      do i = 1,nx1
        ipos = i + (myidx * imax) - 1
        ipatch = int((ipos-0.5)*dx/dxint)+1
        do k = 1,nx2
          kpatch = k
          uwallcurrent = fpc*boundary(4)%v(i,k,itpc)+fmc*boundary(4)%v(i,k,itmc)
          if(uwallcurrent>=0.) then ! Outflow (Radiation)
            vp(i+1,j2,k) = -min(max(boundary(4)%uphase(ipatch,kpatch),uwallcurrent),dy/rdt) * &
              (v0(i+1,j2,k)-v0(i+1,j1,k))*dyi
          else ! Inflow (Dirichlet)
            unext = fpn*boundary(4)%v(i,k,itpn)+fmn*boundary(4)%v(i,k,itmn)
            vp(i+1,j2,k) = ((unext+turb(i,k)) - v0(i+1,j2,k))/tau
          endif
        end do
      end do
    case(5) ! Top
      tau = max(taum,rdt)
      do i = 1,nx1
        ipos = i + (myidx * imax) - 1
        ipatch = int((ipos-0.5)*dx/dxint)+1
        do j = 1,nx2
          jpos = j + (myidy * jmax) - 1
          jpatch = int((jpos-0.5)*dy/dyint)+1
          uwallcurrent = fpc*boundary(5)%w(i,j,itpc)+fmc*boundary(5)%w(i,j,itmc)
          if(uwallcurrent>=0.) then ! Outflow (Radiation)
            wp(i+1,j+1,k1) = -min(max(boundary(5)%uphase(ipatch,jpatch),uwallcurrent),dzf(kmax)/rdt) * &
              (rhobh(k1)*w0(i+1,j+1,k1)-rhobh(kmax)*w0(i+1,j+1,kmax))/(dzf(kmax)*rhobh(k1)) + &
              ibuoy*(grav*(thv0h(i+1,j+1,k1)-thvh(k1))/thvh(k1))
          else ! Inflow (Dirichlet)
            unext = fpn*boundary(5)%w(i,j,itpn)+fmn*boundary(5)%w(i,j,itmn)
            wp(i+1,j+1,k1) = ((unext+turb(i,j)) - w0(i+1,j+1,k1))/tau
          endif
        end do
      end do
    end select
  end subroutine applyboundaryh

  subroutine radcorrection(ib)
    ! Calculates the integrated mass correction term for the boundary normal
    ! velocity components
    use modmpi, only : comm3d,commrow,commcol,myidx,myidy,mpierr, D_MPI_ALLREDUCE, MPI_SUM
    use modglobal, only : jmax,imax,kmax,i1,j1,dx,dy,dzf,i2,j2,k1,dxi,dyi,rtimee,rdt
    use modfields, only : rhobf, up, vp, wp
    implicit none
    integer, intent(in) :: ib
    integer :: ipos,jpos,kpos,ipatch,jpatch,kpatch,i,j,k,itp,itm
    real :: sum,tp,tm,idtb,dubdt

    itm = 1
    if(ntboundary>1) then
      do while(rtimee-rdt>tboundary(itm))
        itm=itm+1
      end do
      if (rtimee-rdt>tboundary(1)) then
        itm=itm-1
      end if
      itp = itm+1
    else
      itp = 1
    endif
    tm = tboundary(itm)
    tp = tboundary(itp)
    idtb = 1./max(1e-6,tp-tm)
    select case(ib) ! Select boundary
    case(1) ! West
      ! Calculate correction term for each patch
      boundary(1)%radcorrsingle = 0.
      do j = 2,j1
        jpos = j + (myidy * jmax) - 1
        jpatch = int((jpos-0.5)*dy/dyint)+1
        do k = 1,kmax
          kpatch = k
          dubdt = (boundary(1)%u(j-1,k,itp)-boundary(1)%u(j-1,k,itm))*idtb
          boundary(1)%radcorrsingle(jpatch,kpatch) = boundary(1)%radcorrsingle(jpatch,kpatch) + &
            rhobf(k)*(-up(2,j,k)+dubdt)*dzf(k)*dy*rhointi(kpatch)/dyint
        end do
      end do
      ! Communicate integration between processes
      call D_MPI_ALLREDUCE(boundary(1)%radcorrsingle,boundary(1)%radcorr,nypatch*nzpatch, &
                         MPI_SUM, commcol,mpierr)
      ! Apply correction term
      do j = 2,j1
        jpos = j + (myidy * jmax) - 1
        jpatch = int((jpos-0.5)*dy/dyint)+1
        do k = 1,kmax
          kpatch = k
          up(2,j,k) = up(2,j,k) + boundary(1)%radcorr(jpatch,kpatch)
        end do
      end do
    case(2) ! East
      ! Calculate correction term for each patch
      boundary(2)%radcorrsingle = 0.
      do j = 2,j1
        jpos = j + (myidy * jmax) - 1
        jpatch = int((jpos-0.5)*dy/dyint)+1
        do k = 1,kmax
          kpatch = k
          dubdt = (boundary(2)%u(j-1,k,itp)-boundary(2)%u(j-1,k,itm))*idtb
          boundary(2)%radcorrsingle(jpatch,kpatch) = boundary(2)%radcorrsingle(jpatch,kpatch) + &
            rhobf(k)*(-up(i2,j,k)+dubdt)*dzf(k)*dy*rhointi(kpatch)/dyint
        end do
      end do
      ! Communicate integration between processes
      call D_MPI_ALLREDUCE(boundary(2)%radcorrsingle,boundary(2)%radcorr,nypatch*nzpatch, &
                         MPI_SUM, commcol,mpierr)
      ! Apply correction term
      do j = 2,j1
        jpos = j + (myidy * jmax) - 1
        jpatch = int((jpos-0.5)*dy/dyint)+1
        do k = 1,kmax
          kpatch = k
          up(i2,j,k) = up(i2,j,k) + boundary(2)%radcorr(jpatch,kpatch)
        end do
      end do
    case(3) ! South
      ! Calculate correction term for each patch
      boundary(3)%radcorrsingle= 0.
      do i = 2,i1
        ipos = i + (myidx * imax) - 1
        ipatch = int((ipos-0.5)*dx/dxint)+1
        do k = 1,kmax
          kpatch = k
          dubdt = (boundary(3)%v(i-1,k,itp)-boundary(3)%v(i-1,k,itm))*idtb
          boundary(3)%radcorrsingle(ipatch,kpatch) = boundary(3)%radcorrsingle(ipatch,kpatch) + &
            rhobf(k)*(-vp(i,2,k)+dubdt)*dzf(k)*dx*rhointi(kpatch)/dxint
        end do
      end do
      ! Communicate integration between processes
      call D_MPI_ALLREDUCE(boundary(3)%radcorrsingle,boundary(3)%radcorr,nxpatch*nzpatch, &
                         MPI_SUM, commrow,mpierr)
      ! Apply correction term
      do i = 2,i1
        ipos = i + (myidx * imax) - 1
        ipatch = int((ipos-0.5)*dx/dxint)+1
        do k = 1,kmax
          kpatch = k
          vp(i,2,k) = vp(i,2,k) + boundary(3)%radcorr(ipatch,kpatch)
        end do
      end do
    case(4) ! North
      ! Calculate correction term for each patch
      boundary(4)%radcorrsingle = 0.
      do i = 2,i1
        ipos = i + (myidx * imax) - 1
        ipatch = int((ipos-0.5)*dx/dxint)+1
        do k = 1,kmax
          kpatch = k
          dubdt = (boundary(4)%v(i-1,k,itp)-boundary(4)%v(i-1,k,itm))*idtb
          boundary(4)%radcorrsingle(ipatch,kpatch) = boundary(4)%radcorrsingle(ipatch,kpatch) + &
            rhobf(k)*(-vp(i,j2,k)+dubdt)*dzf(k)*dx*rhointi(kpatch)/dxint
        end do
      end do
      ! Communicate integration between processes
      call D_MPI_ALLREDUCE(boundary(4)%radcorrsingle,boundary(4)%radcorr,nxpatch*nzpatch, &
                         MPI_SUM, commrow,mpierr)
      ! Apply correction term
      do i = 2,i1
        ipos = i + (myidx * imax) - 1
        ipatch = int((ipos-0.5)*dx/dxint)+1
        do k = 1,kmax
          kpatch = k
          vp(i,j2,k) = vp(i,j2,k) + boundary(4)%radcorr(ipatch,kpatch)
        end do
      end do
    case(5) ! Top
      ! Calculate correction term for each patch
      boundary(5)%radcorrsingle = 0.
      do i = 2,i1
        ipos = i + (myidx * imax) - 1
        ipatch = int((ipos-0.5)*dx/dxint)+1
        do j = 2,j1
          jpos = j + (myidy * jmax) - 1
          jpatch = int((jpos-0.5)*dy/dyint)+1
          dubdt = (boundary(5)%w(i-1,j-1,itp)-boundary(5)%w(i-1,j-1,itm))*idtb
          boundary(5)%radcorrsingle(ipatch,jpatch) = boundary(5)%radcorrsingle(ipatch,jpatch) + &
            (-wp(i,j,k1)+dubdt)*dx*dy/(dxint*dyint)
        end do
      end do
      ! Communicate integration between processes
      call D_MPI_ALLREDUCE(boundary(5)%radcorrsingle,boundary(5)%radcorr,nxpatch*nypatch, &
                         MPI_SUM, comm3d,mpierr)
      ! Apply correction term
      do i = 2,i1
        ipos = i + (myidx * imax) - 1
        ipatch = int((ipos-0.5)*dx/dxint)+1
        do j = 2,j1
          jpos = j + (myidy * jmax) - 1
          jpatch = int((jpos-0.5)*dy/dyint)+1
          wp(i,j,k1) = wp(i,j,k1) + boundary(5)%radcorr(ipatch,jpatch)
        end do
      end do
    end select
  end subroutine radcorrection

  subroutine take_prof(field0,fieldm,prof)
    use modglobal, only : i1,j1,k1,ih,jh,kmax
    implicit none
    real(field_r), intent(inout), dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: field0,fieldm
    real(field_r), intent(in), dimension(k1) :: prof
    integer :: i,j,k
    do k=1,kmax
      do j=2,j1
        do i=2,i1
          field0(i,j,k) = prof(k)
          fieldm(i,j,k) = prof(k)
        end do
      end do
    end do
  end subroutine take_prof

  subroutine handle_err(errcode)

  implicit none

  integer errcode

  write(6,*) 'Error: ', nf90_strerror(errcode)
  stop 2

  end subroutine handle_err

end module modopenboundary
