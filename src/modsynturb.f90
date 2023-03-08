!> \file modopenboundary.f90
!!  Creates synthetic turbulence for the open boundary implementation
!>
!!  Creates synthetic turbulence for the open boundary implementation
!>
!!  \author Frans Liqui Lung
!  This file is part of DALES.
!  To do:
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
module modsynturb
use netcdf
use modglobal,only: lsynturb,iturb,lboundary,lperiodic,boundary,nmodes,lambda,tau,dxturb,dyturb,itot,jtot,dx,dy,kmax
use RandomNumbers, only : getRandomReal,randomNumberSequence,new_RandomNumberSequence
implicit none
real, allocatable, dimension(:,:) :: kn,p,q,vturb,wturb,k_thl,k_qt
real, allocatable, dimension(:) :: omega,omega_thl,omega_qt,p_thl,p_qt,q_thl,q_qt
real, allocatable, dimension(:) :: xf,xh,yf,yh
real :: nisqrt,ctot,nisqrt2
real, dimension(3) :: lambdasxyz
integer :: nxturb,nyturb,nzturb
integer, parameter :: isepsim_mom = 10,isepsim_all=11, isynturb_mom = 0, isynturb_all = 1
integer :: ntturb,itimestep=1
real, allocatable, dimension(:) :: tturb
real, allocatable, dimension(:,:,:) :: uturbin,vturbin,wturbin,thlturbin,qtturbin
character(len = nf90_max_name) :: RecordDimName
integer :: VARID,STATUS,NCID,mpierr,timeID
! ! Uncommend for netcdf output turbulent pertubations west boundary
! character (80) :: fname = 'turbOut.xxx.xxx.nc'
! integer :: ncid
! integer, parameter :: NDIMS = 3
! character (len = *), parameter :: Z_NAME = "z"
! character (len = *), parameter :: Y_NAME = "y"
! character (len = *), parameter :: t_NAME = "time"
! integer :: z_dimid, y_dimid, t_dimid
! integer :: z_varid, y_varid, t_varid
! integer :: uturb_varid, vturb_varid, wturb_varid, thlturb_varid, qtturb_varid
! integer :: dimids(NDIMS)
! integer :: start(NDIMS), count(NDIMS)

type(randomNumberSequence) :: noise
contains
  subroutine initsynturb
    use netcdf
    use modglobal, only : dx,dy,imax,jmax,i1,j1,zf,lambdas,lambdas_x,lambdas_y,lambdas_z,kmax,k1,cexpnr,lmoist
    use modmpi, only : myidx, myidy
    implicit none
    integer :: i,j,ib
    if(.not.lsynturb) return
    if(any(lboundary.and..not.(lperiodic))) then
      if(iturb == isynturb_all .or. iturb == isynturb_mom) then
        ! Constants
        nisqrt  = sqrt(2./nmodes)
        nisqrt2 = sqrt(1./nmodes)
        ctot = lambda/tau
        noise = new_RandomNumberSequence(seed = 100)
        nxturb = int(dx/dxturb*real(itot));
        nyturb = int(dy/dyturb*real(jtot));
        nzturb = kmax
        lambdas = merge(lambda,lambdas,lambdas==-1.)
        lambdasxyz = (/merge(lambdas,lambdas_x,lambdas_x==-1.), &
                    & merge(lambdas,lambdas_y,lambdas_y==-1.), &
                    & merge(lambdas,lambdas_z,lambdas_z==-1.)/)
        ! Allocate variables
        allocate(kn(nmodes,3),q(nmodes,3),p(nmodes,3),omega(nmodes),&
          xf(imax),xh(i1),yf(jmax),yh(j1),k_thl(nmodes,3),k_qt(nmodes,3), &
          omega_thl(nmodes),omega_qt(nmodes), &
          p_thl(nmodes),p_qt(nmodes),q_thl(nmodes),q_qt(nmodes))
        ! allocate(vturb(jmax,kmax),wturb(jmax,kmax))
        ! Calculate coordinates
        xf = (/((i-0.5)*dx,i=1,imax,1)/)+imax*myidx*dx
        xh = (/(i*dx,i=1,i1,1)/)+imax*myidx*dx
        yf = (/((i-0.5)*dy,i=1,jmax,1)/)+jmax*myidy*dy
        yh = (/(i*dy,i=1,j1,1)/)+jmax*myidy*dy
        ! Fill random distributed variables
        do i = 1,3
          do j = 1,nmodes
            kn(j,i) = gaussrand(0.,0.5)
            k_thl(j,i) = gaussrand(0.,0.5)
            k_qt(j,i) = gaussrand(0.,0.5)
            if(i==1) then
              omega(j) = gaussrand(0.,1.)
              omega_thl(j) = gaussrand(0.,1.)
              omega_qt(j) = gaussrand(0.,1.)
              p_thl(j) = gaussrand(0.,1.)
              q_thl(j) = gaussrand(0.,1.)
              p_qt(j)  = gaussrand(0.,1.)
              q_qt(j)  = gaussrand(0.,1.)
            endif
          end do
        end do
        ! Obtain p and q with cross product
        do i = 1,nmodes
          p(i,:) = cross((/gaussrand(0.,1.),gaussrand(0.,1.),gaussrand(0.,1.)/),kn(i,:))
          q(i,:) = cross((/gaussrand(0.,1.),gaussrand(0.,1.),gaussrand(0.,1.)/),kn(i,:))
        end do
        do ib = 1,5
          if(.not.lboundary(ib)) cycle
          allocate(boundary(ib)%eigvec(boundary(ib)%nx1patch,boundary(ib)%nx2patch,3,3), &
                 & boundary(ib)%ci(boundary(ib)%nx1patch,boundary(ib)%nx2patch,3))!, &
                 ! & boundary(ib)%randthl(boundary(ib)%nx1,boundary(ib)%nx2), &
                 ! & boundary(ib)%randqt(boundary(ib)%nx1,boundary(ib)%nx2))
        end do
      elseif((iturb == isepsim_mom .or. iturb == isepsim_all) .and. lboundary(1)) then
        !--- open nc file ---
        STATUS = NF90_OPEN('turb.inp.'//cexpnr//'.nc', nf90_nowrite, NCID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        !--- get time dimensions
        status = nf90_inq_dimid(ncid, "time", timeID)
        if (status /= nf90_noerr) call handle_err(status)
        status = nf90_inquire_dimension(NCID, timeID, len=ntturb, name=RecordDimName)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        !--- read time
        allocate(tturb(ntturb))
        STATUS = NF90_INQ_VARID(NCID, 'time', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tturb, start=(/1/), count=(/ntturb/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! --- allocate fields
        allocate(uturbin(jmax,kmax,ntturb))
        allocate(vturbin(j1,kmax,ntturb))
        allocate(wturbin(jmax,k1,ntturb))
        ! --- read fields
        ! Read u
        STATUS = NF90_INQ_VARID(NCID,'uturbwest', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, uturbin, start=(/myidy*jmax+1,1,1/), &
          & count=(/jmax,kmax,ntturb/))
        ! Read v
        STATUS = NF90_INQ_VARID(NCID,'vturbwest', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, vturbin, start=(/myidy*jmax+1,1,1/), &
          & count=(/j1,kmax,ntturb/))
        ! Read w
        STATUS = NF90_INQ_VARID(NCID,'wturbwest', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, wturbin, start=(/myidy*jmax+1,1,1/), &
          & count=(/jmax,k1,ntturb/))
	if(iturb==isepsim_all) then
	  allocate(thlturbin(jmax,kmax,ntturb))
          ! Read thl
          STATUS = NF90_INQ_VARID(NCID,'thlturbwest', VARID)
          if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          STATUS = NF90_GET_VAR (NCID, VARID, thlturbin, start=(/myidy*jmax+1,1,1/), &
            & count=(/jmax,kmax,ntturb/))
	  if(lmoist) then
	    allocate(qtturbin(jmax,kmax,ntturb))
            ! qt
            STATUS = NF90_INQ_VARID(NCID,'wturbwest', VARID)
            if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
            STATUS = NF90_GET_VAR (NCID, VARID, wturbin, start=(/myidy*jmax+1,1,1/), &
              & count=(/jmax,k1,ntturb/))
          endif
        endif
      endif
      ! ! Uncommend for netcdf output turbulent pertubations west boundary
      ! if(lboundary(1)) then
      !   write(fname,'(A,i3.3,A,i3.3,A)') 'turbOut.', myidx, '.', myidy, '.nc'
      !   call check( nf90_create(fname, nf90_clobber, ncid) )
      !   call check( nf90_def_dim(ncid, z_NAME, kmax, z_dimid) )
      !   call check( nf90_def_dim(ncid, y_NAME, jmax, y_dimid) )
      !   call check( nf90_def_dim(ncid, t_NAME, NF90_UNLIMITED, t_dimid) )
      !   call check( nf90_def_var(ncid, z_NAME, NF90_REAL, z_dimid, z_varid) )
      !   call check( nf90_def_var(ncid, y_NAME, NF90_REAL, y_dimid, y_varid) )
      !   call check( nf90_def_var(ncid, t_NAME, NF90_REAL, t_dimid, t_varid) )
      !   dimids = (/y_dimid, z_dimid, t_dimid/)
      !   call check( nf90_def_var(ncid, 'u', NF90_REAL, dimids, uturb_varid) )
      !   call check( nf90_def_var(ncid, 'v', NF90_REAL, dimids, vturb_varid) )
      !   call check( nf90_def_var(ncid, 'w', NF90_REAL, dimids, wturb_varid) )
      !   call check( nf90_def_var(ncid, 'thl', NF90_REAL, dimids, thlturb_varid) )
      !   call check( nf90_def_var(ncid, 'qt', NF90_REAL, dimids, qtturb_varid) )
      !   call check( nf90_enddef(ncid) )
      !   call check( nf90_put_var(ncid, z_varid, zf(1:kmax)) )
      !   call check( nf90_put_var(ncid, y_varid, yf) )
      !   count = (/ jmax, kmax, 1 /)
      !   start = (/ 1, 1, 1 /)
      ! endif
    endif
  end subroutine initsynturb

  subroutine handle_err(errcode)

  implicit none

  integer errcode

  write(6,*) 'Error: ', nf90_strerror(errcode)
  stop 2

  end subroutine handle_err

  subroutine exitsynturb
    use modglobal, only : lmoist
    implicit none
    integer :: ib
    if(.not.lsynturb) return
    ! ! Uncommend for netcdf output turbulent pertubations
    ! if(lboundary(1)) call check( nf90_close(ncid) )
    if(any(lboundary.and..not.(lperiodic))) then
      if(iturb == isynturb_mom .or. iturb == isynturb_all) then
        do ib = 1,5
          if(.not.lboundary(ib)) cycle
          deallocate(boundary(ib)%u2,boundary(ib)%v2,boundary(ib)%w2,&
           & boundary(ib)%uv,boundary(ib)%uw,boundary(ib)%vw,&
           & boundary(ib)%thl2,boundary(ib)%qt2,boundary(ib)%wthl, &
           & boundary(ib)%wqt,boundary(ib)%eigvec,boundary(ib)%ci)
        end do
        !deallocate(vturb,wturb)
        deallocate(kn,q,p,omega,xf,xh,yf,yh,k_thl,k_qt,omega_thl,omega_qt,p_thl,p_qt,q_thl,q_qt)
      elseif(iturb == isepsim_mom .and. lboundary(1)) then
        deallocate(uturbin,vturbin,wturbin)
      elseif(iturb == isepsim_all .and. lboundary(1)) then
	deallocate(uturbin,vturbin,wturbin,thlturbin)
	if(lmoist) deallocate(qtturbin)
      endif
    endif
  end subroutine exitsynturb

  ! ! Uncommend for netcdf output turbulent pertubations
  ! subroutine check(status)
  !   integer, intent ( in) :: status
  !
  !   if(status /= nf90_noerr) then
  !     print *, trim(nf90_strerror(status))
  !     stop "Stopped"
  !   end if
  ! end subroutine check

  subroutine synturb()
    implicit none
    if(.not.lsynturb) return
    select case(iturb)
    case(isynturb_mom)
      call synturb_mom()
    case(isynturb_all)
      call synturb_all()
    case(isepsim_mom, isepsim_all)
      call sepsim
    end select
  end subroutine synturb

  subroutine synturb_mom()
    use modglobal, only : dx,dy,itot,jtot,zh,zf,imax,jmax,i1,j1,kmax,k1
    implicit none
    integer :: ib
    do ib = 1,5
      if(.not.lboundary(ib).or.lperiodic(ib)) cycle
      call calc_eigdec(ib)
      select case(ib)
      case(1)
        call calc_pert(ib,(/0./),yf,zf,1,jmax,kmax,boundary(ib)%uturb,1)
        call calc_pert(ib,(/0./),yh,zf,1,j1,kmax,boundary(ib)%vturb,2)
        call calc_pert(ib,(/0./),yf,zh,1,jmax,k1,boundary(ib)%wturb,3)
      case(2)
        call calc_pert(ib,(/(itot+1)*dx/),yf,zf,1,jmax,kmax,boundary(ib)%uturb,1)
        call calc_pert(ib,(/(itot+1)*dx/),yh,zf,1,j1,kmax,boundary(ib)%vturb,2)
        call calc_pert(ib,(/(itot+1)*dx/),yf,zh,1,jmax,k1,boundary(ib)%wturb,3)
      case(3)
        call calc_pert(ib,xh,(/0./),zf,i1,1,kmax,boundary(ib)%uturb,1)
        call calc_pert(ib,xf,(/0./),zf,imax,1,kmax,boundary(ib)%vturb,2)
        call calc_pert(ib,xh,(/0./),zh,imax,1,k1,boundary(ib)%wturb,1)
      case(4)
        call calc_pert(ib,xh,(/(jtot+1)*dy/),zf,i1,1,kmax,boundary(ib)%uturb,1)
        call calc_pert(ib,xf,(/(jtot+1)*dy/),zf,imax,1,kmax,boundary(ib)%vturb,2)
        call calc_pert(ib,xh,(/(jtot+1)*dy/),zh,imax,1,k1,boundary(ib)%wturb,1)
      case(5)
        call calc_pert(ib,xh,yf,(/zh(k1)/),i1,jmax,1,boundary(ib)%uturb,1)
        call calc_pert(ib,xf,yh,(/zh(k1)/),imax,j1,1,boundary(ib)%vturb,2)
        call calc_pert(ib,xf,yf,(/zh(k1)/),imax,jmax,1,boundary(ib)%wturb,3)
      end select
    end do
  end subroutine synturb_mom

  subroutine synturb_all()
    use modglobal, only : dx,dy,itot,jtot,zh,zf,imax,jmax,i1,j1,kmax,k1
    implicit none
    integer :: ib
    do ib = 1,5
      if(.not.lboundary(ib).or.lperiodic(ib)) cycle
      call calc_eigdec(ib)
      select case(ib)
      case(1)
        call calc_pert2(ib,(/0./),yf,zf,1,jmax,kmax,boundary(ib)%uturb,1, &
          & boundary(ib)%thlturb,boundary(ib)%qtturb)
        call calc_pert(ib,(/0./),yh,zf,1,j1,kmax,boundary(ib)%vturb,2)
        call calc_pert(ib,(/0./),yf,zh,1,jmax,k1,boundary(ib)%wturb,3)
      case(2)
        call calc_pert2(ib,(/(itot+1)*dx/),yf,zf,1,jmax,kmax,boundary(ib)%uturb,1, &
          & boundary(ib)%thlturb,boundary(ib)%qtturb)
        call calc_pert(ib,(/(itot+1)*dx/),yh,zf,1,j1,kmax,boundary(ib)%vturb,2)
        call calc_pert(ib,(/(itot+1)*dx/),yf,zh,1,jmax,k1,boundary(ib)%wturb,3)
      case(3)
        call calc_pert(ib,xh,(/0./),zf,i1,1,kmax,boundary(ib)%uturb,1)
        call calc_pert2(ib,xf,(/0./),zf,imax,1,kmax,boundary(ib)%vturb,2, &
          & boundary(ib)%thlturb,boundary(ib)%qtturb)
        call calc_pert(ib,xh,(/0./),zh,imax,1,k1,boundary(ib)%wturb,1)
      case(4)
        call calc_pert(ib,xh,(/(jtot+1)*dy/),zf,i1,1,kmax,boundary(ib)%uturb,1)
        call calc_pert2(ib,xf,(/(jtot+1)*dy/),zf,imax,1,kmax,boundary(ib)%vturb,2, &
          & boundary(ib)%thlturb,boundary(ib)%qtturb)
        call calc_pert(ib,xh,(/(jtot+1)*dy/),zh,imax,1,k1,boundary(ib)%wturb,1)
      case(5)
        call calc_pert(ib,xh,yf,(/zh(k1)/),i1,jmax,1,boundary(ib)%uturb,1)
        call calc_pert(ib,xf,yh,(/zh(k1)/),imax,j1,1,boundary(ib)%vturb,2)
        call calc_pert2(ib,xf,yf,(/zh(k1)/),imax,jmax,1,boundary(ib)%wturb,3, &
          boundary(ib)%thlturb,boundary(ib)%qtturb)
      end select
    end do
  end subroutine synturb_all

  subroutine sepsim()
    use modglobal, only : rtimee,lmoist
    use modmpi, only : myid
    implicit none
    integer :: ib
    do ib = 1,5
      if(.not.lboundary(ib).or.lperiodic(ib)) cycle
      select case(ib)
      case(1)
        do while(tturb(itimestep)<real(rtimee))
          itimestep = itimestep+1
        end do
        if(abs(tturb(itimestep)-real(rtimee))>0.01 .and. real(rtimee)>5.) then
          print *, 'mistake in time', itimestep,tturb(itimestep),real(rtimee)
          stop
        endif
        boundary(ib)%uturb = uturbin(:,:,itimestep)
        boundary(ib)%vturb = vturbin(:,:,itimestep)
        boundary(ib)%wturb = wturbin(:,:,itimestep)
	if(iturb==11) then
	  boundary(ib)%thlturb = thlturbin(:,:,itimestep)
	  if(lmoist) boundary(ib)%qtturb = qtturbin(:,:,itimestep)
	endif
      case(2)
        ! Do nothing (needs to be added)
      case(3)
        ! Do nothing (needs to be added)
      case(4)
        ! Do nothing (needs to be added)
      case(5)
        ! Do nothing (needs to be added)
      end select
    end do
  end subroutine sepsim

  subroutine calc_eigdec(ib)
    use modglobal, only : rtimee,tboundary,ntboundary
    implicit none
    integer, intent(in) :: ib
    real*8,dimension(3,3) :: r,eigvec
    real*8,dimension(3) :: eigval
    integer :: i,j,itp,itm,ii
    real :: fm,fp
    ! Interpolate covariance to current time
    itm = 1
    if(ntboundary>1) then
      do while(rtimee>tboundary(itm))
        itm = itm+1
      end do
      if (rtimee>tboundary(1)) then
        itm = itm-1
      end if
      itp = itm+1
      fm = (tboundary(itp)-rtimee)/(tboundary(itp)-tboundary(itm))
      fp = (rtimee-tboundary(itm))/(tboundary(itp)-tboundary(itm))
    else
      itm = 1; itp = 1; fp  = 1.; fm  = 1.
    endif
    ! Calculate eigenvalues and eigenvectors for each patch
    do i = 1,boundary(ib)%nx1patch
      do j = 1,boundary(ib)%nx2patch
        r(1,1) = fp*boundary(ib)%u2(i,j,itp)+fm*boundary(ib)%u2(i,j,itm)
        r(2,2) = fp*boundary(ib)%v2(i,j,itp)+fm*boundary(ib)%v2(i,j,itm)
        r(3,3) = fp*boundary(ib)%w2(i,j,itp)+fm*boundary(ib)%w2(i,j,itm)
        r(1,2) = fp*boundary(ib)%uv(i,j,itp)+fm*boundary(ib)%uv(i,j,itm)
        r(1,3) = fp*boundary(ib)%uw(i,j,itp)+fm*boundary(ib)%uw(i,j,itp)
        r(2,3) = fp*boundary(ib)%vw(i,j,itp)+fm*boundary(ib)%vw(i,j,itm)
        r(2,1) = r(1,2); r(3,1) = r(1,3); r(3,2) = r(2,3)
        call DSYEVJ3(r,eigvec,eigval)
        do ii = 1,3
          !if(eigval(ii)<-1.e-8) print *,"warning negative eigenvalue ",eigval(ii), " value set to 1e-8"
          eigval(ii) = max(eigval(ii),1.e-8)
        end do
        boundary(ib)%eigvec(i,j,:,:) = real(eigvec)
        boundary(ib)%ci(i,j,:) = real(sqrt(eigval))
      end do
    end do
  end subroutine calc_eigdec

  subroutine calc_pert(ib,x,y,z,nx,ny,nz,uturb,iuturb)
    use modglobal, only : rtimee,dxturb,dyturb
    use modmpi, only : myidx,myidy
    implicit none
    real, dimension(:), intent(in) :: x,y,z
    integer, intent(in) :: nx,ny,nz,iuturb,ib
    real, dimension(:,:), intent(out) :: uturb
    integer,target :: i,j,k,ipatch,jpatch,kpatch
    integer :: ii
    integer, pointer :: pi1, pi2,pi1patch,pi2patch
    real, dimension(3) :: xx,ci
    real, dimension(3,3) :: eigvec
    real :: t,utemp,vtemp,wtemp
    t = rtimee
    select case(ib)
    case(1,2)
      pi1 => j; pi2 => k; pi1patch => jpatch; pi2patch => kpatch
    case(3,4)
      pi1 => i; pi2 => k; pi1patch => ipatch; pi2patch => kpatch
    case(5)
      pi1 => i; pi2 => j; pi1patch => ipatch; pi2patch => jpatch
    end select
    do i = 1,nx
      xx(1) = x(i); ipatch = min(int(x(i)/dxturb)+1,nxturb)
    do j = 1,ny
      xx(2) = y(j); jpatch = min(int(y(j)/dyturb)+1,nyturb)
    do k = 1,nz
      xx(3) = z(k); kpatch = min(k,nzturb)
      ci     = boundary(ib)%ci(pi1patch,pi2patch,:)
      eigvec = boundary(ib)%eigvec(pi1patch,pi2patch,:,:)
      utemp = 0.; vtemp = 0.; wtemp = 0.
      do ii = 1,nmodes
        utemp = utemp + p(ii,1)*cos(dot(kn(ii,:)/ci(1)*ctot,xx/lambda)+omega(ii)*t/tau) + &
                      & q(ii,1)*sin(dot(kn(ii,:)/ci(1)*ctot,xx/lambda)+omega(ii)*t/tau)
        vtemp = vtemp + p(ii,2)*cos(dot(kn(ii,:)/ci(2)*ctot,xx/lambda)+omega(ii)*t/tau) + &
                      & q(ii,2)*sin(dot(kn(ii,:)/ci(2)*ctot,xx/lambda)+omega(ii)*t/tau)
        wtemp = wtemp + p(ii,3)*cos(dot(kn(ii,:)/ci(3)*ctot,xx/lambda)+omega(ii)*t/tau) + &
                      & q(ii,3)*sin(dot(kn(ii,:)/ci(3)*ctot,xx/lambda)+omega(ii)*t/tau)
      end do
      ! Scale velocity fields
      utemp = utemp*ci(1)
      vtemp = vtemp*ci(2)
      wtemp = wtemp*ci(3)
      ! Reproject to cartesian velocity pertubations
      uturb(pi1,pi2)  = nisqrt*dot(eigvec(iuturb,:),(/utemp,vtemp,wtemp/))
    end do
    end do
    end do
    nullify(pi1,pi2,pi1patch,pi2patch)
  end subroutine calc_pert

  subroutine calc_pert2(ib,x,y,z,nx,ny,nz,uturb,iuturb,thlturb,qtturb)
    use modglobal, only : rtimee,tboundary,ntboundary
    use modmpi, only : myidx,myidy
    implicit none
    real, dimension(:), intent(in) :: x,y,z
    integer, intent(in) :: nx,ny,nz,iuturb,ib
    real, dimension(:,:), intent(out) :: uturb,thlturb,qtturb
    integer,target :: i,j,k,ipatch,jpatch,kpatch
    integer :: ii,itp,itm
    integer, pointer :: pi1, pi2,pi1patch,pi2patch
    real, dimension(3) :: xx,ci
    real, dimension(3,3) :: eigvec
    real :: wthl,wqt,w2,thl2,qt2,utemp,vtemp,wtemp,wturbf,thltemp,qttemp
    real :: t,fp,fm,rho
    t = rtimee
    itm = 1
    ! Interpolate covariance to current time
    if(ntboundary>1) then
      do while(t>tboundary(itm))
        itm = itm+1
      end do
      if(t>tboundary(1)) then
        itm = itm-1
      end if
      itp = itm+1
      fm = (tboundary(itp)-t)/(tboundary(itp)-tboundary(itm))
      fp = (t-tboundary(itm))/(tboundary(itp)-tboundary(itm))
    else
      itm = 1; itp = 1; fp  = 1.; fm  = 1.
    endif
    select case(ib)
    case(1,2)
      pi1 => j; pi2 => k; pi1patch => jpatch; pi2patch => kpatch
    case(3,4)
      pi1 => i; pi2 => k; pi1patch => ipatch; pi2patch => kpatch
    case(5)
      pi1 => i; pi2 => j; pi1patch => ipatch; pi2patch => jpatch
    end select
    do i = 1,nx
      xx(1) = x(i); ipatch = min(int(x(i)/dxturb)+1,nxturb)
    do j = 1,ny
      xx(2) = y(j); jpatch = min(int(y(j)/dyturb)+1,nyturb)
    do k = 1,nz
      xx(3) = z(k); kpatch = min(k,nzturb)
      ci     = boundary(ib)%ci(pi1patch,pi2patch,:)
      eigvec = boundary(ib)%eigvec(pi1patch,pi2patch,:,:)
      utemp = 0.; vtemp = 0.; wtemp = 0.; thltemp = 0.; qttemp = 0.
      do ii = 1,nmodes
        utemp = utemp + p(ii,1)*cos(dot(kn(ii,:)/ci(1)*ctot,xx/lambda)+omega(ii)*t/tau) + &
                      & q(ii,1)*sin(dot(kn(ii,:)/ci(1)*ctot,xx/lambda)+omega(ii)*t/tau)
        vtemp = vtemp + p(ii,2)*cos(dot(kn(ii,:)/ci(2)*ctot,xx/lambda)+omega(ii)*t/tau) + &
                      & q(ii,2)*sin(dot(kn(ii,:)/ci(2)*ctot,xx/lambda)+omega(ii)*t/tau)
        wtemp = wtemp + p(ii,3)*cos(dot(kn(ii,:)/ci(3)*ctot,xx/lambda)+omega(ii)*t/tau) + &
                      & q(ii,3)*sin(dot(kn(ii,:)/ci(3)*ctot,xx/lambda)+omega(ii)*t/tau)
        thltemp = thltemp + p_thl(ii)*cos(dot(k_thl(ii,:),xx/lambdasxyz)+omega_thl(ii)*t/tau) + &
                          & q_thl(ii)*sin(dot(k_thl(ii,:),xx/lambdasxyz)+omega_thl(ii)*t/tau)
        qttemp = qttemp + p_qt(ii)*cos(dot(k_qt(ii,:),xx/lambdasxyz)+omega_qt(ii)*t/tau) + &
                        & q_qt(ii)*sin(dot(k_qt(ii,:),xx/lambdasxyz)+omega_qt(ii)*t/tau)
      end do
      ! Scale velocity fields
      utemp = utemp*ci(1)
      vtemp = vtemp*ci(2)
      wtemp = wtemp*ci(3)
      ! Reproject to cartesian velocity pertubations
      uturb(pi1,pi2)  = nisqrt*dot(eigvec(iuturb,:),(/utemp,vtemp,wtemp/))
      ! Calculate thlturb and qtturb
      thltemp = thltemp * nisqrt2
      qttemp  = qttemp * nisqrt2
      wturbf = nisqrt*dot(eigvec(3,:),(/utemp,vtemp,wtemp/))
      wthl   = fp*boundary(ib)%wthl(pi1patch,pi2patch,itp) + fm*boundary(ib)%wthl(pi1patch,pi2patch,itm)
      wqt    = fp*boundary(ib)%wqt(pi1patch,pi2patch,itp)  + fm*boundary(ib)%wqt(pi1patch,pi2patch,itm)
      w2     = fp*boundary(ib)%w2(pi1patch,pi2patch,itp)   + fm*boundary(ib)%w2(pi1patch,pi2patch,itm)
      thl2   = fp*boundary(ib)%thl2(pi1patch,pi2patch,itp) + fm*boundary(ib)%thl2(pi1patch,pi2patch,itm)
      qt2    = fp*boundary(ib)%qt2(pi1patch,pi2patch,itp)  + fm*boundary(ib)%qt2(pi1patch,pi2patch,itm)
      if(thl2==0. .or. w2==0.) then
        thlturb(pi1,pi2) = thltemp*sqrt(thl2)
      else
        rho = min(max(wthl/sqrt(thl2*w2),-1.),1.)
        thlturb(pi1,pi2) = (rho*wturbf/sqrt(w2) + sqrt(1-rho**2)*thltemp)*sqrt(thl2)
      endif
      if(qt2==0. .or. w2==0.) then
        qtturb(pi1,pi2)  = qttemp*sqrt(qt2)
      else
        rho = min(max(wqt/sqrt(qt2*w2),-1.),1.)
        qtturb(pi1,pi2)  = (rho*wturbf/sqrt(w2) + sqrt(1-rho**2)*qttemp)*sqrt(qt2)
      endif
      !if(ib==1) then
      !  vturb(pi1,pi2) = nisqrt*dot(eigvec(2,:),(/utemp,vtemp,wtemp/))
      !  wturb(pi1,pi2) = nisqrt*dot(eigvec(3,:),(/utemp,vtemp,wtemp/))
      !endif
    end do
    end do
    end do
    ! ! Uncommend for netcdf output turbulent pertubations
    ! if(ib==1) then
    !   call check( nf90_put_var(ncid, uturb_varid, uturb, start = start, &
    !                            count = count) )
    !   call check( nf90_put_var(ncid, vturb_varid, vturb, start = start, &
    !                            count = count) )
    !   call check( nf90_put_var(ncid, wturb_varid, wturb, start = start, &
    !                            count = count) )
    !   call check( nf90_put_var(ncid, thlturb_varid, thlturb, start = start, &
    !                            count = count) )
    !   call check( nf90_put_var(ncid, qtturb_varid, qtturb, start = start, &
    !                            count = count) )
    !   call check( nf90_put_var(ncid, t_varid, (/t/), start = (/start(3)/), &
    !                            count = (/count(3)/)) )
    !   start(3) = start(3)+1
    ! endif
    nullify(pi1,pi2,pi1patch,pi2patch)
  end subroutine calc_pert2

  function gaussrand(mu,sigma)
    use modglobal, only : pi,eps1
    implicit none
    real :: gaussrand
    real, intent(in) :: mu, sigma
    real :: temp1,temp2
    temp1 = 0.
    do while(temp1<eps1)
      temp1 = real(getRandomReal(noise))
    end do
    temp2 = 0.
    do while(temp2<eps1)
      temp2 = real(getRandomReal(noise))
    end do
    gaussrand   = sqrt(-2.*log(temp1))*cos(2.*pi*temp2)*sigma+mu
    !val(i*2) = sqrt(-2.*log(temp1))*sin(2.*pi*temp2)*sigma+mu
  end function gaussrand

  function cross(a, b)
    implicit none
    real, dimension(3) :: cross
    real, dimension(3), intent(in) :: a, b

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  end function cross

  function dot(a, b)
    implicit none
    real :: dot
    real, dimension(3), intent(in) :: a, b

    dot = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
  end function dot

  ! Required routines for eigenvalues and eigenvectors
  ! ----------------------------------------------------------------------------
  ! Numerical diagonalization of 3x3 matrcies
  ! Copyright (C) 2006  Joachim Kopp
  ! ----------------------------------------------------------------------------
  ! This library is free software; you can redistribute it and/or
  ! modify it under the terms of the GNU Lesser General Public
  ! License as published by the Free Software Foundation; either
  ! version 2.1 of the License, or (at your option) any later version.
  !
  ! This library is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  ! Lesser General Public License for more details.
  !
  ! You should have received a copy of the GNU Lesser General Public
  ! License along with this library; if not, write to the Free Software
  ! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
  ! ----------------------------------------------------------------------------


  ! ----------------------------------------------------------------------------
        SUBROUTINE DSYEVJ3(A, Q, W)
  ! ----------------------------------------------------------------------------
  ! Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
  ! matrix A using the Jacobi algorithm.
  ! The upper triangular part of A is destroyed during the calculation,
  ! the diagonal elements are read but not destroyed, and the lower
  ! triangular elements are not referenced at all.
  ! ----------------------------------------------------------------------------
  ! Parameters:
  !   A: The symmetric input matrix
  !   Q: Storage buffer for eigenvectors
  !   W: Storage buffer for eigenvalues
  ! ----------------------------------------------------------------------------
  !     .. Arguments ..
        DOUBLE PRECISION A(3,3)
        DOUBLE PRECISION Q(3,3)
        DOUBLE PRECISION W(3)

  !     .. Parameters ..
        INTEGER          N
        PARAMETER        ( N = 3 )

  !     .. Local Variables ..
        DOUBLE PRECISION SD, SO
        DOUBLE PRECISION S, C, T
        DOUBLE PRECISION G, H, Z, THETA
        DOUBLE PRECISION THRESH
        INTEGER          I, X, Y, R

  !     Initialize Q to the identitity matrix
  !     --- This loop can be omitted if only the eigenvalues are desired ---
        DO 10 X = 1, N
          Q(X,X) = 1.0D0
          DO 11, Y = 1, X-1
            Q(X, Y) = 0.0D0
            Q(Y, X) = 0.0D0
     11   CONTINUE
     10 CONTINUE

  !     Initialize W to diag(A)
        DO 20 X = 1, N
          W(X) = A(X, X)
     20 CONTINUE

  !     Calculate SQR(tr(A))
        SD = 0.0D0
        DO 30 X = 1, N
          SD = SD + ABS(W(X))
     30 CONTINUE
        SD = SD**2

  !     Main iteration loop
        DO 40 I = 1, 50
  !       Test for convergence
          SO = 0.0D0
          DO 50 X = 1, N
            DO 51 Y = X+1, N
              SO = SO + ABS(A(X, Y))
     51     CONTINUE
     50   CONTINUE
          IF (SO .EQ. 0.0D0) THEN
            RETURN
          END IF

          IF (I .LT. 4) THEN
            THRESH = 0.2D0 * SO / N**2
          ELSE
            THRESH = 0.0D0
          END IF

  !       Do sweep
          DO 60 X = 1, N
            DO 61 Y = X+1, N
              G = 100.0D0 * ( ABS(A(X, Y)) )
              IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X)) &
       &                    .AND. ABS(W(Y)) + G .EQ. ABS(W(Y)) ) THEN
                A(X, Y) = 0.0D0
              ELSE IF (ABS(A(X, Y)) .GT. THRESH) THEN
  !             Calculate Jacobi transformation
                H = W(Y) - W(X)
                IF ( ABS(H) + G .EQ. ABS(H) ) THEN
                  T = A(X, Y) / H
                ELSE
                  THETA = 0.5D0 * H / A(X, Y)
                  IF (THETA .LT. 0.0D0) THEN
                    T = -1.0D0 / (SQRT(1.0D0 + THETA**2) - THETA)
                  ELSE
                    T = 1.0D0 / (SQRT(1.0D0 + THETA**2) + THETA)
                  END IF
                END IF

                C = 1.0D0 / SQRT( 1.0D0 + T**2 )
                S = T * C
                Z = T * A(X, Y)

  !             Apply Jacobi transformation
                A(X, Y) = 0.0D0
                W(X)    = W(X) - Z
                W(Y)    = W(Y) + Z
                DO 70 R = 1, X-1
                  T       = A(R, X)
                  A(R, X) = C * T - S * A(R, Y)
                  A(R, Y) = S * T + C * A(R, Y)
     70         CONTINUE
                DO 80, R = X+1, Y-1
                  T       = A(X, R)
                  A(X, R) = C * T - S * A(R, Y)
                  A(R, Y) = S * T + C * A(R, Y)
     80         CONTINUE
                DO 90, R = Y+1, N
                  T       = A(X, R)
                  A(X, R) = C * T - S * A(Y, R)
                  A(Y, R) = S * T + C * A(Y, R)
     90         CONTINUE

  !             Update eigenvectors
  !             --- This loop can be omitted if only the eigenvalues are desired ---
                DO 100, R = 1, N
                  T       = Q(R, X)
                  Q(R, X) = C * T - S * Q(R, Y)
                  Q(R, Y) = S * T + C * Q(R, Y)
    100         CONTINUE
              END IF
     61     CONTINUE
     60   CONTINUE
     40 CONTINUE

        PRINT *, "DSYEVJ3: No convergence."

        END SUBROUTINE
  ! End of subroutine DSYEVJ3
end module modsynturb
