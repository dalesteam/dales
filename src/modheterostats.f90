!----------------------------------------------------------------------------
! This file is part of DALES.
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
! Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!----------------------------------------------------------------------------
!
!
module modheterostats

use modglobal, only: nsv, kmax,longint

implicit none
private
public :: initheterostats, heterostats, exitheterostats

save
  !namelist variables
  real                 :: dtav
  integer(kind=longint):: idtav,tnext
  logical              :: lheterostats  = .false.
  integer              :: ncklimit

  !VARIABLES FOR STATISTICS
  !id of netcdf file
  integer :: ncid

  !id of dimensions
  integer :: xid, yid, zid, tid

  !id of variables (means)
  integer :: uavgid, vavgid, wavgid, thlavgid, thvavgid, qtavgid, qlavgid, eavgid
  integer, allocatable :: svavgid(:)
  integer :: lwpid

  !id of variables (variances)
  integer :: uvarid, vvarid, wvarid, thlvarid, thvvarid, qtvarid, qlvarid
  integer, allocatable :: svvarid(:)

  !id of variables (covariances)
  integer :: uwcovid, vwcovid
  integer :: uwcovsid, vwcovsid
  integer :: wthlcovid
  integer :: wthlcovsid
  integer :: wthvcovid
  integer :: wthvcovsid
  integer :: wqtcovid
  integer :: wqtcovsid
  integer :: wqlcovid
  integer :: wqlcovsid
  integer :: thlqcovid
  integer, allocatable :: usvcovid(:), vsvcovid(:), wsvcovid(:)
  integer, allocatable :: wsvcovsid(:)

  !Only used in chemistry cases: integer :: OHISOcovid, O3NOcovid


  !COMMON VARIABLES
  !constants
  !integer ncklimit
  !parameter (ncklimit = 50)
  integer :: nccall = 1


contains
  subroutine initheterostats

    use typeSizes
    use netcdf
    use modmpi
    use modglobal

    implicit none

    integer             :: status
    character(len = 20) :: ncfile
    integer             :: n, ierr
    character(20)       :: filename

    namelist/NAMHETEROSTATS/ &
    dtav,lheterostats,ncklimit

    ncklimit = kmax
    dtav = dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMHETEROSTATS,iostat=ierr)
      write(6, NAMHETEROSTATS)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(dtav        ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lheterostats,1,0,comm3d,mpierr)
    call D_MPI_BCAST(ncklimit    ,1,0,comm3d,mpierr)

    if(.not.(lheterostats)) return
    idtav = dtav/tres
    tnext = idtav+btime
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'HETEROSTATS: dtav should be a integer multiple of dtmax'
    end if

    if(myidx .ne. 0) return

    allocate(svavgid(nsv))
    allocate(svvarid(nsv))
    allocate(usvcovid(nsv), vsvcovid(nsv), wsvcovid(nsv))
    allocate(wsvcovsid(nsv))

    ncfile = 'heterostats123.nc'
    write(ncfile(12:14),'(i3.3)') myidy
!    write(6,*) "HETEROSTATS: Creating: ", ncfile

    !create file
    status = nf90_create(ncfile, nf90_clobber, ncid)
    if (status /= nf90_noerr) call nchandle_error(status)

    !create dimensions
    status = nf90_def_dim(ncid, "x", itot, xid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_dim(ncid, "y", jmax, yid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_dim(ncid, "z", ncklimit, zid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_dim(ncid, "t", nf90_unlimited, tid)
    if (status /= nf90_noerr) call nchandle_error(status)

    !create variables
    !means
    status = nf90_def_var(ncid, "uavg", nf90_float, (/yid, zid, tid/), uavgid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "vavg", nf90_float, (/yid, zid, tid/), vavgid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "wavg", nf90_float, (/yid, zid, tid/), wavgid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "thlavg", nf90_float, (/yid, zid, tid/), thlavgid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "thvavg", nf90_float, (/yid, zid, tid/), thvavgid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "qtavg", nf90_float, (/yid, zid, tid/), qtavgid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "qlavg", nf90_float, (/yid, zid, tid/), qlavgid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "eavg", nf90_float, (/yid, zid, tid/), eavgid)
    if (status /= nf90_noerr) call nchandle_error(status)
    do n=1,nsv
      filename = "svnnnavg"
      write (filename(3:5),'(i3.3)') n
      status = nf90_def_var(ncid, filename, nf90_float, (/yid, zid, tid/), svavgid(n))
      if (status /= nf90_noerr) call nchandle_error(status)
    enddo

    status = nf90_def_var(ncid, "lwp", nf90_float, (/xid,yid, tid/), lwpid)
    if (status /= nf90_noerr) call nchandle_error(status)

    !variances
    status = nf90_def_var(ncid, "uvar", nf90_float, (/yid, zid, tid/), uvarid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "vvar", nf90_float, (/yid, zid, tid/), vvarid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "wvar", nf90_float, (/yid, zid, tid/), wvarid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "thlvar", nf90_float, (/yid, zid, tid/), thlvarid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "thvvar", nf90_float, (/yid, zid, tid/), thvvarid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "qtvar", nf90_float, (/yid, zid, tid/), qtvarid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "qlvar", nf90_float, (/yid, zid, tid/), qlvarid)
    if (status /= nf90_noerr) call nchandle_error(status)
    do n=1,nsv
      filename = "svnnnvar"
      write (filename(3:5),'(i3.3)') n
      status = nf90_def_var(ncid, filename, nf90_float, (/yid, zid, tid/), svvarid(n))
      if (status /= nf90_noerr) call nchandle_error(status)
    enddo

    !covariances
    status = nf90_def_var(ncid, "uwcov", nf90_float, (/yid, zid, tid/), uwcovid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "uwcovs", nf90_float, (/yid, zid, tid/), uwcovsid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "vwcov", nf90_float, (/yid, zid, tid/), vwcovid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "vwcovs", nf90_float, (/yid, zid, tid/), vwcovsid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "wthlcov", nf90_float, (/yid, zid, tid/), wthlcovid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "wthlcovs", nf90_float, (/yid, zid, tid/), wthlcovsid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "wthvcov", nf90_float, (/yid, zid, tid/), wthvcovid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "wthvcovs", nf90_float, (/yid, zid, tid/), wthvcovsid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "wqtcov", nf90_float, (/yid, zid, tid/), wqtcovid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "wqlcov", nf90_float, (/yid, zid, tid/), wqlcovid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "wqtcovs", nf90_float, (/yid, zid, tid/), wqtcovsid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "wqlcovs", nf90_float, (/yid, zid, tid/), wqlcovsid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "thlqcov", nf90_float, (/yid, zid, tid/), thlqcovid)
    if (status /= nf90_noerr) call nchandle_error(status)
    do n=1,nsv
      filename = "wsvnnncov"
      write (filename(4:6),'(i3.3)') n
      status = nf90_def_var(ncid, filename, nf90_float, (/yid, zid, tid/), wsvcovid(n))
      if (status /= nf90_noerr) call nchandle_error(status)
      filename = "wsvnnncovs"
      write (filename(4:6),'(i3.3)') n
      status = nf90_def_var(ncid, filename, nf90_float, (/yid, zid, tid/), wsvcovsid(n))
      if (status /= nf90_noerr) call nchandle_error(status)
    enddo

    !turn off define mode
    status = nf90_enddef(ncid)
    if (status /= nf90_noerr) call nchandle_error(status)
  end subroutine initheterostats


  subroutine heterostats

    use modglobal, only : rk3step,dt_lim,timee
    implicit none

    if (.not. lheterostats) return
    if (rk3step/=3) return
    if (timee==0) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+idtav
    dt_lim = min(dt_lim,tnext-timee)

    call do_heterostats
    nccall = nccall + 1

  end subroutine heterostats


  subroutine do_heterostats

    use typeSizes
    use netcdf
    use modfields
    use modmpi,     only : myidx,commrow,mpierr,gatherrow,MPI_SUM,D_MPI_ALLREDUCE &
                         , D_MPI_REDUCE
    use modsurfdata
    use modsubgrid, only : ekm, ekh
    use modglobal,  only : iadv_sv, iadv_kappa, dzf, dzh, rlv, cp, rv, &
                           rd, imax, jmax, i1, j1, k1, ih, jh, itot
    use advec_kappa,only : halflev_kappa
    implicit none

    integer n,i,j,k
    integer status

    real, dimension(jmax,ncklimit)     :: uavg, vavg, wavg, thlavg, thvavg, qtavg, qlavg, eavg,&
    thlhavg, thvhavg, qthavg, qlhavg, vonwavg, uonwavg
    real, dimension(imax,jmax)         :: lwpavg
    real, dimension(jmax,ncklimit)     :: uvar, vvar, wvar, thlvar, thvvar, qtvar, qlvar
    real, dimension(jmax,ncklimit)     :: uwcov, uwcovs, vwcov, vwcovs
    real, dimension(jmax,ncklimit)     :: wthlcov, wthlcovs, wthvcov, wthvcovs, wqtcov, wqtcovs, thlqcov, wqlcov, wqlcovs
    real, dimension(jmax,ncklimit,nsv) :: svavg, svhavg, svvar,  wsvcov, wsvcovs

    real  vonw(2-ih:i1+ih,2-jh:j1+jh,k1),uonw(2-ih:i1+ih,2-jh:j1+jh,k1)
    real(field_r) sv0h(2-ih:i1+ih,2-jh:j1+jh,k1,nsv),thv0(2-ih:i1+ih,2-jh:j1+jh,k1)

    real  buffer(jmax,ncklimit)
    real  bufferij(itot,jmax)

    real  qs0h, t0h, den, c1, c2

    uavg(:,:)     = 0.0
    vavg(:,:)     = 0.0
    wavg(:,:)     = 0.0
    thlavg(:,:)   = 0.0
    thvavg(:,:)   = 0.0
    eavg(:,:)     = 0.0
    qtavg(:,:)    = 0.0
    qlavg(:,:)    = 0.0

    thlhavg(:,:)  = 0.0
    thvhavg(:,:)  = 0.0
    qthavg(:,:)   = 0.0
    qlhavg(:,:)   = 0.0
    vonwavg(:,:)  = 0.0
    uonwavg(:,:)  = 0.0

    lwpavg(:,:)   = 0.0

    uvar(:,:)     = 0.0
    vvar(:,:)     = 0.0
    wvar(:,:)     = 0.0
    thlvar(:,:)   = 0.0
    thvvar(:,:)   = 0.0
    qtvar(:,:)    = 0.0
    qlvar(:,:)    = 0.0

    uwcov(:,:)    = 0.0
    uwcovs(:,:)   = 0.0
    vwcov(:,:)    = 0.0
    vwcovs(:,:)   = 0.0
    wthlcov(:,:)  = 0.0
    wthlcovs(:,:) = 0.0
    wthvcov(:,:)  = 0.0
    wthvcovs(:,:) = 0.0
    wqtcov(:,:)   = 0.0
    wqlcov(:,:)   = 0.0
    wqtcovs(:,:)  = 0.0
    wqlcovs(:,:)  = 0.0
    thlqcov(:,:)  = 0.0

    svavg(:,:,:)  = 0.0
    svhavg(:,:,:) = 0.0
    svvar(:,:,:)  = 0.0
    wsvcov(:,:,:) = 0.0
    wsvcovs(:,:,:)= 0.0

    buffer(:,:)   = 0.0
    bufferij(:,:) = 0.0

    !calculate averages and store them

    !Prepare data
    do k = 2,k1
      do j = 2,j1
        do i = 2,i1
          vonw(i,j,k) = 0.25*(v0(i,j,k) + v0(i,j+1,k) + v0(i,j,k-1) + v0(i,j+1,k-1))
        end do
      end do
    end do

    do k = 2,k1
      do j = 2,j1
        do i = 2,i1
          uonw(i,j,k) = 0.25*(u0(i,j,k) + u0(i+1,j,k) + u0(i,j,k-1) + u0(i+1,j,k-1))
        end do
      end do
    end do

    do n=1,nsv
      if (iadv_sv(n)==iadv_kappa) then
         call halflev_kappa(sv0(2-ih:i1+ih,2-jh:j1+jh,1:k1,n),sv0h(:,:,:,n))
      else
        do  k=2,k1
          do  j=2,j1
            do  i=2,i1
              sv0h(i,j,k,n) = (sv0(i,j,k,n)*dzf(k-1)+sv0(i,j,k-1,n)*dzf(k))/(2*dzh(k))
            enddo
          enddo
        enddo
      end if
    enddo

    do  k=1,k1
      do  j=2,j1
        do  i=2,i1
          thv0(i,j,k) = (thl0(i,j,k)+rlv*ql0(i,j,k)/(cp*exnf(k))) &
                        *(1+(rv/rd-1)*qt0(i,j,k)-rv/rd*ql0(i,j,k))
        enddo
      enddo
    enddo


    !LOOPS ARE NOT PUT IN FUNCTION BECAUSE OF ARRAY DEFINITIONS WHICH DIFFER AMONG VARIABLES!
    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          uavg(j,k) = uavg(j,k) + u0(i+1,j+1,k)
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          vavg(j,k) = vavg(j,k) + v0(i+1,j+1,k)
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          wavg(j,k) = wavg(j,k) + w0(i+1,j+1,k)
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          thlavg(j,k) = thlavg(j,k) + thl0(i+1,j+1,k)
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          thvavg(j,k) = thvavg(j,k) + thv0(i+1,j+1,k)
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          qtavg(j,k) = qtavg(j,k) + qt0(i+1,j+1,k)
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          qlavg(j,k) = qlavg(j,k) + ql0(i+1,j+1,k)
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          eavg(j,k) = eavg(j,k) + (e120(i+1,j+1,k))**2.0
        end do
      end do
    end do

    !create average for thl0h - only used for cov
    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          thlhavg(j,k) = thlhavg(j,k) + thl0h(i+1,j+1,k)
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          thvhavg(j,k) = thvhavg(j,k) + thv0h(i+1,j+1,k)
        end do
      end do
    end do

    !create average for qt0h - only used for cov
    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          qthavg(j,k) = qthavg(j,k) + qt0h(i+1,j+1,k)
        end do
      end do
    end do

    !create average for ql0h - only used for cov
    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          qlhavg(j,k) = qlhavg(j,k) + ql0h(i+1,j+1,k)
        end do
      end do
    end do

    !create average for v projected on w - only used for cov, start at level 2
    do k = 2,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          vonwavg(j,k) = vonwavg(j,k) + vonw(i+1,j+1,k)
        end do
      end do
    end do

    !create average for u projected on w - only used for cov, start at level 2
    do k = 2,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          uonwavg(j,k) = uonwavg(j,k) + uonw(i+1,j+1,k)
        end do
      end do
    end do

    do n = 1,nsv
      do k = 1,ncklimit
        do j = 1,jmax
          do i = 1,imax
            !shift prognostic fields one step as 1st column
            !is dummy column because of MPI and periodicity
            svavg(j,k,n) = svavg(j,k,n) + sv0(i+1,j+1,k,n)
          end do
        end do
      end do

      do k = 2,ncklimit
        do j = 1,jmax
          do i = 1,imax
            !shift prognostic fields one step as 1st column
            !is dummy column because of MPI and periodicity
            svhavg(j,k,n) = svhavg(j,k,n) + sv0h(i+1,j+1,k,n)
          end do
        end do
      end do
    enddo


    uavg   = uavg / itot
    vavg   = vavg / itot
    wavg   = wavg / itot
    thlavg = thlavg / itot
    thvavg = thvavg / itot
    qtavg  = qtavg / itot
    qlavg  = qlavg / itot
    eavg   = eavg / itot
    svavg  = svavg / itot

    thlhavg = thlhavg / itot
    thvhavg = thvhavg / itot
    qthavg  = qthavg / itot
    qlhavg  = qlhavg / itot
    vonwavg = vonwavg / itot
    uonwavg = uonwavg / itot
    svhavg  = svhavg / itot

    ! Communicate all the averages to the processors in the row
    ! copy the results back from the mpi buffer because they are needed
    ! to calculate the variances
    ! TODO: 1. create fast path: the mpi calls are slow and are not needed for nprocx = 1
    !       2. this function could be cleaned up a lot

    call D_MPI_ALLREDUCE(uavg,buffer,ncklimit*jmax,MPI_SUM,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, uavgid, buffer, (/1,1,nccall/), (/jmax, ncklimit , 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif
    uavg = buffer

    call D_MPI_ALLREDUCE(vavg,buffer,ncklimit*jmax,MPI_SUM,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, vavgid, buffer, (/1,1,nccall/), (/jmax, ncklimit , 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif
    vavg = buffer

    call D_MPI_ALLREDUCE(wavg,buffer,ncklimit*jmax,MPI_SUM,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, wavgid, buffer, (/1,1,nccall/), (/jmax, ncklimit , 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif
    wavg = buffer

    call D_MPI_ALLREDUCE(thlavg,buffer,ncklimit*jmax,MPI_SUM,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, thlavgid, buffer, (/1,1,nccall/), (/jmax, ncklimit , 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif
    thlavg = buffer

    call D_MPI_ALLREDUCE(thvavg,buffer,ncklimit*jmax,MPI_SUM,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, thvavgid, buffer, (/1,1,nccall/), (/jmax, ncklimit , 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif
    thvavg = buffer

    call D_MPI_ALLREDUCE(qtavg,buffer,ncklimit*jmax,MPI_SUM,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, qtavgid, buffer, (/1,1,nccall/), (/jmax, ncklimit , 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif
    qtavg = buffer

    call D_MPI_ALLREDUCE(qlavg,buffer,ncklimit*jmax,MPI_SUM,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, qlavgid, buffer, (/1,1,nccall/), (/jmax, ncklimit , 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif
    qlavg = buffer

    call D_MPI_ALLREDUCE(eavg,buffer,ncklimit*jmax,MPI_SUM,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, eavgid, buffer, (/1,1,nccall/), (/jmax, ncklimit , 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif
    eavg = buffer

    do n=1,nsv
      call D_MPI_ALLREDUCE(svavg(:,:,n),buffer,ncklimit*jmax,MPI_SUM,commrow,mpierr)
      if(myidx == 0) then
         status = nf90_put_var(ncid, svavgid(n), buffer, (/1,1,nccall/), (/jmax, ncklimit , 1/))
         if(status /= nf90_noerr) call nchandle_error(status)
      endif
      svavg(:,:,n) = buffer
    enddo

    call D_MPI_ALLREDUCE(thlhavg,buffer,ncklimit*jmax,MPI_SUM,commrow,mpierr)
    thlhavg = buffer

    call D_MPI_ALLREDUCE(thvhavg,buffer,ncklimit*jmax,MPI_SUM,commrow,mpierr)
    thvhavg = buffer

    call D_MPI_ALLREDUCE(qthavg,buffer,ncklimit*jmax,MPI_SUM,commrow,mpierr)
    qthavg = buffer

    call D_MPI_ALLREDUCE(qlhavg,buffer,ncklimit*jmax,MPI_SUM,commrow,mpierr)
    qlhavg = buffer

    call D_MPI_ALLREDUCE(vonwavg,buffer,ncklimit*jmax,MPI_SUM,commrow,mpierr)
    vonwavg = buffer

    call D_MPI_ALLREDUCE(uonwavg,buffer,ncklimit*jmax,MPI_SUM,commrow,mpierr)
    uonwavg = buffer

    do n = 1,nsv
      call D_MPI_ALLREDUCE(svhavg(:,:,n),buffer,ncklimit*jmax,MPI_SUM,commrow,mpierr)
      svhavg(:,:,n) = buffer
    enddo

    !calculate liquid water path and store it
    do j = 1,jmax
      do i = 1,imax
        do k = 1,kmax
          lwpavg(i,j) = lwpavg(i,j) + ql0(i+1,j+1,k)*rhof(k)*dzf(k)
        end do
      end do
    end do

    call gatherrow(lwpavg,bufferij,imax,jmax,itot)
    if(myidx == 0) then
       status = nf90_put_var(ncid, lwpid, bufferij, (/1,1,nccall/), (/itot, jmax, 1/))
       if(status /= nf90_noerr) call nchandle_error(status)
    endif

    !calculate variances and store them
    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          uvar(j,k) = uvar(j,k) + (u0(i+1,j+1,k)-uavg(j,k))**2.
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          vvar(j,k) = vvar(j,k) + (v0(i+1,j+1,k)-vavg(j,k))**2.
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          wvar(j,k) = wvar(j,k) + (w0(i+1,j+1,k)-wavg(j,k))**2.
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          thlvar(j,k) = thlvar(j,k) + (thl0(i+1,j+1,k)-thlavg(j,k))**2.
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          thvvar(j,k) = thvvar(j,k) + (thv0(i+1,j+1,k)-thvavg(j,k))**2.
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          qtvar(j,k) = qtvar(j,k) + (qt0(i+1,j+1,k)-qtavg(j,k))**2.
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          qlvar(j,k) = qlvar(j,k) + (ql0(i+1,j+1,k)-qlavg(j,k))**2.
        end do
      end do
    end do

    do n = 1,nsv
      do k = 1,ncklimit
        do j = 1,jmax
          do i = 1,imax
            svvar(j,k,n) = svvar(j,k,n) + (sv0(i+1,j+1,k,n)-svavg(j,k,n))*(sv0(i+1,j+1,k,n)-svavg(j,k,n))
          end do
        end do
      end do
    end do

    uvar = uvar / itot
    vvar = vvar / itot
    wvar = wvar / itot
    thlvar = thlvar / itot
    thvvar = thvvar / itot
    qtvar = qtvar / itot
    qlvar = qlvar / itot
    svvar = svvar / itot

    call D_MPI_REDUCE(uvar,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, uvarid, buffer, (/1,1,nccall/), (/jmax, ncklimit , 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    call D_MPI_REDUCE(vvar,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, vvarid, buffer, (/1,1,nccall/), (/jmax, ncklimit , 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    call D_MPI_REDUCE(wvar,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, wvarid, buffer, (/1,1,nccall/), (/jmax, ncklimit , 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    call D_MPI_REDUCE(thlvar,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, thlvarid, buffer, (/1,1,nccall/), (/jmax, ncklimit, 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    call D_MPI_REDUCE(thvvar,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, thvvarid, buffer, (/1,1,nccall/), (/jmax, ncklimit, 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    call D_MPI_REDUCE(qtvar,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, qtvarid, buffer, (/1,1,nccall/), (/jmax, ncklimit, 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    call D_MPI_REDUCE(qlvar,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, qlvarid, buffer, (/1,1,nccall/), (/jmax, ncklimit, 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    do n=1,nsv
      call D_MPI_REDUCE(svvar(:,:,n),buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
      if(myidx == 0) then
        status = nf90_put_var(ncid, svvarid(n), buffer, (/1,1,nccall/), (/jmax, ncklimit , 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
      endif
    enddo

    !calculate covariances and store them
    do k = 2,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          uwcov(j,k) = uwcov(j,k) + (uonw(i+1,j+1,k)-uonwavg(j,k))*(w0(i+1,j+1,k)-wavg(j,k))
        end do
      end do
    end do


    do k = 2,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          vwcov(j,k) = vwcov(j,k) + (vonw(i+1,j+1,k)-vonwavg(j,k))*(w0(i+1,j+1,k)-wavg(j,k))
        end do
      end do
    end do


    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          if(k==1) then
            uwcovs(j,k) = uwcovs(j,k) - max(ustar(i,j)**2.0, 1.e-10)
          else
            uwcovs(j,k) = uwcovs(j,k) - 0.5*(dzf(k-1)*ekm(i+1,j+1,k)+dzf(k)*ekm(i+1,j+1,k-1))/dzh(k) *&
                          (0.5*(u0(i+1,j+1,k)+u0(i+2,j+1,k)) - 0.5*(u0(i+1,j+1,k-1)+u0(i+2,j+1,k-1))) / dzh(k)
          endif
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          if(k==1) then
            vwcovs(j,k) = vwcovs(j,k) - max(ustar(i,j)**2.0, 1.e-10)
          else
            vwcovs(j,k) = vwcovs(j,k) - 0.5*(dzf(k-1)*ekm(i+1,j+1,k)+dzf(k)*ekm(i+1,j+1,k-1))/dzh(k) *&
            (0.5*(v0(i+1,j+1,k)+v0(i+1,j+2,k)) - 0.5*(v0(i+1,j+1,k-1)+v0(i+1,j+2,k-1))) / dzh(k)
          endif
        end do
      end do
    end do

    do k = 2,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          wthlcov(j,k) = wthlcov(j,k) + (w0(i+1,j+1,k)-wavg(j,k)) * (thl0h(i+1,j+1,k)-thlhavg(j,k))
        end do
      end do
    end do


    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          if(k==1) then
            wthlcovs(j,k) = wthlcovs(j,k) + thlflux(i+1,j+1)
          else
            wthlcovs(j,k) = wthlcovs(j,k) - 0.5*(dzf(k-1)*ekh(i+1,j+1,k)+dzf(k)*ekh(i+1,j+1,k-1))/dzh(k) &
            * (thl0(i+1,j+1,k) - thl0(i+1,j+1,k-1)) / dzh(k)
          endif
        end do
      end do
    end do

    do k = 2,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          wqtcov(j,k) = wqtcov(j,k) + (w0(i+1,j+1,k)-wavg(j,k)) * (qt0h(i+1,j+1,k)-qthavg(j,k))
        end do
      end do
    end do

    do k = 2,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          wqlcov(j,k) = wqlcov(j,k) + (w0(i+1,j+1,k)-wavg(j,k)) * (ql0h(i+1,j+1,k)-qlhavg(j,k))
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          if(k==1) then
            wqtcovs(j,k) = wqtcovs(j,k) + qtflux(i+1,j+1)
          else
            wqtcovs(j,k) = wqtcovs(j,k) - 0.5*(dzf(k-1)*ekh(i+1,j+1,k)+dzf(k)*ekh(i+1,j+1,k-1))/dzh(k) * &
            (qt0(i+1,j+1,k) - qt0(i+1,j+1,k-1)) / dzh(k)
          endif
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          if(k==1) then
            wqlcovs(j,k) = 0.0
          else
            wqlcovs(j,k) = wqlcovs(j,k) - 0.5*(dzf(k-1)*ekh(i+1,j+1,k)+dzf(k)*ekh(i+1,j+1,k-1))/dzh(k) * &
            (ql0(i+1,j+1,k) - ql0(i+1,j+1,k-1)) / dzh(k)
          endif
        end do
      end do
    end do

    do k = 2,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          wthvcov(j,k) = wthvcov(j,k) + (w0(i+1,j+1,k)-wavg(j,k)) * (thv0h(i+1,j+1,k)-thvhavg(j,k))
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 2,(jmax+1)
        do i = 2,(imax+1)
          if(k==1) then
            c1  = 1.+(rv/rd-1)*qts
            c2  = (rv/rd-1)

            wthvcovs(j-1,k) = wthvcovs(j-1,k) + c1 * thlflux(i,j) + c2 * thls * qtflux(i,j)

          else
            qs0h  =  (qt0h(i,j,k) - ql0h(i,j,k))
            t0h   =  exnh(k)*thl0h(i,j,k) + (rlv/cp)*ql0h(i,j,k)
            den   = 1. + (rlv**2)*qs0h/(rv*cp*(t0h**2))
            if (ql0h(i,j,k)>0) then
              c1    = (1.-qt0h(i,j,k)+rv/rd*qs0h * (1.+rd/rv*rlv/(rd*t0h)))/den
              c2    =  c1*rlv/(t0h*cp)-1.
            else
              c1 = 1. + (rv/rd-1)*qt0h(i,j,k)
              c2 = (rv/rd-1)
            end if

            wthvcovs(j-1,k) = wthvcovs(j-1,k) - &
                              c1 * 0.5*(dzf(k-1)*ekh(i,j,k)+dzf(k)*ekh(i,j,k-1))/dzh(k) * &
                                 (thl0(i,j,k) - thl0(i,j,k-1)) / dzh(k) - &
                              c2 * thl0h(i,j,k) * 0.5*(dzf(k-1)*ekh(i,j,k)+dzf(k)*ekh(i,j,k-1))/dzh(k) * &
                                 (qt0(i,j,k) - qt0(i,j,k-1)) / dzh(k)

          endif
        end do
      end do
    end do

    do n=1,nsv
      do k = 2,ncklimit
        do j = 1,jmax
          do i = 1,imax
            !shift prognostic fields one step as 1st column
            !is dummy column because of MPI and periodicity
            wsvcov(j,k,n) = wsvcov(j,k,n) + (w0(i+1,j+1,k)-wavg(j,k)) * (sv0h(i+1,j+1,k,n)-svhavg(j,k,n))
          end do
        end do
      end do

      do k = 1,ncklimit
        do j = 1,jmax
          do i = 1,imax
            !shift prognostic fields one step as 1st column
            !is dummy column because of MPI and periodicity
            if(k==1) then
              wsvcovs(j,k,n) = wsvcovs(j,k,n) + svflux(i+1,j+1,n)
            else
              wsvcovs(j,k,n) = wsvcovs(j,k,n) - 0.5*(dzf(k-1)*ekh(i+1,j+1,k)+dzf(k)*ekh(i+1,j+1,k-1))/dzh(k) &
              * (sv0(i+1,j+1,k,n) - sv0(i+1,j+1,k-1,n)) / dzh(k)
            endif
          end do
        end do
      end do

    enddo

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          thlqcov(j,k) = thlqcov(j,k) + (thl0(i+1,j+1,k)-thlavg(j,k))*(qt0(i+1,j+1,k)-qtavg(j,k))
        end do
      end do
    end do

    vwcov    = vwcov    / itot
    vwcovs   = vwcovs   / itot
    uwcov    = uwcov    / itot
    uwcovs   = uwcovs   / itot
    wthlcov  = wthlcov  / itot
    wthlcovs = wthlcovs / itot
    wthvcov  = wthvcov  / itot
    wthvcovs = wthvcovs / itot
    wqtcov   = wqtcov   / itot
    wqtcovs  = wqtcovs  / itot
    wqlcov   = wqlcov   / itot
    wqlcovs  = wqlcovs  / itot
    wsvcov   = wsvcov   / itot
    wsvcovs  = wsvcovs  / itot
    thlqcov  = thlqcov  / itot

    call D_MPI_REDUCE(uwcov,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, uwcovid, buffer, (/1,1,nccall/), (/jmax, ncklimit, 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    call D_MPI_REDUCE(uwcovs,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, uwcovsid, buffer, (/1,1,nccall/), (/jmax, ncklimit, 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    call D_MPI_REDUCE(vwcov,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, vwcovid, buffer, (/1,1,nccall/), (/jmax, ncklimit, 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    call D_MPI_REDUCE(vwcovs,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, vwcovsid, buffer, (/1,1,nccall/), (/jmax, ncklimit, 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    call D_MPI_REDUCE(wthlcov,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, wthlcovid, buffer, (/1,1,nccall/), (/jmax, ncklimit, 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    call D_MPI_REDUCE(wthlcovs,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, wthlcovsid, buffer, (/1,1,nccall/), (/jmax, ncklimit, 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    call D_MPI_REDUCE(wthvcov,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, wthvcovid, buffer, (/1,1,nccall/), (/jmax, ncklimit, 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    call D_MPI_REDUCE(wthvcovs,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, wthvcovsid, buffer, (/1,1,nccall/), (/jmax, ncklimit, 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    call D_MPI_REDUCE(wqtcov,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, wqtcovid, buffer, (/1,1,nccall/), (/jmax, ncklimit, 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    call D_MPI_REDUCE(wqtcovs,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, wqtcovsid, buffer, (/1,1,nccall/), (/jmax, ncklimit, 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    call D_MPI_REDUCE(wqlcov,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, wqlcovid, buffer, (/1,1,nccall/), (/jmax, ncklimit, 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    call D_MPI_REDUCE(wqlcovs,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, wqlcovsid, buffer, (/1,1,nccall/), (/jmax, ncklimit, 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    call D_MPI_REDUCE(thlqcov,buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
    if(myidx == 0) then
        status = nf90_put_var(ncid, thlqcovid, buffer, (/1,1,nccall/), (/jmax, ncklimit, 1/))
        if(status /= nf90_noerr) call nchandle_error(status)
    endif

    do n=1,nsv
       call D_MPI_REDUCE(wsvcov(:,:,n),buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
       if(myidx == 0) then
          status = nf90_put_var(ncid, wsvcovid(n), buffer, (/1,1,nccall/), (/jmax, ncklimit, 1/))
          if(status /= nf90_noerr) call nchandle_error(status)
       endif
       call D_MPI_REDUCE(wsvcovs(1,:,:),buffer,ncklimit*jmax,MPI_SUM,0,commrow,mpierr)
       if(myidx == 0) then
          status = nf90_put_var(ncid, wsvcovsid(n), buffer, (/1,1,nccall/), (/jmax, ncklimit, 1/))
          if(status /= nf90_noerr) call nchandle_error(status)
       endif
    enddo

  end subroutine do_heterostats

  subroutine exitheterostats

    use typeSizes
    use netcdf
    use modmpi, only : myidx

    implicit none

    integer status

    if(.not.(lheterostats)) return

    if(myidx .ne. 0) return

    deallocate(svavgid)
    deallocate(svvarid)
    deallocate(usvcovid, vsvcovid, wsvcovid)
    deallocate(wsvcovsid)

!    write(6,*) "HETEROSTATS: Closing: ", myid

    status = nf90_close(ncid)
    if (status /= nf90_noerr) call nchandle_error(status)

  end subroutine exitheterostats

  subroutine nchandle_error(status)

    use typeSizes
    use netcdf
    implicit none

    integer, intent(in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if

  end subroutine nchandle_error

end module modheterostats
