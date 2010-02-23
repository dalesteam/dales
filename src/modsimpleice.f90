!> \file modsimpleice.f90

!>
!!  Ice microphysics.
!>
!! Calculates ice microphysics in a cheap scheme without prognostic Nr
!  simpleice is called from *modmicrophysics*
!! \see  Grabowski, 1998, JAS
!!  \author Steef B\"oing, TU Delft
!!  \par Revision list
!! \todo test
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


module modsimpleice
  use modmicrodata

  implicit none
  save
  contains

!> Initializes and allocates the arrays
  subroutine initsimpleice
    use modglobal, only : ih,i1,jh,j1,k1
    use modmpi,    only : myid
    implicit none

    allocate (qr(2-ih:i1+ih,2-jh:j1+jh,k1)        & ! qr is converted from a scalar variable
             ,qrp(2-ih:i1+ih,2-jh:j1+jh,k1)       & ! qr tendency due to microphysics only, for statistics
             ,thlpmcr(2-ih:i1+ih,2-jh:j1+jh,k1)   & ! thl tendency due to microphysics only, for statistics
             ,qtpmcr(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! qt tendency due to microphysics only, for statistics
             ,sed_qr(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! sedimentation rain droplets mixing ratio
             ,qr_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! time-splitting substep qr
             ,ilratio(2-ih:i1+ih,2-jh:j1+jh,k1)     & ! partition ratio liquid/solid
             ,lambdal(2-ih:i1+ih,2-jh:j1+jh,k1)   & ! slope parameter for liquid phase
             ,lambdai(2-ih:i1+ih,2-jh:j1+jh,k1))    ! slope parameter for ice phase

    allocate (Nr(2-ih:i1+ih,2-jh:j1+jh,k1)        & ! Nrain, only for statistics
             ,qc(2-ih:i1+ih,2-jh:j1+jh,k1))         ! = ql, also for statistics

    allocate (qrmask(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! mask for rain water
             ,qcmask(2-ih:i1+ih,2-jh:j1+jh,k1))     ! mask for cloud water

  end subroutine initsimpleice

!> Cleaning up after the run
  subroutine exitsimpleice
    implicit none
    deallocate(qr,qrp,thlpmcr,qtpmcr,sed_qr,qr_spl,ilratio,lambdal,lambdai,Nr,qc)
    deallocate(qrmask,qcmask) 

  end subroutine exitsimpleice

!> Calculates the microphysical source term.
  subroutine simpleice
    use modglobal, only : ih,jh,i1,j1,k1,dt,rk3step,timee,kmax,rlv,cp
    use modfields, only : sv0,svm,svp,qtp,thlp,qt0,ql0,exnf
    use modbulkmicrostat, only : bulkmicrotend
    use modmpi,    only : myid
    implicit none
    integer:: i,j,k 
    real:: qrsmall, qrsum

    delt = dt/ (4. - dble(rk3step))

    qrsum=0.
    qrsmall=0.
    do k=1,k1
    do i=2,i1
    do j=2,j1
      ! initialise qr and reset microphysics tendencies 
      qr(i,j,k)= sv0(i,j,k,iqr)
      qrp(i,j,k)= 0.0
      thlpmcr(i,j,k)=0.0
      qtpmcr(i,j,k)=0.0
      ! initialise qc and Nr, only for 
      qc(i,j,k)= ql0(i,j,k)
      Nr(i,j,k)= 0
      ! initialise qc mask
      if (ql0(i,j,k) > qcmin) then
        qcmask(i,j,k) = .true.
      else
        qcmask = .false.
      end if
      ! initialise qr mask
      if (l_rain) then
        qrsum = qrsum+ qr(i,j,k)
        if (qr(i,j,k) <= qrmin) then
          qrmask(i,j,k) = .false.
          if(qr(i,j,k)<0) then
          qr(i,j,k)=0
          qrsmall = qrsmall-qr(i,j,k)
          end if
        else
          qrmask(i,j,k)=.true.
        endif
      endif
    enddo
    enddo
    enddo

    if (qrsmall > 0.000001*qrsum) then
      write(*,*)'amount of neg. qr thrown away is too high  ',timee,' sec'
    end if

    if (l_rain) then       
!     call bulkmicrotend   
      call autoconvert
!     call bulkmicrotend   
      call accrete         
!     call bulkmicrotend
      call evaposite       
!     call bulkmicrotend   
      call precipitate
!     call bulkmicrotend
    endif

    do k=1,k1
    do i=2,i1
    do j=2,j1
      svp(i,j,k,iqr)=svp(i,j,k,iqr)+qrp(i,j,k)
      thlp(i,j,k)=thlp(i,j,k)+thlpmcr(i,j,k)
      qtp(i,j,k)=qtp(i,j,k)+qtpmcr(i,j,k)
      ! adjust negative qr tendencies at the end of the time-step
      if (svp(i,j,k,iqr)+svm(i,j,k,iqr)/delt < qrmin) then
        svp(i,j,k,iqr) = - svm(i,j,k,iqr)/delt
        qtp(i,j,k) = qtp(i,j,k) + svm(i,j,k,iqr)/delt
        thlp(i,j,k) = thlp(i,j,k) - (rlv/(cp*exnf(k)))*svm(i,j,k,iqr)/delt
      endif
    enddo
    enddo
    enddo

  end subroutine simpleice

  subroutine autoconvert
    implicit none
  end subroutine autoconvert

  subroutine accrete
    implicit none
  end subroutine accrete

  subroutine precipitate
    implicit none
  end subroutine precipitate

  subroutine evaposite
    implicit none
  end subroutine evaposite

end module modsimpleice


