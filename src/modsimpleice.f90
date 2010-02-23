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
    use modglobal, only : ih,jh,i1,j1,k1,rdt,rk3step,timee,kmax,rlv,cp,tup,tdn
    use modfields, only : sv0,svm,svp,qtp,thlp,qt0,ql0,exnf,rhobf,tmp0
    use modbulkmicrostat, only : bulkmicrotend
    use modmpi,    only : myid
    implicit none
    integer:: i,j,k 
    real:: qrsmall, qrsum

    delt = rdt/ (4. - dble(rk3step))

    qrsum=0.
    qrsmall=0.
    ! reset microphysics tendencies 
    qrp=0.
    thlpmcr=0.
    qtpmcr=0.
    ! Nr=0, only for statistics module
    Nr= 0.

    do k=1,k1
    do i=2,i1
    do j=2,j1
      ! initialise qr
      qr(i,j,k)= sv0(i,j,k,iqr)
      ! initialise qc only for statistics 
      qc(i,j,k)= ql0(i,j,k)
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
          if(qr(i,j,k)<0.) then
          qr(i,j,k)=0.
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

    do k=1,k1
    do i=2,i1
    do j=2,j1
      ilratio(i,j,k)=amax1(0.,amin1(1.,(tmp0(i,j,k)-tdn)/(tup-tdn)))   ! liquid contribution
      lambdal(i,j,k)=(aal*n0rl*gamb1l/rhobf(k)/(ql0(i,j,k)*ilratio(i,j,k)+1.e-6))**(1./(1.+bbl)) ! lambda rain
      lambdai(i,j,k)=(aai*n0ri*gamb1i/rhobf(k)/(ql0(i,j,k)*(1.-ilratio(i,j,k))+1.e-6))**(1./(1.+bbi)) ! lambda ice
    enddo
    enddo
    enddo

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
    use modglobal, only : ih,jh,i1,j1,k1,rlv,cp, tmelt
    use modfields, only : ql0,exnf,rhobf,tmp0
    use modmpi,    only : myid
    implicit none
    real :: qll,qli,ddisp,del2,autl,tc,times,auti,aut
    integer:: i,j,k

    do k=1,k1
    do i=2,i1
    do j=2,j1
      if (qcmask(i,j,k)==.true.) then
        qll=ql0(i,j,k)*ilratio(i,j,k)
        qli=ql0(i,j,k)-qll
        ddisp=0.146-5.964e-2*alog(Nc_0/2000.) ! relative dispersion for Barry autoconversion
        del2=1.e3*rhobf(k)*qll ! liquid water content
        autl=1./rhobf(k)*1.67e-5*del2*del2/(5. + .0366*Nc_0/(ddisp*(del2+1.E-6)))
        tc=tmp0(i,j,k)-tmelt
        times=amin1(1.e3,(3.56*tc+106.7)*tc+1.e3) ! time scale for ice autoconversion
        auti=qli/times
        aut = autl + auti
        qrp(i,j,k) = qrp(i,j,k)+aut
        qtpmcr(i,j,k) = qtpmcr(i,j,k)-aut
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlv/(cp*exnf(k)))*aut
      endif
    enddo
    enddo
    enddo

  end subroutine autoconvert

  subroutine accrete
    use modglobal, only : ih,jh,i1,j1,k1,rlv,cp,pi
    use modfields, only : ql0,exnf,rhobf
    use modmpi,    only : myid
    implicit none
    real :: qrl,qri,conl,coni,massl,massi,diaml,diami,g_acc_l,g_acc_i,acc_l,acc_i,acc
    integer:: i,j,k

    do k=1,k1
    do i=2,i1
    do j=2,j1
      if (qrmask(i,j,k)==.true. .and. qcmask(i,j,k)==.true.) then
        ! ql partitioning
        qrl=qr(i,j,k)*ilratio(i,j,k)
        qri=qr(i,j,k)-qrl
        conl=n0rl/lambdal(i,j,k) ! liquid concentration
        coni=n0ri/lambdai(i,j,k) ! ice concentration
        massl=rhobf(k)*(qrl+1.e-7) / conl  ! mass liquid particles
        massi=rhobf(k)*(qri+1.e-7) / coni  ! mass ice particles
        diaml=(massl/aal)**(1./bbl) ! diameter liquid particles
        diami=(massi/aal)**(1./bbi) ! diameter ice particles
        g_acc_l=pi/4.*ccl*diaml**(2.+ddl)*ceffl*alphai*rhobf(k)*ql0(i,j,k)
        g_acc_i=pi/4.*cci*diami**(2.+ddi)*ceffi*alphal*rhobf(k)*ql0(i,j,k)
        acc_l=conl*g_acc_l* qrl/(qrl+1.e-9)
        acc_i=coni*g_acc_i* qri/(qri+1.e-9)
        acc= acc_l + acc_i  ! growth by accretion
        qrp(i,j,k) = qrp(i,j,k)+acc
        qtpmcr(i,j,k) = qtpmcr(i,j,k)-acc
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlv/(cp*exnf(k)))*acc
      end if
    enddo
    enddo
    enddo

  end subroutine accrete

  subroutine evaposite
    use modglobal, only : ih,jh,i1,j1,k1,rlv,riv,cp,rv,rd,tmelt,es0,pi
    use modfields, only : qt0,ql0,exnf,rhobf,tmp0,presf
    use modmpi,    only : myid
    implicit none
    real :: qrl,qri,esl,qvsl,ssl,esi,qvsi,ssi, conl,coni,massl,massi,diaml,diami,rel,rei,ventl,venti,thfun,g_devap_l,g_devap_i,devap_l,devap_i,devap
    integer:: i,j,k

    do k=1,k1
    do i=2,i1
    do j=2,j1
      if (qrmask(i,j,k)==.true.) then
        ! qr partitioning
        qrl=qr(i,j,k)*ilratio(i,j,k)
        qri=qr(i,j,k)-qrl
        ! saturation ratio wrt liquid phase
        esl=es0*exp(rlv/rv*(1./tmelt-1./tmp0(i,j,k)))
        qvsl=rd/rv*esl/(presf(k)-(1.-rd/rv)*esl)
        ssl=(qt0(i,j,k)-ql0(i,j,k))/qvsl
        ! saturation ratio wrt ice phase
        esi= es0*exp(riv/rv*(1./tmelt-1./tmp0(i,j,k)))
        qvsi=rd/rv*esi/(presf(k)-(1.-rd/rv)*esi)
        ssi=(qt0(i,j,k)-ql0(i,j,k))/qvsi
        conl=n0rl/lambdal(i,j,k) ! liquid concentration
        coni=n0ri/lambdai(i,j,k) ! ice concentration
        massl=rhobf(k)*(qrl+1.e-7) / conl  ! mass liquid particles
        massi=rhobf(k)*(qri+1.e-7) / coni  ! mass ice particles
        diaml=(massl/aal)**(1./bbl) ! diameter liquid particles
        diami=(massi/aal)**(1./bbi) ! diameter ice particles
        rel=ccl*diaml**(ddl+1.)/2.e-5  ! Reynolds number liquid
        rei=cci*diami**(ddi+1.)/2.e-5  ! Reynolds number ice
        ventl=amax1(1.,.78+.27*sqrt(rel))  ! ventilation factor liquid
        venti=amax1(1.,.65+.39*sqrt(rei))  ! ventilation factor ice
        thfun=1.e-7/(2.2*tmp0(i,j,k)/ssl+2.2e-2/tmp0(i,j,k))  ! thermodynamic fun.
        g_devap_l=4.*pi*diaml/betal*(ssl-1.)*ventl*thfun   ! growth/evap
        g_devap_i=4.*pi*diami/betai*(ssi-1.)*venti*thfun   ! growth/evap
        devap_l=conl * g_devap_l * qrl / (qrl + 1.e-9)
        devap_i=coni * g_devap_i * qri / (qri + 1.e-9)
        devap= devap_l + devap_i  ! growth by deposition
        qrp(i,j,k) = qrp(i,j,k)+devap
        qtpmcr(i,j,k) = qtpmcr(i,j,k)-devap
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlv/(cp*exnf(k)))*devap
      end if
    enddo
    enddo
    enddo

  end subroutine evaposite

  subroutine precipitate
    implicit none

  end subroutine precipitate

end module modsimpleice


