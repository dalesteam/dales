!> \file modbulkmicro.f90

!>
!!  Bulk microphysics.
!>
!! Calculates bulk microphysics using a two moment scheme.
!! \see   Seifert and Beheng (Atm. Res., 2001)
!! \see  Seifert and Beheng (Met Atm Phys, 2006)
!! \see  Stevens and Seifert (J. Meteorol. Soc. Japan, 2008)  (rain sedim, mur param)
!! \see  Seifert (J. Atm Sc., 2008) (rain evap)
!! \see  Khairoutdinov and Kogan (2000) (drizzle param : auto, accr, sedim, evap)
!!  \author Olivier Geoffroy, K.N.M.I.
!!  \author Margreet van Zanten, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \par Revision list
!! \todo documentation
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


module modbulkmicro

!*********************************************************************
!
!   Amount of liquid water is splitted into cloud water and precipitable
!   water (as such it is a two moment scheme). Cloud droplet number conc. is
!   fixed in place and time.
!
!   same rhof value used for some diagnostics calculation (in modbulkmicrostat, modtimestat)
!
!   Cond. sampled timeav averaged profiles are weighted with fraction of condition,
!   similarly as is done in sampling.f90
!
!   bulkmicro is called from *modmicrophysics*
!
!*********************************************************************
  use modmicrodata

  implicit none
  real :: gamma25
  real :: gamma3
  real :: gamma35
  real :: wMax
  contains

!> Initializes and allocates the arrays
  subroutine initbulkmicro
    use modglobal, only : ih,i1,jh,j1,k1,dzf
    use modmpi,    only : myid
    implicit none

    allocate( Nr       (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,Nrp      (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,qr       (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,qrp      (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,Nc       (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,sedc     (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,Dvr      (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,xr       (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,mur      (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,lbdr     (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,tau      (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,evap     (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,Nevap    (2-ih:i1+ih,2-jh:j1+jh,k1)  )

    allocate(precep    (2-ih:i1+ih,2-jh:j1+jh,k1)  )
    allocate(qrmask    (2-ih:i1+ih,2-jh:j1+jh,k1)  &
            ,qcmask    (2-ih:i1+ih,2-jh:j1+jh,k1)  )

    allocate(sed_qr    (2-ih:i1+ih,2-jh:j1+jh,k1)  &
            ,sed_Nr    (2-ih:i1+ih,2-jh:j1+jh,k1)  &
            ,qr_spl    (2-ih:i1+ih,2-jh:j1+jh,k1)  &
            ,Nr_spl    (2-ih:i1+ih,2-jh:j1+jh,k1)  )

    gamma25=lacz_gamma(2.5)
    gamma3=2.
    gamma35=lacz_gamma(3.5)

    if (ibulk==0) then ! ibulk has not been set
      print *, "MODBULKMICRO: The use of l_sb is deprecated. Instead use ibulk={1,2,3} to select bulk micro model"
      if (l_sb) then
        print *, "MODBULKMICRO: set ibulk=1 (Seifert and Beheng 2001)"
        ibulk = ibulk_sb01
      else ! l_sb = .false.
        print *, "MODBULKMICRO: set ibulk=3 (Kogan 2013)"
        ibulk = ibulk_k13
      end if
    end if

    select case (ibulk)
      case (ibulk_sb01)
        wMax = 9.9 ! Please check this value, it was set for SB01 in the original code!
      case (ibulk_kk00)
        wMax = wMaxkk ! Approx. 5.8 m/s, corresponding to xrmaxkk in modmicrodata
      case (ibulk_k13)
        wMax = 7.     ! Based on Kogan (2013) Fig. 4
      case default
        print *,"MODBULKMICRO: No valid bulkmicrophysics scheme selected! (ibulk)"
        stop "MODBULKMICRO: No valid bulkmicrophysics scheme selected! (ibulk)"
    end select

  end subroutine initbulkmicro

!> Cleaning up after the run
  subroutine exitbulkmicro
    implicit none

    deallocate(Nr,Nrp,qr,qrp,Nc)

    deallocate(sedc,Dvr,xr,mur,lbdr,tau,evap,Nevap)
    deallocate(sed_qr,sed_Nr,qr_spl,Nr_spl)

    deallocate(precep)

  end subroutine exitbulkmicro

!> Calculates the microphysical source term.
  subroutine bulkmicro
    use modglobal, only : ih,jh,i1,j1,k1,rdt,iTimeInt,iTimeWicker,iTimeTVD,iTimeLowStor,rkStep,timee,kmax,rlv,cp,dzf
    use modfields, only : sv0,svm,svp,qtp,thlp,qt0,ql0,presf,exnf,rhof
    use modbulkmicrostat, only : bulkmicrotend
    use modtstep,  only : rkb
    use modmpi,    only : myid
    implicit none
    integer :: i,j,k
    real,dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: qrtest,nrtest

    ! Initialize variables
    Nrp=0.; qrp=0.; qr=0.; Nr=0.; qrtest=0.; nrtest=0.
    Nc=0.; qrmask=.false.; qcmask=.false.

    ! Set the temporary number concentrations and rain water specific humidities
    do k=1,k1; do j=2,j1; do i=2,i1
      Nr(i,j,k) = sv0(i,j,k,inr)
      qr(i,j,k) = sv0(i,j,k,iqr)
    end do; end do; end do

    if (l_rain) then
      if (sum(qr, qr<0.) > 0.000001*sum(qr)) then
        write(*,*)'amount of neg. qr and Nr is too high  ',timee,' sec'
      end if
      if (sum(Nr, Nr<0.) > 0.000001*sum(Nr)) then
         write(*,*)'amount of neg. qr and Nr is too high  ',timee,' sec'
      end if
    end if   ! l_rain

    do k=1,k1; do j=2,j1; do i=2,i1
      ! determine the rain mask
      if (qr(i,j,k) > qrmin .and. Nr(i,j,k)>0.) then
        qrmask(i,j,k) = .true.
      end if
      ! initialize cloud droplet number Nc and determine the mask
      if (ql0(i,j,k) > qcmin) then
        qcmask(i,j,k) = .true.
        Nc(i,j,k) = Nc_0
      end if
    end do; end do; end do

    ! calculate Rain DSD integral properties & parameters xr, Dvr, lbdr, mur
    if (l_rain) then

      xr   (2:i1,2:j1,1:k1) = 0.
      Dvr  (2:i1,2:j1,1:k1) = 0.
      mur  (2:i1,2:j1,1:k1) = 30.
      lbdr (2:i1,2:j1,1:k1) = 0.

      select case (ibulk)
        case (ibulk_sb01)
        !=== Seifert and Beheng (2001)
        
        do k=1,k1; do j=2,j1; do i=2,i1
          if (qrmask(i,j,k)) then
            xr (i,j,k) = rhof(k)*qr(i,j,k)/Nr(i,j,k) 
            xr (i,j,k) = min(max(xr(i,j,k),xrmin),xrmax) ! to ensure xr is within borders
            Dvr(i,j,k) = (xr(i,j,k)/pirhow)**(1./3.)
          endif
        end do; end do; end do

        if (l_mur_cst) then
        ! mur = cst
          do k=1,k1; do j=2,j1; do i=2,i1
            if (qrmask(i,j,k)) then
              mur(i,j,k) = mur_cst
              lbdr(i,j,k) = ((mur(i,j,k)+3.)*(mur(i,j,k)+2.)*(mur(i,j,k)+1.))**(1./3.)/Dvr(i,j,k)
            endif
          end do; end do; end do
        else
        ! mur = f(Dv)
          do k=1,k1; do j=2,j1; do i=2,i1
            if (qrmask(i,j,k)) then
!
!             mur(2:i1,2:j1,1:k1) = 10. * (1+tanh(1200.*(Dvr(2:i1,2:j1,1:k1)-0.0014))) 
!             Stevens & Seifert (2008) param
!
              mur(i,j,k) = min(30.,- 1. + 0.008/ (qr(i,j,k)*rhof(k))**0.6)  ! G09b
              lbdr(i,j,k) = ((mur(i,j,k)+3.)*(mur(i,j,k)+2.)*(mur(i,j,k)+1.))**(1./3.)/Dvr(i,j,k)
            endif
          end do; end do; end do

        endif

      case (ibulk_kk00)
        !=== Khairoutdinov and Kogan (2000)
        do k=1,k1; do j=2,j1; do i=2,i1
          if (qrmask(i,j,k)) then
            xr  (i,j,k) = rhof(k)*qr(i,j,k)/Nr(i,j,k) ! average mass of a droplet
            xr  (i,j,k) = min(xr(i,j,k),xrmaxkk) ! to ensure x_pw is within borders
            Dvr (i,j,k) = (xr(i,j,k)/pirhow)**(1./3.)
          endif
        end do; end do; end do

      case (ibulk_k13)
        !=== Kogan (2013)
        do k=1,k1; do j=2,j1; do i=2,i1
          if (qrmask(i,j,k)) then
            xr (i,j,k) = rhof(k)*qr(i,j,k)/Nr(i,j,k) ! average mass of a droplet
            Dvr(i,j,k) = (xr(i,j,k)/pirhow)**(1./3.)
            Dvr(i,j,k) = min(Dvr(i,j,k),0.3e-3) ! Limit Dvr to .3mm. Not sure what it is based on.
                                                ! Should check influence later.
          endif
        end do; end do; end do

      end select
    endif   ! l_rain

  !*********************************************************************
  ! call microphysical processes subroutines
  !*********************************************************************
    if (l_sedc)  call sedimentation_cloud

    if (l_rain) then
      call bulkmicrotend
      call autoconversion
      call bulkmicrotend
      call accretion
      call bulkmicrotend
      call evaporation
      call bulkmicrotend
      call sedimentation_rain
      call bulkmicrotend
    endif

    !*********************************************************************
    ! remove negative values and non physical low values
    !*********************************************************************
    ! Add tendency due to cloud droplet sedimentation (done here for bulkmicrotend)
    do k=1,kmax; do j=2,j1; do i=2,i1
      qtp(i,j,k)  = qtp(i,j,k)  + (sedc(i,j,k+1)-sedc(i,j,k))/(dzf(k)*rhof(k))
      thlp(i,j,k) = thlp(i,j,k) - (rlv/(cp*exnf(k))) &
                       *(sedc(i,j,k+1)-sedc(i,j,k))/(dzf(k)*rhof(k))
    end do; end do; end do

    ! Evaporate all qr (and Nr of course) that is smaller than qrmin
    select case (iTimeInt)
      case (iTimeWicker,iTimeTVD)
        ! Calculate a test qr and Nr
        qrtest = svm(:,:,:,iqr)+(svp(:,:,:,iqr)+qrp)*rdt*rkb(rkStep)
        Nrtest = svm(:,:,:,iNr)+(svp(:,:,:,iNr)+Nrp)*rdt*rkb(rkStep)
        where ( qrtest<qrmin               .or. &
                Nrtest<0.                  .or. &
               (Nrtest==0..and.qrtest/=0.) .or. &
               (qrtest==0..and.Nrtest/=0.))
          qrp = -svp(:,:,:,iqr)-svm(:,:,:,iqr)/(rdt*rkb(rkStep))
          Nrp = -svp(:,:,:,iNr)-svm(:,:,:,iNr)/(rdt*rkb(rkStep))
        end where
      case (iTimeLowStor)
        qrtest = sv0(:,:,:,iqr)+(svp(:,:,:,iqr)+qrp)*rdt*rkb(rkStep)
        Nrtest = sv0(:,:,:,iNr)+(svp(:,:,:,iNr)+Nrp)*rdt*rkb(rkStep)

        where ( qrtest<qrmin               .or. &
                Nrtest<0.                  .or. &
               (Nrtest==0..and.qrtest/=0.) .or. &
               (qrtest==0..and.Nrtest/=0.))
          qrp = -svp(:,:,:,iqr)-sv0(:,:,:,iqr)/(rdt*rkb(rkStep))
          Nrp = -svp(:,:,:,iNr)-sv0(:,:,:,iNr)/(rdt*rkb(rkStep))
        end where
      case default
        write(*,*) 'No valid time integration method selected'
    end select

    svp(2:i1,2:j1,1:kmax,iqr)= svp(2:i1,2:j1,1:kmax,iqr) + qrp(2:i1,2:j1,1:kmax)
    svp(2:i1,2:j1,1:kmax,inr)= svp(2:i1,2:j1,1:kmax,inr) + Nrp(2:i1,2:j1,1:kmax)
    qtp(2:i1,2:j1,1:kmax)    = qtp(2:i1,2:j1,1:kmax) - qrp(2:i1,2:j1,1:kmax)
    do j=2,j1; do i=2,i1
      thlp(i,j,1:kmax)   = thlp(i,j,1:kmax) + (rlv/(cp*exnf(1:kmax)))*qrp(i,j,1:kmax)
    end do; end do

  end subroutine bulkmicro

  !> Determine autoconversion rate and adjust qrp and Nrp accordingly
  !!
  !!   The autoconversion rate is formulated for f(x)=A*x**(nuc)*exp(-Bx),
  !!   decaying exponentially for droplet mass x.
  !!   It can easily be reformulated for f(x)=A*x**(nuc)*exp(-Bx**(mu)) and
  !!   by chosing mu=1/3 one would get a gamma distribution in drop diameter
  !!   -> faster rain formation. (Seifert)
  subroutine autoconversion
    use modglobal, only : ih,i1,jh,j1,k1,kmax,eps1,rlv,cp
    use modmpi,    only : myid
    use modfields, only : exnf,rhof,ql0
    implicit none
    integer i,j,k
    real :: xc,nuc,au,phi,k_au

    au=0.

    select case (ibulk)
      case (ibulk_sb01)
        !=== SB autoconversion
        tau(2:i1,2:j1,1:k1) = 0.
        !phi(2:i1,2:j1,1:k1) = 0.
        k_au = k_c/(20*x_s)

        do k=1,k1; do j=2,j1; do i=2,i1
          if (qcmask(i,j,k)) then 
            nuc            = 1.58*(rhof(k)*ql0(i,j,k)*1000.) +0.72-1. !G09a
            xc             = rhof(k) * ql0(i,j,k) / Nc(i,j,k) ! No eps0 necessary
            au             = k_au * (nuc       +2.) * (nuc       +4.) / (nuc       +1.)**2.    &
                    * (ql0(i,j,k) * xc)**2. * 1.225 ! *rho**2/rho/rho (= 1)
            tau    (i,j,k) = 1.0 - ql0(i,j,k) / (ql0(i,j,k)+qr(i,j,k))
            phi            = k_1 * tau(i,j,k)**k_2 * (1.0 -tau(i,j,k)**k_2)**3
            au             = au        * (1.0 + phi       /(1.0 -tau(i,j,k))**2)

            qrp    (i,j,k) = qrp    (i,j,k) + au        
            Nrp    (i,j,k) = Nrp    (i,j,k) + au /x_s
          endif
        end do; end do; end do

      case (ibulk_kk00)
        !=== KK00 autoconversion
        do k=1,k1; do j=2,j1; do i=2,i1
          if (qcmask(i,j,k)) then
            ! au      = 1e7 * ql0(i,j,k)**(3) * (Nc(i,j,k)/1.0E6)**(-1.79)  ! ECMWF-like
            ! au      = .001*(ql0(i,j,k)-qlKess)  ! Kessler-type (SAM, Khairoutdinov and Randall 2003)
            au = 1350.0 * ql0(i,j,k)**(2.47) * (Nc(i,j,k)/1.e6)**(-1.79)
  
            qrp(i,j,k) = qrp(i,j,k) + au       
            Nrp(i,j,k) = Nrp(i,j,k) + au*rhof(k)/(pirhow*D0_kk**3.)   ! D0 is the separation drop diameter (=50 um)
          endif
        end do; end do; end do

      case (ibulk_k13)
        !=== Kogan (2013) autoconversion
        do k=1,k1; do j=2,j1; do i=2,i1
          if (qcmask(i,j,k)) then
            au = 7.98e10 * ql0(i,j,k)**(4.22) * (Nc(i,j,k)/1.e6)**(-3.01)

            ! Add the tendencies
            qrp(i,j,k) = qrp(i,j,k) + au
            Nrp(i,j,k) = Nrp(i,j,k) + au*rhof(k)/(pirhow*D0_k13**3.)  ! D0 is the separation drop diameter (=80 um)
            ! NOTE: In the article, it seems that a factor 3 is missing in Eq. (27) (3*rho_a). I did include it here.
          endif
        end do; end do; end do
      
      end select

  end subroutine autoconversion

  subroutine accretion
  !*********************************************************************
  ! determine accr. + self coll. + br-up rate and adjust qrp and Nrp
  ! accordingly. Break-up : Seifert (2007)
  !*********************************************************************
    use modglobal, only : ih,i1,jh,j1,k1,kmax,eps1,rlv,cp,dzf
    use modfields, only : exnf,rhof,ql0
    use modmpi,    only : myid
    implicit none
!    real , allocatable :: phi_br(:,:,:)
    real :: ac,sc,br,phi,phi_br
    integer :: i,j,k
!   allocate (phi_br(2-ih:i1+ih,2-jh:j1+jh,k1))

    select case (ibulk)
      case (ibulk_sb01)
        !=== SB accretion
        do j=2,j1
        do i=2,i1
        do k=1,k1
           if (qrmask(i,j,k) .and. qcmask(i,j,k)) then
              tau    (i,j,k) = 1.0 - ql0(i,j,k)/(ql0(i,j,k)+qr(i,j,k))
              phi            = (tau(i,j,k)/(tau(i,j,k) + k_l))**4.
              ac             = k_r *rhof(k)*ql0(i,j,k) * qr(i,j,k) * phi * &
                               (1.225/rhof(k))**0.5
              qrp    (i,j,k) = qrp    (i,j,k) + ac       
           endif
        enddo
        enddo
        enddo

        !=== SB self-collection & Break-up
        sc=0.
        br=0.
   
        do j=2,j1
        do i=2,i1
        do k=1,k1
           if (qrmask(i,j,k)) then
              sc = k_rr *rhof(k)* qr(i,j,k) * Nr(i,j,k)  &
                          * (1 + kappa_r/lbdr(i,j,k)*pirhow**(1./3.))**(-9.)* (1.225/rhof(k))**0.5
           endif
           if (Dvr(i,j,k) .gt. 0.30E-3 .and. qrmask(i,j,k)) then
              phi_br = k_br * (Dvr(i,j,k)-D_eq)
              br     = (phi_br + 1.) * sc
           else
              br     = 0. ! (phi_br = -1)
           endif
   
           Nrp(i,j,k) = Nrp(i,j,k) - sc + br
   
        enddo
        enddo
        enddo

      case (ibulk_kk00)
        !=== KK00 ccretion
        do k=1,k1; do j=2,j1; do i=2,i1
          if (qrmask(i,j,k) .and. qcmask(i,j,k)) then
            ac = 67.0 * ( ql0(i,j,k)*qr(i,j,k) )**(1.15)
            ! Add the tendencies
            qrp(i,j,k) = qrp(i,j,k) + ac       
          endif
        end do; end do; end do

      case (ibulk_k13)
        !=== Kogan (2013) accretion and self-collection
        do k=1,k1; do j=2,j1; do i=2,i1
          ! accretion
          if (qrmask(i,j,k) .and. qcmask(i,j,k)) then
            ac = 8.53 * ql0(i,j,k)**(1.05) * qr(i,j,k)**(0.98)
            ! Add the tendencies
            qrp(i,j,k) = qrp(i,j,k) + ac       
          endif
          ! self-collection
          if (qrmask(i,j,k)) then
            sc = 205. * qr(i,j,k)**(1.55) * Nr(i,j,k)**(0.60)
            Nrp(i,j,k) = Nrp(i,j,k) + sc
          end if
        end do; end do; end do
      
    end select


  end subroutine accretion

!> Sedimentation of cloud water
!!
!!   The sedimentation of cloud droplets assumes a lognormal DSD in which the
!!   geometric std dev. is assumed to be fixed at 1.3.
!! sedimentation of cloud droplets
!! lognormal CDSD is assumed (1 free parameter : sig_g)
!! terminal velocity : Stokes velocity is assumed (v(D) ~ D^2)
!! flux is calc. anal.
  subroutine sedimentation_cloud
    use modglobal, only : i1,j1,k1,kmax,eps1,rlv,cp,dzf,pi,rhow
    use modmpi,    only : myid
    use modfields, only : rhof,exnf,ql0
    implicit none
    integer :: i,j,k

    sedc = 0.
    csed = c_St*(3./(4.*pi*rhow))**(2./3.)*exp(5.*log(sig_g)**2)

    do k=1,k1; do j=2,j1; do i=2,i1
      if (qcmask(i,j,k)) then
        sedc(i,j,k) = csed*(Nc(i,j,k))**(-2./3.)*(ql0(i,j,k)*rhof(k))**(5./3.)
      endif
    end do; end do; end do

  end subroutine sedimentation_cloud

!> Sedimentaion of rain
!! sedimentation of drizzle water
!! - gen. gamma distr is assumed. Terminal velocities param according to
!!   Stevens & Seifert. Flux are calc. anal.
!! - l_lognormal =T : lognormal DSD is assumed with D_g and N known and
!!   sig_g assumed. Flux are calc. numerically with help of a
!!   polynomial function
  subroutine sedimentation_rain
    use modglobal, only : ih,i1,jh,j1,k1,kmax,eps1,dzf,pi,rdt
    use modfields, only : rhof
    use modmpi,    only : myid,mpi_max,mpi_integer,mpierr,comm3d
    implicit none
    integer :: i,j,k,jn
    integer :: n_spl      !<  sedimentation time splitting loop
    real    :: pwcont
    real, allocatable :: xr_spl(:,:,:),Dvr_spl(:,:,:),&
                        mur_spl(:,:,:),lbdr_spl(:,:,:),Dgr(:,:,:)
    real,save :: dt_spl
    real,parameter :: Cstab=0.7 ! Presumably stable CFL number
    real :: vqr,vNr ! Fall velocities of qr and Nr respectively

    allocate( xr_spl(2-ih:i1+ih,2-jh:j1+jh,k1)      &!<  for time splitting
              ,Dvr_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    &!<     -
              ,mur_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    &!<     -
              ,lbdr_spl(2-ih:i1+ih,2-jh:j1+jh,k1)   &!<     -
              ,Dgr(2-ih:i1+ih,2-jh:j1+jh,k1))        !<  lognormal geometric diameter

    qr_spl(2:i1,2:j1,1:k1) = qr(2:i1,2:j1,1:k1)
    Nr_spl(2:i1,2:j1,1:k1) = Nr(2:i1,2:j1,1:k1)

    !n_spl = ceiling(wfallmax*subDt/(minval(dzf)))
    n_spl = ceiling(wMax*rdt/(Cstab*minval(dzf))) ! JvdD changed the Runge-Kutta subtimestep to the full timestep
                                            ! Also changed the presumed stable CFL criterion from 1. to 0.7
    dt_spl = rdt/real(n_spl)

    ! Reset precep
    precep = 0.

    do jn = 1 , n_spl ! time splitting loop

      sed_qr(2:i1,2:j1,1:k1) = 0.
      sed_Nr(2:i1,2:j1,1:k1) = 0.

      select case (ibulk)
        case (ibulk_sb01)
          !=== SB sedimentation

          do k=1,k1; do j=2,j1; do i=2,i1
            if (qr_spl(i,j,k)>qrmin .and. Nr_spl(i,j,k)>0.) then
              xr_spl (i,j,k) = rhof(k)*qr_spl(i,j,k)/Nr_spl(i,j,k)
              xr_spl (i,j,k) = min(max(xr_spl(i,j,k),xrmin),xrmax) ! to ensure xr is within borders
              Dvr_spl(i,j,k) = (xr_spl(i,j,k)/pirhow)**(1./3.)
            endif
          end do; end do; end do
  
  
          if (l_lognormal) then
            do k = 1,kmax
            do j = 2,j1
            do i = 2,i1
              if (qr_spl(i,j,k)>qrmin .and. Nr_spl(i,j,k)>0.) then
                Dgr(i,j,k) = (exp(4.5*(log(sig_gr))**2))**(-1./3.)*Dvr_spl(i,j,k) ! correction for width of DSD
                sed_qr(i,j,k) = -1.*sed_flux(Nr_spl(i,j,k),Dgr(i,j,k),log(sig_gr)**2,D_s,3)
                sed_Nr(i,j,k) = -1./pirhow*sed_flux(Nr_spl(i,j,k),Dgr(i,j,k) ,log(sig_gr)**2,D_s,0)
    !        correction for the fact that pwcont .ne. qr_spl
    !        actually in this way for every grid box a fall velocity is determined
                pwcont = liq_cont(Nr_spl(i,j,k),Dgr(i,j,k),log(sig_gr)**2,D_s,3)         ! note : kg m-3
                if (pwcont > eps1) then
                  sed_qr(i,j,k) = (qr_spl(i,j,k)*rhof(k)/pwcont)*sed_qr(i,j,k)  ! or qr_spl*(sed_qr/pwcont) = qr_spl*fallvel.
                end if
              end if ! qr_spl threshold statement
            end do
            end do
            end do
    
          else
            if (l_mur_cst) then
              mur_spl(2:i1,2:j1,1:k1) = mur_cst
            else
              do j=2,j1
              do i=2,i1
              do k=1,k1
                if (qr_spl(i,j,k)>qrmin .and. Nr_spl(i,j,k)>0.) then
    !             mur_spl(i,j,k) = 10. * (1+tanh(1200.*(Dvr_spl(i,j,k)-0.0014))) ! SS08
                  mur_spl(i,j,k) = min(30.,- 1. + 0.008/ (qr_spl(i,j,k)*rhof(k))**0.6)  ! G09b
                endif
              enddo
              enddo
              enddo
    
            endif
    
            do j=2,j1
            do i=2,i1
            do k=1,k1
              if (qr_spl(i,j,k)>qrmin .and. Nr_spl(i,j,k)>0.) then
                  lbdr_spl(i,j,k) = ((mur_spl(i,j,k)+3.)*(mur_spl(i,j,k)+2.)* &
                                     (mur_spl(i,j,k)+1.))**(1./3.)/Dvr_spl(i,j,k)
                  vqr          = min(wMax,max(1.e-3,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr_spl(i,j,k))**(-1.*(mur_spl(i,j,k)+4.)))))
                  vNr          = min(wMax,max(1.e-3,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr_spl(i,j,k))**(-1.*(mur_spl(i,j,k)+1.)))))
                  sed_qr  (i,j,k) = -vqr*qr_spl(i,j,k)*rhof(k)
                  sed_Nr  (i,j,k) = -vNr*Nr_spl(i,j,k)
              endif
            enddo
            enddo
            enddo
    
          endif !l_lognormal

        case (ibulk_kk00)
          !=== KK00 sedimentation
          do k=1,k1; do j=2,j1; do i=2,i1
            if (qr_spl(i,j,k)>qrmin .and. Nr_spl(i,j,k)>0.) then
              xr_spl(i,j,k) = rhof(k)*qr_spl(i,j,k)/(Nr_spl(i,j,k)+eps0)
              xr_spl(i,j,k) = min(xr_spl(i,j,k),xrmaxkk) ! to ensure xr is within borders 
              Dvr_spl(i,j,k) = (xr_spl(i,j,k)/pirhow)**(1./3.)
              vqr        = 0.006*1.0E6*Dvr_spl(i,j,k)- 0.2  ! Original fall speed DALES
              vqr        = min(max(vqr,1e-3),wMax)          ! Limiter on fall speed based on figure 6 in KK00 article
              ! In this setting, this limiter seems to be extremely important;
              ! simulations crash without it!
              ! Sensitivity to value of lower limit seems to be negligible (range 1-1e-3)
              vNr        = 0.0035*1.0E6*Dvr_spl(i,j,k)- 0.1 ! Original fall speed DALES
              vNr        = min(max(vNr,1e-3),wMax)          ! Limiter from UCLA microphys
              sed_Nr(i,j,k) = -vNr*Nr_spl(i,j,k)*rhof(k)
              sed_qr(i,j,k) = -vqr*qr_spl(i,j,k)*rhof(k)
            end if
          end do; end do; end do

        case (ibulk_k13)
          !=== Kogan (2013) sedimentation
          do k=1,k1; do j=2,j1; do i=2,i1
            if (qr_spl(i,j,k)>qrmin .and. Nr_spl(i,j,k)>0.) then
              xr_spl(i,j,k) = rhof(k)*qr_spl(i,j,k)/Nr_spl(i,j,k)
              Dvr_spl(i,j,k) = (xr_spl(i,j,k)/pirhow)**(1./3.)
              Dvr_spl(i,j,k) = max(0.,min(Dvr_spl(i,j,k),0.3e-3)) ! Limit Dvr to 0-3mm. Not sure what it is based on.
              ! NOTE: the article uses the mean volume radius in um, hence the
              ! factor .5e6. Velocities are in cm/s, hence the factor 1/100.
              vqr = 2.4  *(.5e6*Dvr_spl(i,j,k)) - 62.
              vqr = vqr/100.
              vqr        = min(max(vqr,1e-3),wMax)          ! Limiter on fall speed based on figure 6 in KK00 article
              ! These velocities are defined <0 downwards!!

              vNr = 0.385*(.5e6*Dvr_spl(i,j,k)) + 5.76
              vNr = vNr/100.
              vNr = min(max(vNr,1e-3),wMax)
            
              sed_Nr(i,j,k) = -vNr*Nr_spl(i,j,k)*rhof(k)
              sed_qr(i,j,k) = -vqr*qr_spl(i,j,k)*rhof(k)
            end if
          end do;end do; end do
          
      end select

      ! Determine the new Nr and qr within after a split time-step
      
      do k=1,kmax; do j=2,j1; do i=2,i1
        Nr_spl(i,j,k) = Nr_spl(i,j,k) - &
                (sed_Nr(i,j,k+1) - sed_Nr(i,j,k))*dt_spl/dzf(k)
        qr_spl(i,j,k) = qr_spl(i,j,k) - &
                (sed_qr(i,j,k+1) - sed_qr(i,j,k))*dt_spl/(dzf(k)*rhof(k))

        precep(i,j,k) =  precep(i,j,k) - sed_qr(i,j,k)/rhof(k)   ! kg kg-1 m s-1
      end do; end do; end do

    ! end of time splitting loop
    end do 

    ! Average the precipitation rate over all subtimesteps for consistency
    precep = precep/n_spl

!   JvdD replaced by rdt (=total timestep, not the RK subtimestep) tendency
    Nrp(2:i1,2:j1,1:k1)= Nrp(2:i1,2:j1,1:k1) + (Nr_spl(2:i1,2:j1,1:k1) - Nr(2:i1,2:j1,1:k1))/rdt
    qrp(2:i1,2:j1,1:k1)= qrp(2:i1,2:j1,1:k1) + (qr_spl(2:i1,2:j1,1:k1) - qr(2:i1,2:j1,1:k1))/rdt

    deallocate (xr_spl,Dvr_spl,mur_spl,lbdr_spl,Dgr) 
  end subroutine sedimentation_rain

  !*********************************************************************
  !*********************************************************************

  subroutine evaporation
  !*********************************************************************
  ! Evaporation of prec. : Seifert (2008)
  ! Cond. (S>0.) neglected (all water is condensed on cloud droplets)
  !*********************************************************************

    use modglobal, only : ih,i1,jh,j1,k1,kmax,eps1,es0,rd,rv,tmelt,rlv,cp,at,bt,pi,ep,mygamma251,mygamma21,lacz_gamma,rdt,rhow
    use modfields, only : exnf,thl0,qt0,sv0,qvsl,tmp0,ql0,esl,qvsl,rhof,exnf
    use modmpi,    only : myid
    implicit none
    integer :: i,j,k
    real, allocatable :: S(:,:,:),G(:,:,:)
    real :: F
    real :: Cr
    integer :: numel

    allocate( S(2-ih:i1+ih,2-jh:j1+jh,k1), & ! super or undersaturation
              G(2-ih:i1+ih,2-jh:j1+jh,k1)  & ! cond/evap rate of a drop
             )

    evap=0.; Nevap=0.

    do k=1,k1; do j=2,j1; do i=2,i1
      if (qrmask(i,j,k)) then
        S (i,j,k) = min(0.,(qt0(i,j,k)-ql0(i,j,k))/qvsl(i,j,k)- 1.)
        G (i,j,k) = (Rv * tmp0(i,j,k)) / (Dv*esl(i,j,k)) + rlv/(Kt*tmp0(i,j,k))*(rlv/(Rv*tmp0(i,j,k)) -1.)
        G (i,j,k) = G(i,j,k)**(-1) ! NOTE: [G]=kg/m-1 s-1, KK00 assumes units of m^2/s
        ! G is the same as coefficient A3 (Eq. 13-32) in Pullacher and Klett (1978) and is given by the denominator in Eq. (13-28)
        ! This term differs by a factor rho_w from the equation used here! This is corrected later
      end if
    end do; end do; end do
                
    select case (ibulk)
      case (ibulk_sb01)
        !=== SB evaporation
        do k=1,k1; do j=2,j1; do i=2,i1
          if (qrmask(i,j,k)) then
            numel=nint(mur(i,j,k)*100.)
            F = avf * mygamma21(numel)*Dvr(i,j,k) +  &
               bvf*Sc_num**(1./3.)*(a_tvsb/nu_a)**0.5*mygamma251(numel)*Dvr(i,j,k)**(3./2.) * &
               (1.-(1./2.)  *(b_tvsb/a_tvsb)    *(lbdr(i,j,k)/(   c_tvsb+lbdr(i,j,k)))**(mur(i,j,k)+2.5)  &
                  -(1./8.)  *(b_tvsb/a_tvsb)**2.*(lbdr(i,j,k)/(2.*c_tvsb+lbdr(i,j,k)))**(mur(i,j,k)+2.5)  &
                  -(1./16.) *(b_tvsb/a_tvsb)**3.*(lbdr(i,j,k)/(3.*c_tvsb+lbdr(i,j,k)))**(mur(i,j,k)+2.5) &
                  -(5./128.)*(b_tvsb/a_tvsb)**4.*(lbdr(i,j,k)/(4.*c_tvsb+lbdr(i,j,k)))**(mur(i,j,k)+2.5)  )
  ! *lbd(i,j,k)**(mur(i,j,k)+1.)/f_gamma_1(i,j,k) factor moved to F
             evap(i,j,k) = 2*pi*Nr(i,j,k)*G(i,j,k)*F*S(i,j,k)/rhof(k) ! This equation is similar to Eq. (21) of Seifert (2008) ! but not identical ! Where does it come from?
             Nevap(i,j,k) = c_Nevap*evap(i,j,k)*rhof(k)/xr(i,j,k)
          endif
        end do; end do; end do

      case (ibulk_kk00)
        !=== KK00 evaporation
        do k=1,k1; do j=2,j1; do i=2,i1
          if (qrmask(i,j,k)) then
            evap(i,j,k)  = c_evapkk*2*pi*Dvr(i,j,k)*G(i,j,k)*S(i,j,k)*Nr(i,j,k)/rhof(k)
            Nevap(i,j,k) = evap(i,j,k)*rhof(k)/xr(i,j,k)
          end if
        end do; end do; end do

      case (ibulk_k13)
        !=== Kogan (2013) evaporation
        do k=1,k1; do j=2,j1; do i=2,i1
          if (qrmask(i,j,k)) then
            Cr = .45 + 23./(.5e6*Dvr(i,j,k))
            evap(i,j,k) = 3.*Cr*G(i,j,k)/rhow*S(i,j,k)       & ! Note the division of G by rho_w in accord with Pruppacher and Klett (1978)
                         *(4.*pi*rhow/(3.*rhof(k)))**(2./3.) &
                         *qr(i,j,k)**(1./3.)                 &
                         *Nr(i,j,k)**(2./3.)
            ! NOTE: condensation is neglected (S>0.). For condensation, Nevap would be 0. 
            ! For evaporation, Nevap should scale as follows:
            ! (See Kogan, 2013 Eq. 19)
            Nevap(i,j,k) = evap(i,j,k)/qr(i,j,k)*Nr(i,j,k)
          end if
        end do; end do; end do
    end select

    qrp(2:i1,2:j1,1:k1) = qrp(2:i1,2:j1,1:k1) +  evap(2:i1,2:j1,1:k1)
    Nrp(2:i1,2:j1,1:k1) = Nrp(2:i1,2:j1,1:k1) + Nevap(2:i1,2:j1,1:k1)

    deallocate (S,G)

  end subroutine evaporation

  !*********************************************************************

  real function sed_flux(Nin,Din,sig2,Ddiv,nnn)

  !*********************************************************************
  ! Function to calculate numerically the analytical solution of the
  ! sedimentation flux between Dmin and Dmax based on
  ! Feingold et al 1986 eq 17 -20.
  ! fall velocity is determined by alfa* D^beta with alfa+ beta taken as
  ! specified in Rogers and Yau 1989 Note here we work in D and in SI
  ! (in Roger+Yau in cm units + radius)
  ! flux is multiplied outside sed_flux with 1/rho_air to get proper
  ! kg/kg m/s units
  !
  ! M.C. van Zanten    August 2005
  !*********************************************************************
    use modglobal, only : pi,rhow
    implicit none

    real, intent(in) :: Nin, Din, sig2, Ddiv
    integer, intent(in) :: nnn
    !para. def. lognormal DSD (sig2 = ln^2 sigma_g), D sep. droplets from drops
    !,power of of D in integral

    real, parameter ::   C = rhow*pi/6.     &
                        ,D_intmin = 1e-6    &
                        ,D_intmax = 4.3e-3

    real ::  alfa         & ! constant in fall velocity relation
            ,beta         & ! power in fall vel. rel.
            ,D_min        & ! min integration limit
            ,D_max        & ! max integration limit
            ,flux           ![kg m^-2 s^-1]

    integer :: k

    flux = 0.0

    if (Din < Ddiv) then
      alfa = 3.e5*100  ![1/ms]
      beta = 2
      D_min = D_intmin
      D_max = Ddiv
      flux = C*Nin*alfa*erfint(beta,Din,D_min,D_max,sig2,nnn)
    else
      do k = 1,3
        select case(k)
        case(1)        ! fall speed ~ D^2
          alfa = 3.e5*100 ![1/m 1/s]
          beta = 2
          D_min = Ddiv
          D_max = 133e-6
        case(2)        ! fall speed ~ D
          alfa = 4e3     ![1/s]
          beta = 1
          D_min = 133e-6
          D_max = 1.25e-3
        case default         ! fall speed ~ sqrt(D)
          alfa = 1.4e3 *0.1  ![m^.5 1/s]
          beta = .5
          D_min = 1.25e-3
          D_max = D_intmax
        end select
        flux = flux + C*Nin*alfa*erfint(beta,Din,D_min,D_max,sig2,nnn)
      end do
    end if
      sed_flux = flux
  end function sed_flux

  !*********************************************************************
  !*********************************************************************

  real function liq_cont(Nin,Din,sig2,Ddiv,nnn)

  !*********************************************************************
  ! Function to calculate numerically the analytical solution of the
  ! liq. water content between Dmin and Dmax based on
  ! Feingold et al 1986 eq 17 -20.
  !
  ! M.C. van Zanten    September 2005
  !*********************************************************************
    use modglobal, only : pi,rhow
    implicit none

    real, intent(in) :: Nin, Din, sig2, Ddiv
    integer, intent(in) :: nnn
    !para. def. lognormal DSD (sig2 = ln^2 sigma_g), D sep. droplets from drops
    !,power of of D in integral

    real, parameter :: beta = 0           &
                      ,C = pi/6.*rhow     &
                      ,D_intmin = 80e-6    &   ! value of start of rain D
                      ,D_intmax = 3e-3         !4.3e-3    !  value is now max value for sqrt fall speed rel.

    real ::  D_min        & ! min integration limit
            ,D_max          ! max integration limit

    if (Din < Ddiv) then
    D_min = D_intmin
    D_max = Ddiv
    else
    D_min = Ddiv
    D_max = D_intmax
    end if

    liq_cont = C*Nin*erfint(beta,Din,D_min,D_max,sig2,nnn)

  end function liq_cont

  !*********************************************************************
  !*********************************************************************

  real function erfint(beta, D, D_min, D_max, sig2,nnn )

  !*********************************************************************
  ! Function to calculate erf(x) approximated by a polynomial as
  ! specified in 7.1.27 in Abramowitz and Stegun
  ! NB phi(x) = 0.5(erf(0.707107*x)+1) but 1 disappears by substraction
  !
  !*********************************************************************
    implicit none
    real, intent(in) :: beta, D, D_min, D_max, sig2
    integer, intent(in) :: nnn

    real, parameter :: eps = 1e-10       &
                      ,a1 = 0.278393    & !a1 till a4 constants in polynomial fit to the error
                      ,a2 = 0.230389    & !function 7.1.27 in Abramowitz and Stegun
                      ,a3 = 0.000972    &
                      ,a4 = 0.078108
    real :: nn, ymin, ymax, erfymin, erfymax, D_inv

    D_inv = 1./(eps + D)
    nn = beta + nnn

    ymin = 0.707107*(log(D_min*D_inv) - nn*sig2)/(sqrt(sig2))
    ymax = 0.707107*(log(D_max*D_inv) - nn*sig2)/(sqrt(sig2))

    erfymin = 1.-1./((1.+a1*abs(ymin) + a2*abs(ymin)**2 + a3*abs(ymin)**3 +a4*abs(ymin)**4)**4)
    erfymax = 1.-1./((1.+a1*abs(ymax) + a2*abs(ymax)**2 + a3*abs(ymax)**3 +a4*abs(ymax)**4)**4)
    if (ymin < 0.) then
      erfymin = -1.*erfymin
    end if
    if (ymax < 0.) then
      erfymax = -1.*erfymax
    end if
    erfint = D**nn*exp(0.5*nn**2*sig2)*0.5*(erfymax-erfymin)
  !  if (erfint < 0.) write(*,*)'erfint neg'
    if (erfint < 0.) erfint = 0.
  end function erfint



end module modbulkmicro


