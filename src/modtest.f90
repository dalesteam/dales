module modtest
  use modmicrodata
  implicit none
!  real :: wMax

  contains

  subroutine inittest
    use modglobal, only : k1,ih,i1,jh,j1
    use modmpi,    only : myid
    implicit none

    if (imicro==0) then
      l_rain=.false.
      l_sedc=.false.
      return
    end if

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

!    select case (ibulk)
!      case (ibulk_sb01,ibulk_k13)
!        wMax = 9.9 ! Please check this value, it was set for SB01 in the original code!
!      case (ibulk_kk00)
!        wMax = wMaxkk ! Approx. 5.8 m/s, corresponding to xrmaxkk in modmicrodata
!      case default
!        print *,"MODBULKMICRO: No valid bulkmicrophysics scheme selected! (ibulk)"
!        stop "MODBULKMICRO: No valid bulkmicrophysics scheme selected! (ibulk)"
!    end select

    allocate( Nr(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate( qr(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(Dvr(2-ih:i1+ih,2-jh:j1+jh,k1))
    Dvr=0.; Nr=0.; qr=0.

    allocate(qrevap(k1), qrsed(k1), qrsrc(k1))
    qrevap=0.; qrsed=0.; qrsrc=0.
   
  end subroutine inittest

  subroutine exittest
    implicit none
    deallocate(Dvr,qr,Nr)
    deallocate(qrevap, qrsed, qrsrc)
  end subroutine exittest

  subroutine cloud_sedimentation
    use modglobal, only : i1,j1,kmax,eps1,rlv,cp,dzf,pi,rhow,rkStep,rkMaxStep,rdt
    use modmpi,    only : myid
    use modfields, only : qt0,thl0,rhof,exnf,ql0
    implicit none
    integer :: i,j,k

    if (rkStep/=rkMaxStep) return
    if (.not. l_sedc) return

    sedc = 0.
    csed = c_St*(3./(4.*pi*rhow))**(2./3.)*exp(5.*log(sig_g)**2)

    do k=1,kmax
      do j=2,j1 
        do i=2,i1
          if (ql0(i,j,k)>0.) then
            sedc(i,j,k) = csed*(Nc_0)**(-2./3.)*(ql0(i,j,k)*rhof(k))**(5./3.)

            qt0(i,j,k)  = qt0(i,j,k)  + (sedc(i,j,k+1)-sedc(i,j,k))/(dzf(k)*rhof(k))*rdt
            thl0(i,j,k) = thl0(i,j,k) - (rlv/(cp*exnf(k))) &
                            *(sedc(i,j,k+1)-sedc(i,j,k))/(dzf(k)*rhof(k))*rdt
          endif
        end do
      end do
    end do

  end subroutine cloud_sedimentation

!! Semi-Lagrangian description of the rain sedimentation process,
!! based on Juang and Hong (MWR 2010)
!! Using this method, the time-splitting will hopefully not be necessary
  subroutine rain_sedimentation_sl
    use modglobal, only : i1,j1,kmax,k1,rdt,dzf,dzh,zh,zf,rkStep,rkMaxStep
    use modfields, only : rhof,sv0
    use modmpi,    only : myid
    implicit none
    real, dimension(1:k1)           :: qr_sedim,Nr_sedim  ! Variables after sedimentation
    real, dimension(1:k1)           :: wqr,wNr            ! Sedimentation velocities
    real    :: dummy,xrr,Dvrr
    integer :: i,j,k,kk

    ! Only perform this routine at the end of the last subtimestep
    if (rkStep/=rkMaxStep) return
    if (.not. l_rain) return

    qr_sedim=0.; Nr_sedim=0.; wqr=0.; wNr=0.
    Dvr(:,:,:)=0.

    select case (ibulk)
      case (ibulk_sb01)
      !===================================================================================!
      !=========== Seifert and Beheng (2001) sedimentation velocities ====================!
      !===================================================================================!
        if (l_lognormal) then
          stop 'STOP: Rain sedimentation for this set of parameters not implemented.'
        else
          do j=2,j1; do i=2,i1
            do k=1,k1
              if (sv0(i,j,k,iqr) > qrmin) then
                wqr(k) = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr(i,j,k))**(-1.*(mur(i,j,k)+4.))))
                wNr(k) = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr(i,j,k))**(-1.*(mur(i,j,k)+1.))))
              endif
            end do
            ! qr and Nr after sedimentation (qr_sedim and Nr_sedim) are determined
            ! in a separate subroutine
            call advec_rain_sl(wqr,sv0(i,j,:,iqr)*rhof,qr_sedim(:))
            call advec_rain_sl(wNr,sv0(i,j,:,iNr),Nr_sedim(:))
            sv0(i,j,:,iqr) = qr_sedim(:)/rhof(:)
            sv0(i,j,:,iNr) = Nr_sedim(:)
          end do; end do

        endif !l_lognormal

      case (ibulk_kk00)
      !===================================================================================!
      !=========== Khairoutdinov and Kogan (2000) sedimentation velocities ===============!
      !===================================================================================!

        ! Sedimentation is only in the vertical direction, so only one large loop over j and i
        do j=2,j1; do i=2,i1
          wqr(:)=0.; wNr(:)=0.

          ! First, the terminal velocity is calculated.
          do k=1,kmax
            xrr  = rhof(k)*sv0(i,j,k,iqr)/(sv0(i,j,k,iNr)+1e-8) ! average mass of a droplet
!            xrr  = min(xrr,xrmaxkk) ! to ensure x_pw is within borders
            Dvrr = min((xrr/pirhow)**(1./3.),Dvrmaxkk) ! Ensure that the droplets don't become too large, max 3mm arbitrarily chosen
            Dvr(i,j,k) = Dvrr ! Easy fix for statistics
  
            ! Calculate terminal velocities based on the description by KK00
            ! Note: value limited by the maximum drop diameter Dvrmax
            wqr(k) = max(0.1  ,6000.*Dvrr - 0.2) ! 1.e-3 should be 0.1 m/s (=6000*D0_kk-0.2)
            wNr(k) = max(0.075,3500.*Dvrr - 0.1) ! 1.e-3 should be 0.075 m/s (=3500*D0_kk-0.1)
          end do

          ! qr and Nr after sedimentation (qr_sedim and Nr_sedim) are determined
          ! in a separate subroutine
          call advec_rain_sl(wqr,sv0(i,j,:,iqr)*rhof,qr_sedim(:))
          call advec_rain_sl(wNr,sv0(i,j,:,iNr),Nr_sedim(:))

          qrsed(:) = qrsed(:) + (sv0(i,j,:,iqr)-qr_sedim/rhof(:))

          do k=kmax,1,-1
            precep(i,j,k) = precep(i,j,k+1) - (sv0(i,j,k,iqr)*rhof(k)-qr_sedim(k))*dzf(k)/rdt
          end do

          sv0(i,j,:,iqr) = qr_sedim(:)/rhof(:)
          sv0(i,j,:,iNr) = Nr_sedim(:)

          ! Make sure that the droplet concentration is zero when there is no
          ! rain water.
          where (sv0(i,j,:,iqr)<=0. .or. sv0(i,j,:,iNr)<=0.)
            sv0(i,j,:,iqr) = 0.
            sv0(i,j,:,iNr) = 0.
          end where

        end do; end do

      case (ibulk_k13)
        stop 'STOP: Rain sedimentation for this set of parameters not implemented.'
    end select

    ! Statistics assume downward flux of precipitation > 0.
    precep = -precep

  end subroutine rain_sedimentation_sl

  subroutine micro_kk00
    use modglobal, only : i1,j1,kmax,k1,rdt,Rv,Rd,rlv,cp,pi,rkStep,rkMaxStep, &
                          iTimeInt,iTimeWicker,iTimeTVD
    use modfields, only : thl0,thlm,qt0,qtm,ql0,sv0,svm,qvsl,esl,tmp0,rhof,exnf
    implicit none

    integer :: i,j,k
    real :: qcc,dq,auto,accr,S,G,Dvrr,xrr

    if (rkStep/=rkMaxStep) return
    if (.not. l_rain) return

    qr(:,:,:) = sv0(:,:,:,iqr)
    Nr(:,:,:) = sv0(:,:,:,iNr)

!    if (any(qr<0.)) write(*,*) 'Negative qr!!'
    where (qr<0.) qr=0. ! Remove any negative qr that is probably the result of advection
    
    tint = tint+rdt
    do k=1,kmax
      do j=2,j1
        do i=2,i1
          ! Find out if there is cloud or precipitation
          if ((ql0(i,j,k)+qr(i,j,k)) > 0.) then
            if (ql0(i,j,k) > 0.) then
              ! Do autoconversion and accretion
              qcc = ql0(i,j,k)
              ! Autoconversion and accretion as in Khairoutdinov and Kogan (2000)
              auto = 1350.*qcc**1.47 * (Nc_0/1.e6)**(-1.79)
              accr =   67.*qcc**0.15 * qr(i,j,k)**1.15
  
              ! Method implemented by Marat Khairoutdinov in SAM
              qcc  = qcc/(1.+rdt*(auto+accr))
              auto = auto*qcc
              accr = accr*qcc
  
              ! Limit the tendency so no negative concentrations are produced
              dq = min( ql0(i,j,k), rdt*(auto+accr) )

              ! Apply the tendencies
              sv0(i,j,k,iqr) = sv0(i,j,k,iqr) + dq
              sv0(i,j,k,iNr) = sv0(i,j,k,iNr) + max(0.,dq-rdt*accr)*rhof(k)/(pirhow*D0_kk**3)
              qt0(i,j,k)     = qt0(i,j,k) - dq
              thl0(i,j,k)    = thl0(i,j,k) + (rlv/(cp*exnf(k)))*dq
              ! Store the tendency for the statistics
              qrsrc(k) = qrsrc(k) + dq
              
            elseif (qr(i,j,k)>qrmin .and. ql0(i,j,k)==0.) then
!              xrr  = rhof(k)*qr(i,j,k)/Nr(i,j,k) ! average mass of a droplet
!              xrr  = min(xrr,xrmaxkk) ! to ensure x_pw is within borders
!              Dvrr = min((xrr/pirhow)**(1./3.),Dvrmax)

              S = min(0.,qt0(i,j,k)/qvsl(i,j,k)- 1.)
              G = (Rv * tmp0(i,j,k)) / (Dv*esl(i,j,k)) + rlv/(Kt*tmp0(i,j,k))*(rlv/(Rv*tmp0(i,j,k)) -1.)
              G = G**(-1)/rhow ! NOTE: [G]=kg/m-1 s-1, KK00 assumes units of m^2/s
              ! G is the same as coefficient A3 (Eq. 13-32) in Pruppacher and Klett (1978) and is given by the denominator in Eq. (13-28)
              ! This term differs by a factor rho_w from the equation used here! This is corrected later
!              dq = c_evapkk*2*pi*Dvrr*G*S*Nr(i,j,k)/rhof(k) * rdt
              dq = 3.*c_evapkk*G*(4.*pi*rhow/(3.*rhof(k)))**(2./3.) * &
                      (qr(i,j,k)*Nr(i,j,k)**2)**(1./3.)*S*rdt

              ! Limit the tendency
              dq = max(-.5*qr(i,j,k),dq) ! Note: dq.le.0
              
              ! Apply the tendencies
              sv0(i,j,k,iNr) = sv0(i,j,k,iNr) + dq/sv0(i,j,k,iqr)*sv0(i,j,k,iNr)
              sv0(i,j,k,iqr) = sv0(i,j,k,iqr) + dq
              qt0(i,j,k)     = qt0(i,j,k) - dq
              thl0(i,j,k)    = thl0(i,j,k) + (rlv/(cp*exnf(k)))*dq
              ! Store the tendency for the statistics
              qrevap(k) = qrevap(k) + dq
            else
              ! Evaporate the small amount of precip water immediately
              qt0(i,j,k) = qt0(i,j,k) + sv0(i,j,k,iqr)
              thl0(i,j,k) = thl0(i,j,k) - (rlv/(cp*exnf(k)))*sv0(i,j,k,iqr)

              ! Store the tendency for the statistics
              qrevap(k) = qrevap(k) - sv0(i,j,k,iqr)

              sv0(i,j,k,iqr) = 0.
              sv0(i,j,k,iNr) = 0.
              
            end if
          end if

          if (sv0(i,j,k,iqr)<=0. .or. sv0(i,j,k,iNr)<=0.) then
            sv0(i,j,k,:) = 0. ! Remove negative concentrations
          end if
          if (sv0(i,j,k,iNr)>0.) then
            ! Set a minimum to the number of rain drops
            ! If you don't, qr will somehow accumulate in the subcloud layer.
!            sv0(i,j,k,iNr) = max(10.,sv0(i,j,k,iNr)) ! Without this limiter, qr starts to accumulate in the subcloud layer
                                                    ! SAM also has this limiter but uses max(qr*coef,Nr) where coef is unknown
            sv0(i,j,k,iNr) = max(sv0(i,j,k,iqr)*rhof(k)/(pirhow*Dvrmaxkk**3),sv0(i,j,k,iNr))
          end if
        end do
      end do
    end do

    ! Depending on the time integration method, the previous fields have to be
    ! set equal to the new fields.
    if (iTimeInt==iTimeWicker .or. iTimeInt==iTimeTVD) then
      svm (2:i1,2:j1,1:k1,:) = sv0 (2:i1,2:j1,1:k1,:)
      qtm (2:i1,2:j1,1:k1)   = qt0 (2:i1,2:j1,1:k1)
      thlm(2:i1,2:j1,1:k1)   = thl0(2:i1,2:j1,1:k1)
    end if
    
  end subroutine micro_kk00

  !=========================================================================!
  !======== This routine handles the advection of qr and Nr using the ======!
  !======== sedimentation velocities calculated according to KK00 or SB01 ==!
  ! Original routines: WRF module_mp_wsm3.F                                 !
  ! After: Juang H.-M. H. and S.-Y. Hong, Monthly Weather Review (2010)     !
  ! Johan van der Dussen 2013                                               !
  !=========================================================================!
  subroutine advec_rain_sl(wSed,varf,varOut)
    use modglobal,   only : kmax,k1,zh,dzf,dzh,rdt,rtimee, &
                            rkStep,rkMaxStep, & ! Information about the RK scheme is required to evaluate the amount of surface prec
                            sedimMethod,sedimPCM,sedimPLM,sedimPPM
    use modfields,   only : rhof
    use modmpi,      only : myid
    implicit none
    ! In- and output
    real,dimension(1:k1),intent(in)  :: wSed,   & ! Sedimentation velocity at full level (positive downwards!)
                                        varf      ! Variable to be advected at full level (qr*rhof or Nr)
    real,dimension(1:k1),intent(out) :: varOut    ! Variable after sedimentation

    ! Additional variables
    real,dimension(1:k1) :: wSedh,zAh,dzAf,delth,deltf,Qplus,Qminus,QA
    real,parameter       :: cLimit=0.05 ! Coefficient used to limit the velocies to avoid negative concentrations
    real                 :: qsumb,qsum,qsumt,zsumb,zsum,zsumt,qrsedimsum,znorm,tl,th,dql,qql,dip,dim
    real                 :: dq,dqMax,dqMin
    real,dimension(1:k1) :: dqMono
    integer              :: k,ki,kb,kt
   
    varOut=0. 
    wSedh=0.; delth=0.; QA=0.; Qplus=0.; Qminus=0.

    do k=1,kmax
      ! These velocities are defined at the full gridlevels
      ! At half levels, the velocities can be determined using 2nd order interpolation
      wSedh(k+1) = (wSed(k+1)*dzf(k) + wSed(k)*dzf(k+1))/(2*dzh(k+1))
    end do
    wSedh(1) = wSed(1)

    ! Limit the top of the rain group?
    do k=2,k1
      if ( wSed(k).eq.0. ) wSedh(k)=wSed(k-1)
    end do

    ! Limit the velocities to avoid negative arrival cell depths
    ! (diffusivity)
    do k=kmax,1,-1
      if ((wSedh(k+1)-wSedh(k))*rdt/dzf(k) > cLimit ) then
        wSedh(k) = wSedh(k+1)-cLimit*dzf(k)/rdt
      end if
    end do

    ! Calculate the dz due to sedimentation
    ! Determine the arrival heights zA, using half gridlevels as departure
    ! heights
    do k=1,k1
      zAh(k) = zh(k)-wSedh(k)*rdt
    end do
      
    ! First calculate the arrival height difference
    dzAf(1:kmax) = (zAh(2:k1)-zAh(1:kmax))
    dzAf(k1) = zh(k1)-zAh(k1)

    ! Using the arrival heights, the arrival mass can be calculated:
    do k=1,kmax
      QA(k) = varf(k)*dzf(k)/(zAh(k+1)-zAh(k))
    end do
    QA(k1) = 0.

    !===== Selection for the sedimentation method
    select case (sedimMethod)
    case (sedimPCM) ! Piecewise constant method
      kb=1
      kt=1
      do k=1,kmax
        ! Find the rain water content at full gridlevels by finding the
        ! closest neighbours
        kb=max(kb-1,1)
        do while (zh(k) > zAh(kb) .and. kb < k1)
          kb=kb+1
        end do
        kb=max(kb-1,1) ! Note that the max is necessary for the lowest level

        ! The index of the arrival level just above the next full grid level
        kt=max(kt-1,1)
        do while (zh(k+1) > zAh(kt) .and. kt<k1)
          kt=kt+1
        end do
        kt=kt-1

        if ( kt-kb == 0 ) then
          varOut(k) = QA(kb)
        else ! if ( kt-kb >= 1 )
          zsumb = zAh(kb+1)-zh(k)
          qsumb = QA(kb)*zsumb    ! The quanity of rho*qr that is now in the bottom part of the gridcell
          zsumt = zh(k+1)-zAh(kt)
          qsumt = QA(kt)*zsumt    ! The quanity of rho*qr that is now in the top part of the gridcell
          zsum = 0.               ! If the fall velocities diverge strongly there is possibly qr from
          qsum = 0.               !  several complete grid cells compressed into one level...
          if (kt-kb > 1) then     ! ...that quantity is calculated here
            zsum = zsum + sum(dzAf(kb+1:kt-1))
            qsum = qsum + sum(QA(kb+1:kt-1)*dzAf(kb+1:kt-1))
          end if
          varOut(k) = (qsumb+qsum+qsumt)/(zsumb+zsum+zsumt)
        end if
      end do

    case (sedimPLM) ! Piecewise linear method

      do k=2,kmax
        dip=(QA(k+1)-QA(k))/(dzAf(k+1)+dzAf(k))
        dim=(QA(k)-QA(k-1))/(dzAf(k-1)+dzAf(k))
        if (dip*dim<=0.) then
          Qminus(k)=QA(k)
          Qplus(k) =QA(k)
        else
          Qplus(k) = QA(k)+.5*(dip+dim)*dzAf(k)
          Qminus(k) = 2.*QA(k)-Qplus(k)
          if (Qplus(k)<0..or.Qminus(k)<0.) then
            Qplus(k)=QA(k)
            Qminus(k)=QA(k)
          end if
        end if
      end do
      Qplus(1) = QA(1)
      Qminus(1) = QA(1)
      Qplus(k1)  = QA(k1)
      Qminus(k1) = QA(k1)

      ! Find the rain water content at full gridlevels by finding the
      ! closest neighbours
      kb=1
      kt=1
      do k=1,kmax
        ! Find the rain water content at full gridlevels by finding the
        ! closest neighbours
        kb=max(kb-1,1)
        do while (zh(k) > zAh(kb) .and. kb < k1)
          kb=kb+1
        end do
        kb=max(kb-1,1) ! Note that the max is necessary for the lowest level

        ! The index of the arrival level just above the next full grid level
        kt=max(kt-1,1)
        do while (zh(k+1) > zAh(kt) .and. kt<k1)
          kt=kt+1
        end do
        kt=kt-1

        ! Determine the normalized level height between arrival interface heights
        tl=(zh(k)   - zAh(kb))/dzAf(kb)  ! Lower grid cell face
        th=(zh(k+1) - zAh(kt))/dzAf(kt)  ! Upper grid cell face

        if ( kt-kb==0 ) then ! Multiple grid levels between 2 arrival levels
          ! Next, determine the integral of qr between levels zh(k) and zh(k+1)
          ! (integrate Eq. A5 Xiao et al. 2010)
          varOut(k) = .5*(Qplus(kb)-Qminus(kb))*(th**2-tl**2) &
                          + Qminus(kb)*(th-tl)
          varOut(k) = varOut(k)/(th-tl)

        else ! if kt-kb>0 at least one arrival level between two grid faces
          ! First determine the bottom part 
          dql=QA(kb)-.5*(Qplus(kb)-Qminus(kb))*tl**2-Qminus(kb)*tl
          qsum  = dql*dzAf(kb)
          zsum  = (1.-tl)*dzAf(kb)
          ! Sum the content of any intermediate levels
          if (kt-kb > 1) then
            zsum = zsum + sum(dzAf(kb+1:kt-1))
            qsum = qsum + sum(QA(kb+1:kt-1)*dzAf(kb+1:kt-1))
          end if
          ! Now the top part
          qsumt = .5*(Qplus(kt)-Qminus(kt))*th**2 &
                    + Qminus(kt)*th
          qsum  = qsum + qsumt*dzAf(kt)
          zsum  = zsum + th*dzAf(kt)

          varOut(k) = qsum/zsum
        end if

      end do

    case (sedimPPM) ! Piecewise polynomial method, most accurate and most expensive

      ! Lower boundary condition for dqMono
      dq=.5*(QA(2)-QA(1))
      dqMax=max(QA(2),QA(1))-QA(1)
      dqMin=QA(1)-min(QA(2),QA(1))
      dqMono(1) = sign(min(abs(dq),dqMax,dqMin),dq)
      ! Loop through the rest of the levels
      do k=2,kmax
        dq = .25*(QA(k+1)-QA(k-1))
        dqMax = max(QA(k+1),QA(k),QA(k-1)) - QA(k)
        dqMin = QA(k) - min(QA(k+1),QA(k),QA(k-1))
        dqMono(k) = sign(min(abs(dq),dqMax,dqMin),dq)
      end do

      ! Determine cell edge values for the polynomial interpolation
      do k=2,kmax
        Qminus(k)   = (QA(k-1)*(zAh(k+1)-zAh(k)) + QA(k)*(zAh(k)-zAh(k-1)))/(zAh(k+1)-zAh(k-1))-(dqMono(k)-dqMono(k-1))/3.
        Qplus (k-1) = Qminus(k)
      end do
      ! Set lower and upper boundary conditions
      Qminus(1)    = QA(1)
      Qplus (kmax:k1) = QA(kmax:k1)
      Qminus(k1) = QA(k1)
      ! Apply monotonic limiters
      do k=1,k1
        Qminus(k) = QA(k) - sign(min(abs(2*dqMono(k)),abs(Qminus(k)-QA(k))),dqMono(k))
        Qplus (k) = QA(k) + sign(min(abs(2*dqMono(k)),abs(Qplus (k)-QA(k))),dqMono(k))
      end do
  
      ! Find the rain water content at full gridlevels by finding the
      ! closest neighbours
      kb=1
      kt=1
      do k=1,kmax
        ! Find the rain water content at full gridlevels by finding the
        ! closest neighbours
        kb=max(kb-1,1) ! Start slightly below previous kb->more efficient than starting each loop at 1
        do while (zh(k) > zAh(kb) .and. kb < k1)
          kb=kb+1
        end do
        kb=max(kb-1,1) ! Note that the max is necessary for the lowest level

        ! The index of the arrival level just above the next full grid level
        kt=max(kt-1,1)
        do while (zh(k+1) > zAh(kt) .and. kt<k1)
          kt=kt+1
        end do
        kt=kt-1

        ! Determine the normalized level height between arrival interface heights
        tl=(zh(k)   - zAh(kb))/dzAf(kb)  ! Lower grid cell face
        th=(zh(k+1) - zAh(kt))/dzAf(kt)  ! Upper grid cell face

        if ( kt-kb==0 ) then ! Multiple grid levels between 2 arrival levels
          varOut(k) = (Qplus(kb)+Qminus(kb)-2*QA(kb))*(th**3-tl**3)   - &
                        (Qplus(kb)+2*Qminus(kb)-3*QA(kb))*(th**2-tl**2) + &
                         Qminus(kb)*(th-tl)
          varOut(k) = varOut(k)/(th-tl)
        else
          ! Bottom part
          qql = (Qplus(kb)+Qminus(kb)-2*QA(kb))*tl**3 - &
                (Qplus(kb)+2*Qminus(kb)-3*QA(kb))*tl**2 + &
                 Qminus(kb)*tl
          dql = QA(kb)-qql
          qsum = dql*dzAf(kb)
          zsum = (1.-tl)*dzAf(kb)
          ! Middle
          if (kt-kb > 1) then
            zsum = zsum + sum(dzAf(kb+1:kt-1))
            qsum = qsum + sum(QA(kb+1:kt-1)*dzAf(kb+1:kt-1))
          end if
          ! Top part
          qsumt = (Qplus(kt)+Qminus(kt)-2*QA(kt))*th**3 - &
                  (Qplus(kt)+2*Qminus(kt)-3*QA(kt))*th**2 + &
                   Qminus(kt)*th
          qsum  = qsum + qsumt*dzAf(kt)
          zsum  = zsum + th*dzAf(kt)

          varOut(k) = qsum/zsum
        end if

      end do

    end select

    psAccumLast = psAccum
    ! Note: [psAccum]=kg/m^3
    do k=1,kmax
      if (zAh(k) < 0.) then
        if (zAh(k+1)<0.) then
          psAccum = psAccum + QA(k)*dzAf(k)
        else
          psAccum = psAccum + QA(k)*(0.-zAh(k))
        end if
      end if
    end do

    !precs = precLast + precp*ticb(substep)
    !precp = tica(substep)*precp
    !if (substep==nSubStep) precLast = precs
    !precs = precLast + precp*rkb(rkStep)
    !precp = 0.
    !if (rkStep==rkMaxStep) precLast = precs

  end subroutine advec_rain_sl

end module modtest
