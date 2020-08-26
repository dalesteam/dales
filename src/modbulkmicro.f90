!> \file modbulkmicro.f90

!>
!!  Bulk microphysics.
!>
!! Calculates bulk microphysics using a two moment scheme.
!! \see  Seifert and Beheng (Atm. Res., 2001)
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
!*********************************************************************
  implicit none
  real :: gamma25
  real :: gamma3
  real :: gamma35
  contains

!> Initializes and allocates the arrays
  subroutine initbulkmicro
    use modglobal, only : i1,j1,k1
    use modmicrodata, only : lacz_gamma, Nr, Nrp, qltot, qr, qrp, thlpmcr, &
                             qtpmcr, sedc, Dvr, xr, mur, &
                             lbdr, au, ac, sc, br, evap, Nevap, &
                             precep, qrmask, qcmask
    implicit none

                                        ! Fields accessed by:
    allocate(Nr       (2:i1,2:j1,k1)  & ! dobulkmicrostat, dosimpleicestat
            ,qr       (2:i1,2:j1,k1)  & ! dobulkmicrostat, dosimpleicestat
            ,Nrp      (2:i1,2:j1,k1)  & ! bulkmicrotend, simpleicetend
            ,qrp      (2:i1,2:j1,k1)  & ! bulkmicrotend, simpleicetend
            ,Dvr      (2:i1,2:j1,k1)  & ! dobulkmicrostat
            ,precep   (2:i1,2:j1,k1)  ) ! dobulkmicrostat, dosimpleicestat, docape

    allocate(qltot    (2:i1,2:j1,k1)  & !
            ,thlpmcr  (2:i1,2:j1,k1)  & !
            ,qtpmcr   (2:i1,2:j1,k1)  & !
            ,sedc     (2:i1,2:j1,k1)  & !
            ,xr       (2:i1,2:j1,k1)  & !
            ,mur      (2:i1,2:j1,k1)  & !
            ,lbdr     (2:i1,2:j1,k1)  & !
            ,au       (2:i1,2:j1,k1)  & !
            ,ac       (2:i1,2:j1,k1)  & !
            ,sc       (2:i1,2:j1,k1)  & !
            ,br       (2:i1,2:j1,k1)  & !
            ,evap     (2:i1,2:j1,k1)  & !
            ,Nevap    (2:i1,2:j1,k1)  & !
            ,qrmask   (2:i1,2:j1,k1)  & !
            ,qcmask   (2:i1,2:j1,k1)  )

    gamma25=lacz_gamma(2.5)
    gamma3=2.
    gamma35=lacz_gamma(3.5)
  end subroutine initbulkmicro

!> Cleaning up after the run
  subroutine exitbulkmicro
  !*********************************************************************
  ! subroutine exitbulkmicro
  !*********************************************************************
    use modmicrodata, only : Nr,Nrp,qltot,qr,qrp,thlpmcr,qtpmcr, &
                             sedc,Dvr,xr,mur,lbdr, &
                             au,ac,sc,br,evap,Nevap, &
                             precep,qrmask,qcmask
    implicit none

    deallocate(Nr,Nrp,qltot,qr,qrp,thlpmcr,qtpmcr)

    deallocate(sedc,Dvr,xr,mur,lbdr, &
               au,ac,sc,br,evap,Nevap)

    deallocate(precep,qrmask,qcmask)

  end subroutine exitbulkmicro

!> Calculates the microphysical source term.
  subroutine bulkmicro
    use modglobal, only : i1,j1,k1,rdt,rk3step,timee,rlv,cp
    use modfields, only : sv0,svm,svp,qtp,thlp,ql0,exnf,rhof
    use modbulkmicrostat, only : bulkmicrotend
    use modmpi,    only : myid
    use modmicrodata, only : Nr, qr, Nrp, qrp, thlpmcr, qtpmcr, delt, &
                             l_sb, l_sedc, l_mur_cst, l_lognormal, l_mur_cst, l_rain, &
                             qrmask, qrmin, qcmask, qcmin, &
                             xr, Dvr, mur, lbdr, qltot, &
                             pirhow, inr, iqr, &
                             xrmin, xrmax, xrmaxkk, mur_cst
    implicit none
    integer :: i,j,k
    real :: qrtest,nrtest

    Nr = sv0(2:i1,2:j1,1:k1,inr)
    qr = sv0(2:i1,2:j1,1:k1,iqr)

    Nrp    = 0.0
    qrp    = 0.0
    thlpmcr = 0.0
    qtpmcr  = 0.0

    delt = rdt/ (4. - dble(rk3step))

    if (timee.eq.0 .and. rk3step.eq.1 .and. myid.eq.0) then
      write(*,*) 'l_lognormal',l_lognormal
      write(*,*) 'rhof(1)', rhof(1),' rhof(10)', rhof(10)
      write(*,*) 'l_mur_cst',l_mur_cst,' mur_cst',mur_cst
      write(*,*) 'nuc = param'
    endif

    !*********************************************************************
    ! remove neg. values of Nr and qr
    !*********************************************************************
    if (l_rain) then
       if (sum(qr, qr<0.) > 0.000001*sum(qr)) then
         write(*,*)'amount of neg. qr and Nr thrown away is too high  ',timee,' sec'
       end if
       if (sum(Nr, Nr<0.) > 0.000001*sum(Nr)) then
          write(*,*)'amount of neg. qr and Nr thrown away is too high  ',timee,' sec'
       end if

       Nr = max(0.,Nr)
       qr = max(0.,qr)
    end if   ! l_rain

    qrmask = qr.gt.qrmin.and.Nr.gt.0

    !*********************************************************************
    ! calculate qltot
    !*********************************************************************

    qltot = ql0(2:i1,2:j1,1:k1) + qr
    qcmask = ql0(2:i1,2:j1,1:k1).gt.qcmin

    !*********************************************************************
    ! calculate Rain DSD integral properties & parameters xr, Dvr, lbdr, mur
    !*********************************************************************
    if (l_rain) then

      xr = 0.
      Dvr = 0.
      mur = 30.
      lbdr = 0.

      if (l_sb) then
        do k=1,k1
        do j=2,j1
        do i=2,i1
           if (qrmask(i,j,k)) then
             ! REMOVED: JvdD Added eps0 to avoid floating point exception
             ! NOTE: Nr.gt.0 because of definition of qrmask above
             xr (i,j,k) = rhof(k)*qr(i,j,k)/Nr(i,j,k)

             ! to ensure xr is within borders
             xr (i,j,k) = min(max(xr(i,j,k),xrmin),xrmax)
             Dvr(i,j,k) = (xr(i,j,k)/pirhow)**(1./3.)
           endif
        enddo
        enddo
        enddo

        if (l_mur_cst) then
          ! mur = cst
          do k=1,k1
          do j=2,j1
          do i=2,i1
            if (qrmask(i,j,k)) then
               mur(i,j,k) = mur_cst
               lbdr(i,j,k) = ((mur(i,j,k)+3.)*(mur(i,j,k)+2.)*(mur(i,j,k)+1.))**(1./3.)/Dvr(i,j,k)
            endif
          enddo
          enddo
          enddo
        else
          ! mur = f(Dv)
          do k=1,k1
          do j=2,j1
          do i=2,i1
            if (qrmask(i,j,k)) then
              mur(i,j,k) = min(30.,- 1. + 0.008/ (qr(i,j,k)*rhof(k))**0.6)  ! G09b
              lbdr(i,j,k) = ((mur(i,j,k)+3.)*(mur(i,j,k)+2.)*(mur(i,j,k)+1.))**(1./3.)/Dvr(i,j,k)
            endif
          enddo
          enddo
          enddo
        endif
      else ! l_sb
         do k=1,k1
         do j=2,j1
         do i=2,i1
            if (qrmask(i,j,k)) then
              ! REMOVED: JvdD Added eps0 to avoid floating point exception
              ! NOTE: Nr.gt.0 because of definition of qrmask above
              xr(i,j,k) = rhof(k)*qr(i,j,k)/Nr(i,j,k)

              ! to ensure x_pw is within borders
              xr(i,j,k) = min(xr(i,j,k),xrmaxkk)
              Dvr(i,j,k) = (xr(i,j,k)/pirhow)**(1./3.)
            endif
         enddo
         enddo
         enddo
      endif ! l_sb
    endif ! l_rain

    !*********************************************************************
    ! call microphysical processes subroutines
    !*********************************************************************
    if (l_sedc) then
      call sedimentation_cloud
    endif

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

    sv0(2:i1,2:j1,1:k1,inr) = Nr
    sv0(2:i1,2:j1,1:k1,iqr) = qr

    !*********************************************************************
    ! remove negative values and non physical low values
    !*********************************************************************
    do k=1,k1
    do j=2,j1
    do i=2,i1
      qrtest=svm(i,j,k,iqr)+(svp(i,j,k,iqr)+qrp(i,j,k))*delt
      nrtest=svm(i,j,k,inr)+(svp(i,j,k,inr)+Nrp(i,j,k))*delt
      if ((qrtest < qrmin) .or. (nrtest < 0.)) then
        ! correction, after Jerome's implementation in Gales
        qtp(i,j,k) = qtp(i,j,k) + qtpmcr(i,j,k) + svm(i,j,k,iqr)/delt + svp(i,j,k,iqr) + qrp(i,j,k)
        thlp(i,j,k) = thlp(i,j,k) +thlpmcr(i,j,k) - (rlv/(cp*exnf(k)))*(svm(i,j,k,iqr)/delt + svp(i,j,k,iqr) + qrp(i,j,k))
        svp(i,j,k,iqr) = - svm(i,j,k,iqr)/delt
        svp(i,j,k,inr) = - svm(i,j,k,inr)/delt
      else
        svp(i,j,k,iqr) = svp(i,j,k,iqr)+qrp(i,j,k)
        svp(i,j,k,inr) = svp(i,j,k,inr)+Nrp(i,j,k)
        thlp(i,j,k)=thlp(i,j,k)+thlpmcr(i,j,k)
        qtp(i,j,k)=qtp(i,j,k)+qtpmcr(i,j,k)
        ! adjust negative qr tendencies at the end of the time-step
      endif
    enddo
    enddo
    enddo
  end subroutine bulkmicro

  !> Determine autoconversion rate and adjust qrp and Nrp accordingly
  !!
  !!   The autoconversion rate is formulated for f(x)=A*x**(nuc)*exp(-Bx),
  !!   decaying exponentially for droplet mass x.
  !!   It can easily be reformulated for f(x)=A*x**(nuc)*exp(-Bx**(mu)) and
  !!   by chosing mu=1/3 one would get a gamma distribution in drop diameter
  !!   -> faster rain formation. (Seifert)
  subroutine autoconversion
    use modglobal, only : i1,j1,k1,kmax,rlv,cp
    use modmpi,    only : myid
    use modfields, only : exnf,rhof,ql0
    use modmicrodata, only : au, qrp, Nrp, qtpmcr, thlpmcr, &
                             qltot, qcmask, pirhow, x_s, &
                             d0_kk, delt, k_1, k_2, k_au, k_c, Nc_0, &
                             l_sb
    implicit none
    integer i,j,k
    real :: tau  !  internal time scale
    real :: phi  !  correction function (see SB2001)
    real :: xc   !  mean mass of cloud water droplets
    real :: nuc  !  width parameter of cloud DSD

    au = 0.

    if (l_sb) then
      !
      ! SB autoconversion
      !
      k_au = k_c/(20*x_s)

      do k=1,k1
      do j=2,j1
      do i=2,i1
         if (qcmask(i,j,k)) then
            nuc = 1.58*(rhof(k)*ql0(i,j,k)*1000.) +0.72-1. !G09a
            xc  = rhof(k) * ql0(i,j,k) / Nc_0 ! No eps0 necessary
            au(i,j,k) = k_au * (nuc+2.) * (nuc+4.) / (nuc+1.)**2.    &
                        * (ql0(i,j,k) * xc)**2. * 1.225 ! *rho**2/rho/rho (= 1)

            tau = 1.0 - ql0(i,j,k) / qltot(i,j,k)
            phi = k_1 * tau**k_2 * (1.0 - tau**k_2)**3
            au(i,j,k) = au(i,j,k) * (1.0 + phi/(1.0 - tau)**2)
         endif
      enddo
      enddo
      enddo

      qrp = qrp + au
      Nrp = Nrp + au/x_s
      qtpmcr = qtpmcr - au
      do k=1,k1
        thlpmcr(:,:,k) = thlpmcr(:,:,k) + (rlv/(cp*exnf(k)))*au(:,:,k)
      enddo
    else
      !
      ! KK00 autoconversion
      !
      do k=1,k1
      do j=2,j1
      do i=2,i1
         if (qcmask(i,j,k)) then
            au     (i,j,k) = 1350.0 * ql0(i,j,k)**(2.47) * (Nc_0/1.0E6)**(-1.79)
         endif
      enddo
      enddo
      enddo

      qrp = qrp + au
      qtpmcr = qtpmcr - au
      do k=1,k1
        thlpmcr(:,:,k) = thlpmcr(:,:,k) + (rlv/(cp*exnf(k)))*au(:,:,k)
        Nrp(:,:,k) = Nrp(:,:,k) + au(:,:,k) * rhof(k)/(pirhow*D0_kk**3.)
      enddo
    end if !l_sb

    if (any(ql0(2:i1,2:j1,1:kmax) .lt. au(:,:,1:kmax)*delt)) then
      write(6,*)'au too large', count(ql0(2:i1,2:j1,1:kmax)/delt - au(2:i1,2:j1,1:kmax) .lt. 0.),myid
    end if
  end subroutine autoconversion

  subroutine accretion
  !*********************************************************************
  ! determine accr. + self coll. + br-up rate and adjust qrp and Nrp
  ! accordingly. Break-up : Seifert (2007)
  !*********************************************************************
    use modglobal, only : i1,j1,k1,kmax,rlv,cp
    use modfields, only : exnf,rhof,ql0
    use modmpi,    only : myid
    use modmicrodata, only : ac, sc, br, Nrp, lbdr, &
                             qrmask, qcmask, qr, Nr, qrp, qtpmcr, thlpmcr, qltot, &
                             D_eq, delt, Dvr, &
                             k_br, k_l, k_r, k_rr, kappa_r, pirhow, &
                             l_sb
    implicit none
    integer :: i,j,k

    real :: phi     !  correction function (see SB2001)
    real :: phi_br
    real :: tau     !  internal time scale

    ac = 0.
    sc = 0.
    br = 0.

    if (l_sb) then
      !
      ! SB accretion
      !
      do k=1,k1
      do j=2,j1
      do i=2,i1
        if (qrmask(i,j,k) .and. qcmask(i,j,k)) then
           tau            = 1.0 - ql0(i,j,k)/(qltot(i,j,k))
           phi            = (tau/(tau + k_l))**4.
           ac(i,j,k) = k_r * rhof(k) * ql0(i,j,k) * qr(i,j,k) * phi * (1.225/rhof(k))**0.5
        endif
      enddo
      enddo
      enddo

      qrp = qrp + ac
      qtpmcr = qtpmcr - ac
      do k=1,k1
        thlpmcr(:,:,k) = thlpmcr(:,:,k) + (rlv/(cp*exnf(k)))*ac(:,:,k)
      enddo

      !
      ! SB self-collection & Break-up
      !
      do k=1,k1
      do j=2,j1
      do i=2,i1
        if (qrmask(i,j,k)) then
           sc(i,j,k) = k_rr *rhof(k)* qr(i,j,k) * Nr(i,j,k)  &
                       * (1 + kappa_r/lbdr(i,j,k)*pirhow**(1./3.))**(-9.)* (1.225/rhof(k))**0.5
           if (Dvr(i,j,k) .gt. 0.30E-3) then
              phi_br = k_br * (Dvr(i,j,k)-D_eq)
              br(i,j,k) = (phi_br + 1.) * sc(i,j,k)
           endif
        endif
      enddo
      enddo
      enddo

      Nrp = Nrp - sc + br
    else
      !
      ! KK00 accretion
      !
      do k=1,k1
      do j=2,j1
      do i=2,i1
        if (qrmask(i,j,k) .and. qcmask(i,j,k)) then
          ac(i,j,k) = 67.0 * (ql0(i,j,k) * qr(i,j,k))**1.15
        endif
      enddo
      enddo
      enddo

      qrp = qrp + ac
      qtpmcr = qtpmcr - ac
      do k=1,k1
        thlpmcr(:,:,k) = thlpmcr(:,:,k) + (rlv/(cp*exnf(k)))*ac(:,:,k)
      enddo
    end if !l_sb

    if (any(ql0(2:i1,2:j1,1:kmax) .lt. ac(2:i1,2:j1,1:kmax) * delt)) then
      write(6,*)'ac too large', count(ql0(2:i1,2:j1,1:kmax)/delt - ac(2:i1,2:j1,1:kmax) .lt. 0.),myid
    end if
  end subroutine accretion

!> Sedimentation of cloud water ((Bretherton et al,GRL 2007))
!!
!!   The sedimentation of cloud droplets assumes a lognormal DSD in which the
!!   geometric std dev. is assumed to be fixed at 1.3.
!! sedimentation of cloud droplets
!! lognormal CDSD is assumed (1 free parameter : sig_g)
!! terminal velocity : Stokes velocity is assumed (v(D) ~ D^2)
!! flux is calc. anal.
  subroutine sedimentation_cloud
    use modglobal, only : i1,j1,k1,kmax,rlv,cp,dzf,pi
    use modfields, only : rhof,exnf,ql0
    use modmicrodata, only : sedc,csed,c_St,rhow,sig_g,Nc_0, &
                             qtpmcr,thlpmcr,qcmask
    implicit none
    integer :: i,j,k

    sedc = 0.
    csed = c_St*(3./(4.*pi*rhow))**(2./3.)*exp(5.*log(sig_g)**2.)

    do k=1,k1
    do j=2,j1
    do i=2,i1
      if (qcmask(i,j,k)) then
        sedc(i,j,k) = csed*Nc_0**(-2./3.)*(ql0(i,j,k)*rhof(k))**(5./3.)
      endif
    enddo
    enddo
    enddo

    do k=1,kmax
      qtpmcr(:,:,k)  = qtpmcr(:,:,k)  + (sedc(:,:,k+1)-sedc(:,:,k))/(dzf(k)*rhof(k))
      thlpmcr(:,:,k) = thlpmcr(:,:,k) - (rlv/(cp*exnf(k))) &
                       *(sedc(:,:,k+1)-sedc(:,:,k))/(dzf(k)*rhof(k))
    enddo
  end subroutine sedimentation_cloud


!> Sedimentaion of rain
!! sedimentation of drizzle water
!! - gen. gamma distr is assumed. Terminal velocities param according to
!!   Stevens & Seifert. Flux are calc. anal.
!! - l_lognormal =T : lognormal DSD is assumed with D_g and N known and
!!   sig_g assumed. Flux are calc. numerically with help of a
!!   polynomial function
  subroutine sedimentation_rain
    use modglobal, only : i1,j1,k1,kmax,eps1,dzf
    use modfields, only : rhof
    use modmpi,    only : myid
    use modmicrodata, only : Nr, Nrp, qr, qrp, precep, &
                             l_sb, l_lognormal, l_mur_cst, delt, &
                             qrmin, pirhow, sig_gr, &
                             xrmax, xrmin, xrmaxkk, &
                             D_s, a_tvsb, b_tvsb, c_tvsb, mur_cst, eps0

    implicit none
    integer :: i,k,jn
    integer :: n_spl      !<  sedimentation time splitting loop

    real,save :: dt_spl,wfallmax

    logical,allocatable,dimension(:,:,:) :: mask
    integer :: nmasked
    real :: pwcont
    real, allocatable :: qr_spl(:) &
                        ,Nr_spl(:) &
                        ,sed_qr(:) &
                        ,sed_Nr(:) &
                        ,rhof_spl(:)


    real, allocatable :: rhof_3d(:,:,:) &
                        ,qr_spl_3d(:,:,:) &
                        ,Nr_spl_3d(:,:,:) &
                        ,sed_qr_3d(:,:,:) &
                        ,sed_Nr_3d(:,:,:)

    real :: xr, Dvr, mur, lbdr, Dgr, wfall_qr, wfall_Nr

    allocate(qr_spl_3d(2:i1,2:j1,1:k1))
    allocate(Nr_spl_3d(2:i1,2:j1,1:k1))
    allocate(sed_qr_3d(2:i1,2:j1,1:k1))
    allocate(sed_Nr_3d(2:i1,2:j1,1:k1))
    allocate(rhof_3d(2:i1,2:j1,1:k1))
    do k=1,k1
      rhof_3d(:,:,k) = rhof(k)
    enddo

    nmasked = i1 * j1 * k1 / 5
    ! Allocate work variables
    allocate(sed_qr(nmasked))
    allocate(sed_Nr(nmasked))

    wfallmax = 9.9
    n_spl = ceiling(wfallmax*delt/(minval(dzf)))
    dt_spl = delt/real(n_spl)

    qr_spl_3d = qr
    Nr_spl_3d = Nr

    do jn=1,n_spl ! time splitting loop
      mask = qr_spl_3d .gt. qrmin
      nmasked = count(mask)
      if (nmasked > i1 * j1 * k1 / 5) then
        write(*,*) 'Too many masked points'
      endif

      ! Pack the input variables to a 1D array
      qr_spl = pack(qr_spl_3d, mask)
      Nr_spl = pack(Nr_spl_3d, mask)
      rhof_spl = pack(rhof_3d, mask)

      if (l_sb) then
        if (l_lognormal) then
          do i=1,nmasked
            ! JvdD Added eps0 to avoid division by zero
            ! to ensure xr is within borders
            xr = rhof_spl(i)*qr_spl(i)/(Nr_spl(i)+eps0)
            xr = min(max(xr,xrmin),xrmax)
            Dvr = (xr/pirhow)**(1./3.)

            ! correction for width of DSD
            Dgr = (exp(4.5*(log(sig_gr))**2))**(-1./3.)*Dvr
            sed_Nr(i) = 1./pirhow*sed_flux(Nr_spl(i),Dgr,log(sig_gr)**2,D_s,0)

            ! correction for the fact that pwcont .ne. qr_spl
            ! actually in this way for every grid box a fall velocity is determined

            sed_qr(i) = 1.*sed_flux(Nr_spl(i),Dgr,log(sig_gr)**2,D_s,3)
            pwcont = liq_cont(Nr_spl(i),Dgr,log(sig_gr)**2,D_s,3)       ! note : kg m-3
            if (pwcont > eps1) then
              ! or qr_spl*(sed_qr/pwcont) = qr_spl*fallvel.
              sed_qr(i) = (qr_spl(i)*rhof_spl(i)/pwcont)*sed_qr(i)
            endif
          enddo
        else
          !
          ! SB rain sedimentation
          !
          do i=1,nmasked
            ! TODO move out of loop
            if (l_mur_cst) then
              mur = mur_cst
            else
              mur = min(30.,- 1. + 0.008/ (qr_spl(i)*rhof_spl(i))**0.6)  ! G09b
            endif

            ! JvdD Added eps0 to avoid division by zero
            xr = rhof_spl(i)*qr_spl(i)/(Nr_spl(i)+eps0)

            ! to ensure xr is within borders
            xr = min(max(xr,xrmin),xrmax)

            Dvr = (xr/pirhow)**(1./3.)
            lbdr = ((mur+3.)*(mur+2.)*(mur+1.))**(1./3.)/Dvr
            wfall_qr = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr)**(-1.*(mur+4.))))
            wfall_Nr = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr)**(-1.*(mur+1.))))

            sed_qr(i) = wfall_qr*qr_spl(i)*rhof_spl(i)
            sed_Nr(i) = wfall_Nr*Nr_spl(i)
          enddo
        endif ! l_lognormal
      else
        !
        ! KK00 rain sedimentation
        !
        do i=1,nmasked
          ! JvdD added eps0 to avoid division by zero
          xr = rhof_spl(i)*qr_spl(i)/(Nr_spl(i)+eps0)

          ! to ensure xr is within borders
          xr = min(xr,xrmaxkk)

          Dvr = (xr/pirhow)**(1./3.)
          sed_qr(i) = max(0., 0.006*1.0E6*Dvr - 0.2) * qr_spl(i)*rhof_spl(i)
          sed_Nr(i) = max(0.,0.0035*1.0E6*Dvr - 0.1) * Nr_spl(i)
        enddo
      endif ! l_sb

      sed_Nr_3d = unpack(sed_Nr, mask, sed_Nr_3d)
      sed_qr_3d = unpack(sed_qr, mask, sed_qr_3d)
      do k=1,kmax
        Nr_spl_3d(:,:,k) = Nr_spl_3d(:,:,k) + &
                 (sed_Nr_3d(:,:,k+1) - sed_Nr_3d(:,:,k))*dt_spl/dzf(k)
        qr_spl_3d(:,:,k) = qr_spl_3d(:,:,k) + &
                 (sed_qr_3d(:,:,k+1) - sed_qr_3d(:,:,k))*dt_spl/(dzf(k)*rhof(k))
      enddo

      if (any(qr_spl_3d(:,:,1:kmax) .lt. 0.)) then
        write(6,*)'sed_qr too large', count(qr_spl_3d(:,:,1:kmax) .lt. 0.),myid
      endif

      if (jn==1) then
        do k=1,kmax
          precep(:,:,k) = sed_qr_3d(:,:,k)/rhof(k)   ! kg kg-1 m s-1
        enddo
      endif

      ! Clean up 1D vars
      deallocate(qr_spl,Nr_spl,rhof_spl)

    enddo ! time splitting loop

    Nrp = Nrp + (Nr_spl_3d - Nr)/delt
    qrp = qrp + (qr_spl_3d - qr)/delt

    ! Clean up 1D vars
    deallocate(sed_qr,sed_Nr)

    ! Clean up 3d vars
    deallocate(qr_spl_3d, Nr_spl_3d, sed_qr_3d, sed_Nr_3d, rhof_3d)
  end subroutine sedimentation_rain

  !*********************************************************************
  !*********************************************************************

  subroutine evaporation
  !*********************************************************************
  ! Evaporation of prec. : Seifert (2008)
  ! Cond. (S>0.) neglected (all water is condensed on cloud droplets)
  !*********************************************************************

    use modglobal, only : i1,j1,k1,Rv,rlv,cp,pi,mygamma251,mygamma21,lacz_gamma
    use modfields, only : exnf,qt0,svm,qvsl,tmp0,ql0,esl,rhof
    use modmicrodata, only : evap, Nevap, Nr, mur, Dv, &
                             inr, iqr, Kt, &
                             l_sb, &
                             a_tvsb, b_tvsb, c_tvsb, &
                             nu_a, Sc_num, avf, bvf, &
                             c_Nevap, c_evapkk, delt, &
                             qrmask, lbdr, xr, dvr, qrp, Nrp, &
                             qtpmcr, thlpmcr
    implicit none
    integer :: i,j,k
    integer :: numel

    real :: F !< ventilation factor
    real :: S !< super or undersaturation
    real :: G !< cond/evap rate of a drop

    evap = 0.
    Nevap = 0.

    if (l_sb) then
       do k=1,k1
       do j=2,j1
       do i=2,i1
         if (qrmask(i,j,k)) then
           numel=nint(mur(i,j,k)*100.)
           F = avf * mygamma21(numel)*Dvr(i,j,k) +  &
              bvf*Sc_num**(1./3.)*(a_tvsb/nu_a)**0.5*mygamma251(numel)*Dvr(i,j,k)**(3./2.) * &
              (1.-(1./2.)  *(b_tvsb/a_tvsb)    *(lbdr(i,j,k)/(   c_tvsb+lbdr(i,j,k)))**(mur(i,j,k)+2.5)  &
                 -(1./8.)  *(b_tvsb/a_tvsb)**2.*(lbdr(i,j,k)/(2.*c_tvsb+lbdr(i,j,k)))**(mur(i,j,k)+2.5)  &
                 -(1./16.) *(b_tvsb/a_tvsb)**3.*(lbdr(i,j,k)/(3.*c_tvsb+lbdr(i,j,k)))**(mur(i,j,k)+2.5) &
                 -(5./128.)*(b_tvsb/a_tvsb)**4.*(lbdr(i,j,k)/(4.*c_tvsb+lbdr(i,j,k)))**(mur(i,j,k)+2.5)  )
           S = min(0.,(qt0(i,j,k)-ql0(i,j,k))/qvsl(i,j,k)- 1.)
           G = (Rv * tmp0(i,j,k)) / (Dv*esl(i,j,k)) + rlv/(Kt*tmp0(i,j,k))*(rlv/(Rv*tmp0(i,j,k)) -1.)
           G = 1./G

           evap(i,j,k) = 2*pi*Nr(i,j,k)*G*F*S/rhof(k)
           Nevap(i,j,k) = c_Nevap*evap(i,j,k)*rhof(k)/xr(i,j,k)
         endif
       enddo
       enddo
       enddo
    else
       do k=1,k1
       do j=2,j1
       do i=2,i1
         if (qrmask(i,j,k)) then
           S = min(0.,(qt0(i,j,k)-ql0(i,j,k))/qvsl(i,j,k)- 1.)
           G = (Rv * tmp0(i,j,k)) / (Dv*esl(i,j,k)) + rlv/(Kt*tmp0(i,j,k))*(rlv/(Rv*tmp0(i,j,k)) -1.)
           G = 1./G

           evap(i,j,k) = c_evapkk*2*pi*Dvr(i,j,k)*G*S*Nr(i,j,k)/rhof(k)
           Nevap(i,j,k) = evap(i,j,k)*rhof(k)/xr(i,j,k)
         endif
       enddo
       enddo
       enddo
    endif

    do k=1,k1
    do j=2,j1
    do i=2,i1
      if (evap(i,j,k) < -svm(i,j,k,iqr)/delt .and. qrmask(i,j,k)) then
        Nevap(i,j,k) = - svm(i,j,k,inr)/delt
        evap (i,j,k) = - svm(i,j,k,iqr)/delt
      endif
    enddo
    enddo
    enddo

    qrp = qrp + evap
    Nrp = Nrp + Nevap
    qtpmcr = qtpmcr - evap
    do k=1,k1
      thlpmcr(:,:,k) = thlpmcr(:,:,k) + (rlv/(cp*exnf(k)))*evap(:,:,k)
    enddo
  end subroutine evaporation

  !*********************************************************************
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

    flux = 0.0

    if (Din < Ddiv) then
      alfa = 3.e5*100  ![1/ms]
      beta = 2
      D_min = D_intmin
      D_max = Ddiv
      flux = C*Nin*alfa*erfint(beta,Din,D_min,D_max,sig2,nnn)
    else
      ! fall speed ~ D^2
      alfa = 3.e5*100 ![1/m 1/s]
      beta = 2
      D_min = Ddiv
      D_max = 133e-6
      flux = flux + C*Nin*alfa*erfint(beta,Din,D_min,D_max,sig2,nnn)

      ! fall speed ~ D
      alfa = 4e3     ![1/s]
      beta = 1
      D_min = 133e-6
      D_max = 1.25e-3
      flux = flux + C*Nin*alfa*erfint(beta,Din,D_min,D_max,sig2,nnn)

      ! fall speed ~ sqrt(D)
      alfa = 1.4e3 *0.1  ![m^.5 1/s]
      beta = .5
      D_min = 1.25e-3
      D_max = D_intmax
      flux = flux + C*Nin*alfa*erfint(beta,Din,D_min,D_max,sig2,nnn)
    endif
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
