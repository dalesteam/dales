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
  use modprecision, only : field_r
  use modtimer
  implicit none
  private
  public initbulkmicro, exitbulkmicro, bulkmicro

  real :: gamma25
  real :: gamma3
  real :: gamma35
  integer :: qrbase, qrroof, qcbase, qcroof
  contains

!> Initializes and allocates the arrays
  subroutine initbulkmicro
    use modglobal, only : i1,j1,k1,ih,jh
    use modmicrodata, only : lacz_gamma, Nr, Nrp, qr, qrp, thlpmcr, &
                             qtpmcr, Dvr, xr, mur, &
                             lbdr, &
                             precep, qrmask, qcmask
    implicit none

                                        ! Fields accessed by:
    allocate(Nr       (2:i1,2:j1,k1)  & ! dobulkmicrostat, dosimpleicestat
            ,qr       (2:i1,2:j1,k1)  & ! dobulkmicrostat, dosimpleicestat
            ,Nrp      (2:i1,2:j1,k1)  & ! bulkmicrotend, simpleicetend
            ,qrp      (2:i1,2:j1,k1)  & ! bulkmicrotend, simpleicetend
            ,Dvr      (2:i1,2:j1,k1)  & ! dobulkmicrostat
            ,precep   (2:i1,2:j1,k1)  ) ! dobulkmicrostat, dosimpleicestat, docape

    allocate(thlpmcr  (2:i1,2:j1,k1)  & !
            ,qtpmcr(2-ih:i1+ih,2-jh:j1+jh,k1) & ! ghost cells added here for modvarbudget
            ,xr       (2:i1,2:j1,k1)  & !
            ,mur      (2:i1,2:j1,k1)  & !
            ,lbdr     (2:i1,2:j1,k1)  & !
            ,qrmask   (2:i1,2:j1,k1)  & !
            ,qcmask   (2:i1,2:j1,k1)  )

    gamma25=lacz_gamma(2.5)
    gamma3=2.
    gamma35=lacz_gamma(3.5)

    !$acc enter data copyin(Nr, qr, Nrp, qrp, Dvr, precep, &
    !$acc&                  thlpmcr, qtpmcr, xr, mur, lbdr, qrmask, qcmask)

  end subroutine initbulkmicro

!> Cleaning up after the run
  subroutine exitbulkmicro
  !*********************************************************************
  ! subroutine exitbulkmicro
  !*********************************************************************
    use modmicrodata, only : Nr,Nrp,qr,qrp,thlpmcr,qtpmcr, &
                             Dvr,xr,mur,lbdr, &
                             precep,qrmask,qcmask
    implicit none

    !$acc exit data delete(Nr, qr, Nrp, qrp, Dvr, precep, &
    !$acc&                 thlpmcr, qtpmcr, xr, mur, lbdr, qrmask, qcmask)

    deallocate(Nr,Nrp,qr,qrp,thlpmcr,qtpmcr)
    deallocate(Dvr,xr,mur,lbdr)
    deallocate(precep,qrmask,qcmask)

  end subroutine exitbulkmicro

!> Calculates rain DSD integral properties & parameters xr, Dvr, lbdr, mur
  subroutine calculate_rain_parameters(Nr, qr)
    use modmicrodata, only : xr, Dvr, lbdr, mur, &
                             l_sb, l_mur_cst, mur_cst, pirhow, &
                             qrmask, xrmin, xrmax, xrmaxkk
    use modglobal, only : i1,j1,k1
    use modfields, only : rhof

    implicit none

    real(field_r), intent(in)    :: Nr  (2:i1, 2:j1, 1:k1), &
                                    qr  (2:i1, 2:j1, 1:k1)
    integer :: i,j,k

    call timer_tic('modbulkmicro/calculate_rain_parameters', 1)

    if (qrbase .gt. qrroof) return

    if (l_sb) then
      !$acc parallel loop collapse(3) default(present)
      do k = qrbase, qrroof
        do j = 2, j1
          do i = 2, i1
            if (qrmask(i,j,k)) then
              xr(i,j,k) = rhof(k) * qr(i,j,k) / Nr(i,j,k)

              ! to ensure xr is within bounds
              xr (i,j,k) = min(max(xr(i,j,k), xrmin), xrmax)
              Dvr(i,j,k) = (xr(i,j,k) / pirhow)**(1./3.)
            endif
          enddo
        enddo
      enddo

      if (l_mur_cst) then
        !$acc parallel loop collapse(3) default(present)
        do k = qrbase, qrroof
          do j = 2, j1
            do i = 2, i1
              mur(i,j,k) = mur_cst
            enddo
          enddo
        enddo

        !$acc parallel loop collapse(3) default(present)
        do k = qrbase, qrroof
          do j = 2, j1
            do i = 2, i1
              if (qrmask(i,j,k)) then
                lbdr(i,j,k) = ((mur_cst+3.)*(mur_cst+2.)*(mur_cst+1.))**(1./3.)/Dvr(i,j,k)
              endif
            enddo
          enddo
        enddo
      else
        ! mur = f(Dv)
        !$acc parallel loop collapse(3) default(present)
        do k = qrbase, qrroof
          do j = 2, j1
            do i = 2, i1
              if (qrmask(i,j,k)) then
                mur(i,j,k) = min(30.,- 1. + 0.008/ (qr(i,j,k)*rhof(k))**0.6)  ! G09b
                lbdr(i,j,k) = ((mur(i,j,k)+3.)*(mur(i,j,k)+2.)*(mur(i,j,k)+1.))**(1./3.)/Dvr(i,j,k)
              endif
            enddo
          enddo
        enddo
      endif
    else ! l_sb
      !$acc parallel loop collapse(3) default(present)
      do k = qrbase, qrroof
        do j = 2, j1
          do i = 2, i1
            if (qrmask(i,j,k)) then
              xr(i,j,k) = rhof(k) * qr(i,j,k) / Nr(i,j,k)

              ! to ensure x_pw is within bounds
              xr(i,j,k) = min(xr(i,j,k),xrmaxkk)
              Dvr(i,j,k) = (xr(i,j,k)/pirhow)**(1./3.)
            endif
          enddo
        enddo
      enddo
    endif ! l_sb

    call timer_toc('modbulkmicro/calculate_rain_parameters')
  end subroutine calculate_rain_parameters

!> Calculates the microphysical source term.
  subroutine bulkmicro
    use modglobal, only : i1,j1,k1,rdt,rk3step,timee,rlv,cp
    use modfields, only : sv0,svm,svp,qtp,thlp,ql0,exnf,rhof
    use modbulkmicrostat, only : bulkmicrotend
    use modmpi,    only : myid
    use modmicrodata, only : Nr, qr, Nrp, qrp, thlpmcr, qtpmcr, delt, &
                             l_sedc, l_mur_cst, l_lognormal, l_rain, &
                             qrmask, qrmin, qcmask, qcmin, &
                             mur_cst, inr, iqr
    implicit none
    integer :: i, j, k
    real :: qrtest,nr_cor,qr_cor
    real :: qrsum_neg, qrsum, Nrsum_neg, Nrsum

    !$acc parallel loop collapse(3) default(present)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          Nr(i,j,k) = sv0(i,j,k,inr)
          qr(i,j,k) = sv0(i,j,k,iqr)
          Nrp(i,j,k)     = 0.0
          qrp(i,j,k)     = 0.0
          thlpmcr(i,j,k) = 0.0
          qtpmcr(i,j,k)  = 0.0
        enddo
      enddo
    enddo

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
      qrsum_neg = 0.0
      qrsum = 0.0
      Nrsum_neg = 0.0
      Nrsum = 0.00
      !$acc parallel loop collapse(3) default(present) reduction(+: qrsum_neg, qrsum, Nrsum_neg, Nrsum)
      do k = 1, k1
        do j = 2, j1
          do i = 2, i1
            qrsum = qrsum + qr(i,j,k)
            Nrsum = Nrsum + Nr(i,j,k)
            if (qr(i,j,k) < 0.0) then
              qrsum_neg = qrsum_neg + qr(i,j,k)
              qr(i,j,k) = 0.0
            end if
            if (Nr(i,j,k) < 0.0) then
              Nrsum_neg = Nrsum_neg + Nr(i,j,k)
              Nr(i,j,k) = 0.0
            end if
          enddo
        enddo
      enddo

      ! LE: Commenting those out for now, popping up too often.
      !if ( -qrsum_neg > 0.000001*qrsum) then
      !  write(*,*)'amount of neg. qr thrown away is too high  ',timee,' sec'
      !end if
      !if ( -Nrsum_neg > 0.000001*Nrsum) then
      !   write(*,*)'amount of neg. Nr thrown away is too high  ',timee,' sec'
      !end if
    end if   ! l_rain

    !*********************************************************************
    ! Find gridpoints where the microphysics scheme should run
    !*********************************************************************

#if defined(DALES_GPU)
    ! Faster with OpenACC acceleration as it enables collapse(3)
    qrbase = k1 + 1
    qrroof = 1 - 1
    qcbase = k1 + 1
    qcroof = 1 - 1
    !$acc parallel loop collapse(3) default(present) reduction(min:qrbase,qcbase)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          ! Update mask prior to using it
          qrmask(i,j,k) = (qr(i,j,k) > qrmin .and. Nr(i,j,k) > 0.0)
          qcmask(i,j,k) = ql0(i,j,k) > qcmin
          if (qrmask(i,j,k)) then
            qrbase = min(k, qrbase)
          endif
          if (qcmask(i,j,k)) then
            qcbase = min(k, qcbase)
          endif
        enddo
      enddo
    enddo
    qrbase = max(1, qrbase)
    qcbase = max(1, qcbase)

    if (qrbase.le.k1 .or. qcbase.le.k1) then
      !$acc parallel loop collapse(3) default(present) reduction(max:qrroof,qcroof)
      do k = min(qrbase,qcbase), k1
        do j = 2, j1
          do i = 2, i1
            if (qrmask(i,j,k)) then
              qrroof = max(k, qrroof)
            endif
            if (qcmask(i,j,k)) then
              qcroof = max(k, qcroof)
            endif
          enddo
        enddo
      enddo
      qrroof = min(k1, qrroof)
      qcroof = min(k1, qcroof)
    endif
#else
    qrmask = qr.gt.qrmin.and.Nr.gt.0
    qrbase = k1 + 1
    qrroof = 1 - 1
    do k=1,k1
      if (any(qrmask(:,:,k))) then
        qrbase = max(1, k)
        exit
      endif
    enddo
    if (qrbase.le.k1) then
      do k=k1,qrbase,-1
        if (any(qrmask(:,:,k))) then
          qrroof = min(k1, k)
          exit
        endif
      enddo
    endif

    qcmask = ql0(2:i1,2:j1,1:k1).gt.qcmin
    qcbase = k1 + 1
    qcroof = 1 - 1
    do k=1,k1
      if (any(qcmask(:,:,k))) then
        qcbase = max(1, k)
        exit
      endif
    enddo
    if (qcbase.le.k1) then
      do k=k1,qcbase,-1
        if (any(qcmask(:,:,k))) then
          qcroof = min(k1, k)
          exit
        endif
      enddo
    endif
#endif

    ! if there is nothing to do, we can return at this point
    ! if (min(qrbase,qcbase).gt.max(qrroof,qcroof)) return

    !*********************************************************************
    ! calculate Rain DSD integral properties & parameters xr, Dvr, lbdr, mur
    !*********************************************************************
    if (l_rain) then
      call calculate_rain_parameters(Nr, qr)
    end if   ! l_rain

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
#if defined(DALES_GPU)
      call sedimentation_rain_gpu
#else
      call sedimentation_rain
#endif
      call bulkmicrotend
    endif

    !*********************************************************************
    ! remove negative values and non physical low values
    !*********************************************************************
    ! qcbase/qcroof are based on ql0.gt.qcmin and
    ! qrbase/qrroof are based on qr.gt.qrmin
    ! but we need boundaries to update qtp/thlp.
    !
    ! The difference between them comes from:
    !  * sedimentation_cloud updated qtpmcr/thlpmcr at qcbase-1
    !  * sedimentation_rain updated qrbase/qrroof,
    !    but at those levels qtmpcr/thlpmcr are either zero or
    !    already in the qc boundaries.
    if (qcbase.le.k1) qcbase = max(1, qcbase - 1)

    if (min(qrbase,qcbase) .gt. max(qrroof, qcroof)) return

    !$acc parallel loop collapse(3) default(present)
    do k = min(qrbase,qcbase), max(qrroof, qcroof)
      do j = 2, j1
        do i = 2, i1
          qrtest=svm(i,j,k,iqr)/delt+svp(i,j,k,iqr)+qrp(i,j,k)
          qr_cor = (0.5 - sign(0.5, qrtest*delt - qrmin)) * qrtest
          nr_cor = min(svm(i,j,k,inr)/delt+svp(i,j,k,inr)+Nrp(i,j,k), 0.)

          ! correction, after Jerome's implementation in Gales
          qtp (i,j,k) = qtp (i,j,k) + qtpmcr (i,j,k) + qr_cor
          thlp(i,j,k) = thlp(i,j,k) + thlpmcr(i,j,k) - qr_cor * (rlv/(cp*exnf(k)))

          svp(i,j,k,iqr) = svp(i,j,k,iqr) + qrp(i,j,k) - qr_cor
          svp(i,j,k,inr) = svp(i,j,k,inr) + Nrp(i,j,k) - nr_cor
          ! adjust negative qr tendencies at the end of the time-step
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
    use modglobal, only : i1,j1,rlv,cp
    use modmpi,    only : myid
    use modfields, only : exnf,rhof,ql0
    use modmicrodata, only : qrp, Nrp, qtpmcr, thlpmcr, &
                             qr, qcmask, pirhow, x_s, &
                             D0_kk, delt, k_1, k_2, k_au, k_c, Nc_0, &
                             l_sb
    implicit none
    integer i,j,k
    real :: au
    real :: tau  !  internal time scale
    real :: phi  !  correction function (see SB2001)
    real :: xc   !  mean mass of cloud water droplets
    real :: nuc  !  width parameter of cloud DSD

    call timer_tic('modbulkmicro/autoconversion', 1)

    if (qcbase.gt.qcroof) return

    if (l_sb) then
      !
      ! SB autoconversion
      !
      k_au = k_c/(20*x_s)

      !$acc parallel loop collapse(3) default(present)
      do k = qcbase, qcroof
        do j = 2, j1
          do i = 2, i1
             if (qcmask(i,j,k)) then
                nuc = 1.58*(rhof(k)*ql0(i,j,k)*1000.) +0.72-1. !G09a
                xc  = rhof(k) * ql0(i,j,k) / Nc_0 ! No eps0 necessary
                au = k_au * (nuc+2.) * (nuc+4.) / (nuc+1.)**2.    &
                          * (ql0(i,j,k) * xc)**2. * 1.225 ! *rho**2/rho/rho (= 1)

                tau = qr(i,j,k)/(ql0(i,j,k)+qr(i,j,k))
                phi = k_1 * tau**k_2 * (1.0 - tau**k_2)**3
                au = au * (1.0 + phi/(1.0 - tau)**2)

                qrp(i,j,k) = qrp(i,j,k) + au
                Nrp(i,j,k) = Nrp(i,j,k) + au / x_s
                qtpmcr(i,j,k) = qtpmcr(i,j,k) - au
                thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k))) * au

                !if (ql0(i,j,k) .lt. au*delt) then
                !  write(6,*)'au too large', au*delt, ql0(i,j,k), i, j, k, myid
                !end if
             endif
          enddo
        enddo
      enddo
    else
      !
      ! KK00 autoconversion
      !
      !$acc parallel loop collapse(3) default(present)
      do k = qcbase, qcroof
        do j = 2, j1
          do i = 2, i1
             if (qcmask(i,j,k)) then
                au = 1350.0 * ql0(i,j,k)**(2.47) * (Nc_0/1.0E6)**(-1.79)
                qrp(i,j,k) = qrp(i,j,k) + au
                qtpmcr(i,j,k) = qtpmcr(i,j,k) - au
                thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv / (cp*exnf(k)))*au
                Nrp(i,j,k) = Nrp(i,j,k) + au * rhof(k) / (pirhow*D0_kk**3.)

                !if (ql0(i,j,k) .lt. au*delt) then
                !  write(6,*)'au too large', au*delt, ql0(i,j,k), i, j, k, myid
                !end if
             endif
          enddo
        enddo
      enddo
    end if !l_sb

    call timer_toc('modbulkmicro/autoconversion')

  end subroutine autoconversion

  subroutine accretion
  !*********************************************************************
  ! determine accr. + self coll. + br-up rate and adjust qrp and Nrp
  ! accordingly. Break-up : Seifert (2007)
  !*********************************************************************
    use modglobal, only : i1,j1,rlv,cp
    use modfields, only : exnf,rhof,ql0
    use modmpi,    only : myid
    use modmicrodata, only : qr, Nr, qrp, Nrp, &
                             qrmask, qcmask, qtpmcr, thlpmcr, &
                             D_eq, delt, lbdr, Dvr, &
                             k_br, k_l, k_r, k_rr, kappa_r, pirhow, &
                             l_sb
    implicit none
    integer :: i,j,k

    real :: ac, sc, br
    real :: phi     !  correction function (see SB2001)
    real :: phi_br
    real :: tau     !  internal time scale

    call timer_tic('modbulkmicro/accretion', 1)

    if (l_sb) then
      !
      ! SB accretion
      !

      if (max(qrbase,qcbase).gt.min(qrroof,qcroof)) return

      !$acc parallel loop collapse(3) default(present)
      do k = max(qrbase,qcbase), min(qrroof, qcroof)
        do j = 2, j1
          do i = 2, i1
            if (qrmask(i,j,k) .and. qcmask(i,j,k)) then
               tau = qr(i,j,k)/(ql0(i,j,k)+qr(i,j,k))
               phi = (tau/(tau + k_l))**4.
               ac  = k_r * rhof(k) * ql0(i,j,k) * qr(i,j,k) * phi * (1.225/rhof(k))**0.5

               qrp(i,j,k) = qrp(i,j,k) + ac
               qtpmcr(i,j,k) = qtpmcr(i,j,k) - ac
               thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*ac

               !if (ql0(i,j,k) .lt. ac * delt) then
               !  write(6,*)'ac too large', ac*delt, ql0(i,j,k), i, j, k, myid
               !end if
            endif
          enddo
        enddo
      enddo

      !
      ! SB self-collection & Break-up
      !
      if (qrbase.gt.qrroof) return

      !$acc parallel loop collapse(3) default(present)
      do k = qrbase, qrroof
        do j = 2, j1
          do i = 2, i1
            if (qrmask(i,j,k)) then
               sc = k_rr *rhof(k)* qr(i,j,k) * Nr(i,j,k)  &
                    * (1 + kappa_r/lbdr(i,j,k)*pirhow**(1./3.))**(-9.)* (1.225/rhof(k))**0.5
               if (Dvr(i,j,k) .gt. 0.30E-3) then
                 phi_br = k_br * (Dvr(i,j,k)-D_eq)
                 br = (phi_br + 1.) * sc
               else
                 br = 0.
               endif

               Nrp(i,j,k) = Nrp(i,j,k) - sc + br
            endif
          enddo
        enddo
      enddo
    else
      !
      ! KK00 accretion
      !
      if (max(qrbase,qcbase).gt.min(qcroof,qcroof)) return

      !$acc parallel loop collapse(3) default(present)
      do k = max(qrbase,qcbase), min(qcroof,qrroof)
        do j = 2, j1
          do i = 2, i1
            if (qrmask(i,j,k) .and. qcmask(i,j,k)) then
              ac = 67.0 * (ql0(i,j,k) * qr(i,j,k))**1.15
              qrp(i,j,k) = qrp(i,j,k) + ac
              qtpmcr(i,j,k) = qtpmcr(i,j,k) - ac
              thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*ac

              !if (ql0(i,j,k) .lt. ac * delt) then
              !  write(6,*)'ac too large', ac*delt, ql0(i,j,k), i, j, k, myid
              !end if
            endif
          enddo
        enddo
      enddo
    end if !l_sb

    call timer_toc('modbulkmicro/accretion')
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
    use modglobal, only : i1,j1,rlv,cp,dzf,pi
    use modfields, only : rhof,exnf,ql0
    use modmicrodata, only : csed,c_St,rhow,sig_g,Nc_0, &
                             qtpmcr,thlpmcr,qcmask
    implicit none
    integer :: i, j, k
    real :: sedc

    call timer_tic('modbulkmicro/sedimentation_cloud', 1)

    if (qcbase .gt. qcroof) return

    csed = c_St*(3./(4.*pi*rhow))**(2./3.)*exp(5.*log(sig_g)**2.)

    !$acc parallel loop collapse(3) default(present)
    do k = qcbase, qcroof
      do j = 2, j1
        do i = 2, i1
          if (qcmask(i,j,k)) then
            sedc = csed*Nc_0**(-2./3.)*(ql0(i,j,k)*rhof(k))**(5./3.)

            !$acc atomic update
            qtpmcr(i,j,k)  = qtpmcr (i,j,k) - sedc /(dzf(k)*rhof(k))
            !$acc atomic update
            thlpmcr(i,j,k) = thlpmcr(i,j,k) + sedc * (rlv/(cp*exnf(k)))/(dzf(k)*rhof(k))

            if (k > 1) then
              !$acc atomic update
              qtpmcr(i,j,k-1)  = qtpmcr(i,j,k-1) + sedc / (dzf(k-1)*rhof(k-1))
              !$acc atomic update
              thlpmcr(i,j,k-1) = thlpmcr(i,j,k-1) - sedc * (rlv/(cp*exnf(k-1)))/(dzf(k-1)*rhof(k-1))
            end if
          endif
        enddo
      enddo
    enddo

    call timer_toc('modbulkmicro/sedimentation_cloud')

  end subroutine sedimentation_cloud

!> Sedimentaion of rain
!! sedimentation of drizzle water
!! - gen. gamma distr is assumed. Terminal velocities param according to
!!   Stevens & Seifert. Flux are calc. anal.
!! - l_lognormal =T : lognormal DSD is assumed with D_g and N known and
!!   sig_g assumed. Flux are calc. numerically with help of a
!!   polynomial function
!! - this version is reworked with a temporary holder which enables
!!   more collaspe with OpenACC acceleration, but slower on CPU.
  subroutine sedimentation_rain_gpu
    use modglobal, only : i1,j1,k1,eps1,dzf
    use modfields, only : rhof
    use modmpi,    only : myid
    use modmicrodata, only : Nr, Nrp, qr, qrp, precep, &
                             l_sb, l_lognormal, delt, &
                             qrmask, qrmin, pirhow, sig_gr, &
                             D_s, a_tvsb, b_tvsb, c_tvsb, &
                             Dvr, mur, lbdr

    implicit none
    integer :: i,j,k,jn,sedimbase
    integer :: n_spl      !<  sedimentation time splitting loop
    real    :: pwcont
    real :: delt_inv
    real :: Dgr           !<  lognormal geometric diameter
    real :: wfall_qr      !<  fall velocity for qr
    real :: wfall_Nr      !<  fall velocity for Nr
    real :: sed_qr
    real :: sed_Nr
    real(field_r), allocatable     :: qr_spl(:,:,:), Nr_spl(:,:,:)
    real(field_r), allocatable     :: qr_tmp(:,:,:), Nr_tmp(:,:,:)

    real,save :: dt_spl,wfallmax

    call timer_tic('modbulkmicro/sedimentation_rain', 1)

    ! zero the precipitation flux field
    ! the update below is not always performed

    !$acc parallel loop collapse(3) default(present)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          precep(i,j,k) = 0.0
        enddo
      enddo
    enddo

    if (qrbase .gt. qrroof) return

    allocate(qr_spl(2:i1,2:j1,1:k1))
    allocate(Nr_spl(2:i1,2:j1,1:k1))
    allocate(qr_tmp(2:i1,2:j1,1:k1))
    allocate(Nr_tmp(2:i1,2:j1,1:k1))

    !$acc enter data create(qr_spl, Nr_spl, qr_tmp, Nr_tmp)

    wfallmax = 9.9
    n_spl = ceiling(wfallmax*delt/(minval(dzf)))
    dt_spl = delt/real(n_spl)

    do jn = 1, n_spl ! time splitting loop

      if (jn .eq. 1) then
        !$acc parallel loop collapse(3) default(present)
        do k = 1, k1
          do j = 2, j1
            do i = 2, i1
              qr_spl(i,j,k) = qr(i,j,k)
              Nr_spl(i,j,k) = Nr(i,j,k)
              qr_tmp(i,j,k) = qr(i,j,k)
              Nr_tmp(i,j,k) = Nr(i,j,k)
            enddo
          enddo
        enddo
      else
        !Copy from tmp into spl
        !$acc parallel loop collapse(3) default(present)
        do k = 1, k1
          do j = 2, j1
            do i = 2, i1
              qr_spl(i,j,k) = qr_tmp(i,j,k)
              Nr_spl(i,j,k) = Nr_tmp(i,j,k)

              ! Update mask
              qrmask(i,j,k) = (qr_spl(i,j,k) > qrmin .and. Nr_spl(i,j,k) > 0.0)
            enddo
          enddo
        enddo

        ! lower the rain base by one level to include the rain fall
        ! from the previous step
        qrbase = max(1, qrbase - 1)

        call calculate_rain_parameters(Nr_spl, qr_spl)
      endif

      ! Compute precep
      if (jn == 1) then
        if (l_sb) then
          if (l_lognormal) then
            !$acc parallel loop collapse(3) default(present) private(Dgr)
            do k = qrbase, qrroof
              do j = 2, j1
                do i = 2, i1
                  if (qrmask(i,j,k)) then
                    Dgr = (exp(4.5*(log(sig_gr))**2))**(-1./3.)*Dvr(i,j,k)
                    sed_qr = 1.*sed_flux(Nr_spl(i,j,k),Dgr,log(sig_gr)**2,D_s,3)
                    pwcont = liq_cont(Nr_spl(i,j,k),Dgr,log(sig_gr)**2,D_s,3)
                    if (pwcont > eps1) then
                      sed_qr = (qr_spl(i,j,k)*rhof(k)/pwcont)*sed_qr
                    endif
                    precep(i,j,k) = sed_qr/rhof(k)   ! kg kg-1 m s-1
                  endif
                enddo
              enddo
            enddo
          else ! l_lognormal
            !$acc parallel loop collapse(3) default(present)
            do k = qrbase, qrroof
              do j = 2, j1
                do i = 2, i1
                  if (qrmask(i,j,k)) then
                    wfall_qr = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr(i,j,k))**(-1.*(mur(i,j,k)+4.))))
                    sed_qr  = wfall_qr*qr_spl(i,j,k)*rhof(k)
                    precep(i,j,k) = sed_qr/rhof(k)   ! kg kg-1 m s-1
                  endif
                enddo
              enddo
            enddo
          endif ! l_lognormal
        else ! l_sb
          !$acc parallel loop collapse(3) default(present)
          do k = qrbase, qrroof
            do j = 2, j1
              do i = 2, i1
                if (qrmask(i,j,k)) then
                  precep(i,j,k) = max(0., 0.006*1.0E6*Dvr(i,j,k) - 0.2) * qr_spl(i,j,k)
                endif
              enddo
            enddo
          enddo
        endif ! l_sb
      endif ! jn == 1

      sedimbase = qrbase

      ! k qrbase if == 1
      if (qrbase == 1) then
        sedimbase = sedimbase + 1
        k = 1
        if (l_sb) then
          if (l_lognormal) then
            !$acc parallel loop collapse(2) default(present) private(Dgr)
            do j = 2, j1
              do i = 2, i1
                if (qrmask(i,j,k)) then
                  ! correction for width of DSD
                  Dgr = (exp(4.5*(log(sig_gr))**2))**(-1./3.)*Dvr(i,j,k)
                  sed_qr = 1.*sed_flux(Nr_spl(i,j,k),Dgr,log(sig_gr)**2,D_s,3)
                  sed_Nr = 1./pirhow*sed_flux(Nr_spl(i,j,k),Dgr,log(sig_gr)**2,D_s,0)

                  ! correction for the fact that pwcont .ne. qr_spl
                  ! actually in this way for every grid box a fall velocity is determined
                  pwcont = liq_cont(Nr_spl(i,j,k),Dgr,log(sig_gr)**2,D_s,3)       ! note : kg m-3
                  if (pwcont > eps1) then
                    sed_qr = (qr_spl(i,j,k)*rhof(k)/pwcont)*sed_qr
                    ! or:
                    ! qr_spl*(sed_qr/pwcont) = qr_spl*fallvel.
                  endif

                  qr_tmp(i,j,k) = qr_tmp(i,j,k) - sed_qr*dt_spl/(dzf(k)*rhof(k))
                  Nr_tmp(i,j,k) = Nr_tmp(i,j,k) - sed_Nr*dt_spl/dzf(k)
                endif
              enddo
            enddo
          else ! l_lognormal
            !$acc parallel loop collapse(2) default(present)
            do j = 2, j1
              do i = 2, i1
                if (qrmask(i,j,k)) then
                  wfall_qr = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr(i,j,k))**(-1.*(mur(i,j,k)+4.))))
                  wfall_Nr = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr(i,j,k))**(-1.*(mur(i,j,k)+1.))))

                  sed_qr  = wfall_qr*qr_spl(i,j,k)*rhof(k)
                  sed_Nr  = wfall_Nr*Nr_spl(i,j,k)

                  qr_tmp(i,j,k) = qr_tmp(i,j,k) - sed_qr*dt_spl/(dzf(k)*rhof(k))
                  Nr_tmp(i,j,k) = Nr_tmp(i,j,k) - sed_Nr*dt_spl/dzf(k)
                endif
              enddo
            enddo
          endif ! l_lognormal
        else
          !
          ! KK00 rain sedimentation
          !
          !$acc parallel loop collapse(2) default(present)
          do j = 2, j1
            do i = 2, i1
              if (qrmask(i,j,k)) then
                sed_qr = max(0., 0.006*1.0E6*Dvr(i,j,k) - 0.2) * qr_spl(i,j,k)*rhof(k)
                sed_Nr = max(0.,0.0035*1.0E6*Dvr(i,j,k) - 0.1) * Nr_spl(i,j,k)

                qr_tmp(i,j,k) = qr_tmp(i,j,k) - sed_qr*dt_spl/(dzf(k)*rhof(k))
                Nr_tmp(i,j,k) = Nr_tmp(i,j,k) - sed_Nr*dt_spl/dzf(k)
              endif
            enddo
          enddo
        endif ! l_sb
      endif ! qrbase == 1

      if (l_sb) then
        !
        ! SB rain sedimentation
        !
        if (l_lognormal) then
          !$acc parallel loop collapse(3) default(present) private(Dgr)
          do k = sedimbase, qrroof
            do j = 2, j1
              do i = 2, i1
                if (qrmask(i,j,k)) then
                  ! correction for width of DSD
                  Dgr = (exp(4.5*(log(sig_gr))**2))**(-1./3.)*Dvr(i,j,k)
                  sed_qr = 1.*sed_flux(Nr_spl(i,j,k),Dgr,log(sig_gr)**2,D_s,3)
                  sed_Nr = 1./pirhow*sed_flux(Nr_spl(i,j,k),Dgr,log(sig_gr)**2,D_s,0)

                  ! correction for the fact that pwcont .ne. qr_spl
                  ! actually in this way for every grid box a fall velocity is determined
                  pwcont = liq_cont(Nr_spl(i,j,k),Dgr,log(sig_gr)**2,D_s,3)       ! note : kg m-3
                  if (pwcont > eps1) then
                    sed_qr = (qr_spl(i,j,k)*rhof(k)/pwcont)*sed_qr
                    ! or:
                    ! qr_spl*(sed_qr/pwcont) = qr_spl*fallvel.
                  endif

                  !$acc atomic update
                  qr_tmp(i,j,k) = qr_tmp(i,j,k) - sed_qr*dt_spl/(dzf(k)*rhof(k))
                  !$acc atomic update
                  Nr_tmp(i,j,k) = Nr_tmp(i,j,k) - sed_Nr*dt_spl/dzf(k)

                  !$acc atomic update
                  qr_tmp(i,j,k-1) = qr_tmp(i,j,k-1) + sed_qr*dt_spl/(dzf(k-1)*rhof(k-1))
                  !$acc atomic update
                  Nr_tmp(i,j,k-1) = Nr_tmp(i,j,k-1) + sed_Nr*dt_spl/dzf(k-1)
                endif
              enddo
            enddo
          enddo
        else
          !$acc parallel loop collapse(3) default(present)
          do k = sedimbase, qrroof
            do j = 2, j1
              do i = 2, i1
                if (qrmask(i,j,k)) then
                  wfall_qr = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr(i,j,k))**(-1.*(mur(i,j,k)+4.))))
                  wfall_Nr = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr(i,j,k))**(-1.*(mur(i,j,k)+1.))))

                  sed_qr  = wfall_qr*qr_spl(i,j,k)*rhof(k)
                  sed_Nr  = wfall_Nr*Nr_spl(i,j,k)

                  !$acc atomic update
                  qr_tmp(i,j,k) = qr_tmp(i,j,k) - sed_qr*dt_spl/(dzf(k)*rhof(k))
                  !$acc atomic update
                  Nr_tmp(i,j,k) = Nr_tmp(i,j,k) - sed_Nr*dt_spl/dzf(k)

                  !$acc atomic update
                  qr_tmp(i,j,k-1) = qr_tmp(i,j,k-1) + sed_qr*dt_spl/(dzf(k-1)*rhof(k-1))
                  !$acc atomic update
                  Nr_tmp(i,j,k-1) = Nr_tmp(i,j,k-1) + sed_Nr*dt_spl/dzf(k-1)
                endif
              enddo
            enddo
          enddo
        endif ! l_lognormal
      else
        !
        ! KK00 rain sedimentation
        !
        !$acc parallel loop collapse(3) default(present)
        do k = sedimbase, qrroof
          do j = 2, j1
            do i = 2, i1
              if (qrmask(i,j,k)) then
                sed_qr = max(0., 0.006*1.0E6*Dvr(i,j,k) - 0.2) * qr_spl(i,j,k)*rhof(k)
                sed_Nr = max(0.,0.0035*1.0E6*Dvr(i,j,k) - 0.1) * Nr_spl(i,j,k)

                !$acc atomic update
                qr_tmp(i,j,k) = qr_tmp(i,j,k) - sed_qr*dt_spl/(dzf(k)*rhof(k))
                !$acc atomic update
                Nr_tmp(i,j,k) = Nr_tmp(i,j,k) - sed_Nr*dt_spl/dzf(k)

                !$acc atomic update
                qr_tmp(i,j,k-1) = qr_tmp(i,j,k-1) + sed_qr*dt_spl/(dzf(k-1)*rhof(k-1))
                !$acc atomic update
                Nr_tmp(i,j,k-1) = Nr_tmp(i,j,k-1) + sed_Nr*dt_spl/dzf(k-1)
              endif
            enddo
          enddo
        enddo
      endif ! l_sb

    enddo ! time splitting loop

    ! the last time splitting step lowered the base level
    ! and we still need to adjust for it
    qrbase = max(1,qrbase-1)

    delt_inv = 1.0 / delt
    !$acc parallel loop collapse(3) default(present)
    do k = qrbase, qrroof
      do j = 2, j1
        do i = 2, i1
          Nrp(i,j,k) = Nrp(i,j,k) + (Nr_tmp(i,j,k) - Nr(i,j,k)) * delt_inv
          qrp(i,j,k) = qrp(i,j,k) + (qr_tmp(i,j,k) - qr(i,j,k)) * delt_inv
        enddo
      enddo
    enddo

    !$acc exit data delete(qr_spl, Nr_spl, qr_tmp, Nr_tmp)

    deallocate(qr_spl, Nr_spl, qr_tmp, Nr_tmp)

    call timer_toc('modbulkmicro/sedimentation_rain')

  end subroutine sedimentation_rain_gpu

  subroutine sedimentation_rain
    use modglobal, only : i1,j1,k1,eps1,dzf
    use modfields, only : rhof
    use modmpi,    only : myid
    use modmicrodata, only : Nr, Nrp, qr, qrp, precep, &
                             l_sb, l_lognormal, delt, &
                             qrmask, qrmin, pirhow, sig_gr, &
                             D_s, a_tvsb, b_tvsb, c_tvsb, &
                             Dvr, mur, lbdr

    implicit none
    integer :: i,j,k,jn
    integer :: n_spl      !<  sedimentation time splitting loop
    real    :: pwcont
    real :: Dgr           !<  lognormal geometric diameter
    real :: wfall_qr      !<  fall velocity for qr
    real :: wfall_Nr      !<  fall velocity for Nr
    real :: sed_qr
    real :: sed_Nr
    real(field_r), allocatable     :: qr_spl(:,:,:), Nr_spl(:,:,:)

    real,save :: dt_spl,wfallmax

    call timer_tic('modbulkmicro/sedimentation_rain', 1)

    precep = 0 ! zero the precipitation flux field
               ! the update below is not always performed

    if (qrbase.gt.qrroof) return

    allocate(qr_spl(2:i1,2:j1,1:k1))
    allocate(Nr_spl(2:i1,2:j1,1:k1))

    wfallmax = 9.9
    n_spl = ceiling(wfallmax*delt/(minval(dzf)))
    dt_spl = delt/real(n_spl)

    do jn=1,n_spl ! time splitting loop

      if (jn .eq. 1) then
        qr_spl = qr
        Nr_spl = Nr
      else
        ! update parameters after the first iteration

        ! a new mask
        qrmask = (qr_spl .gt. qrmin).and.(Nr_spl .gt. 0) ! BUG: added Nr_spl

        ! lower the rain base by one level to include the rain fall
        ! from the previous step
        qrbase = max(1, qrbase - 1)

        call calculate_rain_parameters(Nr_spl, qr_spl)
      endif

      if (l_sb) then
        if (l_lognormal) then
          do k = qrbase,qrroof
          do j = 2,j1
          do i = 2,i1
            if (qrmask(i,j,k)) then
              ! correction for width of DSD
              Dgr = (exp(4.5*(log(sig_gr))**2))**(-1./3.)*Dvr(i,j,k)
              sed_qr = 1.*sed_flux(Nr_spl(i,j,k),Dgr,log(sig_gr)**2,D_s,3)
              sed_Nr = 1./pirhow*sed_flux(Nr_spl(i,j,k),Dgr,log(sig_gr)**2,D_s,0)

              ! correction for the fact that pwcont .ne. qr_spl
              ! actually in this way for every grid box a fall velocity is determined
              pwcont = liq_cont(Nr_spl(i,j,k),Dgr,log(sig_gr)**2,D_s,3)       ! note : kg m-3
              if (pwcont > eps1) then
                sed_qr = (qr_spl(i,j,k)*rhof(k)/pwcont)*sed_qr
                ! or:
                ! qr_spl*(sed_qr/pwcont) = qr_spl*fallvel.
              endif

              qr_spl(i,j,k) = qr_spl(i,j,k) - sed_qr*dt_spl/(dzf(k)*rhof(k))
              Nr_spl(i,j,k) = Nr_spl(i,j,k) - sed_Nr*dt_spl/dzf(k)

              if (k .gt. 1) then
                qr_spl(i,j,k-1) = qr_spl(i,j,k-1) + sed_qr*dt_spl/(dzf(k-1)*rhof(k-1))
                Nr_spl(i,j,k-1) = Nr_spl(i,j,k-1) + sed_Nr*dt_spl/dzf(k-1)
              endif
              if (jn==1) then
                precep(i,j,k) = sed_qr/rhof(k)   ! kg kg-1 m s-1
              endif
            endif ! qr_spl threshold statement
          enddo
          enddo
          enddo
        else
          !
          ! SB rain sedimentation
          !
          do k=qrbase,qrroof
          do j=2,j1
          do i=2,i1
            if (qrmask(i,j,k)) then
              wfall_qr = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr(i,j,k))**(-1.*(mur(i,j,k)+4.))))
              wfall_Nr = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr(i,j,k))**(-1.*(mur(i,j,k)+1.))))

              sed_qr  = wfall_qr*qr_spl(i,j,k)*rhof(k) ! m/s * kg/m3
              sed_Nr  = wfall_Nr*Nr_spl(i,j,k)

              qr_spl(i,j,k) = qr_spl(i,j,k) - sed_qr*dt_spl/(dzf(k)*rhof(k))
              Nr_spl(i,j,k) = Nr_spl(i,j,k) - sed_Nr*dt_spl/dzf(k)

              if (k .gt. 1) then
                qr_spl(i,j,k-1) = qr_spl(i,j,k-1) + sed_qr*dt_spl/(dzf(k-1)*rhof(k-1))
                Nr_spl(i,j,k-1) = Nr_spl(i,j,k-1) + sed_Nr*dt_spl/dzf(k-1)
              endif
              if (jn==1) then
                precep(i,j,k) = sed_qr/rhof(k)   ! kg kg-1 m s-1
              endif
            endif
          enddo
          enddo
          enddo
        endif ! l_lognormal
      else
        !
        ! KK00 rain sedimentation
        !
        do k=qrbase,qrroof
        do j=2,j1
        do i=2,i1
          if (qrmask(i,j,k)) then
            sed_qr = max(0., 0.006*1.0E6*Dvr(i,j,k) - 0.2) * qr_spl(i,j,k)*rhof(k)
            sed_Nr = max(0.,0.0035*1.0E6*Dvr(i,j,k) - 0.1) * Nr_spl(i,j,k)

            qr_spl(i,j,k) = qr_spl(i,j,k) - sed_qr*dt_spl/(dzf(k)*rhof(k))
            Nr_spl(i,j,k) = Nr_spl(i,j,k) - sed_Nr*dt_spl/dzf(k)

            if (k .gt. 1) then
              qr_spl(i,j,k-1) = qr_spl(i,j,k-1) + sed_qr*dt_spl/(dzf(k-1)*rhof(k-1))
              Nr_spl(i,j,k-1) = Nr_spl(i,j,k-1) + sed_Nr*dt_spl/dzf(k-1)
            endif
            if (jn==1) then
              precep(i,j,k) = sed_qr/rhof(k)   ! kg kg-1 m s-1
            endif
          endif
        enddo
        enddo
        enddo
      endif ! l_sb

    enddo ! time splitting loop

    ! the last time splitting step lowered the base level
    ! and we still need to adjust for it
    qrbase = max(1,qrbase-1)

    Nrp(:,:,qrbase:qrroof) = Nrp(:,:,qrbase:qrroof) + &
      (Nr_spl(:,:,qrbase:qrroof) - Nr(:,:,qrbase:qrroof))/delt

    qrp(:,:,qrbase:qrroof) = qrp(:,:,qrbase:qrroof) + &
      (qr_spl(:,:,qrbase:qrroof) - qr(:,:,qrbase:qrroof))/delt

    deallocate(qr_spl, Nr_spl)

    call timer_toc('modbulkmicro/sedimentation_rain')

  end subroutine sedimentation_rain

  !*********************************************************************
  !*********************************************************************

  subroutine evaporation
  !*********************************************************************
  ! Evaporation of prec. : Seifert (2008)
  ! Cond. (S>0.) neglected (all water is condensed on cloud droplets)
  !*********************************************************************

    use modglobal, only : i1,j1,Rv,rlv,cp,pi,mygamma251,mygamma21
    use modfields, only : exnf,qt0,svm,qvsl,tmp0,ql0,esl,rhof
    use modmicrodata, only : Nr, mur, Dv, &
                             inr, iqr, Kt, &
                             l_sb, &
                             a_tvsb, b_tvsb, c_tvsb, &
                             nu_a, Sc_num, avf, bvf, &
                             c_Nevap, c_evapkk, delt, &
                             qrmask, lbdr, xr, Dvr, qrp, Nrp, &
                             qtpmcr, thlpmcr
    implicit none
    integer :: i,j,k
    integer :: numel

    real :: F !< ventilation factor
    real :: S !< super or undersaturation
    real :: G !< cond/evap rate of a drop

    real :: evap, Nevap

    call timer_tic('modbulkmicro/evaporation', 1)

    if (qrbase.gt.qrroof) return

    if (l_sb) then

      !$acc parallel loop collapse(3) default(present)
      do k = qrbase, qrroof
        do j = 2, j1
          do i = 2, i1
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

              evap = 2*pi*Nr(i,j,k)*G*F*S/rhof(k)
              Nevap = c_Nevap*evap*rhof(k)/xr(i,j,k)

              if (evap < -svm(i,j,k,iqr)/delt) then
                Nevap = - svm(i,j,k,inr)/delt
                evap  = - svm(i,j,k,iqr)/delt
              endif

              qrp(i,j,k) = qrp(i,j,k) + evap
              Nrp(i,j,k) = Nrp(i,j,k) + Nevap

              qtpmcr(i,j,k) = qtpmcr(i,j,k) - evap
              thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*evap
            endif
          enddo
        enddo
      enddo
    else ! l_sb
      !$acc parallel loop collapse(3) default(present)
      do k = qrbase, qrroof
        do j = 2, j1
          do i = 2, i1
            if (qrmask(i,j,k)) then
              S = min(0.,(qt0(i,j,k)-ql0(i,j,k))/qvsl(i,j,k)- 1.)
              G = (Rv * tmp0(i,j,k)) / (Dv*esl(i,j,k)) + rlv/(Kt*tmp0(i,j,k))*(rlv/(Rv*tmp0(i,j,k)) -1.)
              G = 1./G

              evap = c_evapkk*2*pi*Dvr(i,j,k)*G*S*Nr(i,j,k)/rhof(k)
              Nevap = evap*rhof(k)/xr(i,j,k)

              if (evap < -svm(i,j,k,iqr)/delt) then
                Nevap = - svm(i,j,k,inr)/delt
                evap  = - svm(i,j,k,iqr)/delt
              endif

              qrp(i,j,k) = qrp(i,j,k) + evap
              Nrp(i,j,k) = Nrp(i,j,k) + Nevap

              qtpmcr(i,j,k) = qtpmcr(i,j,k) - evap
              thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*evap
            endif
          enddo
        enddo
      enddo
    endif

    call timer_toc('modbulkmicro/evaporation')

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

    real(field_r), intent(in) :: Nin
    real         , intent(in) :: Din, sig2, Ddiv
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

    real(field_r), intent(in) :: Nin
    real         , intent(in) :: Din, sig2, Ddiv
    integer, intent(in) :: nnn
    !para. def. lognormal DSD (sig2 = ln^2 sigma_g), D sep. droplets from drops
    !,power of of D in integral

    real, parameter :: beta = 0           &
                      ,C = pi/6.*rhow     &
                      ,D_intmin = 80e-6    &   ! value of start of rain D
                      ,D_intmax = 3e-3         !4.3e-3    !  value is now max value for sqrt fall speed rel.

    real ::  D_min        & ! min integration limit
            ,D_max        & ! max integration limit
            ,sn

    sn = sign(0.5, Din - Ddiv)
    D_min = (0.5 - sn) * D_intmin + (0.5 + sn) * Ddiv
    D_max = (0.5 - sn) * Ddiv     + (0.5 + sn) * D_intmax

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

    erfymin = sign(erfymin, ymin)
    erfymax = sign(erfymax, ymax)

    erfint = max(0., D**nn*exp(0.5*nn**2*sig2)*0.5*(erfymax-erfymin))
  end function erfint
end module modbulkmicro
