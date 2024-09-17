!> \file modbulkmicro_aerosol.f90
!!  Bulk microphysics with explicit aerosols.

!>
!! Calculates bulk microphysics using a two moment scheme and aerosol microphysics
!! using the M7 scheme.
!! \see  Seifert and Beheng (Atm. Res., 2001)
!! \see  Seifert and Beheng (Met Atm Phys, 2006)
!! \see  Stevens and Seifert (J. Meteorol. Soc. Japan, 2008)  (rain sedim, mur param)
!! \see  Seifert (J. Atm Sc., 2008) (rain evap)
!! \see  Khairoutdinov and Kogan (2000) (drizzle param : auto, accr, sedim, evap)
!! \see  Vignati, Wilson and Stier (2004) (M7 aerosol microphysics)
!!  \author Olivier Geoffroy, K.N.M.I.
!!  \author Margreet van Zanten, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Caspar Jungbacker, TU Delft
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
!  Copyright 1993-2024 Delft University of Technology, Wageningen University, Utrecht University, KNMI


module modbulkmicro_aerosol
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

  use modprecision, only: field_r
  use modtimer,     only: timer_tic, timer_toc
  use modmicrodata_aerosol
  use modglobal,    only: i1, j1, k1

  implicit none

  private

  public :: initaerosol, exitaerosol
  public :: bulkmicro_aerosol

  real :: gamma25
  real :: gamma3
  real :: gamma35
  integer :: qrbase, qrroof, qcbase, qcroof

contains

  !> Setup the modes
  subroutine initaerosol

    use modtracers, only: defined_tracers, tracer_prop, add_tracer

    integer :: itrac, imod

    ! First, set all static properties
    do imod = 1, maxmodes
      call modes(imod) % construct(name=modenames(imod), long_name=longnames(imod), sigma_g=sigma_g(imod))
    end do

    ! Ok, now add the species to the modes
    do itrac = 1, defined_tracers
      if (.not. tracer_prop(itrac) % laero) cycle

      ! Find the mode index
      imod = findloc(modenames, tracer_prop(itrac) % mode, dim=1)

      ! Add aerosol to mode
      call modes(imod) % add(tracer_prop(itrac))
    end do

    ! Finally, allocate memory
    do imod = 1, maxmodes
      call modes(imod) % init
      call modes(imod) % print
    end do

    ! To prevent this file from getting too big, just include the source
    ! of the scavenging LUT's
    ! CJ: on GPU, optimize lookup tables with CUDA texture memory?
    include 'scavenging.inc'

  end subroutine initaerosol

  !> Cleaning up after the run
  subroutine exitaerosol
    
  end subroutine exitaerosol

  !> Calculates the microphysical source term.
  subroutine bulkmicro_aerosol

    use modglobal,    only: rdt
    use modfields,    only: sv0
    use modbulkmicro, only: remove_negative_values, calculate_rain_parameters
    use modmicrodata, only: inr, iqr, delt, Nr, qr, Nrp, qrp, thlpmcr, qtpmcr, l_rain

    real(field_r) :: qrsum, qrsum_neg, Nrsum, Nrsum_neg
    integer :: i, j, k, s, imod

    delt = rdt

    ! Populate the fields, reset tendencies
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

    do imod = 1, maxmodes
      if (.not. modes(imod) % enabled) cycle

      do s = 1, modes(imod) % nspecies
        do k = 1, k1
          do j = 2, j1
            do i = 2, i1
              modes(imod) % aer_conc(i,j,k,s) = sv0(i,j,k,modes(imod) % species(s) % trac_idx)
              modes(imod) % aer_tend = 0.0_field_r
            end do
          end do
        end do
      end do
    end do
        
  end subroutine bulkmicro_aerosol

  !> Determine autoconversion rate and adjust qrp and Nrp accordingly
  !!
  !!   The autoconversion rate is formulated for f(x)=A*x**(nuc)*exp(-Bx),
  !!   decaying exponentially for droplet mass x.
  !!   It can easily be reformulated for f(x)=A*x**(nuc)*exp(-Bx**(mu)) and
  !!   by chosing mu=1/3 one would get a gamma distribution in drop diameter
  !!   -> faster rain formation. (Seifert)

  subroutine autoconversion_kk00_aerosol

    use modglobal,            only: i1, j1, rlv,cp
    use modfields,            only: exnf, rhof, ql0
    use modmicrodata,         only: Nc_0, qrp, qtpmcr, thlpmcr, Nrp, pirhow, D0_kk, qcmask
    use modmicrodata_aerosol, only: modes, iINC, iINR

    integer       :: i, j, k, s
    real(field_r) :: au

    call timer_tic('modbulkmicro_aerosol/autoconversion_aerosol', 1)

    if (qcbase.gt.qcroof) return
    
    !$acc parallel loop collapse(3) default(present) private(au)
    do k = qcbase, qcroof
      do j = 2, j1
        do i = 2, i1
           if (qcmask(i,j,k)) then
              au = 1350 * ql0(i,j,k)**(2.47_field_r) * (Nc_0/1.0E6_field_r)**(-1.79_field_r)

              qrp(i,j,k) = qrp(i,j,k) + au
              qtpmcr(i,j,k) = qtpmcr(i,j,k) - au
              thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv / (cp*exnf(k))) * au

              Nrp(i,j,k) = Nrp(i,j,k) + au * rhof(k) / (pirhow * D0_kk**3)

              ! From KK00 paper:
              ! "The sink of the cloud drop concentration due to autoconversion is combined with the sink
              !  due to accretion"
              ! Hence, we don't do anything to Nc here

              do s = 1, modes(iINC) % nspecies
                modes(iINC) % aer_tend(i,j,k,s) = modes(iINC) % aer_tend(i,j,k,s) - au / ql0(i,j,k) * modes(iINC) % aer_conc(i,j,k,s)
                modes(iINR) % aer_tend(i,j,k,s) = modes(iINR) % aer_tend(i,j,k,s) + au / ql0(i,j,k) * modes(iINC) % aer_conc(i,j,k,s)
              end do
           endif
        enddo
      enddo
    enddo

    call timer_toc('modbulkmicro_aerosol/autoconversion/aerosol')

  end subroutine autoconversion_kk00_aerosol

  !> Autoconversion according to Seifert & Beheng (2006), with aerosols
  subroutine autoconversion_sb06_aerosol

    use modfields,    only: rhof, ql0, exnf
    use modglobal,    only: i1, j1, rlv, cp
    use modmicrodata, only: qr, qrp, Nrp, qtpmcr, thlpmcr, k_1, k_2, k_au, x_s, qcmask
    use modmicrodata_aerosol, only: modes, iINC, iINR

    real(field_r) :: nuc, xc, au, tau, phi
    integer       :: i, j, k, s
    
    !$acc parallel loop collapse(3) default(present) private(nuc, xc, au, tau, phi)
    do k = qcbase, qcroof
      do j = 2, j1
        do i = 2, i1
           if (qcmask(i,j,k)) then
              nuc = 1.58_field_r * (rhof(k) * ql0(i,j,k) * 1000) + 0.72_field_r - 1 !G09a
              xc = rhof(k) * ql0(i,j,k) / Nc(i,j,k) ! No eps0 necessary
              au = k_au * (nuc+2) * (nuc+4) / (nuc+1)**2    &
                        * (ql0(i,j,k) * xc)**2 * 1.225_field_r ! *rho**2/rho/rho (= 1)

              tau = qr(i,j,k) / (ql0(i,j,k) + qr(i,j,k))
              phi = k_1 * tau**k_2 * (1 - tau**k_2)**3
              au = au * (1 + phi / (1 - tau)**2)

              qrp(i,j,k) = qrp(i,j,k) + au
              Nrp(i,j,k) = Nrp(i,j,k) + au / x_s
              qtpmcr(i,j,k) = qtpmcr(i,j,k) - au
              thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k))) * au

              modes(iINC) % aer_tend(i,j,k,1) = modes(iINC) % aer_tend(i,j,k,1) - au / xc * rhof(k)
              modes(iINR) % aer_tend(i,j,k,1) = modes(iINR) % aer_tend(i,j,k,1) + au / x_s * rhof(k) 

              do s = 1, modes(iINC) % nspecies
                modes(iINC) % aer_tend(i,j,k,s) = modes(iINC) % aer_tend(i,j,k,s) - au / ql0(i,j,k) * modes(iINC) % aer_conc(i,j,k,s)
                modes(iINR) % aer_tend(i,j,k,s) = modes(iINR) % aer_tend(i,j,k,s) + au / ql0(i,j,k) * modes(iINC) % aer_conc(i,j,k,s)
              end do
           endif
        enddo
      enddo
    enddo
  end subroutine autoconversion_sb06_aerosol

  subroutine accretion_kk00_aerosol

    use modfields,            only: ql0, exnf
    use modglobal,            only: i1, j1, rlv, cp
    use modmicrodata,         only: qr, qrp, qtpmcr, thlpmcr, qrmask, qcmask, rhow
    use modmicrodata_aerosol, only: modes, iINC

    real(field_r) :: ac, au, r_vc, xc
    integer       :: i, j, k, s

    associate(Nc => modes(iINC) % aer_conc(:,:,:,1), Ncp => modes(iINC) % aer_tend(:,:,:,1), &
              qap_inc => modes(iINC) % aer_conc, qap_inr => modes(iINR) % aer_tend)

    !$acc parallel loop collapse(3) default(present)
    do k = max(qrbase,qcbase), min(qcroof,qrroof)
      do j = 2, j1
        do i = 2, i1
          if (qrmask(i,j,k) .and. qcmask(i,j,k)) then
            ac = 67 * (ql0(i,j,k) * qr(i,j,k))**1.15_field_r
            qrp(i,j,k) = qrp(i,j,k) + ac
            qtpmcr(i,j,k) = qtpmcr(i,j,k) - ac
            thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*ac

            xc = ql0(i,j,k) / Nc(i,j,k)
            Ncp(i,j,k) = Ncp(i,j,k) - ac / xc
            
            do s = 2, modes(iINC) % aer_tend(i,j,k,s)
              qap_inc = qap_inc - ac / ql0(i,j,k) * qap
              modes(iINR) % aer_tend(i,j,k,s) = modes(iINR) % aer_tend(i,j,k,s) &
                                                + ac / ql0(i,j,k) * modes(iINR) % aer_conc(i,j,k,s)
            end do
          end if
        end do
      end do
    end do

    end associate

  end subroutine accretion_kk00_aerosol

  subroutine accretion_sb06_aerosol

    use modfields,    only: ql0, rhof, exnf
    use modglobal,    only: i1, j1, rlv, cp
    use modmicrodata, only: qr, qrp, qtpmcr, thlpmcr, Nr, Nrp, lbdr, Dvr, &
      &                     D_eq, k_br, pirhow, kappa_r, qrmask, qcmask, k_rr, &
      &                     k_l, k_r
    
    real(field_r) :: tau, phi, ac, sc, phi_br, br, xc
    integer       :: i, j, k, s

    if (max(qrbase,qcbase).gt.min(qrroof,qcroof)) return

    associate(Nc => modes(iINC) % aer_conc(:,:,:,1), Ncp => modes(iINC) % aer_tend(:,:,:,1), &
              qa_inc => modes(iINC) % aer_conc, qa_inr => modes(iINR) % aer_conc, &
              qap_inc => modes(iINC) % aer_tend, qap_inr => modes(iINR) % aer_tend)

    !$acc parallel loop collapse(3) default(present)
    do k = max(qrbase,qcbase), min(qrroof, qcroof)
      do j = 2, j1
        do i = 2, i1
          if (qrmask(i,j,k) .and. qcmask(i,j,k)) then
            tau = qr(i,j,k)/(ql0(i,j,k)+qr(i,j,k))
            phi = (tau/(tau + k_l))**4.
            ac = k_r * rhof(k) * ql0(i,j,k) * qr(i,j,k) * phi * (1.225/rhof(k))**0.5

            qrp(i,j,k) = qrp(i,j,k) + ac
            qtpmcr(i,j,k) = qtpmcr(i,j,k) - ac
            thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*ac

            xc = ql0(i,j,k) / Nc(i,j,k)
            Ncp(i,j,k) = Ncp(i,j,k) - ac / xc

            do s = 2, modes(iINC) % nspecies 
              qap_inc(i,j,k,s) = qap_inc(i,j,k,s) - ac / ql0(i,j,k) * qa_inc(i,j,k,s)
              qap_inr(i,j,k,s) = qap_inr(i,j,k,s) + ac / ql0(i,j,k) * qa_inc(i,j,k,s)
            end do
          endif
        enddo
      enddo
    enddo

    end associate

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
            br = merge((k_br * (Dvr(i,j,k) - D_eq) + 1) * sc, 0._field_r, Dvr(i,j,k) > 0.3e-3_field_r)
            Nrp(i,j,k) = Nrp(i,j,k) - sc + br

            ! Self-collection and break-up directly act on the rain drop number concentration
            modes(iINR) % aer_tend(i,j,k,1) = modes(iINR) % aer_tend(i,j,k,1) - sc + br
          endif
        enddo
      enddo
    enddo 
  end subroutine accretion_sb06_aerosol

!> Sedimentation of cloud water ((Bretherton et al,GRL 2007))
!!
!!   The sedimentation of cloud droplets assumes a lognormal DSD in which the
!!   geometric std dev. is assumed to be fixed at 1.3.
!! sedimentation of cloud droplets
!! lognormal CDSD is assumed (1 free parameter : sig_g)
!! terminal velocity : Stokes velocity is assumed (v(D) ~ D^2)
!! flux is calc. anal.
  subroutine sedimentation_cloud_aerosol
    use modglobal, only : i1,j1,rlv,cp,dzf,pi
    use modfields, only : rhof,exnf,ql0
    use modmicrodata, only : csed,c_St,rhow,sig_g,Nc_0, &
                             qtpmcr,thlpmcr,qcmask
    use modmicrodata_aerosol, only: modes, iINC
    implicit none
    integer :: i, j, k, s
    real :: sedc, sedc_n

    call timer_tic('modbulkmicro/sedimentation_cloud', 1)

    if (qcbase .gt. qcroof) return

    csed = c_St*(3./(4.*pi*rhow))**(2./3.)*exp(5.*log(sig_g)**2.)

    !$acc parallel loop gang vector collapse(3) default(present)
    do k = qcbase, qcroof
      do j = 2, j1
        do i = 2, i1
          if (qcmask(i,j,k)) then
            sedc = csed * Nc(i,j,k)**(-2._field_r/3._field_r) * (ql0(i,j,k) * rhof(k))**(5._field_r/3._field_r)
            sedc_n = sedc / (rhof(k) * ql0(i,j,k)) * Nc(i,j,k)

            !CJ: this is a lot of atomics. Alternative: recalculate sedc for layer below as well
            !$acc atomic update
            qtpmcr(i,j,k)  = qtpmcr (i,j,k) - sedc /(dzf(k)*rhof(k))
            !$acc atomic update
            thlpmcr(i,j,k) = thlpmcr(i,j,k) + sedc * (rlv/(cp*exnf(k)))/(dzf(k)*rhof(k))
            !$acc atomic update 
            modes(iINC) % aer_tend(i,j,k,1) = modes(iINC) % aer_tend(i,j,k,1) - sedc_n / dzf(k)

            do s = 2, modes(iINC) % nspecies
              !$acc atomic update
              modes(iINC) % aer_tend(i,j,k,s) = modes(iINC) % aer_tend(i,j,k,s) &
                                                - sedc / (rhof(k) * ql0(i,j,k)) * modes(iINC) % aer_conc(i,j,k,s)
            end do

            if (k > 1) then
              !$acc atomic update
              qtpmcr(i,j,k-1)  = qtpmcr(i,j,k-1) + sedc / (dzf(k-1)*rhof(k-1))
              !$acc atomic update
              thlpmcr(i,j,k-1) = thlpmcr(i,j,k-1) - sedc * (rlv/(cp*exnf(k-1)))/(dzf(k-1)*rhof(k-1))
              !$acc atomic update
              modes(iINC) % aer_tend(i,j,k-1,1) = modes(iINC) % aer_tend(i,j,k-1,1) + sedc_n / dzf(k-1)

              do s = 2, modes(iINC) % nspecies
                !$acc atomic update
                modes(iINC) % aer_tend(i,j,k-1,s) = modes(iINC) % aer_tend(i,j,k-1,s) &
                                                    + sedc / (rhof(k-1) * ql0(i,j,k-1)) * modes(iINC) % aer_conc(i,j,k-1,s)
              end do
            end if
          endif
        enddo
      enddo
    enddo

    call timer_toc('modbulkmicro/sedimentation_cloud')

  end subroutine sedimentation_cloud_aerosol

  !> Calculates rain drop sedimentation based on Khairoutdinov and Kogan (2000).
  !! Also takes into acoount in-rain aerosols.
  !! 
  !! @param Nr Rain drop number concentration.
  !! @param Nrp Tendency of the rain drop number concentration.
  !! @param qr Rain water specific humidity.
  !! @param qrp Tendency of the rain water specific humidty
  !! @param qa Aerosol mass mixing ratio.
  !! @param qap Tendency of aerosol mass mixing ratio.
  !! @param precep Precipitation.
  !! @param rhof Slab-averaged density at full model level.
  !! @param dzf Thickness of full model levels.
  !! @param qrbase Lowest model level which contains precipitation.
  !! @param qrroof Highest model level which contains precipitation.
  !! @param qrmask Rain mask.
  subroutine sedimentation_rain_kk00_aerosol(Nr, Nrp, qr, qrp, qa, qap, & 
    &                                        precep, rhof, dzf, qrbase, &
    &                                        qrroof, qrmask, Dvr)
    real(field_r), intent(inout) :: Nr(:,:,:), Nrp(:,:,:)
    real(field_r), intent(inout) :: qr(:,:,:), qrp(:,:,:)
    real(field_r), intent(inout) :: qa(:,:,:), qap(:,:,:)
    real(field_r), intent(out)   :: precep(:,:,:)
    real(field_r), intent(in)    :: rhof(:)
    real(field_r), intent(in)    :: dzf(:)
    integer,       intent(inout) :: qrbase, qrroof
    logical,       intent(in)    :: qrmask(:,:,:)
    real(field_r), intent(in)    :: Dvr(:,:,:) ! TODO: calculate Dvr on the fly?

    real(field_r), allocatable :: qr_spl(:,:,:), qr_tmp(:,:,:)
    real(field_r), allocatable :: Nr_spl(:,:,:), Nr_tmp(:,:,:)
    real(field_r), allocatable :: qa_spl(:,:,:), qa_tmp(:,:,:)
    real(field_r)              :: wfallmax, dt_spl
    integer                    :: i, j, k, s, jn
    integer                    :: sedimbase

    call timer_tic('modbulkmicro/sedimentation_rain_kk00_aerosol', 1)

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
    allocate(qa_spl(2:i1,2:j1,1:k1))
    allocate(qa_tmp(2:i1,2:j1,1:k1))

    !$acc enter data create(qr_spl, Nr_spl, qr_tmp, Nr_tmp)

    wfallmax = 9.9
    n_spl = ceiling(wfallmax*delt/(minval(dzf)))
    dt_spl = delt/real(n_spl)

    do jn = 1, n_spl

      if (jn == 1) then
        !$acc parallel loop collapse(3) default(present)
        do k = 1, k1
          do j = 2, j1
            do i = 2, i1
              qr_spl(i,j,k) = qr(i,j,k)
              Nr_spl(i,j,k) = Nr(i,j,k)
              qr_tmp(i,j,k) = qr(i,j,k)
              Nr_tmp(i,j,k) = Nr(i,j,k)
            end do
          end do
        end do 
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
            end do
          end do
        end do

        ! lower the rain base by one level to include the rain fall
        ! from the previous step
        qrbase = max(1, qrbase - 1)

        call calculate_rain_parameters(Nr_spl, qr_spl)
      end if ! jn == 1
      
      if (jn == 1) then
        !$acc parallel loop collapse(3) default(present)
        do k = qrbase, qrroof
          do j = 2, j1
            do i = 2, i1
              if (qrmask(i,j,k)) then
                precep(i,j,k) = max(0., 0.006*1.0E6*Dvr(i,j,k) - 0.2) * qr_spl(i,j,k)
              end if
            end do
          end do
        end do
      end if ! jn == 1

      sedimbase = qrbase

      if (qrbase == 1) then
        sedimbase = sedimbase + 1
        k = 1
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
      end if ! qrbase == 1

      !$acc parallel loop collapse(3) default(present)
      do k = sedimbase, qrroof
        do j = 2, j1
          do i = 2, i1
            if (qrmask(i,j,k)) then
              sed_qr = max(0._field_r, 0.006*1.0E6*Dvr(i,j,k) - 0.2_field_r) * qr_spl(i,j,k)*rhof(k)
              sed_Nr = max(0._field_r,0.0035*1.0E6*Dvr(i,j,k) - 0.1_field_r) * Nr_spl(i,j,k)

              ! Save sed_qr for scavenging
              sed_qr_dup(i,j,k) = sed_qr

              !$acc atomic update
              qr_tmp(i,j,k) = qr_tmp(i,j,k) - sed_qr*dt_spl/(dzf(k)*rhof(k))
              !$acc atomic update
              Nr_tmp(i,j,k) = Nr_tmp(i,j,k) - sed_Nr*dt_spl/dzf(k)

              !$acc atomic update
              qr_tmp(i,j,k-1) = qr_tmp(i,j,k-1) + sed_qr*dt_spl/(dzf(k-1)*rhof(k-1))
              !$acc atomic update
              Nr_tmp(i,j,k-1) = Nr_tmp(i,j,k-1) + sed_Nr*dt_spl/dzf(k-1)

              do s = 1, nspec
                ! What goes to the level below
                a_out = merge( &
                  min(rhof(k) * dzf(k) / dt_spl, sed_qr(i,j,k) / qr_spl(i,j,k)), &
                  0._field_r, &
                  qr_spl(i,j,k) > 0._field_r
                ) * qa_spl(i,j,k,s)

                !$acc atomic update
                qa_spl(i,j,k,s) = qa_spl(i,j,k,s) - a_out * dt_spl / (dzf(k) * rhof(k))

                !$acc atomic update
                qa_spl(i,j,k-1,s) = qa_spl(i,j,k-1,s) + a_out * dt_spl / (dzf(k-1) * rhof(k-1))
              end do
            end if
          end do
        end do
      end do
    end do

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
        end do
      end do
    end do

    !$acc exit data delete(qr_spl, Nr_spl, qr_tmp, Nr_tmp)

    deallocate(qr_spl, Nr_spl, qr_tmp, Nr_tmp, qa_spl, qa_tmp)

    call timer_toc('modbulkmicro/sedimentation_rain_kk00_aerosol')

  end subroutine sedimentation_rain_kk00_aerosol

  !> Calculates sedimentation of rain drops according to Seifert & Beheng (2006). 
  !! Also takes into account the in-rain aerosols.
  !!
  !! @param Nr Rain drop number concentration.
  !! @param Nrp Tendency of rain drop number concentration.
  !! @param qr Rain water mixing ratio.
  !! @param qrp Tendency of rain water mixing ratio.
  !! @param qa Aerosol mass mixing ratio.
  !! @param qap Tendency of aerosol mass mixing ratio.
  !! @param precip Precipitation. 
  !! @param rhof Density at full model levels.
  !! @param dzf Thickness of full model levels.
  !! @param qrbase Lowest model level which contains precipitation.
  !! @param qrroof Highest model level which contains precipitation.
  !! @param qrmask Rain mask.
  !! @param Dvr
  subroutine sedimentation_rain_sb06_aerosol(Nr, Nrp, qr, qrp, qa, qap, & 
    &                                        precip, rhof, dzf, qrbase, &
    &                                        qrroof, qrmask, Dvr)
    real(field_r), intent(inout) :: Nr, Nrp
    real(field_r), intent(inout) :: qr, qrp
    real(field_r), intent(inout) :: qa, qap
    real(field_r), intent(out)   :: precip
    real(field_r), intent(in)    :: rhof
    real(field_r), intent(in)    :: dzf
    integer,       intent(inout) :: qrbase, qrroof
    logical,       intent(in)    :: qrmask
    real(field_r), intent(in)    :: Dvr

    real(field_r) :: qr_spl(:,:,:), qr_tmp(:,:,:)
    real(field_r) :: Nr_spl(:,:,:), Nr_tmp(:,:,:)
    real(field_r) :: qa_spl(:,:,:,:), qa_tmp(:,:,:,:) 
    real(field_r) :: wfallmax, dt_spl
    integer       :: i, j, k, s, jn
    integer       :: sedimbase

    call timer_tic('modbulkmicro/sedimentation_rain_sb06_aerosol', 1)

    !$acc parallel loop collapse(3) default(present)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          precip(i,j,k) = 0.0
        enddo
      enddo
    enddo

    if (qrbase .gt. qrroof) return

    ! TODO: don't reallocate these temporary arrays every time step
    allocate(qr_spl(2:i1,2:j1,1:k1))
    allocate(Nr_spl(2:i1,2:j1,1:k1))
    allocate(qr_tmp(2:i1,2:j1,1:k1))
    allocate(Nr_tmp(2:i1,2:j1,1:k1))
    allocate(qa_spl(2:i1,2:j1,1:k1,1:nspec))
    allocate(qa_tmp(2:i1,2:j1,1:k1,1:nspec))

    !$acc enter data create(qr_spl, Nr_spl, qr_tmp, Nr_tmp)

    wfallmax = 9.9
    n_spl = ceiling(wfallmax*delt/(minval(dzf)))
    dt_spl = delt/real(n_spl)

    do jn = 1, n_spl ! time splitting loop
      
      if (jn == 1) then
        !$acc parallel loop collapse(3) default(present)
        do k = 1, k1
          do j = 2, j1
            do i = 2, i1
              qr_spl(i,j,k) = qr(i,j,k)
              Nr_spl(i,j,k) = Nr(i,j,k)
              qr_tmp(i,j,k) = qr(i,j,k)
              Nr_tmp(i,j,k) = Nr(i,j,k)
            end do
          end do
        end do 

        !$acc parallel loop collapse(4) default(present)
        do s = 1, nspec
          do k = 1, k1
            do j = 2, j1
              do i = 2, i1
                qa_spl(i,j,k,s) = qa(i,j,k,s)
                qa_tmp(i,j,k,s) = qa(i,j,k,s)
              end do
            end do
          end do
        end do
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
            end do
          end do
        end do

        ! lower the rain base by one level to include the rain fall
        ! from the previous step
        qrbase = max(1, qrbase - 1)

        call calculate_rain_parameters(Nr_spl, qr_spl)
      end if ! jn == 1

      if (jn == 1) then
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
                  end if
                  precip(i,j,k) = sed_qr/rhof(k)   ! kg kg-1 m s-1
                end if
              end do
            end do
          end do 
        else ! l_nognormal
          !$acc parallel loop collapse(3) default(present)
          do k = qrbase, qrroof
            do j = 2, j1
              do i = 2, i1
                if (qrmask(i,j,k)) then
                  wfall_qr = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr(i,j,k))**(-1.*(mur(i,j,k)+4.))))
                  sed_qr  = wfall_qr*qr_spl(i,j,k)*rhof(k)
                  precip(i,j,k) = sed_qr/rhof(k)   ! kg kg-1 m s-1
                end if
              end do
            end do
          end do  
        end if ! l_lognormal
      end if ! jn == 1

      sedimbase = qrbase

      if (qrbase == 1) then
        sedimbase = sedimbase + 1
        k = 1
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
                end if

                qr_tmp(i,j,k) = qr_tmp(i,j,k) - sed_qr*dt_spl/(dzf(k)*rhof(k))
                Nr_tmp(i,j,k) = Nr_tmp(i,j,k) - sed_Nr*dt_spl/dzf(k)
              end if
            end do
          end do
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
              end if
            end do
          end do
        end if ! l_lognormal
      end if ! qrbase == 1
      
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

                do s = 1, nspec
                  ! What goes to the level below
                  a_out = merge( &
                    min(rhof(k) * dzf(k) / dt_spl, sed_qr(i,j,k) / qr_spl(i,j,k)), &
                    0._field_r, &
                    qr_spl(i,j,k) > 0._field_r
                  ) * qa_spl(i,j,k,s)

                  !$acc atomic update
                  qa_spl(i,j,k,s) = qa_spl(i,j,k,s) - a_out * dt_spl / (dzf(k) * rhof(k))

                  !$acc atomic update
                  qa_spl(i,j,k-1,s) = qa_spl(i,j,k-1,s) + a_out * dt_spl / (dzf(k-1) * rhof(k-1))
                end do
              endif
            enddo
          enddo
        enddo
      endif ! l_lognormal
    end do ! time splitting loop

    qrbase = max(1, qrbase-1)
    delt_inv = 1.0_field_r / delt

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

    deallocate(qr_spl, Nr_spl, qr_tmp, Nr_tmp, qa_spl, qa_tmp) 

    call timer_toc('modbulkmicro/sedimentation_rain_sb06_aerosol')

  end subroutine sedimentation_rain_sb06_aerosol

  !> 
  subroutine evaporation_rain_sb06_aerosol(Nr, Nrp, Nrm, qr, qrp, qrm, qa, &
    &                                      qap, qt0, ql0, qvsl, rhof, exnf, & 
    &                                      Tf, esl, Dvr, mur, lbdr, xr, &
    &                                      qrbase, qrroof, qtpmcr, thlpmcr, &
    &                                      mode_inr, mode_cos, mode_acs, delt)
    real(field_r), intent(in)    :: Nr(:,:,:), qr(:,:,:), qa(:,:,:)
    real(field_r), intent(inout) :: Nrp(:,:,:), qrp(:,:,:), qap(:,:,:)
    real(field_r), intent(in)    :: Nrm(:,:,:), qrm(:,:,:)
    real(field_r), intent(in)    :: qt0(:,:,:)
    real(field_r), intent(in)    :: ql0(:,:,:)
    real(field_r), intent(in)    :: qvsl(:,:,:)
    real(field_r), intent(in)    :: rhof(:)
    real(field_r), intent(in)    :: exnf(:)
    real(field_r), intent(in)    :: Tf(:,:,:)
    real(field_r), intent(in)    :: esl(:,:,:)
    real(field_r), intent(in)    :: Dvr(:,:,:)
    real(field_r), intent(in)    :: mur(:,:,:)
    real(field_r), intent(in)    :: lbdr(:,:,:)
    real(field_r), intent(in)    :: xr(:,:,:)
    integer,       intent(in)    :: qrbase, qrroof
    real(field_r), intent(inout) :: qtpmcr(:,:,:)
    real(field_r), intent(inout) :: thlpmcr(:,:,:)
    type(mode_t),  intent(inout) :: mode_inr
    type(mode_t),  intent(inout) :: mode_cos
    type(mode_t),  intent(inout) :: mode_acs
    real(field_r), intent(in)    :: delt

    real(field_r) :: F !< Ventilation factor.
    real(field_r) :: S !< Supersaturation.
    real(field_r) :: G !< Condensation rate.
    real(field_r) :: evap, Nevap
    integer       :: i, j, k
    integer       :: numel
    real(field_r) :: f_evp
    real(field_r) :: eps
    real(field_r) :: evapt !< Tracer evaporation
    real(field_r) :: rho, E, V

    real(field_r), parameter :: sigma = 1.5
    real(field_r), parameter :: Dc
    real(field_r), save      :: Dm_fac = exp(3 * log(sigma)**2)

    !CJ: lbdr, mur, Dvr and xr are all functions of qr, Nr and rho.
    ! We could save some memory by computing these quantities on the fly.
    
    associate(qa_inr => mode_inr % aer_conc, qap_inr => mode_inr % aer_tend)

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
            G = (Rv * Tf(i,j,k)) / (Dv*esl(i,j,k)) + rlv/(Kt*Tf(i,j,k))*(rlv/(Rv*Tf(i,j,k)) -1.)
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

            ! Calculate how much of the in-rain aerosol is resuspended
            ! Based on Gong et al., 2006

            ! Could split this part into a separate kernel, if that helps
            ! with GPU registry usage
            f_evp = min(-evap / qr(i,j,k) * delt, 1._field_r)
            eps = (1 - (1 + 2 * sqrt(f_evp) + 2 * f_evp + &
                  (4._field_r / 3) * f_evp**(3._field_r / 2)) * exp(-2*sqrt(f_evp))) &
                  * (1 - f_evp) + f_evp * f_evp

            do s = 2, mode_inr % nspecies
              evapt = eps * f_evp * qa_inr(i,j,k,s) / delt
              qap_inr(i,j,k,s) = qap_inr(i,j,k,s) - evapt

              E = E + max(0.0_field_r, evapt)
              V = V + E / mode_inr % species(s) % rho
            end do

            rho = E / V
            
            Nevap = max(0._field_r, -1 * Nevap)
            Dn = 1e6 * (6 * E / (pi * Nevap * rho * 1e9))**(1._field_r / 3) &
                 * exp(-(3._field_r / 2) * log(sigma)**2)
            Dm = Dn * Dm_fac

            Fn = 0.5_field_r * erfc(-log(Dc/Dn) / (sqrt(2) * log(sigma)))
            Fm = 0.5_field_r * erfc(-log(Dc/Dm) / (sqrt(2) * log(sigma)))

            do s = 2, mode_inr % nspecies
              ! TODO: redistribute over accumulation and coarse modes
            end do

          endif
        enddo
      enddo
    enddo

    ! TODO: take care of in-rain aerosols in cells where we don't have rain

    end associate

  end subroutine evaporation_rain_sb06_aerosol

  subroutine activation

    use modmicrodata_aerosol, only: lkohler
    
    if (lkohler) then
      call activation_kohler
    else
      call activation_updraft
    end if

  end subroutine activation

  subroutine activation_kohler

    use modfields,            only: thl0, ql0, sv0, exnf
    use modglobal,            only: i1, j1, k1, rlv, cp, rhow, pi
    use modmicrodata,         only: delt
    use modmicrodata_aerosol, only: modes, nmodes, iINC, ssat
    use modmode_type,         only: mode_t

    real(field_r), parameter :: R = 8.1314, & ! Gas constant
      &                         mw = 0.018, & ! Molar mass of water
      &                         sw = 0.072    ! Surface tension of water
    real(field_r)            :: T, A, eps, kappa
    real(field_r)            :: d_crit, dn_mean, dm_mean
    real(field_r)            :: fn, fm
    real(field_r)            :: mass_total, volume_total
    integer                  :: i, j, k, s, ss, sss, imod

    do imod = 1, nmodes
      if (modes(imod) % lactivation) then ! Excluding in-cloud and in-rain
        !$acc parallel loop collapse(4) default(present) private(mass_total, volume_total, T, A, d_crit, dm_mean, fm) async(1)
        do k = 1, k1
          do j = 2, j1
            do i = 2, i1
              ! Calculate total aerosol mass and volume in mode
              do s = 2, modes(imod) % nspecies
                mass_total = mass_total + modes(imod) % aer_conc(i,j,k,s)
                volume_total = volume_total + (modes(imod) % aer_conc(i,j,k,s) / modes(imod) % species(s) % rho)
              end do

              ! Calculate kappa (PK07, eq 7)
              do s = 2, modes(imod) % nspecies
                eps = (modes(imod) % aer_conc(i,j,k,s) / modes(imod) % species(s) % rho) / volume_total
                ! Some species have zero hygroscopicity
                kappa = kappa + eps * modes(imod) % species(s) % kappa
              end do
              
              T = thl0(i,j,k) * exnf(k) + (rlv/cp) * ql0(i,j,k)

              ! Eq 10 of PK07
              A = 4._field_r * sw * mw / (R * T * rhow)

              ! The critical diameter
              ! Instead of using _field_r, use integers?
              d_crit = (4 * A**3 / (27 * kappa * log(ssat / 100 + 1)**2))**(1._field_r/3._field_r)

              ! Calculate the mean mass diameter
              dm_mean = (6 * mass_total / (pi * (modes(imod) % aer_conc(i,j,k,1) * 1.e9_field_r) * volume_total))**(1._field_r/3._field_r) &
                        * exp(3 * log(modes(imod) % sigma_g)**2)

              ! Activated fraction
              fm = 1._field_r - 0.5_field_r * erfc(-log(d_crit / dm_mean) / sqrt(2._field_r) * log(modes(imod) % sigma_g))

              ! Tendency for source aerosol - negative
              modes(imod) % aer_acti(i,j,k,s) = -fm * modes(imod) % aer_conc(i,j,k,s) / delt
              modes(imod) % aer_tend(i,j,k,s) = modes(imod) % aer_tend(i,j,k,s) + modes(imod) % aer_acti(i,j,k,s)

              ! Tendency for in-cloud aerosol - positive
              modes(iINC) % aer_acti(i,j,k,s) = modes(iINC) % aer_acti(i,j,k,s) - modes(imod) % aer_acti(i,j,k,s)
              modes(iINC) % aer_tend(i,j,k,s) = modes(iINC) % aer_tend(i,j,k,s) - modes(imod) % aer_acti(i,j,k,s)
            end do
          end do
        end do

        ! Do the same but for the number concentration
        ! CJ: separated from the mass to keep the code clean, but from a performance perspective it might be
        ! beneficial to merge the two into one kernel. Or just use async.
        !$acc parallel loop collapse(3) default(present) private(mass_total, volume_total, T, A, d_crit, dn_mean, fn) async(2)
        do k = 1, k1
          do j = 2, j1
            do i = 2, i1
              do ss = 2, modes(imod) % nspecies
                mass_total = mass_total + modes(imod) % aer_conc(i,j,k,ss)
                volume_total = volume_total + (modes(imod) % aer_conc(i,j,k,ss) / modes(imod) % species(ss) % rho)
              end do

              do ss = 2, modes(imod) % nspecies
                eps = (modes(imod) % aer_conc(i,j,k,s) / modes(imod) % species(ss) % rho) / volume_total
                if (modes(imod) % species(ss) % kappa > 0.) kappa = kappa + eps * modes(imod) % species(ss) % kappa
              end do

              T = thl0(i,j,k) * exnf(k) + (rlv/cp) * ql0(i,j,k)
              A = 4._field_r * sw * mw / (R * T * rhow)
              d_crit = (4 * A**3 / (27 * kappa * log(ssat / 100 + 1)**2))**(1._field_r/3._field_r)

              ! Number mean diameter 
              dn_mean = (6 * mass_total / (pi * (modes(imod) % aer_conc(i,j,k,1) * 1.e9) * volume_total))**(1._field_r/3._field_r)
              fn = 1._field_r - 0.5_field_r * erfc(-log(d_crit / dn_mean) / sqrt(2._field_r) * log(modes(imod) % sigma_g))

              modes(imod) % aer_tend(i,j,k,1) = modes(imod) % aer_tend(i,j,k,1) - fn * modes(imod) % aer_conc(i,j,k,1) / delt
              modes(iINC) % aer_tend(i,j,k,1) = modes(iINC) % aer_tend(i,j,k,1) - modes(imod) % aer_acti(i,j,k,1)
            end do
          end do
        end do
      end if
    end do

    !$acc wait(1,2)

  end subroutine activation_kohler

  subroutine activation_updraft

  end subroutine activation_updraft

  subroutine scavenging
    
    use modfields,            only: rhof, ql0
    use modglobal,            only: k1, j1, i1, pi, rhow
    use modmicrodata,         only: delt, eps0, qcmask, qrmask
    use modmicrodata_aerosol, only: modes, maxmodes, ncmin, Nc, sed_qr_dup

    real(field_r) :: mass_total, volume_total !< Mode mean quantities
    real(field_r) :: meancld, rainrate
    real(field_r) :: dcdt_inc_n, dcdt_inc_m   !< In-cloud scavenging
    real(field_r) :: dcdt_blc_n, dcdt_blc_m   !< Below-cloud (rain) scavenging
    real(field_r) :: rm, rm_inc, rm_blc       !< Various radii
    
    integer :: i, j, k, s, ss, imod
    integer :: t_inc, t_inr !< Index of target in-cloud and in-rain modes

    do imod = 1, maxmodes - 2 ! Exclude in-cloud and in-rain aerosols
      if (.not. modes(imod) % enabled ) cycle
        do k = min(qcbase, qrbase), max(qcroof, qrroof)
          do j = 2, j1
            do i = 2, i1
              ! TODO: what is this?
              meancld = 1e6_field_r * (3._field_r * ql0(i,j,k) * rhof(k) / (4._field_r * pi * Nc(i,j,k) * rhow + eps0))**(1._field_r/3._field_r)
              meancld = min(max(meancld, 5.001), 49.999)

              do s = 2, modes(imod) % nspecies
                mass_total = mass_total + modes(imod) % aer_conc(i,j,k,s)
                volume_total = volume_total + (modes(imod) % aer_conc(i,j,k,s) / modes(imod) % species(s) % rho)
              end do

              ! Geometric mean aerosol radius
              ! TODO: CJ: check of deze formule klopt
              rm = 0.5_field_r * (6 * mass_total / (pi * modes(imod) % aer_conc(i,j,k,1) * 1.e9_field_r * volume_total))**(1._field_r/3._field_r) &
                   * exp(-(3._field_r/2._field_r))

              if (qcmask(i,j,k) .and. Nc(i,j,k) > ncmin) then
                ! Make sure rm is within bounds of the lookup table
                rm_inc = max(min(100._field_r * rm, 8.0e-3_field_r), 1.0e-8_field_r)
                
                ! In-cloud scavenging
                dcdt_inc_n = 1e-6_field_r * Nc(i,j,k) * &
                             lookup_interpolate(ncld, logcldrad, naer_inc, logaerrad, incnumb, log(meancld), log(rm))
                dcdt_inc_m = 1e-6_field_r * Nc(i,j,k) * & 
                             lookup_interpolate(ncld, logcldrad, naer_inc, logaerrad, incmass, log(meancld), log(rm))

                dcdt_inc_n = merge(1._field_r / delt, dcdt_inc_n, (dcdt_inc_n * delt > 1._field_r) .or. (dcdt_inc_m * delt > 1._field_r))
                dcdt_inc_m = merge(1._field_r / delt, dcdt_inc_m, (dcdt_inc_n * delt > 1._field_r) .or. (dcdt_inc_m * delt > 1._field_r))

                modes(imod) % aer_tend(i,j,k,1) = modes(imod) % aer_tend(i,j,k,1) - dcdt_inc_n * max(0._field_r, modes(imod) % aer_conc(i,j,k,1))
                modes(imod) % aer_tend(i,j,k,1) = modes(imod) % aer_tend(i,j,k,1) + dcdt_inc_n * max(0._field_r, modes(imod) % aer_conc(i,j,k,1))

                do s = 2, modes(imod) % nspecies
                  t_inc = modes(imod) % species(s) % iINC
                  modes(imod) % aer_tend(i,j,k,s) = modes(imod) % aer_tend(i,j,k,s) - dcdt_inc_m * max(0._field_r, modes(imod) % aer_conc(i,j,k,s))
                  modes(iINC) % aer_tend(i,j,k,t_inc) = modes(iINC) % aer_tend(i,j,k,t_inc) + dcdt_inc_m * max(0._field_r, modes(imod) % aer_conc(i,j,k,s))
                end do
              end if
                
              ! GPU: probably too much thread divergence? 
              if (qrmask(i,j,k) .and. sed_qr_dup(i,j,k) * 3600_field_r > 0.01_field_r) then
                ! Below-cloud scavenging
                rm_blc = max(min(.9999e3_field_r, rm * 1.e6_field_r), 1.001e-3_field_r)

                dcdt_blc_n = lookup_interpolate(nrai, lograinrate, naer_blc, logaerrad_blc, blcnumb, log(rainrate), log(rm))
                dcdt_blc_m = lookup_interpolate(nrai, lograinrate, naer_blc, logaerrad_blc, blcmass, log(rainrate), log(rm))

                do s = 2, modes(imod) % nspecies
                  t_inr = modes(imod) % species(s) % iINR
                  modes(imod) % aer_tend(i,j,k,s) = modes(imod) % aer_tend(i,j,k,s) - dcdt_blc_m * max(0._field_r, modes(imod) % aer_conc(i,j,k,s))
                  modes(imod) % aer_tend(i,j,k,t_inr) = modes(imod) % aer_tend(i,j,k,t_inr) - dcdt_blc_m * max(0._field_r, modes(imod) % aer_conc(i,j,k,s))
                end do
              end if 
            end do
          end do
        end do
      end do
  end subroutine scavenging

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

  function binarysearch(length, array, value, delta)
     ! Given an array and a value, returns the index of the element
     ! that
     ! is closest to, but less than, the given value.
     ! Uses a binary search algorithm.
     ! "delta" is the tolerance used to determine if two values are
     ! equal
     ! if ( abs(x1 - x2) <= delta) then
     !    assume x1 = x2
     ! endif

     implicit none
     integer, intent(in) :: length
     real(field_r), dimension(length), intent(in) :: array
     real(field_r), intent(in) :: value
     real(field_r), intent(in), optional :: delta

     integer :: binarysearch

     integer :: left, middle, right
     real(field_r) :: d

     if (present(delta) .eqv. .true.) then
         d = delta
     else
         d = 1e-9
     endif

     left = 1
     right = length
     do
         if (left > right) then
             exit
         endif
         middle = nint((left+right) / 2.0)
         if ( abs(array(middle) - value) <= d) then
             binarySearch = middle
             return
         else if (array(middle) > value) then
             right = middle - 1
         else
             left = middle + 1
         end if
     end do
     binarysearch = right

    end function binarysearch

    real function lookup_interpolate(x_len, x_array, y_len, y_array, f, x, y, delta)
        ! This function uses bilinear interpolation to estimate the
        ! value of a function f at point (x,y)
        ! f is assumed to be sampled on a regular grid, with the grid x
        ! values specified by x_array and the grid y values specified by y_array
        ! Reference: http://en.wikipedia.org/wiki/Bilinear_interpolation
        implicit none
        integer, intent(in) :: x_len, y_len
        real(field_r), dimension(x_len), intent(in) :: x_array
        real(field_r), dimension(y_len), intent(in) :: y_array
        real(field_r), dimension(x_len, y_len), intent(in) :: f
        real(field_r), intent(in) :: x,y
        real(field_r), intent(in), optional :: delta

        real(field_r) :: denom, x1, x2, y1, y2
        integer :: i,j

        i = binarysearch(x_len, x_array, x)
        j = binarysearch(y_len, y_array, y)

        x1 = x_array(i)
        x2 = x_array(i+1)

        y1 = y_array(j)
        y2 = y_array(j+1)

        denom = (x2 - x1)*(y2 - y1)

        lookup_interpolate = (f(i,j)*(x2-x)*(y2-y) + f(i+1,j)*(x-x1)*(y2-y) + &
            f(i,j+1)*(x2-x)*(y-y1) + f(i+1, j+1)*(x-x1)*(y-y1))/denom

    end function lookup_interpolate

end module modbulkmicro_aerosol
