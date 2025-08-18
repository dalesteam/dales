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
  use modglobal,    only: dzf, ih, jh, i1, j1, pi, rhow, rlv, cp, &
                          ifnamopt, checknamelisterror
  use modprecision, only : field_r
  use modtimer,     only: timer_tic, timer_toc
  use modmicrodata, only: Nc_0, sig_g, qtpmcr, thlpmcr, l_rain, l_homogenize
  use modbulkmicro_data, only: qrbase, qrroof, qcbase, qcroof, qcmin, l_sb, &
                          l_sedc, l_mur_cst, l_lognormal, mur_cst, &
                          sig_gr, c_St, mygamma21, mygamma251
  use bulkmicro_sb, only: autoconversion_sb, &
                          accretion_sb, evaporation_sb, sedimentation_rain_sb
  use bulkmicro_kk, only: autoconversion_kk, &
                          accretion_kk, evaporation_kk, sedimentation_rain_kk
  use modbulkmicro_stat, only: init_bulkmicro_stat, bulkmicro_stat
  use modmpi, only: myid, D_MPI_BCAST, comm3d, mpierr, print_info_stderr
  implicit none
  private

  character(len=*), parameter :: modname = 'modbulkmicro'

  public initbulkmicro, exitbulkmicro, bulkmicro

  real :: gamma25
  real :: gamma3
  real :: gamma35
  contains

!> Initializes and allocates the arrays
  subroutine initbulkmicro
    use modglobal, only : i1,j1,k1,ih,jh
    use modmicrodata, only: iqr, inr, lstat, precep
    use modbulkmicro_data, only : Nr, Nrp, qr, qrp 
    use modtracers,   only: add_tracer
    implicit none

    integer :: m

    ! Setup two tracers for precipitation
    call add_tracer("qr", long_name="rain water mixing ratio", &
                    unit="kg/kg", lmicro=.true., isv=iqr)

    call add_tracer("Nr", long_name="rain droplet number concentration", &
                    unit="1/m^3", lmicro=.true., isv=inr)

                                        ! Fields accessed by:
    allocate(Nr       (2:i1,2:j1,k1)  & ! dobulkmicrostat, dosimpleicestat
            ,qr       (2:i1,2:j1,k1)  & ! dobulkmicrostat, dosimpleicestat
            ,Nrp      (2:i1,2:j1,k1)  & ! bulkmicrotend, simpleicetend
            ,qrp      (2:i1,2:j1,k1)  & ! bulkmicrotend, simpleicetend
            ,precep   (2:i1,2:j1,k1)  ) ! dobulkmicrostat, dosimpleicestat, docape

    allocate(thlpmcr  (2:i1,2:j1,k1)  & !
            ,qtpmcr(2-ih:i1+ih,2-jh:j1+jh,k1))  ! ghost cells added here for modvarbudget

    precep = 0
    gamma25=gamma(2.5)
    gamma3=2.
    gamma35=gamma(3.5)

    ! Setup lookup tables for ventilation factor in SB evaporation.
    ! Entries are computed in double precision, and stored in field_r precision.
    ! TODO: on GPU, it might be faster to inline the computation.
    if (l_sb) then
      mygamma21(-100) = 0.0
      mygamma251(-100) = 0.0

      do m = -99, 4000
        mygamma21(m) = max(0.0_field_r, &
          gamma(m/100.0 + 2.0) / gamma(m/100.0 + 1.0) &
          * (((m/100.0 + 3.0)*(m/100.0 + 2.0)*(m/100.0 + 1.0))**(-1/3.0)))
        mygamma251(m) = max(0.0, &
          gamma(m/100.0 + 2.5) / gamma(m/100.0 + 1.0) &
          * (((m/100.0 + 3.0)*(m/100.0 + 2.0)*(m/100.0 + 1.0))**(-1/2.0)))
      end do

      !$acc update device(mygamma21, mygamma251)
    end if

    !$acc enter data copyin(Nr, qr, Nrp, qrp, precep, thlpmcr, qtpmcr)

    if (lstat) call init_bulkmicro_stat

  end subroutine initbulkmicro

!> Cleaning up after the run
  subroutine exitbulkmicro
  !*********************************************************************
  ! subroutine exitbulkmicro
  !*********************************************************************
    use modmicrodata,      only : precep, qtpmcr, thlpmcr
    use modbulkmicro_data, only : Nr,Nrp,qr,qrp
    implicit none

    !$acc exit data delete(Nr, qr, Nrp, qrp, precep, thlpmcr, qtpmcr)

    deallocate(Nr,Nrp,qr,qrp,thlpmcr,qtpmcr)
    deallocate(precep)

  end subroutine exitbulkmicro

!> Calculates the microphysical source term.
  subroutine bulkmicro
    use modglobal, only : i1,j1,kmax,k1,rdt,rk3step,timee,rlv,cp,dzf,ijtot
    use modfields, only : sv0,svm,svp,qtp,thlp,ql0,exnf,rhof, esl, qt0, qvsl, tmp0
    use modbulkmicrostat, only : bulkmicrotend
    use modmpi,    only : myid, slabsum
    use modbulkmicro_data, only : Nr, qr, Nrp, qrp,  &
                             l_sedc, l_mur_cst, l_lognormal, l_rain, &
                             qrmin, qcmin, &
                             mur_cst, l_sb
    use modmicrodata, only: iqr, inr, lstat, precep, delt, qtpmcr, thlpmcr
    use modmicroutil, only: zero_field, sum_fields
    use modstat_profiles, only: sample_field
    implicit none
    integer :: i, j, k
    real :: qrtest,nr_cor,qr_cor
    real :: qrsum_neg, qrsum, Nrsum_neg, Nrsum
    real(field_r), allocatable :: qrp_tmp(:,:,:), nrp_tmp(:,:,:)
    real(field_r), allocatable :: qtp_tmp(:,:,:), thlp_tmp(:,:,:) ! temp fields, for homogenized evaporation
    real(field_r) :: qtpevap(1:k1), thlpevap(1:k1)
    !$acc enter data create(qtpevap, thlpevap)

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
          if (qr(i,j,k) > qrmin) then
            qrbase = min(k, qrbase)
          endif
          if (ql0(i,j,k) > qcmin) then
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
            if (qr(i,j,k) > qrmin) then
              qrroof = max(k, qrroof)
            endif
            if (ql0(i,j,k) > qcmin) then
              qcroof = max(k, qcroof)
            endif
          enddo
        enddo
      enddo
      qrroof = min(k1, qrroof)
      qcroof = min(k1, qcroof)
    endif

    ! if there is nothing to do, we can return at this point
    ! if (min(qrbase,qcbase).gt.max(qrroof,qcroof)) return

    if (l_sedc) then
      call sedimentation_cloud(ql0, rhof, exnf, qcbase, qcroof, qtpmcr, thlpmcr)
      if(lstat) call sample_field('qtpsedc', qtpmcr) ! First process, no need to zero beforehand
    endif

    ! Rain processes
    if (l_rain) then
      allocate(qrp_tmp(2:i1,2:j1,1:k1), nrp_tmp(2:i1,2:j1,1:k1))

      !$acc enter data create(qrp_tmp, nrp_tmp)

      call zero_field(qrp_tmp)
      call zero_field(nrp_tmp)

      ! 1. Autoconversion
      if (l_sb) then
        call autoconversion_sb(ql0, qr, exnf, rhof, qcbase, qcroof, thlpmcr, &
                               qtpmcr, qrp_tmp, Nrp_tmp)
      else
        call autoconversion_kk(ql0, rhof, exnf, qcbase, qcroof, thlpmcr, &
                               qtpmcr, qrp_tmp, Nrp_tmp)
      end if

      if (lstat) then
        call sample_field('qrpauto', qrp_tmp)
        call sample_field('npauto', nrp_tmp)
      end if

      call sum_fields(qrp_tmp, qrp)
      call sum_fields(nrp_tmp, nrp)

      call zero_field(qrp_tmp)
      call zero_field(nrp_tmp)

      ! 2. Accretion
      if (l_sb) then
        call accretion_sb(ql0, qr, Nr, exnf, rhof, qcbase, qcroof, qrbase, qrroof, &
                          thlpmcr, qtpmcr, qrp_tmp, Nrp_tmp)
      else
        call accretion_kk(ql0, qr, exnf, qcbase, qcroof, qrbase, qrroof, &
                          thlpmcr, qtpmcr, qrp_tmp)
      end if

      if (lstat) then
        call sample_field('qrpaccr', qrp_tmp)
        call sample_field('npaccr', nrp_tmp)
      end if

      call sum_fields(qrp_tmp, qrp)
      call sum_fields(nrp_tmp, nrp)

      call zero_field(qrp_tmp)
      call zero_field(nrp_tmp)

      ! 3. Evaporation
      if (l_homogenize) then
         ! qt, thl tendencies from rain evaporation are homogenized horizontally
         allocate(qtp_tmp(2:i1,2:j1,1:k1), thlp_tmp(2:i1,2:j1,1:k1))
         !$acc enter data create(qtp_tmp, thlp_tmp)
         call zero_field(qtp_tmp)
         call zero_field(thlp_tmp)
         if(l_sb) then
            call evaporation_sb(ql0, qt0, svm(:,:,:,iqr), svm(:,:,:,inr), qvsl, tmp0, &
                                esl, exnf, rhof, Nr, qr, qrbase, qrroof, &
                                qrp_tmp, Nrp_tmp, delt, qtp_tmp, thlp_tmp)
         else
            call evaporation_kk(ql0, qt0, qvsl, esl, tmp0, svm(:,:,:,iqr), svm(:,:,:,iNr), &
                                Nr, qr, rhof, exnf, qrbase, qrroof, delt, &
                                thlp_tmp, qtp_tmp, qrp_tmp, Nrp_tmp)
         end if
         qtpevap = 0
         thlpevap = 0
         call slabsum(qtpevap,1,k1,qtp_tmp,2,i1,2,j1,1,k1,2,i1,2,j1,1,k1)
         call slabsum(thlpevap,1,k1,thlp_tmp,2,i1,2,j1,1,k1,2,i1,2,j1,1,k1)
         qtpevap = qtpevap/ijtot
         thlpevap = thlpevap/ijtot
         !$acc exit data delete(qtp_tmp, thlp_tmp)
         deallocate(qtp_tmp, thlp_tmp)

         !$acc parallel loop collapse(3) default(present)
         do k=1,k1   ! note: must use full k-range here, other MPI tiles may have a wider/different range in qrbase...qrroof
            do j=2,j1
               do i=2,i1
                  qtpmcr(i,j,k) = qtpmcr(i,j,k) + qtpevap(k)
                  thlpmcr(i,j,k) = thlpmcr(i,j,k) + thlpevap(k)
               end do
            end do
         end do
      else ! l_homogenize is false
         if(l_sb) then
            call evaporation_sb(ql0, qt0, svm(:,:,:,iqr), svm(:,:,:,inr), qvsl, tmp0, &
                                esl, exnf, rhof, Nr, qr, qrbase, qrroof, &
                                qrp_tmp, Nrp_tmp, delt, qtpmcr, thlpmcr)
         else
            call evaporation_kk(ql0, qt0, qvsl, esl, tmp0, svm(:,:,:,iqr), svm(:,:,:,iNr), &
                                Nr, qr, rhof, exnf, qrbase, qrroof, delt, &
                                thlpmcr, qtpmcr, qrp_tmp, Nrp_tmp)
         end if
      end if

      if (lstat) then
        call sample_field('qrpevap', qrp_tmp)
        call sample_field('npevap', nrp_tmp)
      end if

      call sum_fields(qrp_tmp, qrp)
      call sum_fields(nrp_tmp, nrp)

      call zero_field(qrp_tmp)
      call zero_field(nrp_tmp)

      ! 4. Sedimentation
      if (l_sb) then
        call sedimentation_rain_sb(qr, Nr, rhof, dzf, qrbase, qrroof, &
                                   l_lognormal, delt, qrp_tmp, Nrp_tmp, precep)
      else
        call sedimentation_rain_kk(qr, Nr, rhof, dzf, qrbase, qrroof, delt, &
                                   qrp_tmp, Nrp_tmp, precep)
      end if

      if (lstat) then
        call sample_field('qrpsed', qrp_tmp)
        call sample_field('npsed', nrp_tmp)
      end if

      call sum_fields(qrp_tmp, qrp)
      call sum_fields(nrp_tmp, nrp)

      call zero_field(qrp_tmp)
      call zero_field(nrp_tmp)

    end if

    !*********************************************************************
    ! remove negative values and non physical low values
    !*********************************************************************
    !$acc parallel loop collapse(3) default(present) private(qr_cor, Nr_cor)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          qr_cor = min(svp(i,j,k,iqr) + qrp(i,j,k) + (svm(i,j,k,iqr) / delt), &
                       0.0_field_r)
          Nr_cor = min(svp(i,j,k,iNr) + Nrp(i,j,k) + (svm(i,j,k,iNr) / delt), &
                       0.0_field_r)

          qrp_tmp(i,j,k) = - qr_cor
          Nrp_tmp(i,j,k) = - Nr_cor
        end do
      end do
    end do

    if (lstat) then
      call sample_field('qrpclip', qrp_tmp)
      call sample_field('npclip', nrp_tmp)
    end if

    call sum_fields(qrp_tmp, qrp)
    call sum_fields(nrp_tmp, nrp)

    if (lstat) then
      call sample_field('qrptot', qrp)
      call sample_field('nptot', nrp)
    end if

    !$acc parallel loop collapse(3) default(present)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          qtp (i,j,k) = qtp (i,j,k) + qtpmcr (i,j,k)
          thlp(i,j,k) = thlp(i,j,k) + thlpmcr(i,j,k)

          svp(i,j,k,iqr) = svp(i,j,k,iqr) + qrp(i,j,k)
          svp(i,j,k,inr) = svp(i,j,k,inr) + Nrp(i,j,k)
        enddo
      enddo
    enddo

    !$acc exit data delete(qrp_tmp, nrp_tmp)
    deallocate(qrp_tmp, nrp_tmp)

    if (lstat) call bulkmicro_stat

    !$acc exit data delete(qtpevap, thlpevap)
  end subroutine bulkmicro

  !> Sedimentation of cloud water ((Bretherton et al,GRL 2007))
  !!
  !!   The sedimentation of cloud droplets assumes a lognormal DSD in which the
  !!   geometric std dev. is assumed to be fixed at 1.3.
  !! sedimentation of cloud droplets
  !! lognormal CDSD is assumed (1 free parameter : sig_g)
  !! terminal velocity : Stokes velocity is assumed (v(D) ~ D^2)
  !! flux is calc. anal.
  subroutine sedimentation_cloud(ql, rhof, exnf, qcbase, qcroof, qtpmcr, thlpmcr)

    real(field_r), intent(in)    :: ql(2:,2:,:)
    real(field_r), intent(in)    :: rhof(:)
    real(field_r), intent(in)    :: exnf(:)
    integer,       intent(in)    :: qcbase, qcroof

    real(field_r), intent(inout) :: qtpmcr(2-ih:,2-jh:,:)
    real(field_r), intent(inout) :: thlpmcr(2:,2:,:)

    character(len=*), parameter :: routine = modname//'/sedimentation_cloud'

    integer       :: i, j, k
    real(field_r) :: csed
    real(field_r) :: sedc

    call timer_tic(routine, 1)

    if (qcbase > qcroof) return

    csed = c_St*(3./(4.*pi*rhow))**(2./3.)*exp(5.*log(sig_g)**2.)

    !$acc parallel loop collapse(3) default(present)
    do k = qcbase, qcroof
      do j = 2, j1
        do i = 2, i1
          if (ql(i,j,k) > qcmin) then
            sedc = csed*Nc_0**(-2./3.)*(ql(i,j,k)*rhof(k))**(5./3.)

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

    call timer_toc(routine)

  end subroutine sedimentation_cloud

end module modbulkmicro
