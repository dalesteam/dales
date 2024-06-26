module modbulkmicro3_column
  use modprecision, only: field_r
  use modmicrodata
  use modmicrodata3
  implicit none
  private
  public column_processes, nucleation3

contains

! Notes on sedimentation calculation:
!
! Some parts of the sedimentation calculations cannot be vectorized (fi, x^y),
! and the many different fields, parameters, etc. make the code cache-unfriendly.
! To prevent calculating sedimentation when there is nothing happening,
! we adjust the boundaries of the two loops over the vertical.
!
! Sketch of the sedimentation algorithm:
!
! set the boundaries of loop1 to the full vertical k_low1=1,k_high1=k1
! Time splitting loop to make sure things fall only one level per timestep
!   k loop1: boundaries k_low1, k_high1
!     calculate the sedimentation (amount, velocity, ...)
!     extend the boundaries of loop 2 to only include the active part
!
!   short-cut: if there is nothing to sediment, immediately return
!
!   k loop2: boundaries k_low2, k_high2
!     have the sedimentation actually fall down a level, also k_low2 => k_low2 - 1
!     set boundaries of loop 1 to cover the updated levels with activity
! update tendencies

subroutine column_processes(sv0, svp, thlpmcr, qtpmcr, &
                            precep_hr,precep_ci,precep_hs,precep_hg,&
                            tend)
  use modglobal, only : k1
  implicit none
  real, intent(in),    dimension(ncols,k1)     :: sv0
  real, intent(inout), dimension(ncols,k1)     :: svp

  real, intent(inout), dimension(k1)        :: thlpmcr, qtpmcr
  real, intent(out)                         :: precep_hr, precep_ci, precep_hs, precep_hg
  real, intent(out),   dimension(ntends,k1) :: tend

  ! sedimentation
  ! -----------------------------------------------------------------
  call sedim_rain3(sv0(iq_hr,1:k1), sv0(in_hr,1:k1) &
                  ,svp(iq_hr,1:k1), svp(in_hr,1:k1) &
                  ,precep_hr,tend)

  call sedim_cl3(sv0(iq_cl,1:k1), sv0(in_cl,1:k1) &
                ,svp(iq_cl,1:k1), svp(in_cl,1:k1) &
                ,svp(in_cc,1:k1)                  &
                ,qtpmcr,thlpmcr,tend)

  call sedim_ice3(sv0(iq_ci,1:k1), sv0(in_ci,1:k1) &
                 ,svp(iq_ci,1:k1), svp(in_ci,1:k1) &
                 ,precep_ci,tend)

  call sedim_snow3(sv0(iq_hs,1:k1), sv0(in_hs,1:k1) &
                  ,svp(iq_hs,1:k1), svp(in_hs,1:k1) &
                  ,precep_hs,tend)

  call sedim_graupel3(sv0(iq_hg,1:k1), sv0(in_hg,1:k1) &
                     ,svp(iq_hg,1:k1), svp(in_hg,1:k1) &
                     ,precep_hg,tend)

end subroutine column_processes


!> Cloud nucleation
!! Written to prognostically evaluate the cloud water number content [ kg^{-1}]
!! directly follows Seifert&Beheng scheme
subroutine nucleation3(qt0, qvsl, w0, q_cl, q_clp, n_cl, n_clp, n_cc, statistics,tend)
  use modglobal, only : dzf,k1
  implicit none
  real, intent(in)    :: qt0(k1), qvsl(k1), w0(k1)
  real, intent(in)    :: q_cl(k1), n_cl(k1), n_cc(k1)
  real, intent(inout) :: q_clp(k1), n_clp(k1)
  real, intent(out)   :: tend(ntends, k1)
  real, intent(inout) :: statistics(nmphys,k1)

  real    :: dn_cl_nu(k1)  !< droplet nucleation rate
  integer :: k
  real    :: coef_ccn, n_act

  ! note that supersaturation is
  ! not always how supersaturated is water vapour,
  ! depending on a flag, it can also include water already in droplets
  real :: ssat_u      & ! supersaturation at (...,k+1)
         ,ssat(k1)    & !                 at (...,k)
         ,ssat_d      & !                 at (...,k-1)
         ,wdssatdz(k1)  ! derivation of supersaturation

  dn_cl_nu = 0.

  coef_ccn  = 1.0/sat_max**kappa_ccn ! 1.0
  ! allows to keep both definitions consistent

  ! calculating supersaturation
  if (l_sb_nuc_sat) then
    ! calculating supersaturation of water vapour only
    ssat = (100./qvsl)*(qt0-q_cl-qvsl)
  else ! l_sb_nuc_sat
    ! ie. cloud liquid water is also included supersaturation
    ssat = (100./qvsl)*(qt0-qvsl)
  endif ! l_sb_nuc_sat

  ! calculating the derivation - second order estimation?
  ! ? add switches for different derivation calculation? - so first the first order only
  ! option A: central differences lax

  ! just an approximation - same as for the second level,
  ! or 0 to prevent condensation there
  wdssatdz(1) = w0(2)*(ssat(2)-ssat(1))/dzf(2)
  do k=2,k1 - 1
    wdssatdz = 0.5*(w0(k+1)+w0(k))*(ssat(k+1) - ssat(k-1))/(dzf(k)+dzf(k-1))
  enddo
  ! first order approximation of the difference
  wdssatdz(k1) = w0(k1)*(ssat(k1) - ssat(k1-1)) / dzf(k1-1)

  do k=1,k1
    ! NOTE: original code went out of bounds for k=1 and k=k1
    !       for now, take d ssat / dz = 0 at boundaries
    !       Altough we dont expect nucleation at top and or lowest level,
    !       to be consistent, just process them anyways.
    if (k.gt.1) then
      ssat_d = ssat(k-1)
    else
      ssat_d = ssat(1)
    endif
    if (k.eq.k1) then
      ssat_u = ssat(k1)
    else
      ssat_u = ssat(k+1)
    endif

    if (l_sb_nuc_expl) then ! calculation of explicit nucleation

      if(ssat(k).gt.sat_min) then ! of course only in cloud

        if (l_c_ccn) then ! c_ccn is constant

          if (l_sb_nuc_diff) then

            if (l_sb_sat_max) then

              if ((ssat(k).lt.sat_max).and.(wdssatdz(k).gt.0.0)) then ! condition for nucleation
                dn_cl_nu(k) = c_ccn*kappa_ccn*wdssatdz(k)*ssat(k)**(kappa_ccn-1.0) ! (1/rhof) *rhof = 1
              endif

            else ! NOT l_sb_sat_max  AND l_sb_nuc_diff

              if (wdssatdz(k).gt.0.0) then ! condition for nucleation
                dn_cl_nu(k) = c_ccn*kappa_ccn*wdssatdz(k)*ssat(k)**(kappa_ccn-1.0)  ! (1/rhof) *rhof = 1
              endif

            endif ! l_sb_sat_max

          else ! NOT l_sb_nuc_diff

            if (l_sb_sat_max) then ! l_sb_sat_max AND NOT l_sb_nuc_diff

              if ((ssat(k).lt.sat_max).and.(w0(k).gt.0.0)) then ! condition for nucleation
                dn_cl_nu(k) = (c_ccn/dzf(k-1))*max(0.0,w0(k)*      &
                   (ssat(k)**kappa_ccn-ssat_d**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif
              if ((ssat(k).lt.sat_max).and.(w0(k+1).lt.0.0)) then ! condition for nucleation
                dn_cl_nu(k) = (c_ccn/dzf(k-1))*max(0.0,w0(k+1)*    &
                   (ssat(k)**kappa_ccn-ssat_u**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif

            else ! NOT l_sb_sat_max AND NOT l_sb_nuc_diff

              if (w0(k).gt.0.0) then ! condition for nucleation
                dn_cl_nu(k) = (c_ccn/dzf(k-1))*max(0.0,w0(k)*      &
                   (ssat(k)**kappa_ccn-ssat_d**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif
              if (w0(k+1).lt.0.0) then ! condition for nucleation
                dn_cl_nu(k) = (c_ccn/dzf(k-1))*max(0.0,w0(k+1)*    &
                   (ssat(k)**kappa_ccn-ssat_u**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif

            endif  ! NOT l_sb_sat_max  AND NOT l_sb_nuc_diff

          endif ! NOT l_sb_nuc_diff

          ! basic limiting
          dn_cl_nu(k) = max(0.0,min(dn_cl_nu(k),(n_clmax-n_cl(k))/delt))

        else ! c_ccn is not constant, but dependent on n_cc

          if (l_sb_nuc_diff) then

            if (l_sb_sat_max) then ! l_sb_sat_max

              if ((ssat(k).lt.sat_max).and.(wdssatdz(k).gt.0.0)) then ! condition for nucleation
                dn_cl_nu(k) = coef_ccn*n_cc(k)*kappa_ccn*wdssatdz(k)*ssat(k)**(kappa_ccn-1.0)  ! (1/rhof) *rhof = 1
              endif

            else  ! l_sb_sat_max

              if (wdssatdz(k).gt.0.0) then ! condition for nucleation
                dn_cl_nu(k) = coef_ccn*n_cc(k)*kappa_ccn*wdssatdz(k)*ssat(k)**(kappa_ccn-1.0)  ! (1/rhof) *rhof = 1
              endif

            endif  ! l_sb_sat_max

          else  ! l_sb_nuc_diff

            if (l_sb_sat_max)  then ! l_sb_sat_max  AND NOT l_sb_nuc_diff

              if ((ssat(k).lt.sat_max).and.(w0(k).gt.0.0)) then ! condition for nucleation
                dn_cl_nu(k) = (coef_ccn*n_cc(k)/dzf(k-1))*max(0.0,w0(k)*   &
                  (ssat(k)**kappa_ccn-ssat_d**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif
              if ((ssat(k).lt.sat_max).and.(w0(k+1).lt.0.0)) then ! condition for nucleation
                dn_cl_nu(k) = (coef_ccn*n_cc(k)/dzf(k-1))*max(0.0,w0(k+1)* &
                   (ssat(k)**kappa_ccn-ssat_u**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif

            else ! NOT l_sb_sat_max  AND NOT l_sb_nuc_diff

              if (w0(k).gt.0.0) then ! condition for nucleation
                dn_cl_nu(k) = (coef_ccn*n_cc(k)/dzf(k-1))*max(0.0,w0(k)*   &
                    (ssat(k)**kappa_ccn-ssat_d**kappa_ccn))   ! (1/rhof) *rhof = 1
              endif
              if (w0(k+1).lt.0.0) then ! condition for nucleation
                dn_cl_nu(k) = (coef_ccn*n_cc(k)/dzf(k-1))*max(0.0,w0(k+1)* &
                   (ssat(k)**kappa_ccn-ssat_u**kappa_ccn))    ! (1/rhof) *rhof = 1
              endif

            endif  ! l_sb_sat_max

          endif  ! l_sb_nuc_diff

          ! basic limiting
          dn_cl_nu(k) = max(0.0,min(dn_cl_nu(k),(n_cc(k)-n_cl(k))/delt))

        endif ! l_c_ccn

      endif ! ssat.gt.sat_min

    else ! l_sb_nuc_expl

      if(ssat(k).gt.0.0) then

        ! calculate number of activated n_ccn
        n_act = coef_ccn*n_cc(k)*min(sat_max,ssat(k))**kappa_ccn
        n_act = max(n_cc(k),n_act)

        if (n_act.gt.n_cl(k)) then
          dn_cl_nu(k) = (n_act-n_cl(k))/delt
        endif

        ! basic limiting - not needed in this case
        ! dn_cl_nu = min(dn_cl_nu, (n_cc - n_cl))/delt)

      endif ! ssat.gt.0.0

    endif ! l_sb_nuc_expl

    !if (l_sb_dbg) then
    !  ! warning if too high liquid water
    !  if (ssat.gt.0.0 .and.( &
    !  (qt0-qvsl-svm(i,j,k,iq_cl)-svm(i,j,k,iq_ci))/delt-x_cnuc*dn_cl_nu.lt. 0.)
    !  ) then
    !    write(6,*) 'WARNING: cloud nucleation too high'
    !    write(6,*) ' removing too much water'
    !    ! count(ssat.gt.0.0 .and.( &
    !    ! (qt0-qvsl-svm(i,j,k,iq_cl)-svm(i,j,k,iq_ci))/delt-x_cnuc*dn_cl_nu.lt.0.))
    !  end if
    !endif ! l_sb_dbg
  enddo

  ! increase in cloud water number
  n_clp = n_clp + dn_cl_nu

  ! update water density [kg kg^{-1}]
  q_clp = q_clp + x_cnuc * dn_cl_nu

  if (l_tendencies) then
    tend(idn_cl_nu,:) = dn_cl_nu(:)
  endif
  if (l_statistics) then
    statistics(imphys_cond,:) = statistics(imphys_cond,:) + x_cnuc * dn_cl_nu(:)
  endif
end subroutine  nucleation3


!> Sedimentaion of rain
!! sedimentation of drizzle water
!! - gen. gamma distr is assumed. Terminal velocities param according to
!!   Stevens & Seifert. Flux are calc. anal.
!! - l_lognormal =T : lognormal DSD is assumed with D_g and N known and
!!   sig_g assumed. Flux are calc. numerically with help of a
!!   polynomial function
subroutine sedim_rain3(q_hr, n_hr, q_hrp, n_hrp, precep_hr, tend)
  use modglobal, only : k1,kmax,eps1,dzf
  use modfields, only : rhof
  implicit none
  real, intent(in)     :: q_hr(k1), n_hr(k1)
  real, intent(inout)  :: q_hrp(k1), n_hrp(k1)
  real, intent(out)    :: precep_hr
  real, intent(out)    :: tend(ntends, k1)

  integer :: k,jn,n_spl
  integer :: k_low1,k_high1  & ! Boundaries for the first k loop
            ,k_low2,k_high2    ! Boundaries for the second k loop

  real :: wvar        &!< work variable
         ,dt_spl      &!<
         ,xr_spl      &!< for time splitting
         ,Dvr_spl     &!<     -
         ,mur_spl     &!<     -
         ,lbdr_spl    &!<     -
         ,Dgr         &!< lognormal geometric diameter
         ,N_r0        &!< rain integral stuff
         ,pwcont       !<

  real :: qr_spl(k1), Nr_spl(k1)
  real :: sed_qr(k1), sed_Nr(k1)

  real :: wfall

  n_spl = ceiling(split_factor*wfallmax_hr*delt/(minval(dzf)))
  dt_spl = delt/real(n_spl)

  qr_spl = q_hr
  Nr_spl = n_hr

  ! First time for the first k loop, go over the full column
  k_low1 = 1
  k_high1 = k1

  ! The boundaries for the second k loop are extended in the first loop
  k_low2 = k1
  k_high2 = 0

  do jn = 1 , n_spl ! time splitting loop

    sed_qr = 0.
    sed_Nr = 0.

    do k=k_low1,k_high1
      ! NOTE: code is a duplicate from subroutine integrals_bulk3
      if (qr_spl(k) > q_hr_min.and.(Nr_spl(k) > 0.0)) then
        if (l_sb) then
          xr_spl   = qr_spl(k)/Nr_spl(k)
          xr_spl   = max(1.0*xrmin,min(1.0*xrmax,xr_spl)) ! 1.0 to convert to double and keep nvfortran happy
          Dvr_spl  = (xr_spl/pirhow)**(1./3.)

          if (l_sb_classic) then
            ! limiting procedure (as per S&B)
            N_r0     = rhof(k)*n_hr(k)/Dvr_spl ! rhof(k)*n_hr(k)/Dvr_spl(k)
            N_r0     = max(N_0min,min(N_0max,N_r0))

            lbdr_spl = (pirhow*N_r0/(rhof(k)*qr_spl(k)))**0.25
            lbdr_spl = max(lbdr_min, min(lbdr_max,lbdr_spl))

            ! calculation of velocities
            wfall = max(0.,((rho0s/rhof(k))**0.5)*(a_tvsbc  &
                      -b_tvsbc*(1.+c_tvsbc/lbdr_spl)**(-4.0))) ! k=1
            sed_qr(k) = wfall*qr_spl(k)*rhof(k)

            wfall = max(0.,((rho0s/rhof(k))**0.5)*(a_tvsbc  &
                      -b_tvsbc*(1.+c_tvsbc/lbdr_spl)**(-1.0))) ! k=0
            sed_Nr(k) = wfall*Nr_spl(k)*rhof(k)
          else  ! l_sb_classic
            if (l_lognormal) then
              ! correction for width of DSD
              ! BUG: Dvr_spl unset? reusing the one from l_sb_classic
              Dgr = (exp(4.5*(log(sig_gr))**2))**(-1./3.)*Dvr_spl

              sed_qr(k) = 1.*sed_flux3(Nr_spl(k),Dgr,log(sig_gr)**2,D_s,3)
              sed_Nr(k) = 1./pirhow*sed_flux3(Nr_spl(k),Dgr,log(sig_gr)**2,D_s,0)

              ! correction for the fact that pwcont .ne. qr_spl
              ! actually in this way for every grid box a fall velocity is determined
              pwcont = liq_cont3(Nr_spl(k),Dgr,log(sig_gr)**2,D_s,3)         ! note : kg m-3
              if (pwcont > eps1) then
                sed_qr(k) = (qr_spl(k)*rhof(k)/pwcont)*sed_qr(k)
              end if
            else ! l_lognormal
              !
              ! SB rain sedimentation
              !
              if (l_mur_cst) then
                mur_spl = mur_cst
              else
                ! SS08
                ! NOTE: Dvr_spl is unset
                ! mur_spl = 10. * (1+tanh(1200.*(Dvr_spl(k)-0.0014)))

                ! G09b
                mur_spl = min(30.,- 1. + 0.008/ (qr_spl(k)*rhof(k))**0.6)
              endif ! l_mur_cst

              lbdr_spl  = ((mur_spl+3.)*(mur_spl+2.)*(mur_spl+1.))**(1./3.)/Dvr_spl ! BUG: Dvr_spl is unset like above
              wfall = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr_spl)**(-1.*(mur_spl+4.))))
              sed_qr(k) = wfall * qr_spl(k) * rhof(k)

              wfall = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr_spl)**(-1.*(mur_spl+1.))))
              sed_Nr(k) = wfall * Nr_spl(k) * rhof(k)
            endif ! l_lognormal
          endif ! l_sb_classic
        else ! l_sb
          !
          ! KK00 rain sedimentation
          !
          xr_spl = rhof(k)*qr_spl(k)/Nr_spl(k)

          ! to ensure xr is within borders
          xr_spl = min(xr_spl,1.0*xrmaxkk) ! 1.0 to convert to double and keep nvfortran happy

          Dvr_spl = (xr_spl/pirhow)**(1./3.)
          sed_qr(k) = max(0., 0.006*1.0E6*Dvr_spl - 0.2) * qr_spl(k)*rhof(k)
          sed_Nr(k) = max(0.,0.0035*1.0E6*Dvr_spl - 0.1) * Nr_spl(k)
        end if ! l_sb

        ! Adjust boundaries for the second k loop
        k_low2 = min(k_low2, k)
        k_high2 = max(k_high2, k)
      endif ! qr_spl
    enddo ! first k loop

    if (k_high2 == 0 .and. jn == 1) return ! no sedimentation and no updates

    ! As the rain falls down, we need to adjust the lower boundary
    k_low2 = max(1,k_low2 - 1)

    do k = k_low2,min(kmax,k_high2) ! second k loop
      qr_spl(k) = max(0.0, qr_spl(k) + (sed_qr(k+1) - sed_qr(k))*dt_spl/(dzf(k)*rhof(k)))
      Nr_spl(k) = max(0.0, Nr_spl(k) + (sed_Nr(k+1) - sed_Nr(k))*dt_spl/(dzf(k)*rhof(k)))
    enddo  ! second k loop

    ! BUG: check this part properly later
    if (jn == 1) then
      precep_hr = sed_qr(1)/rhof(1)          ! kg kg-1 m s-1
    endif

    ! Adjust boundaries for the first k loop
    k_low1 = k_low2
    k_high1 = k_high2

  enddo ! time splitting loop

  ! updates
  do k=k_low2,k_high2
    n_hrp(k) = n_hrp(k) + (Nr_spl(k) - n_hr(k))/delt
    q_hrp(k) = q_hrp(k) + (qr_spl(k) - q_hr(k))/delt

    if (l_tendencies) then
      tend(idn_hr_se,k) = (Nr_spl(k) - n_hr(k))/delt
      tend(idq_hr_se,k) = (qr_spl(k) - q_hr(k))/delt
    endif
  enddo
end subroutine sedim_rain3


! sedimentation of snow
! ---------------------
subroutine sedim_snow3(q_hs, n_hs, q_hsp, n_hsp, precep_hs, tend)
  use modglobal, only : k1,kmax,dzf
  use modfields, only : rhof
  implicit none
  real, intent(in)    :: q_hs(k1), n_hs(k1)
  real, intent(inout) :: q_hsp(k1), n_hsp(k1)
  real, intent(out)   :: precep_hs
  real, intent(out)   :: tend(ntends, k1)

  integer :: k,jn,n_spl
  integer :: k_low1,k_high1  & ! Boundaries for the first k loop
            ,k_low2,k_high2    ! Boundaries for the second k loop

  real  :: qip_spl(k1), nip_spl(k1)
  real  :: sed_qip(k1), sed_nip(k1)

  real :: wfall, xip_spl, wvar
  real :: dt_spl

  qip_spl = q_hs
  nip_spl = n_hs

  n_spl = ceiling(split_factor*wfallmax_hs*delt/(minval(dzf)))
  dt_spl = delt/real(n_spl)

  ! First time for the first k loop, go over the full column
  k_low1 = 1
  k_high1 = k1

  ! The boundaries for the second k loop are extended in the first loop
  k_low2 = k1
  k_high2 = 0

  do jn = 1 , n_spl ! time splitting loop
    sed_qip = 0.
    sed_nip = 0.

    do k=k_low1,k_high1
      ! terminal fall velocity
      if ((qip_spl(k) > qsnowmin).and.(nip_spl(k) > 0.0)) then
        xip_spl = qip_spl(k)/nip_spl(k)
        xip_spl = min(max(xip_spl,x_hs_bmin),x_hs_bmax) ! to ensure xr is within borders

        wfall = max(0.0, c_v_s1 * xip_spl**be_hs)
        sed_qip(k) = wfall*qip_spl(k)*rhof(k)

        sed_nip(k) = wfall*nip_spl(k)*rhof(k)
        wfall = max(0.0, c_v_s0 * xip_spl**be_hs)

        ! Adjust boundaries for the second k loop
        k_low2 = min(k_low2, k)
        k_high2 = max(k_high2, k)
      endif ! qs_spl
    enddo ! first k loop

    if (k_high2 == 0 .and. jn == 1) return ! no sedimentation and no updates

    ! As the snow falls down, we need to adjust the lower boundary
    k_low2 = max(1,k_low2 - 1)

    do k = k_low2,min(kmax,k_high2) ! second k loop
      qip_spl(k) = max(0.0, qip_spl(k) + (sed_qip(k+1) - sed_qip(k))*dt_spl/(dzf(k)*rhof(k)))
      nip_spl(k) = max(0.0, nip_spl(k) + (sed_nip(k+1) - sed_nip(k))*dt_spl/(dzf(k)*rhof(k)))
    enddo  ! second k loop

    ! BUG: check this part properly later
    if (jn == 1) then
      precep_hs = sed_qip(1)/rhof(1)        ! kg kg-1 m s-1
    endif

    ! Adjust boundaries for the first k loop
    k_low1 = k_low2
    k_high1 = k_high2

  enddo ! time splitting loop

  do k=k_low2,k_high2
    ! updates
    n_hsp(k) = n_hsp(k) + (nip_spl(k) - n_hs(k))/delt
    q_hsp(k) = q_hsp(k) + (qip_spl(k) - q_hs(k))/delt

    if (l_tendencies) then
      tend(idn_hs_se,k) = (nip_spl(k) - n_hs(k))/delt
      tend(idq_hs_se,k) = (qip_spl(k) - q_hs(k))/delt
    endif
  enddo
end subroutine sedim_snow3


! sedimentation of graupel
! ------------------------
subroutine sedim_graupel3(q_hg, n_hg, q_hgp, n_hgp, precep_hg, tend)
  use modglobal, only : k1,kmax,dzf
  use modfields, only : rhof
  implicit none
  real, intent(in)    :: q_hg(k1), n_hg(k1)
  real, intent(inout) :: q_hgp(k1), n_hgp(k1)
  real, intent(out)   :: precep_hg
  real, intent(out)   :: tend(ntends, k1)

  integer :: k,jn,n_spl
  integer :: k_low1,k_high1  & ! Boundaries for the first k loop
            ,k_low2,k_high2    ! Boundaries for the second k loop

  real  :: qip_spl(k1), nip_spl(k1)
  real  :: sed_qip(k1), sed_nip(k1)

  real :: wvar,xip_spl
  real :: dt_spl,wfall

  qip_spl = q_hg
  nip_spl = n_hg

  n_spl = ceiling(split_factor*wfallmax_hg*delt/(minval(dzf)))
  dt_spl = delt/real(n_spl)

  ! First time for the first k loop, go over the full column
  k_low1 = 1
  k_high1 = k1

  ! The boundaries for the second k loop are extended in the first loop
  k_low2 = k1
  k_high2 = 0

  do jn = 1 , n_spl ! time splitting loop
    sed_qip = 0.
    sed_nip = 0.

    do k=k_low1,k_high1
      if ((qip_spl(k) > qgrmin).and.(nip_spl(k) > 0.0)) then
        xip_spl = qip_spl(k)/nip_spl(k)
        xip_spl = min(max(xip_spl,x_hg_bmin),x_hg_bmax) ! to ensure xr is within borders

        wfall = max(0.0,c_v_g1 * xip_spl**be_hg)
        sed_qip(k) = wfall*qip_spl(k)*rhof(k)

        wfall = max(0.0,c_v_g0 * xip_spl**be_hg)
        sed_nip(k) = wfall*nip_spl(k)*rhof(k)

        ! Adjust boundaries for the second k loop
        k_low2 = min(k_low2, k)
        k_high2 = max(k_high2, k)
      endif ! qip_spl
    enddo ! first k loop

    if (k_high2 == 0 .and. jn == 1) return ! no sedimentation and no updates

    ! As the graupel falls down, we need to adjust the lower boundary
    k_low2 = max(1,k_low2 - 1)

    do k = k_low2,min(kmax,k_high2)
      qip_spl(k) = max(0.0, qip_spl(k) + (sed_qip(k+1) - sed_qip(k))*dt_spl/(dzf(k)*rhof(k)))
      nip_spl(k) = max(0.0, nip_spl(k) + (sed_nip(k+1) - sed_nip(k))*dt_spl/(dzf(k)*rhof(k)))
    enddo  ! second k loop

    ! BUG: check this part properly later
    if (jn == 1) then
      precep_hg = sed_qip(1)/rhof(1)          ! kg kg-1 m s-1
    endif

    ! Adjust boundaries for the first k loop
    k_low1 = k_low2
    k_high1 = k_high2

  enddo ! time splitting loop

  do k=k_low2,k_high2
    ! updates
    n_hgp(k) = n_hgp(k) + (nip_spl(k) - n_hg(k))/delt
    q_hgp(k) = q_hgp(k) + (qip_spl(k) - q_hg(k))/delt

    if (l_tendencies) then
      tend(idn_hg_se,k) = (nip_spl(k) - n_hg(k))/delt
      tend(idq_hg_se,k) = (qip_spl(k) - q_hg(k))/delt
    endif
  enddo
end subroutine sedim_graupel3


! sedimentation of cloud ice
! --------------------------
subroutine sedim_ice3(q_ci, n_ci, q_cip, n_cip, precep_ci, tend)
  use modglobal, only : k1,kmax,dzf
  use modfields, only : rhof
  implicit none
  real, intent(in)    :: q_ci(k1), n_ci(k1)
  real, intent(inout) :: q_cip(k1), n_cip(k1)
  real, intent(out)   :: precep_ci
  real, intent(out)   :: tend(ntends, k1)

  integer :: k,jn,n_spl
  integer :: k_low1,k_high1  & ! Boundaries for the first k loop
            ,k_low2,k_high2    ! Boundaries for the second k loop

  real :: qip_spl(k1), nip_spl(k1)
  real :: sed_qip(k1), sed_nip(k1)

  real :: dt_spl, xip_spl, wvar, wfall

  n_spl = ceiling(split_factor*wfallmax_ci*delt/(minval(dzf)))
  dt_spl = delt/real(n_spl)

  qip_spl = q_ci
  nip_spl = n_ci

  ! First time for the first k loop, go over the full column
  k_low1 = 1
  k_high1 = k1

  ! The boundaries for the second k loop are extended in the first loop
  k_low2 = k1
  k_high2 = 0

  do jn = 1 , n_spl ! time splitting loop
    sed_qip = 0.
    sed_nip = 0.

    do k=k_low1,k_high1
      if ((qip_spl(k) > qicemin).and.(nip_spl(k) > 0.0)) then
        xip_spl = qip_spl(k)/nip_spl(k)
        xip_spl = min(max(xip_spl,x_ci_bmin),x_ci_bmax) ! to ensure xr is within borders

        ! terminal fall velocity
        wfall = max(0.0,c_v_s1 * xip_spl**be_ci)
        sed_qip(k) = wfall*qip_spl(k)*rhof(k)

        wfall = max(0.0,c_v_s0 * xip_spl**be_ci)
        sed_nip(k) = wfall*nip_spl(k)*rhof(k)

        ! Adjust boundaries for the second k loop
        k_low2 = min(k_low2, k)
        k_high2 = max(k_high2, k)
      endif ! qip_spl
    enddo ! first k loop

    if (k_high2 == 0 .and. jn == 1) return ! no sedimentation and no updates

    ! As the graupel falls down, we need to adjust the lower boundary
    k_low2 = max(1,k_low2 - 1)

    ! segmentation over levels
    do k = k_low2,min(kmax,k_high2)
      qip_spl(k) = max(0.0, qip_spl(k) + (sed_qip(k+1) - sed_qip(k))*dt_spl/(dzf(k)*rhof(k)))
      nip_spl(k) = max(0.0, nip_spl(k) + (sed_nip(k+1) - sed_nip(k))*dt_spl/(dzf(k)*rhof(k)))
    enddo  ! second k loop

    ! BUG: check this part properly later
    if (jn == 1) then
      precep_ci = sed_qip(1)/rhof(1)         ! kg kg-1 m s-1
    endif

    ! Adjust boundaries for the first k loop
    k_low1 = k_low2
    k_high1 = k_high2

  enddo ! time splitting loop

  do k=k_low2,k_high2
    ! updates
    n_cip(k) = n_cip(k) + (nip_spl(k) - n_ci(k))/delt
    q_cip(k) = q_cip(k) + (qip_spl(k) - q_ci(k))/delt

    if (l_tendencies) then
      tend(idn_ci_se,k) = (nip_spl(k) - n_ci(k))/delt
      tend(idq_ci_se,k) = (qip_spl(k) - q_ci(k))/delt
    endif
  enddo
end subroutine sedim_ice3


!*********************************************************************
! sedimentation of cloud water
!*********************************************************************
subroutine sedim_cl3(q_cl, n_cl, q_clp, n_clp, n_ccp, qtpmcr, thlpmcr, tend)
  use modglobal, only : k1,kmax,dzf,rlv,cp
  use modfields, only : rhof, exnf
  implicit none
  real, intent(in)    :: q_cl(k1), n_cl(k1)
  real, intent(inout) :: q_clp(k1), n_clp(k1), n_ccp(k1)
  real, intent(inout) :: qtpmcr(k1), thlpmcr(k1)
  real, intent(out)   :: tend(ntends, k1)

  integer :: k, jn, n_spl
  integer :: k_low1,k_high1  & ! Boundaries for the first k loop
            ,k_low2,k_high2    ! Boundaries for the second k loop

  real :: qip_spl(k1), nip_spl(k1)
  real :: sed_qip(k1), sed_nip(k1)

  real :: dt_spl, xip_spl, wvar, wfall

  qip_spl = q_cl
  nip_spl = n_cl

  n_spl = ceiling(split_factor*wfallmax_cl*delt/(minval(dzf)))
  dt_spl = delt/real(n_spl)

  ! First time for the first k loop, go over the full column
  k_low1 = 1
  k_high1 = k1

  ! The boundaries for the second k loop are extended in the first loop
  k_low2 = k1
  k_high2 = 0

  do jn = 1 , n_spl ! time splitting loop

    sed_qip = 0.
    sed_nip = 0.

    do k=k_low1,k_high1
      if ((qip_spl(k) > qcliqmin).and.(nip_spl(k) > 0.0)) then
        xip_spl = qip_spl(k)/nip_spl(k)
        xip_spl = min(max(xip_spl,x_cl_bmin),x_cl_bmax) ! to ensure xr is within borders

        ! terminal fall velocity
        wfall = max(0.0,c_v_c1 * xip_spl**be_cl)
        sed_qip(k)   = wfall*qip_spl(k)*rhof(k)

        wfall = max(0.0,c_v_c0 * xip_spl**be_cl)
        sed_nip(k) = wfall*nip_spl(k)*rhof(k)

        ! Adjust boundaries for the second k loop
        k_low2 = min(k_low2, k)
        k_high2 = max(k_high2, k)
      endif ! qip_spl
    enddo ! first k loop

    if (k_high2 == 0 .and. jn == 1) return ! no sedimentation and no updates

    ! As the graupel falls down, we need to adjust the lower boundary
    k_low2 = max(1,k_low2 - 1)

    do k = k_low2,min(kmax,k_high2)
      qip_spl(k) = max(0.0, qip_spl(k) + (sed_qip(k+1) - sed_qip(k))*dt_spl/(dzf(k)*rhof(k)))
      nip_spl(k) = max(0.0, nip_spl(k) + (sed_nip(k+1) - sed_nip(k))*dt_spl/(dzf(k)*rhof(k)))
    enddo  ! second k loop

    ! Adjust boundaries for the first k loop
    k_low1 = k_low2
    k_high1 = k_high2

  enddo ! time splitting loop

  do k=k_low2,k_high2
    ! updates
    n_clp(k) = n_clp(k) + (nip_spl(k) - n_cl(k))/delt
    q_clp(k) = q_clp(k) + (qip_spl(k) - q_cl(k))/delt

    ! also qtpmcr and thlpmcr change
    qtpmcr(k)  = qtpmcr(k) + (qip_spl(k) - q_cl(k))/delt
    thlpmcr(k) = thlpmcr(k) - (rlv/(cp*exnf(k)))*(qip_spl(k) - q_cl(k))/delt

    ! NOTE: moved here from recover_cc point process
    ! recovery of ccn
    n_ccp(k) = n_ccp(k) + (nip_spl(k) - n_cl(k))/delt

    if (l_tendencies) then
      tend(idn_cl_se,k) = (nip_spl(k) - n_cl(k))/delt
      tend(idq_cl_se,k) = (qip_spl(k) - q_cl(k))/delt
    endif
  enddo
end subroutine sedim_cl3


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
! ---------------------------------------------------------------------
real function sed_flux3(Nin,Din,sig2,Ddiv,nnn)
  use modglobal, only : pi,rhow
  implicit none

  real, intent(in)    :: Nin, Din, sig2
  real(field_r), intent(in)    :: Ddiv
  integer, intent(in) :: nnn

  ! para. def. lognormal DSD (sig2 = ln^2 sigma_g), D sep. droplets from drops
  !,power of of D in integral
  real, parameter ::   C = rhow*pi/6.     &
                      ,D_intmin = 1e-6    &
                      ,D_intmax = 4.3e-3

  real ::  alfa         & ! constant in fall velocity relation
          ,beta         & ! power in fall vel. rel.
          ,D_min        & ! min integration limit
          ,D_max        & ! max integration limit
          ,flux           ![kg m^-2 s^-1]

  if (Din < Ddiv) then
    alfa = 3.e5*100  ![1/ms]
    beta = 2
    D_min = D_intmin
    D_max = Ddiv
    flux = C*Nin*alfa*erfint3(beta,Din,D_min,D_max,sig2,nnn)
  else
    ! fall speed ~ D^2
    alfa = 3.e5*100 ![1/m 1/s]
    beta = 2
    D_min = Ddiv
    D_max = 133e-6
    flux = flux + C*Nin*alfa*erfint3(beta,Din,D_min,D_max,sig2,nnn)

    ! fall speed ~ D
    alfa = 4e3     ![1/s]
    beta = 1
    D_min = 133e-6
    D_max = 1.25e-3
    flux = flux + C*Nin*alfa*erfint3(beta,Din,D_min,D_max,sig2,nnn)

    ! fall speed ~ sqrt(D)
    alfa = 1.4e3 *0.1  ![m^.5 1/s]
    beta = .5
    D_min = 1.25e-3
    D_max = D_intmax
    flux = flux + C*Nin*alfa*erfint3(beta,Din,D_min,D_max,sig2,nnn)
  end if
  sed_flux3 = flux
end function sed_flux3


! Function to calculate numerically the analytical solution of the
! liq. water content between Dmin and Dmax based on
! Feingold et al 1986 eq 17 -20.
!
! M.C. van Zanten    September 2005
! -----------------------------------------------------------------
real function liq_cont3(Nin,Din,sig2,Ddiv,nnn)
  use modglobal, only : pi,rhow
  implicit none

  real, intent(in)    :: Nin, Din, sig2
  real(field_r), intent(in) :: Ddiv
  integer, intent(in) :: nnn

  ! para. def. lognormal DSD (sig2 = ln^2 sigma_g), D sep. droplets from drops
  ! ,power of of D in integral
  real, parameter :: beta = 0           &
                    ,C = pi/6.*rhow     &
                    ,D_intmin = 80e-6   &   ! value of start of rain D
                    ,D_intmax = 3e-3        !4.3e-3    ! value is now max value for sqrt fall speed rel.

  real ::  D_min        & ! min integration limit
          ,D_max          ! max integration limit

  if (Din < Ddiv) then
    D_min = D_intmin
    D_max = Ddiv
  else
    D_min = Ddiv
    D_max = D_intmax
  end if

  liq_cont3 = C*Nin*erfint3(beta,Din,D_min,D_max,sig2,nnn)
end function liq_cont3


! Function to calculate erf(x) approximated by a polynomial as
! specified in 7.1.27 in Abramowitz and Stegun
! NB phi(x) = 0.5(erf(0.707107*x)+1) but 1 disappears by substraction
! -------------------------------------------------------------------
real function erfint3(beta, D, D_min, D_max, sig2,nnn )
  implicit none
  real, intent(in)    :: beta, D, D_min, D_max, sig2
  integer, intent(in) :: nnn

  real, parameter :: eps = 1e-10      &
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

  erfint3 = D**nn*exp(0.5*nn**2*sig2)*0.5*(erfymax-erfymin)
  if (erfint3 < 0.) erfint3 = 0.
end function erfint3

end module modbulkmicro3_column
