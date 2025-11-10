module modsimpleice_data

  use modprecision, only: field_r

  implicit none

  public

  ! Parameters
  real(field_r), parameter ::     &
    qrmin = 1.0E-13,              &
    qcmin = 1.0E-7,               &
    ! Mass-diameter parameters A and B, terminal velocity parameters C, and D
    ! GRABOWSKI
    aar = 5.2e2,                  &
    bbr = 3.,                     &
    ccr = 130.,                   &
    ddr = 0.5,                    &
    !ccr = 842. & coefficients in Khairoutdinov and Randall
    !ddr = 0.8 & coefficients in Khairoutdinov and Randall
    ! For snow
    ! GRABOWSKI
    aas = 2.5e-2,                 &
    bbs = 2.,                     &
    ccs = 4.,                     &
    dds = 0.25,                   &
    ! For graupel (if present, following Tomita 2008 for terminal velocities and using mass-diameter
    ! relationship as for rain, but with only 40% of density)
    aag = 2.e2,                   &
    bbg = 3.,                     &
    ccg = 82.5,                   &
    ddg = 0.25,                   &
    ceffrl=0.8,                   &
    ceffsl=0.06,                  & ! probably 0.8 is better, wsa exp 156
    ceffgl=0.06,                  & ! probably 0.8 is better
    ceffri=0.8,                   &
    ceffsi=0.06,                  &
    ceffgi=0.06,                  &
    ! Shape factors beta GRABOWSKI
    betar=2.,                     &
    betas=3.,                     &
    betag=2.,                     &
    ! N_0 in Marshall-Palmer Distribution following Grabowski
    n0rr=2.e7,                    &
    n0rs=2.e7,                    &
    n0rg=2.e7,                    &
    ! N_0 in Marshall-Palmer Distribution following Tomita
    ! n0rr=8.e6,                    &
    ! n0rs=4.e6,                    &
    ! n0rg=3.e6,                    &
    ! Gamma distribution parameters, calculated only once
    ! Parameters for Kessler/Lin type autoconversion
    timekessl=0.001,              &
    betakessi=0.001,              &
    qll0=0.001,                   &
    qli0=0.0001,                  &
    ! Diagnostic division between rain, snow and graupel
    tuprsg=268.,                  &
    tdnrsg=253.,                  &
    tupsg=283.,                   & ! Following Khairoutdinov and Randall
    tdnsg=223.                      ! Following Khairoutdinov and Randall

  ! User settings
  logical ::            &
    l_berry = .true.,   & !< Berry-Hsie autoconversion vs Kessler-Lin.
    l_graupel = .true., & !< Switch for graupel.
    l_warm = .false.,   & !< Run ice micro in warm mode, as a check.
    l_mp = .true.         !< Use Marshall-Palmer distribution for rain.

  real(field_r) ::  &
    evapfactor = 1, & !< Prefactor to reduce evaporation.
    courantp = 1      !< CFLmax-criterion for precipitation.

  ! Arrays
  real(field_r), allocatable :: &
    qr(:,:,:),                  & !< Total precipitation specific mixing ratio.
    qrp(:,:,:),                 & !< Tendency of precipitation specific mixing ratio.
    qr_spl(:,:,:),              &
    sed_qr(:,:,:),              &
    ilratio(:,:,:),             &
    rsgratio(:,:,:),            &
    sgratio(:,:,:),             &
    lambdar(:,:,:),             &
    lambdas(:,:,:),             &
    lambdag(:,:,:),             &
    ccrz(:),                    &
    ccsz(:),                    &
    ccgz(:),                    &
    ccrz2(:),                   &
    ccsz2(:),                   &
    ccgz2(:)                      

end module modsimpleice_data
