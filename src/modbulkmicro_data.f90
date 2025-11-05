!> Variables for SB/KK bulk microphysics parameterizations.
module modbulkmicro_data

  use modprecision, only: field_r

  implicit none

  public

  ! Constants
  real(field_r), parameter :: &
    c_st = 1.19E8 ! Stokes fall velocity [m^-1 s^-1]

  ! Threshold values
  real(field_r), parameter :: &
    qcmin = 1.0E-7,           & !< Cloud specific mixing ratio treshold for calculations.
    qrmin = 1.0E-13,          & !< Rain specific mixing ratio treshold for calculations.
    epscloud = 0.01E-3,       & !< Threshold for statistics.
    epsprec = 3.65E-5,        & !< Threshold for statistics.
    epsqr = 1.0E-8              !< Threshold for statistics.

  ! User options
  logical ::              &
    l_sb = .true.,        & !< Switch between SB or KK scheme.
    l_sedc = .true.,      & !< Switch for cloud droplet sedimentation.
    l_rain = .true.,      & !< Switch for rain formation.
    l_mur_cst = .false.,  & !< Use constant value for mu in DSD.
    l_lognormal = .false.   !< Use lognormal distribution for rain terminal velocities.

  !$acc declare create(l_mur_cst)

  real(field_r) :: &
    mur_cst = 5,   & !< Mu value if l_mur_cst = .true.
    sig_gr = 1.5     !< GSD of rain droplet DSD.

  !$acc declare create(mur_cst)

  integer :: &
    qcbase,  & !< Lowest model layer with cloud.
    qcroof,  & !< Highest model layer with cloud.
    qrbase,  & !< Lowest model layer with rain.
    qrroof     !< Highest model layer with rain.

  ! Arrays
  real(field_r), pointer :: &
    Nr(:,:,:),              & !< Rain droplet number concentration.
    qr(:,:,:),              & !< Rain specific mixing ratio.
    Nrp(:,:,:),             & !< Tendency of rain droplet number concentration.
    qrp(:,:,:)                !< Tendency of rain specific mixing ratio.

  ! Gamma function lookup tables
  real(field_r) :: &
    mygamma21(-100:4000), &
    mygamma251(-100:4000)

  !$acc declare create(mygamma21, mygamma251)

end module modbulkmicro_data
