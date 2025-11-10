module modlsmdata
  implicit none

  public

  logical :: llsm            ! On/off switch LSM
  logical :: lfreedrainage   ! Free drainage bottom BC for soil moisture
  logical :: lags            ! Switch for A-Gs scheme

  ! Interpolation types soil from full to half level
  integer :: iinterp_t, iinterp_theta
  integer, parameter :: iinterp_amean = 1  ! val = 0.5*(v1+v2)
  integer, parameter :: iinterp_gmean = 2  ! val = sqrt(v1*v2)
  integer, parameter :: iinterp_hmean = 3  ! val = ((dz1+dz2)*v1*v2)/(dz1*v2+dz2*v1)
  integer, parameter :: iinterp_max   = 4  ! val = max(a,b)

  ! Soil grid
  integer :: kmax_soil = -1
  real :: z_size_soil = -1
  real, allocatable :: z_soil(:), zh_soil(:)
  real, allocatable :: dz_soil(:), dzh_soil(:)
  real, allocatable :: dzi_soil(:), dzhi_soil(:)

  ! Soil index in `van_genuchten_parameters.nc` lookup table:
  integer, allocatable :: soil_index(:,:,:)

  ! Source term in soil water budget
  real, allocatable :: phiw_source(:,:,:)

  ! Precipitation, interception et al.
  real, allocatable :: throughfall(:,:)
  real, allocatable :: interception(:,:)
  real, allocatable :: wl_max(:,:)

  ! Dependency factor canopy resistance on VPD (high veg only)
  real, allocatable :: gD(:,:)

  ! Reduction functions canopy resistance
  real, allocatable :: f1(:,:), f2b(:,:)

  ! Random
  real, allocatable :: du_tot(:,:), thv_1(:,:), land_frac(:,:), cveg(:,:)

  ! A-Gs
  real, allocatable :: an_co2(:,:), resp_co2(:,:)
  integer :: co2_index = -1

  ! Data structure for sub-grid tiles
  type T_lsm_tile
  ! Fixed LU properties
      ! Land use name
      character(len=64) :: luname
      character(len=3)  :: lushort
      ! Check if LU type is vegetation
      logical           :: lveg
      ! Check if LU type is water
      logical           :: laqu
      ! Static properties:
      real, allocatable :: z0m(:,:), z0h(:,:)
      ! Base tile fraction (i.e. without liquid water)
      real, allocatable :: base_frac(:,:)
      ! Conductivity skin layer:
      real, allocatable :: lambda_stable(:,:), lambda_unstable(:,:)
      ! Vegetation properties:
      real, allocatable :: lai(:,:)
      ! Surface properties:
      real, allocatable :: rs_min(:,:)
      ! Root fraction parameters
      real, allocatable :: a_r(:,:), b_r(:,:)

  ! Dynamic land surface properties
      ! Dynamic tile fraction:
      real, allocatable :: frac(:,:)
      ! Monin-obukhov / surface layer:
      real, allocatable :: obuk(:,:), ustar(:,:), ra(:,:)
      ! Surface fluxes:
      real, allocatable :: H(:,:), LE(:,:), G(:,:)
      real, allocatable :: wthl(:,:), wqt(:,:)
      ! Surface (potential) temperature and humidity:
      real, allocatable :: tskin(:,:), thlskin(:,:), qtskin(:,:)
      ! Buoyancy difference surface - atmosphere
      real, allocatable :: db(:,:)
      ! Vegetation properties:
      real, allocatable :: rs(:,:)
      real, allocatable :: f2(:,:), f3(:,:), gD(:,:)
      ! Root fraction in soil
      real, allocatable :: root_frac(:,:,:)
      ! Root fraction weighted mean soil water content
      real, allocatable :: phiw_mean(:,:)

  ! LU dependent deposition parameters
      ! In-canopy resistance parameters
!      real, allocatable :: R_inc_b(:,:), R_inc_h(:,:)
      ! SAI = SAI_a * LAI + SAI_b
      real, allocatable :: SAI_a(:,:), SAI_b(:,:)
      ! Minimum correction factor for stomatal resistance
!      real, allocatable :: fmin(:,:)
      ! Alpha value for light correction of stomatal resistance
!      real, allocatable :: alpha(:,:)
      ! Min, optimum and max emperatures for temperature correction of stomatal resistance
!      real, allocatable :: Tmin(:,:), Topt(:,:), Tmax(:,:)
      ! Maximum leaf stomatal conductance for ozone
!      real, allocatable :: gs_max(:,:)
      ! Minimum and maximum vapour pressure deficit parameters
!      real, allocatable :: vpd_min(:,:), vpd_max(:,:)
      ! Gamma parameter for calculating stomata compensation point
!      real, allocatable :: gamma_stom(:,:)
      ! Gamma correction factor for calculating soil compensation point
!      real, allocatable :: gamma_soil_c_fac(:,:)
      ! Gamma parameter for calculating soil compensation point
!      real, allocatable :: gamma_soil_default(:,:)

  end type T_lsm_tile

  !Tiles for all LU types
  integer :: ilu, nlu, ilu_ws
  type(T_lsm_tile), allocatable :: tile(:)

  ! Land-surface / van Genuchten parameters from NetCDF input table.
  real, allocatable :: &
      theta_res(:), theta_wp(:), theta_fc(:), theta_sat(:), &
      gamma_theta_sat(:), vg_a(:), vg_l(:), vg_n(:)
  ! Derived soil parameters
  real, allocatable :: &
      vg_m(:), &
      lambda_theta_min(:), lambda_theta_max(:), &
      gamma_theta_min(:), gamma_theta_max(:), &
      gamma_t_dry(:), rho_C(:)

end module modlsmdata