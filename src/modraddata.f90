!> \file modraddata.f90
!! Variable definitions and auxilary routines for radiation

!>
!! Variable definitions and auxilary routines for radiation.
!>
!! This routine should have no dependency on any other routine, save perhaps modglobal or modfields.
!!  \author Thijs Heus, MPI-M
!!  \todo Documentation
!!  \par Revision list
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



module modraddata

! implicit none
  use modglobal, only : longint,kind_rb,SHR_KIND_IN,SHR_KIND_R4,kind_im
  use modprecision, only : field_r
SAVE

  integer, parameter :: irad_none  = 0   !< 0=no radiation
  integer, parameter :: irad_full  = 1   !< 1=full radiation
  integer, parameter :: irad_par   = 2   !< 2=parameterized radiation
  integer, parameter :: irad_lsm   = 3   !< 3=simple surface radiation for land surface model
  integer, parameter :: irad_rrtmg = 4   !< 4=radiation using the rapid radiative transfer model
  integer, parameter :: irad_user  = 10  !< 10=user specified radiation
  integer, parameter :: irad_rte_rrtmgp = 5 !< 5=radiation using the rapid radiative transfer model parallel

  logical :: rad_ls      = .true.        !< prescribed radiative forcing
  logical :: rad_longw   = .true.        !< parameterized longwave radiative forcing
  logical :: rad_shortw  = .true.        !< parameterized shortwave radiative forcing
  logical :: rad_smoke   = .false.       !< longwave divergence for smoke cloud
  logical :: useMcICA    = .true.        !< Use the Monte Carlo Independent Column Approach

  logical :: lcloudshading = .false.     !< Let clouds shade the surface for rad_lsm

  real              :: timerad = 0       !<  timescale of the radiation scheme
  integer(kind=longint)  :: itimerad = 0 !<  timescale of the radiation scheme
  integer (kind=longint) :: tnext   = 0  !<  time of the first upcoming call of the radiation scheme
  real :: rka        = 130.              !< extinction coefficient in radpar scheme
  real :: dlwtop     = 74.               !< longwave radiative flux divergence at top of domain
  real :: dlwbot     = 0.                !< longwave radiative flux divergence near the surface
  real :: sw0        = 1368.22           !< Solar constant (in W/m2). SWD at TOA = sw0*cos(mu)
                                         !< NOTE: when using delta-Eddington (iradiation=2) this represents the downwelling solar
                                         !        radiation at the top of the domain/cloud

  real :: gc         = 0.85              !< asymmetry factor of droplet scattering angle distribution
  real :: SSA        = 0.999             !< typical single scattering albedo for clouds
  integer :: iDE     = 1                 !< scalar field to be used as extinction
  logical :: laero   = .false.           !< .true. for aeosols .false. for clouds

  real :: reff       = 1.e-5             !< cloud droplet effective radius (m)
  integer :: isvsmoke = 1                !< number of passive scalar to be used for optical depth calculation
  integer :: iradiation = irad_none      !< Selection parameter for type of radiation scheme
  integer :: irad    = -1                !< Deprecated selection parameter for the type of radiation scheme
  logical :: lCnstZenith = .false.       !< Switch to disable the diurnal cycle and use diurnally averaged SW radiation (e.g. CGILS)
  logical :: lCnstAlbedo = .true.        !< Switch to disable the surface albedo parameterization in RRTMG
  real :: cnstZenith=0.                  !< constant zenith angle, only used when lCnstZenith=.true. (degrees!)

  ! Options in NAMRADIATION that apply to the rrtmg script
  integer(kind=kind_im) :: iaer = 0      ! Aerosol option (SW flag); 0: No aerosol; 6: ECMWF method; 10: Input aerosol optical properties
  integer(kind=kind_im) :: ioverlap = 2  ! Cloud overlap method; 0: Clear only; 1: Random; 2: Maximum/random; 3: Maximum
  integer(kind=kind_im) :: inflglw = 2   ! 0:inp. cld fr and opt. depth; 1:cf and LWP are input; 2:also ice fraction inp.
  integer(kind=kind_im) :: iceflglw = 3  ! 0,1,2,3: ice influence calculations
  integer(kind=kind_im) :: liqflglw = 1  ! 0:optical depths computed; 1:drop eff. rad. is input, opt. depth computed
  integer(kind=kind_im) :: inflgsw = 2   ! 0:inp. cld fr and opt. depth; 1:cf and LWP are input; 2:also ice fraction inp.
  integer(kind=kind_im) :: iceflgsw = 3  ! 0,1,2,3: ice influence calculations
  integer(kind=kind_im) :: liqflgsw = 1  ! 0:optical depths computed; 1:drop eff. rad. is input, opt. depth computed
  logical :: ocean  = .false.            ! if true, run is over ocean.
  logical :: usero3 = .false.            ! if true, the o3 profile is taken from backrad.inp, otherwise from stnd prof RRTMG
  real    :: co2_fraction = -1.          ! If given in namoptions, the CO2 volume fraction is set to this value for RRTMG (CGILS)
  real    :: ch4_fraction = -1.          ! If given in namoptions, the CH4 volume fraction is set to this value for RRTMG (RCEMIP)
  real    :: n2o_fraction = -1.          ! If given in namoptions, the N2O volume fraction is set to this value for RRTMG (RCEMIP)
  logical :: doperpetual = .false.       ! if true, no diurnal cycle is used, but rather a diurnally averaged forcing
  logical :: doseasons = .true.          ! if false, the same day will be repeated, otherwise, next day is taken
  integer(SHR_KIND_IN) :: iyear = 1992   ! The year of the simulation

  ! Options in NAMRTERRTMGP that apply to the RTE-RRTMGP library
  logical :: doclearsky = .false.
  logical :: usepade = .false.
  integer :: nbatch = 0

  ! Logicals and variables that are used in the modradrrtmg module
  logical :: isInitializedRrtmg = .false.           ! used as a initialization check
  logical :: isReadSounding = .false.               ! used as a check for reading the sounding file
  logical :: isAllocated_RadInputsOutputs = .false. ! used as a check for allocating radiation variables
  logical :: isAllocated_TraceGases = .false.       ! check for allocating tracegas variables
  logical :: isReadTraceProfiles =.false.           ! check for reading trace gas profiles for netcdf file

  real(kind=kind_rb),allocatable,dimension(:,:) :: tabs_slice,     &    ! Absolute temperature (2D slice)
                                                   qv_slice,       &    ! Water vapour content (2D slice)
                                                   qcl_slice,      &    ! Liquid water content (2D slice)
                                                   qci_slice,      &    ! Ice content          (2D slice)
                                                   o3_slice,       &    ! Ozon content         (2D slice)
                                                   rho_slice,      &    ! Density              (2D slice)
                                                   lwHR_slice,     &    ! Heating rate due to longwave rad          (2D slice)
                                                   lwHRCS_slice,   &    ! Heating rate due to longwave rad,clear sky value          (2D slice)
                                                   swDownDif_slice,&    ! Downwelling shortwave diffuse rad         (2D slice)
                                                   swHR_slice,     &    ! Heating rate due to shortwave rad         (2D slice)
                                                   swHRCS_slice         ! Heating rate due to shortwave rad,clear sky value         (2D slice)
  real(kind=kind_rb),allocatable,target,dimension(:,:) :: lwUp_slice,     & ! Upwelling longwave rad                    (2D slice)
                                                          lwDown_slice,   & ! Downwelling longwave rad                  (2D slice)
                                                          swUp_slice,     & ! Upwelling shortwave rad                   (2D slice)
                                                          swDown_slice,   & ! Downwelling shortwave rad                 (2D slice)
                                                          swDownDir_slice,& ! Downwelling shortwave direct rad          (2D slice)
                                                          lwUpCS_slice,   &  ! Upwelling longwave rad, clear sky value   (2D slice)
                                                          lwDownCS_slice, &  ! Downwelling longwave rad, clear sky value (2D slice)
                                                          swUpCS_slice,   &  ! Upwelling shortwave rad, clear sky value  (2D slice)
                                                          swDownCS_slice     ! Downwelling shortwave rad, clear sky value(2D slice)

  real(kind=kind_rb),allocatable,dimension(:) :: solarZenithAngleCos  ! The zenith angle of a slice
  real(kind=kind_rb),allocatable,dimension(:) :: asdir,asdif,aldir,aldif                         ! Albedos ...

  real(kind=kind_rb),allocatable,dimension(:,:) :: layerP,    &
                                                   layerT,    &
                                                   interfaceP,&
                                                   interfaceT,&
                                                   h2ovmr,    &
                                                   o3vmr,     &
                                                   co2vmr,    &
                                                   ch4vmr,    &
                                                   n2ovmr,    &
                                                   o2vmr,     &
                                                   cfc11vmr,  &
                                                   cfc12vmr,  &
                                                   cfc22vmr,  &
                                                   ccl4vmr,   &
                                                   tracevmr,  &
                                                   emis
  real(kind=kind_rb),allocatable,dimension(:,:) :: LWP_slice,IWP_slice ,cloudFrac,liquidRe,iceRe
  real(kind=kind_rb),allocatable,dimension(:,:,:) :: taucldlw, &
                                                     tauaerlw, &
                                                     taucldsw, &
                                                     ssacldsw, &
                                                     asmcldsw, &
                                                     fsfcldsw, &
                                                     tauaersw, &
                                                     ssaaersw, &
                                                     asmaersw, &
                                                     ecaersw
  real(SHR_KIND_R4) :: eccen,   &  ! Earth's eccentricity factor (unitless) (typically 0 to 0.1)
                       obliq,   &  ! Earth's obliquity angle (deg) (-90 to +90) (typically 22-26)
                       obliqr,  &  ! Earths obliquity in radians
                       lambm0,  &  ! Mean long of perihelion at the vernal equinox (radians)
                       mvelp,   &  ! Earth's moving vernal equinox at perhelion (deg)(0 to 360.0)
                       mvelpp,  &  ! moving vernal equinox longitude of perihelion plus pi (radians) 
                       delta,   &  ! Solar declination angle in rad
                       eccf        ! Earth-sun distance factor (ie. (1/r)**2)

  real,parameter :: mwdry = 28.966, &
                    mwh2o = 18.016, &
                    mwo3  = 47.998
  !real ,parameter :: scon = 765.84!471.5!1367                          ! Solar constant in W/m^2; code Thijs uses 1365
  real  :: scon
  real  :: mu0_cgils
  integer :: cgils_case_nr
  real, parameter :: tmelt = 273.16
  real,allocatable,dimension(:)   :: presf_input,     &   ! Full-level pressure (sounding patched to domain)
                                     presh_input          ! Halflevel  pressure (sounding patched to domain)
  real,allocatable,dimension(:)   :: tg_slice             ! Sea surface temperature of a 2D slice

  real(kind_rb),allocatable,dimension(:)   :: &
       o3, co2, ch4, n2o, o2, cfc11, cfc12, cfc22, ccl4   ! Profiles of trace gases
  integer :: npatch_start,npatch_end
  integer :: nzrad                                        ! Number of levels in the patched radiation profiles
  integer :: kradmax, krad1, krad2                        ! New variables (stephan), kradmax = nzrad-1, krad1=nzrad, krad2=nzrad+1 (like kmax, k1)

  ! background sounding
  integer, parameter :: nzsnd = 1000
  real,allocatable,dimension(:) :: psnd, & ! pressure sounding read in from SoundingFileName, mb (hPa)
                                   tsnd, & ! temperature sounding read in from SoundingFileName, K
                                   qsnd, & ! water vapor sounding read in from SoundingFileName, kg/kg
                                   o3snd   ! ozon sounding read in from SoundingFileName (if usero3=true)

  real mu                            !< cosine of the solar zenith angle

  real(field_r), allocatable :: thlprad(:,:,:)!<   the radiative tendencies
  real(field_r), allocatable :: swd(:,:,:)    !<   shortwave downward radiative flux
  real(field_r), allocatable :: swdir(:,:,:)  !<   Direct shortwave downward radiative flux
  real(field_r), allocatable :: swdif(:,:,:)  !<   Difuse shortwave downward radiative flux
  real(field_r), allocatable :: lwc(:,:,:)    !<   Liquid water content calculated in rrtmg
  real(field_r), allocatable :: swu(:,:,:)    !<   shortwave upward radiative flux
  real(field_r), allocatable :: lwd(:,:,:)    !<   longwave downward radiative flux
  real(field_r), allocatable :: lwu(:,:,:)    !<   longwave upward radiative flux
!
  real(field_r), allocatable :: swdca(:,:,:)  !<  clear air shortwave downward radiative flux
  real(field_r), allocatable :: swuca(:,:,:)  !<  clear air shortwave upward radiative flux
  real(field_r), allocatable :: lwdca(:,:,:)  !<  clear air longwave downward radiative flux
  real(field_r), allocatable :: lwuca(:,:,:)  !<  clear air longwave upward radiative flux

  real(field_r), allocatable :: SW_up_TOA(:,:), SW_dn_TOA(:,:), LW_up_TOA(:,:), LW_dn_TOA(:,:) !< Top of the atmosphere radiative fluxes
  real(field_r), allocatable :: SW_up_ca_TOA(:,:), SW_dn_ca_TOA(:,:), LW_up_ca_TOA(:,:), LW_dn_ca_TOA(:,:)


contains
!< Calculation of the cosine of the zenith angle
!< \param time UTC Time of the simulation
!< \param xday Day at the start of the simulation
!< \param xlat Latitude of the domain
!< \param xlon Longitude of the domain
  real function zenith(time, xday, xlat,xlon)
    use modglobal, only : pi
!     implicit none
    real, intent(in) :: time, xday, xlat, xlon
    real :: phi,el,obliq,xlam,declin,hora
    real :: day,daytime

    if (.not.lCnstZenith) then
      day    = xday + floor(time/86400.)
      daytime= mod(time,86400.)

      phi    = xlat * pi/180.
      el     = xlon * pi/180.
      obliq  = 23.45 * pi/180.
      xlam   = 4.88 + 0.0172 * day
      declin = asin(sin(obliq)*sin(xlam))
      hora   = el-pi + 2.*pi*(daytime/86400.)
      zenith = max(0.,sin(declin)*sin(phi)+cos(declin)*cos(phi)* &
                                                         cos(hora))
    else
      zenith = cos(cnstZenith*pi/180.)
    end if
  end function zenith

end module modraddata
