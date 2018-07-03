!> \file modmicrodata.f90
!!  Variables necessary for the microphysics

!>
!!  Variables necessary for the microphysics
!>
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

  module modmicrodata

  use modglobal, only : rhow,lacz_gamma
  implicit none
  save
  integer :: imicro = 0

  integer, parameter :: imicro_none    = 0
  integer, parameter :: imicro_drizzle = 1
  integer, parameter :: imicro_bulk    = 2
  integer, parameter :: imicro_bin     = 3
  integer, parameter :: imicro_sice    = 5
  integer, parameter :: imicro_sice2   = 6
  integer, parameter :: imicro_user    = 10
  logical :: l_sb        = .true. , &!< SB scheme (.true.) / KK00 scheme (.false.)   (in namelist NAMMICROPHYSICS)
             l_sedc      = .true. , & !<  cloud droplet sedimentation flag             (in namelist NAMMICROPHYSICS)
             l_rain      = .true. , & !<  rain formation / evolution flag              (in namelist NAMMICROPHYSICS)
             l_mur_cst   = .false. ! false = no constant value of mur (mur=f(Dv)) (in namelist NAMMICROPHYSICS)

  real    :: mur_cst     = 5        & !<  mur value if l_mur_cst=T                     (in namelist NAMMICROPHYSICS)
                 ,Nc_0 = 70e6       & !<  initial cloud droplet number (#/m3)
                 ,sig_g = 1.34      & !<  geom. std dev of cloud droplet DSD
                 ,sig_gr = 1.5        !<  geometric std dev of rain drop DSD

  logical :: l_lognormal = .false.    !<  log param of rain terminal velocities for rain sedim

  integer :: inr = 1, iqr=2

  real, parameter ::  D0_kk = 50e-6     & !<  diameter sep. cloud and prec. in KK00 scheme
                     ,qcmin = 1.0e-7     & !<  Cloud specific mixing ratio treshold for calculations
                     ,qrmin = 1.0e-13    & !<  Rain  specific mixing ratio treshold for calculations
!                     ,nuc = 0           & !< width parameter of cloud DSD
                     ,eps0 = 1e-20      & !< parameter used to avoid division by zero floating point exceptions
                     ,epscloud= 0.01e-3 &
                     ,epsprec = 3.65e-5 & !<  RICO threshold
                     ,epsqr = 1.0e-8     &
!  values picked by Verica Savic-Jovcic to optimize for Sc, note x and D have to be chosen consistently
!                ,xcmin = 4.2e-15       & !<  min mean mass of cw
!                ,xcmax = 6.5e-11       & !<  max mean mass of cw
!                ,xrmin = xcmax         & !<  min mean mass of pw
                     ,xrmax = 5.0e-6       & !<  max mean mass of pw
                     ,xrmaxkk = 5.2e-7     & !<  max mean mass of pw in KK00 scheme
!                ,Dvcmin = 2.0e-6       & !<  min mean diam. of cw
!                     ,Dvcmax = 49.8e-6      & !<  max mean diam. of cw
!                     ,Dvrmin = Dvcmax       & !<  min mean diam. of pw
!                     ,Dvrmax = 1000.0e-6    & !<  max mean diam. of pw
!  values given by SB2001
              ,xcmin = 4.2e-15     & !< \param xcmin  min mean mass of cw (D = 2.0e-6m)
              ,xcmax = 2.6e-10     & !<  max mean mass of cw (D = 80e-6m)
              ,xrmin = xcmax       & !<  min mean mass of pw
!               ,xrmax = 6.0e-07      & !<  max mean mass of pw
               ,Dvcmin = 2.0e-6     & !<  min mean diam. of cw
               ,Dvcmax = 79.2e-6    & !<  max mean diam. of cw
               ,Dvrmin = Dvcmax     & !<  min mean diam. of pw
               ,Dvrmax = 3000.0e-6  & !<  max mean diam. of pw
! NB in table1 in SB2006 komen weer andere getallen voor
! NB x_s is 'scheidingsdrop massa' en die mag dus best groter zijn dan bovengrens
! xcmax omdat die voor mean droplet mass staat!<  -> in gedachten houden
! NB de microphysica is heel gevoelig voor de waarde van de scheidingsdrop massa!<
        ,x_s = xcmax   & !<  drop mass sep. cloud and prec. part of DSD
        ,D_s = Dvcmax  & !<  diameter sep. cloud and prec. part of DSD
!          ,x_s = 2.6e-10  &
!          ,D_s = 79.2e-6  &
!          ,k_c = 9.44e9   & !<  Long Kernel coef. SB2001 [m^3 kg^-2 s^-1]
!          ,k_1 = 6.0e2    & !<  k_1 + k_2: coef. for phi function
!          ,k_2 = 0.68     & !<  in autoconversion rate SB2001
         ,k_c = 10.58e9 & !<  Long Kernel coef. SB2006 (k'cc)
!          ,k_c = 4.44e9   &  !<  Long Kernel coef. SB2006 (k'cc) for test
        ,k_1 = 4.0e2   & !<  k_1 + k_2: coef. for phi function
        ,k_2 = 0.70    & !<  in autoconversion rate SB2006
!
!          ,k_r = 5.78     & !<  Kernel coef. SB2001 [m^3 kg^-1 s^-1]
!          ,k_l = 5.e-4    & !<  coef for phi function in accr. rate
!          ,kappa_r = 0,   &
!          ,k_rr= k_r      &
         ,k_r = 5.25     & !<  Kernel SB2006
         ,k_l = 5.e-5    & !<  coef. for phi function in accr. rate
         ,kappa_r = 60.7 & !<  see eq. 11 SB2006
         ,k_rr = 7.12    & !<  idem dito

         ,Kt    = 2.5e-2  & !<  conductivity of heat [J/(sKm)]
         ,Dv    = 2.4e-5   & !<  diffusivity of water vapor [m2/s]
!  NB (see table 7.1 in Rogers: given Kt is for ~15 C while Dv is for > 30 C  2.4e-5
!  is value for ~ 15C How sensitive is G for this Aug 2006, ~5% -> Dv changed to 15 C value?
        ,c_St  = 1.19e8  & !<  Stokes fall vel. coef. [m^-1 s^-1]
!          ,pirhow = (pi*rhow)/6. & !< used in conversion of mass to diameter
                     ,pirhow = 3.14159*rhow/6.        &
         ,Rv = 461.5       & !<  specific gas constant for water vapor
         ,avf = 0.78       & !<  constants in vent. factor fv   (fv = 1. --> av=1,
         ,bvf = 0.308      & !<                                             bv=0 )
         ,nu_a = 1.41e-5   & !<  kin. viscosity of air [m2s-1]
         ,c_Nevap = 0.7    & !<  coeff for evap
         ,c_evapkk = 0.87  & !<  coeff for evap in KK00 scheme
         ,Sc_num = 0.71    &    !<  Schmidt number
         ,a_tvsb = 9.65    & !<  coeff in terminal velocity param
         ,b_tvsb = 9.8     & !<  coeff in terminal velocity param
         ,c_tvsb = 600.      !<  coeff in terminal velocity param

  real,allocatable, dimension(:,:,:) :: qc  & !<  cloud droplets specific mixing ratio [kg_w/kg_a]
                                       ,Nc  & !<  cloud droplets number conc.  [#/m^3]
                                       ,nuc & !<  width parameter of cloud DSD
                                       ,rhoz  !< slab averaged density in 3 dimensions

  real,allocatable, dimension(:,:,:) :: qr_spl, Nr_spl
                             !< prec. liq. water and conc. for sedim. time splitting
  real,allocatable, dimension(:,:,:) :: sedc,   & !<  sedimentation cloud droplets mix. ratio
                                        sed_qr, & !<  sedimentation rain drops mix. ratio
                                        sed_Nr    !<  sedimentation rain drop number conc.
  real ::  rho_c             &      !<  term to correct for density dep. of fall vel.
    ,k_au                     !<  coeff. for autoconversion rate
  real,allocatable, dimension(:,:,:) ::  &
    presz              &      !<  3D pressure
    ,Dvc               &      !<  cloud water mean diameter
    ,xc                &      !<  mean mass of cloud water droplets
    ,Dvr               &      !<  prec water mean diameter
    ,xr                &      !<  mean mass of prec. water drops
    ,mur               &      !<  mu parameter in rain gamma distribution
    ,lbdr              &      !<  slope parameter (lambda) in rain gamma distribution
    ,au                &      !<  autoconversion rate
    ,phi               &      !<  correction function (see SB2001)
    ,tau               &      !<  internal time scale
    ,ac                &      !<  accretion rate
    ,sc                &      !<  self collection rate
    ,br                &      !<  break-up rate
    ,evap              &      !<  mass tendency due to rain evap/cond
    ,Nevap             &      !<  concentration tendency due to rain evap/cond
    ,wfall_qr          &      !<  fall velocity for qr
    ,wfall_Nr                 !<  fall velocity for Nr

  real :: csed                      !<  parameter in cloud water grav. settling formula

  real, parameter ::  D_eq = 1.1E-3,  & !<  Parameters for break-up
            k_br = 1000       !<

   real,allocatable,dimension(:,:,:) :: Nr,Nrp,qltot,qr,qrp,thlpmcr,qtpmcr
   real,allocatable,dimension(:,:,:) :: precep

  real :: delt

  logical ,allocatable,dimension(:,:,:):: qcmask,qrmask

! Parameters for simple ice microphysics (Grabowski, JAS, 1998)
! With extension to graupel class if l_graupel=.true.
! Latent heats saved in modglobal.f90
! Tup and Tdn for distinction between cloud water and cloud ice  saved in modglobal.f90
! Drop concentration from initial droplet concentration input Nc_0 saved in modglobal.f90


  ! user settings
  logical :: l_berry = .true.   !  Berry-Hsie (Grabowski, 1998) autoconversion vs Kessler-Lin (Khairoutdinov and Randall, 2006)
  logical :: l_graupel = .true. !  Switch for graupel
  logical :: l_warm = .false.   !  Run ice micro in warm mode, as a check
  logical :: l_mp = .true.      !  Use marshall-palmer distribution for rain?
  real :: evapfactor = 1.0      !  Prefactor to reduce evaporation
  real :: courantp = 1.0        !  CFLmax-criterion for precipitation

  real, parameter :: &
     ! Mass-diameter parameters A and B, terminal velocity parameters C, and D
     ! GRABOWSKI
     aar=5.2e2 &
     ,bbr=3. &
     ,ccr=130. &
     ,ddr=0.5 &
!     ,ccr=842. & coefficients in Khairoutdinov and Randall
!     ,ddr=0.8 & coefficients in Khairoutdinov and Randall
     ! For snow
     ! GRABOWSKI
     ,aas=2.5e-2 &
     ,bbs=2. &
     ,ccs=4. &
     ,dds=0.25 &
     ! For graupel (if present, following Tomita 2008 for terminal velocities and using mass-diameter
     ! relationship as for rain, but with only 40% of density)
     ,aag=2.e2 &
     ,bbg=3.&
     ,ccg=82.5 &
     ,ddg=0.25 &
     ! Collection efficiency matrix, alpha factor of Grabowski has been absorbed here GRABOSKWI
     ,ceffrl=0.8 &
     ,ceffsl=0.06 & ! probably 0.8 is better, wsa exp 156
     ,ceffgl=0.06 & ! probably 0.8 is better
     ,ceffri=0.8 &
     ,ceffsi=0.06 &
     ,ceffgi=0.06 &
     ! Shape factors beta GRABOWSKI
     ,betar=2. &
     ,betas=3. &
     ,betag=2. &
     ! N_0 in Marshall-Palmer Distribution following Grabowski
     ,n0rr=2.e7 &
     ,n0rs=2.e7 &
     ,n0rg=2.e7 &
     ! N_0 in Marshall-Palmer Distribution following Tomita
     ! ,n0rr=8.e6 &
     ! ,n0rs=4.e6 &
     ! ,n0rg=3.e6 &
     ! Gamma distribution parameters, calculated only once
     ! Parameters for Kessler/Lin type autoconversion
     ,timekessl=0.001 &
     ,betakessi=0.001 &
     ,qll0=0.001 &
     ,qli0=0.0001 &
     ! Diagnostic division between rain, snow and graupel
     ,tuprsg=268. &
     ,tdnrsg=253. &
     ,tupsg=283. & ! Following Khairoutdinov and Randall
     ,tdnsg=223.   ! Following Khairoutdinov and Randall

   ! Fields related to ice-liquid partitioning and slope of distribution
   real,allocatable,dimension(:,:,:) :: ilratio,rsgratio,sgratio,lambdar,lambdas,lambdag
   ! Density-corrected A coefficients for terminal velocity
   real,allocatable,dimension(:) :: ccrz,ccsz,ccgz
   real,allocatable,dimension(:) :: ccrz2,ccsz2,ccgz2
  end module modmicrodata
