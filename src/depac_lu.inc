!#################################################################
!
! depac_lu.inc file for definition of landuse classes and 
!   land use dependent settings
!
!
! DEPAC uses the "LBG" landuse clasification, 
! while LOTOS-EUROS reads a map following "Smiatek" clasification:
!
! LBG class                       composed from lotos-euros landuse
! ------------------------------  ----------------------------------------
!  1 = grass                           3
!  2 = arable                       2
!  3 = permanent crops              2
!  4 = coniferous forest                     5 6 
!  5 = deciduous forest                   4    6    8         12 13
!  6 = water                                     7  8      11
!  7 = urban                      1
!  8 = other (short grassy area)                         10
!  9 = desert                                         9
!
!
! LBG class                       composed from glc landuse
! ------------------------------  ----------------------------------------
!  1 = grass                                               11 12 13 14 15       18
!  2 = arable                                                             16    18
!  3 = permanent crops                                                    16 17 18 
!  4 = coniferous forest                  4  5 6  7 8 9 10 11 12       15    17
!  5 = deciduous forest           1  2  3      6  7 8 9 10 11 12    14 15    17
!  6 = water                                      7 8                  15             20 21
!  7 = urban                                                                                22
!  8 = other (short grassy area)                                    14 15                      23 
!  9 = desert                                              11                      19
!  10 = ice                       (not used yet)
!  11 = wheat                     (only used for O3 deposition = landuse dependent)
!  12 = beech                     (only used for O3 deposition = landuse dependent)
!  13 = spruce                    (only used for O3 deposition = landuse dependent)
!  14 = semi-natural grass        (only used for O3 deposition = landuse dependent)
!
!
!  References
!  ----------
!
!    Description of the DEPAC module : Dry deposition modelling with DEPAC_GCN2010
!    M.C. van Zanten, F.J. Sauter, R.J. Wichink Kruit, J.A. van Jaarsveld, W.A.J. van Pul
!    2010 ; RIVM rapport 680180001 ; http://www.rivm.nl/bibliotheek/rapporten/680180001.pdf
!
!
!  Proposal for new convention
!  ---------------------------
!
!    2012-03 Roy Wichink Kruit, Arjo Segers (TNO)
!
!    New convention for canopy .
!    Consider a vegetation of trees:
!       vegetation_height = 20  ! m
!    Assume that the 'orography' in the model is at 2/3 of this;
!    the distance from the actual soil to the orogrpahy is called:
!       displacement_height = 2.0/3.0 * vegetation_height
!    The "canopy" (="cover") top is the height of the remaining cover.
!
!    A rule-of-thumb for the roughness length is:
!        z0 = 10% veg_height
!
!
!
!        --------------------------------------------
!        ^             ( )( )( )( )    ^
!        |             ( )( )( )( )    |
!        |             ( )( )( )( )    |
!        |             ( )( )( )( )    |
!        |             ( )( )( )( )    | zcanopy_top = 1/3 veg_height
!        |             ( )( )( )( )    |
!        |             ( )( )( )( )    |
!        |             ( )( )( )( )    |
!        |             ( )( )( )( )    |                                --+---- surface concentration model (z0 + 2.5 m)
!        |             ( )( )( )( )    |                                  |
!        |             ( )( )( )( )    |                                --+------
!        |             ( )( )( )( )    v                                  | z0
!        | veg_height  ( )( )( )( ) ---------------- model surface -------+------         
!        |             ( )( )( )( )   ^                                        
!        |             ( )( )( )( )   | displacement_height = 2/3 veg_height   
!        v             || |} || ||    v                                         
!       =========================================================================  soil                     
!                                                  ^                               
!                                                  v  orography
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  sea level
!
!    All parameters summarized:
!
!        ------------------  ----------------  ----------  ------- ------   --------
!        landuse class       includes           veg.h.     disp.h. zcanopy   z0 2)
!                                                          2/3 vh  1/3 vh
!        ------------------  ----------------  ----------  ------- ------   --------
!        grass                                  0.30   1)   0.20    0.10      0.03   
!        arable              corn               2.500  1)   1.667   0.833     0.25 
!        permanent_crops     fruit trees        2.500  1)   1.667   0.833     0.25 
!        conifereous_forest                    17.4    2)  11.600   5.800     2.00
!        deciduous_forest                      19.4    2)  12.933   6.467     2.00
!        water                                         0)                     0.002
!        urban                                         0)                     3.00
!        other               tundra, savanne?   2.5    1)   1.667   0.833     0.25
!        desert                                        0)                     0.10
!        ice                                           0)                     0.0005
!        wheat                                         0)                     0.25
!        beech                                         0)                     2.00
!        spruce                                        0)                     2.00
!        semi_natural_grass                            0)                     0.25
!        ------------------  ----------------  ----------  ------- ------   --------
!
!        Notes:
!          0) no vegetation
!          1) use veg.height = 10 * z0
!          2) from exisiting code ...
!
!
!  History
!  -------
!
!  2011-02-21 Arjo Segers, TNO
!    Set minimum canopy height to 'z0m_depac', otherwise problems with 
!    negative range [z0,ctop] in stability module.
!    And the idea of trees being higher than z0 is not that stange actually ...
!  2011-09 Arjo Segers, TNO
!    Collected content DEPAC related content from 'landuse' and
!    'depoparameters' module into new module 'landuse_depac' .
!  2015-03 Roy Wichink Kruit, RIVM
!    Created inc file to set DEPAC constants
!
!#################################################################


! depac classes:
integer, parameter  ::  ilu_grass               = 1
integer, parameter  ::  ilu_arable              = 2
integer, parameter  ::  ilu_permanent_crops     = 3
integer, parameter  ::  ilu_coniferous_forest   = 4
integer, parameter  ::  ilu_deciduous_forest    = 5
integer, parameter  ::  ilu_water_sea           = 6
integer, parameter  ::  ilu_urban               = 7
integer, parameter  ::  ilu_other               = 8  ! e.g. short grassy
integer, parameter  ::  ilu_desert              = 9
integer, parameter  ::  ilu_ice                 = 10
integer, parameter  ::  ilu_savanna             = 11
integer, parameter  ::  ilu_tropical_forest     = 12
integer, parameter  ::  ilu_water_inland        = 13
integer, parameter  ::  ilu_mediterrean_scrub   = 14
integer, parameter  ::  ilu_semi_natural_veg    = 15
integer, parameter  ::  ilu_wheat               = 16 ! only for o3fl project
integer, parameter  ::  ilu_beech               = 17 ! only for o3fl project
integer, parameter  ::  ilu_spruce              = 18 ! only for o3fl project

integer, parameter  ::  nlu_normal              = 15 ! number of landuses
integer, parameter  ::  nlu_ozone_specials      = 18

! names:
 character(len=*), parameter  ::  lu_name(nlu_ozone_specials) = &
          (/ 'grass             ', &
             'arable            ', &
             'permanent_crops   ', &
             'conifereous_forest', &
             'deciduous_forest  ', &
             'water_sea         ', &
             'urban             ', &
             'other             ', &
             'desert            ', &
             'ice               ', &
             'savanna           ', &
             'tropical_forest   ', &
             'water_inland      ', &
             'mediterrean_scrub ', &
             'semi_natural_veg  ', &
             'wheat             ', &
             'beech             ', &
             'spruce            ' /)
             
  ! abreviations:
  character(len=*), parameter  ::  lu_name_abbr(nlu_ozone_specials) = &
          (/ 'grs', &  ! grass             
             'ara', &  ! arable            
             'crp', &  ! permanent_crops   
             'cnf', &  ! conifereous_forest
             'dec', &  ! deciduous_forest  
             'wat', &  ! water sea           
             'urb', &  ! urban             
             'oth', &  ! other             
             'dsr', &  ! desert 
             'ice', &  ! ice
             'sav', &  ! savanna
             'trf', &  ! tropical forest
             'wai', &  ! water inland
             'med', &  ! mediterrean scrub
             'sem', &  ! semi natural vegetation
             'whe', &  ! wheat
             'bee', &  ! beeches
             'spr' /)  ! spruces

! set temperatures per land use for F_temp
real, dimension(nlu_ozone_specials), parameter :: Tmin =   (/12.0, 12.0,  12.0,  0.0,  0.0, -999., -999., 12.0, -999., -999., 12.0,  0.0, -999., 12.0,  8.0, 12.0,  0.0,  0.0 /)
real, dimension(nlu_ozone_specials), parameter :: Topt =   (/26.0, 26.0,  26.0, 18.0, 20.0, -999., -999., 26.0, -999., -999., 26.0, 20.0, -999., 26.0, 24.0, 26.0, 16.0, 14.0 /)
real, dimension(nlu_ozone_specials), parameter :: Tmax =   (/40.0, 40.0,  40.0, 36.0, 35.0, -999., -999., 40.0, -999., -999., 40.0, 35.0, -999., 40.0, 39.0, 40.0, 33.0, 35.0 /)
  
! set F_min:
real, dimension(nlu_ozone_specials), parameter :: F_min   =(/0.01, 0.01,  0.01, 0.1,  0.1,   -999., -999.,0.01, -999.,-999.,0.01,0.1,-999.,0.01, 0.04,0.01, 0.13, 0.16/)

! set alpha for f_light (4.57 is conversion factor from 1./(mumol m-2 s-1) naar W m-2
real, dimension(nlu_ozone_specials), parameter :: alpha   =(/0.009,0.009, 0.009,0.006,0.006, -999., -999.,0.009,-999.,-999.,0.009,0.006,-999.,0.009,0.008,0.0105,0.006,0.010/)*4.57
  
! Set maximal conductance (m/s)
! (R T/P) = 1/41000 mmol/m3 is given for 20 deg C to go from  mmol O3/m2/s to m/s
! Could be refined to a function of T and P. in Jones
real, dimension(nlu_ozone_specials), parameter :: g_max   =(/270., 300.,  300., 140., 150.,  -999., -999.,270., -999.,-999.,270., 150.,-999.,300., 422.,500., 162., 135./)/41000

! Set max, min for vapour pressure deficit vpd;
! note: names max and min seemingly the other way around, but vpd is a deficit!
real, dimension(nlu_ozone_specials), parameter :: vpd_max =(/1.3,  0.9,   0.9,  0.5,  1.0,   -999., -999.,1.3,  -999.,-999.,1.3,1.0,  -999.,0.9, 2.8,1.2,  1.0,  0.5/) 
real, dimension(nlu_ozone_specials), parameter :: vpd_min =(/3.0,  2.8,   2.8,  3.0,  3.25,  -999., -999.,3.0,  -999.,-999.,3.0,3.25, -999.,2.8, 4.5,3.2,  3.10, 3.0/) 

! b = empirical constant for computation of rinc (in canopy resistance) (= 14 m-1 or -999 if not applicable)
! h = vegetation height (m)                             grass arabl  crops conif decid   water  urba   othr  desr   ice  sav  trf   wai   med semi wheat beech spruce
real, dimension(nlu_ozone_specials), parameter :: b = (/ -999,  14,    14,   14,   14,    -999,  -999, -999, -999, -999, -999, 14,  -999, 14, 14,   14,   14,   14 /)
real, dimension(nlu_ozone_specials), parameter :: h = (/ -999,  1,     1,    20,   20,    -999,  -999, -999, -999, -999, -999, 20,  -999, 1 ,  1,    1,   20,   20 /)  
             
! Define type for LAI (leaf area index) parameters:
type laitype
   integer     :: sgs50     ! start growing season at 50 degrees latitude (days)
   real        :: dsgs      ! shift in start growing season (days/degree latitude)
                            ! dsgs > 0 -> start-growing-season earlier in the south, later in the north
   integer     :: egs50     ! end growing season at 50 degrees latitude (days)
   real        :: degs      ! shift in end growing season (days/degree latitude)
                            ! degs < 0 -> end-growing-season later in the south, earlier in the north
   real        :: laimin    ! leaf area index at start and end of growing season (m2 leaf/m2 ground surface)
                            ! xxx outside growing season LAI = 0
   real        :: laimax    ! maximal leaf area index (m2 leaf/m2 ground surface)
   integer     :: s_lai_len ! length of starting phase of LAI (days)
   integer     :: e_lai_len ! length of end phase of LAI (days)
end type laitype
! lai_par: Note that for the LAI of arable land and permanent crops the EMEP class "root crop" has been used;
!          for other parameters, "temp_crop" has been used.
!
!          sgs50    dsgs  egs50    degs  laimin  laimax  s_lai_len  e_lai_len
type(laitype), parameter :: lai_par(nlu_ozone_specials) = (/ &                      ! Emberson   DEPAC
    laitype(    0,    0.0,   366,    0.0,    2.0,    3.5,       140,       135), &  ! grass      1. grassland
    laitype(  130,    0.0,   250,    0.0,    0.0,    4.2,        35,        65), &  ! root crop  2. arable land
    laitype(  130,    0.0,   250,    0.0,    0.0,    4.2,        35,        65), &  ! root crop  3. permanent crops
    laitype(    0,    0.0,   366,    0.0,    5.0,    5.0,         1,         1), &  ! temp_conif 4. coniferous forest
    laitype(  100,    1.5,   307,   -2.0,    0.0,    4.0,        20,        30), &  ! temp_decid 5. deciduous forest
    laitype( -999, -999.0,  -999, -999.0, -999.0, -999.0,      -999,      -999), &  ! water      6. water_sea
    laitype( -999, -999.0,  -999, -999.0, -999.0, -999.0,      -999,      -999), &  ! urban      7. urban
    laitype(    0,    0.0,   366,    0.0,    2.0,    3.5,       140,       135), &  ! grass      8. other
    laitype( -999, -999.0,  -999, -999.0, -999.0, -999.0,      -999,      -999), &  ! desert     9. desert
    laitype( -999, -999.0,  -999, -999.0, -999.0, -999.0,      -999,      -999), &  ! water      10.ice 
    laitype(    0,    0.0,   366,    0.0,    2.0,    3.5,       140,       135), &  ! grass      11.Savanna
    laitype(  100,    1.5,   307,   -2.0,    0.0,    4.0,        20,        30), &  ! temp_decid 12.tropical forest
    laitype( -999, -999.0,  -999, -999.0, -999.0, -999.0,      -999,      -999), &  ! water      13.water_inland
    laitype(  130,    0.0,   250,    0.0,    0.0,    4.2,        35,        65), &  ! root crop  14.mediterrean scrub
    laitype(    0,    0.0,   366,    0.0,    2.0,    3.5,       140,       135), &  ! grass      15.semi_natural grassland        
    laitype(  130,    0.0,   250,    0.0,    0.0,    4.2,        35,        65), &  ! root crop  16.wheat 
    laitype(  100,    1.5,   307,   -2.0,    0.0,    4.0,        20,        30), &  ! temp_decid 17.beech          
    laitype(    0,    0.0,   366,    0.0,    5.0,    5.0,         1,         1) /)  ! temp_conif 18.spruce       


  !--end Emberson specific declarations------------


!#################################################################
! component numbers:
integer, parameter  ::  icmp_o3   = 1
integer, parameter  ::  icmp_so2  = 2
integer, parameter  ::  icmp_no2  = 3
integer, parameter  ::  icmp_no   = 4
integer, parameter  ::  icmp_nh3  = 5
integer, parameter  ::  icmp_co   = 6
integer, parameter  ::  icmp_no3  = 7
integer, parameter  ::  icmp_hno3 = 8
integer, parameter  ::  icmp_n2o5 = 9
integer, parameter  ::  icmp_h2o2 = 10
integer, parameter  ::  icmp_vbs_pog1 = 11
integer, parameter  ::  icmp_vbs_pog2 = 12
integer, parameter  ::  icmp_vbs_pog3 = 13
integer, parameter  ::  icmp_vbs_pog4 = 14
integer, parameter  ::  icmp_vbs_pog5 = 15
integer, parameter  ::  icmp_vbs_pog6 = 16
integer, parameter  ::  icmp_vbs_pog7 = 17
integer, parameter  ::  icmp_vbs_pog8 = 18
integer, parameter  ::  icmp_vbs_pog9 = 19
integer, parameter  ::  icmp_vbs_sisog1 = 20
integer, parameter  ::  icmp_vbs_sisog2 = 21
integer, parameter  ::  icmp_vbs_sisog3 = 22
integer, parameter  ::  icmp_vbs_sisog4 = 23
integer, parameter  ::  icmp_vbs_sisog5 = 24
integer, parameter  ::  icmp_vbs_sisog6 = 25
integer, parameter  ::  icmp_vbs_sisog7 = 26
integer, parameter  ::  icmp_vbs_sisog8 = 27
integer, parameter  ::  icmp_vbs_asog1 = 28
integer, parameter  ::  icmp_vbs_asog2 = 29
integer, parameter  ::  icmp_vbs_asog3 = 30
integer, parameter  ::  icmp_vbs_asog4 = 31
integer, parameter  ::  icmp_vbs_asog5 = 32
integer, parameter  ::  icmp_vbs_asog6 = 33
integer, parameter  ::  icmp_vbs_bsog1 = 34
integer, parameter  ::  icmp_vbs_bsog2 = 35
integer, parameter  ::  icmp_vbs_bsog3 = 36
integer, parameter  ::  icmp_vbs_bsog4 = 37
integer, parameter  ::  icmp_vbs_bsog5 = 38
integer, parameter  ::  icmp_vbs_bsog6 = 39

integer, parameter  ::  ncmp      = 39  ! number of components

! diffusion coefficients:
real, parameter     :: dwat = 0.21e-4   ! diffusion coefficient of water vapour (m2/s)
real, parameter     :: dO3  = 0.13e-4   ! diffusion coefficient of ozon (m2/s)

! Soil resistance (numbers matched with lu_classes and component numbers)
real, parameter     :: rsoil(nlu_ozone_specials,ncmp) = reshape( (/ &
                     ! grs    ara    crp    cnf    dec    wat    urb    oth    des    ice    sav    trf    wai    med    sem    whe    bee    spr
                      1000.,  200.,  200.,  200.,  200., 2000.,  400., 1000., 2000., 2000., 1000.,  200., 2000.,  200.,  400.,  200.,  200.,  200., & ! O3
                      1000., 1000., 1000., 1000., 1000.,   10., 1000., 1000., 1000.,  500., 1000., 1000.,   10., 1000., 1000., 1000., 1000., 1000., & ! SO2
                      1000., 1000., 1000., 1000., 1000., 2000., 1000., 1000., 1000., 2000., 1000., 1000., 2000., 1000., 1000., 1000., 1000., 1000., & ! NO2
                      -999., -999., -999., -999., -999., 2000., 1000., -999., 2000., 2000., -999., -999., 2000., -999., -999., -999., -999., -999., & ! NO
                       100.,  100.,  100.,  100.,  100.,   10.,  100.,  100.,  100., 1000.,  100.,  100.,   10.,  100.,  100.,  100.,  100.,  100., & ! NH3
                      -999., -999., -999., -999., -999., 2000., 1000., -999., 2000., 2000., -999., -999., 2000., -999., -999., -999., -999., -999., & ! CO 
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! NO3 
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! HNO3
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! N2O5
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! H2O2   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_POG1   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_POG2   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_POG3   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_POG4   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_POG5   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_POG6   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_POG7   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_POG8   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_POG9   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_SISOG1   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_SISOG2   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_SISOG3   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_SISOG4   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_SISOG5   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_SISOG6   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_SISOG7   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_SISOG8   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_ASOG1   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_ASOG2   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_ASOG3   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_ASOG4   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_ASOG5   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_ASOG6   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_BSOG1   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_BSOG2   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_BSOG3   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_BSOG4   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., & ! VBS_BSOG5   
                      -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999. & ! VBS_BSOG6   
                      /), (/nlu_ozone_specials,ncmp/) )

                                             !  o3       so2       no2        no       nh3        co    no3   hno3   n2o5   h2o2  vbs_pog1  vbs_pog2  vbs_pog3  vbs_pog4  vbs_pog5  vbs_pog6  vbs_pog7  vbs_pog8  vbs_pog9 vbs_sisog1  vbs_sisog2  vbs_sisog3  vbs_sisog4  vbs_sisog5  vbs_sisog6  vbs_sisog7  vbs_sisog8 vbs_asog1 vbs_asog2 vbs_asog3 vbs_asog4 vbs_asog5 vbs_asog6 vbs_bsog1 vbs_bsog2 vbs_bsog3 vbs_bsog4 vbs_bsog5 vbs_bsog6  

real, parameter     :: diffc(ncmp)        = (/ dO3 , 0.11e-4 , 0.13e-4 , 0.16e-4 , 0.21e-4 , 0.16e-4 , -999., -999.,  -999.,  -999.,  -999. , -999.,  -999., &
                                              -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999., &
                                              -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999. /)

integer, parameter  :: ipar_snow(ncmp)    = (/   1 ,       2 ,       1 ,       1 ,       2 ,       1 , -999 , -999 , -999 , -999 ,  -999  , -999 ,    -999 , &
                                              -999 ,  -999 ,  -999 ,  -999 ,  -999 ,  -999 ,  -999 ,  -999 ,  -999 ,  -999 ,  -999 ,  -999 ,  -999 , &
                                              -999 ,  -999 ,  -999 ,  -999 ,  -999 ,  -999 ,  -999 ,  -999 ,  -999 ,  -999 ,  -999 ,  -999 ,  -999  /)

real, parameter     :: rsoil_wet(ncmp)    = (/2000.,    10.,    2000.,    -999.,      10.,    -999., -999., -999., -999., -999.,  -999. , -999.,    -999., &
                                              -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999., &
                                              -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999. /)

real, parameter     :: rsoil_frozen(ncmp) = (/2000.,     500.,    2000.,    -999.,    1000.,    -999., -999., -999., -999., -999.,  -999. , -999.,    -999., &
                                              -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999., &
                                              -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999.,  -999. /)
