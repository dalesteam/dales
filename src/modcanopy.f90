!> \file modcanopy.f90
!!  Canopy parameterization
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
!  Copyright 1993-2014 Delft University of Technology, Wageningen University, Utrecht University, KNMI, MPIC
!
module modcanopy
  implicit none
  save

  ! Namoptions
  logical :: lcanopy   = .false.       !< Switch to enable canopy representation
  integer :: ncanopy   = 10            !< Amount of layers to represent the canopy
  real    :: cd        = 0.15          !< Drag coefficient in the canopy
  real    :: lai_can   = 2             !< Leaf Area Index (or actually plant area index) of the canopy
  logical :: lpaddistr = .false.       !< Switch to customize the general plant area density distribution (at half levels)
  integer :: npaddistr = 11            !< (if lpaddistr): number of half levels for prescribed general plant area density distribution

  logical :: wth_total = .false.       !< Switch: prescribed SH flux is added to surface flux if .false., it contains (the effect of) the surface flux if .true.
  logical :: wqt_total = .false.       !< Switch: prescribed LE flux is added to surface flux if .false., it contains (the effect of) the surface flux if .true.
  logical :: wsv_total(100) = .false.  !< Switch: prescribed sv flux is added to surface flux if .false., it contains (the effect of) the surface flux if .true.

  real    :: wth_can  = 0.0            !< prescribed SH canopy flux
  real    :: wqt_can  = 0.0            !< prescribed LE canopy flux
  real    :: wsv_can(100) = 0.0        !< prescribed scalar canopy flux

  real    :: wth_alph = 0.6            !< Decay constant for flux as function of the vertically integrated PAI (from canopy top)
  real    :: wqt_alph = 0.6            !< Decay constant for flux as function of the vertically integrated PAI (from canopy top)
  real    :: wsv_alph(100) = 0.6       !< Decay constant for flux as function of the vertically integrated PAI (from canopy top)

  ! Fields
  real, allocatable :: padfactor(:)    !< prescribed weighing factor for plant area density
  real, allocatable :: ppad(:)         !< resulting prescribed plant area density
  real, allocatable :: zpad(:)         !< heights of prescribed plant area density
  real, allocatable :: padtemp(:)      !< temporary plant area density used for calculations
  real, allocatable :: padf(:)         !< plant area density field full level
  real, allocatable :: padh(:)         !< plant area density field half level
  real, allocatable :: pai(:)          !< plant area index of the column starting in this grid cell up to the canopy top

  real              :: f_lai_h         !< average plant area density [m2/m2 / m]
  
  ! for canopy energy balance
  
  ! Namoptions
  logical :: lcanopyeb     = .false.   !< Switch to enable canopy surface energy balance per vertical level
  
  real    :: leaf_eps      = 0.95      !< Emissivity of leaves in the LW
  real    :: transpiretype = 1         !< type of transpirer(1=hypostomatous, 2=amphistomatous,
                                       !< 1.25=hypostomatous with some transpiration through cuticle)
  real     :: lwidth       = 0.02      !< !leaf width/shoot diameter [m]                              
  real     :: llength      = 0.1       !< leaf/shoot length [m]                                      
  logical  :: lrelaxgc_can = .false.   !< Switch to delay plant response at canopy.Timescale is equal to 1/kgc_can
  real     :: kgc_can      = 0.00113   !< Standard stomatal response rate at canopy (corresponding to a time scale of 14.75 m min.) [1/s]
  logical  :: gccan_old_set = .false.  !< Only apply relaxing function to canopy after initial gc is calculated once for surface
  logical  :: lrelaxci_can = .false.   !< Switch to delay internal CO2 concentration in canopy leafs.Timescale is equal to 1/kgc_can
  real     :: kci_can      = 0.00113   !< Standard internal CO2 concentration response rate at canopy (corresponding to a time scale of 14.75 m min.) [1/s]
  logical  :: cican_old_set = .false.  !< Only apply relaxing function to canopy after initial ci is calculated once for surface

  real, allocatable :: cPARleaf_shad   (:)     !< PAR at shaded leaves per canopy vertical level [W m-2]
  real, allocatable :: cPARleaf_sun    (:,:)   !< PAR at sunny leaves per vertical level and leaf orientation [W m-2]
  real, allocatable :: cPARleaf_allsun (:)     !< PAR at sunny leaves per vertical level averaged over all leaf orientations [W m-2]
  real, allocatable :: cfSL            (:)     !< Fraction of sunlit leaves per evrtical level [-]
  real, allocatable :: iLAI_can        (:)     !< LAI from that vertical level up to canopy top [m2_leaf m-2_ground]
  real, allocatable :: humidairpa      (:)     !< Atmospheric vapor pressure inside the canopy [Pa]
  real, allocatable :: windsp          (:)     !< Horizontal windspeed inside the canopy  [m s-1]
  real, allocatable :: LWin_leafshad   (:)     !< Inwards LW radiation at shaded leaves [W m-2]
  real, allocatable :: LWout_leafshad  (:)     !< Outwards LW radiation at shaded leaves [W m-2]
  real, allocatable :: LWnet_leafshad  (:)     !< Net LW radiation at shaded leaves [W m-2]
  real, allocatable :: tleaf_shad      (:)     !< Leaf temperature of shaded leaves [K]
  real, allocatable :: rs_leafshad     (:)     !< Stomatal resistance of shaded leaves [s m-1]
  real, allocatable :: rb_leafshad     (:)     !< Leaf boundary layer resistance for shaded leaves[s m-1]
  real, allocatable :: sh_leafshad     (:)     !< Sensible heat flux at shaded leaves [W m-2_leaf]
  real, allocatable :: le_leafshad     (:)     !< Latent heat flux at shaded leaves [W m-2_leaf]
  real, allocatable :: An_leafshad     (:)     !< Net carbon uptake at shaded leaves [mg C s-1 m-2_leaf)]
  real, allocatable :: LWin_leafsun    (:)     !< Inwards LW radiation at sunny leaves [W m-2]   
  real, allocatable :: LWout_leafsun   (:)     !< Outwards LW radiation at sunny leaves [W m-2_leaf]
  real, allocatable :: LWnet_leafsun   (:)     !< Net LW radiation at sunny leaves [W m-2]
  real, allocatable :: tleaf_sun       (:)     !< Leaf temperature of sunny leaves [K]
  real, allocatable :: rs_leafsun      (:)     !< Stomatal resistance of shaded leaves [s m-1]
  real, allocatable :: rb_leafsun      (:)     !< Leaf boundary layer resistance for shaded leaves[s m-1]
  real, allocatable :: sh_leafsun      (:)     !< Sensible heat flux at shaded leaves [W m-2_leaf]
  real, allocatable :: le_leafsun      (:)     !< Latent heat flux at shaded leaves [W m-2_leaf]
  real, allocatable :: An_leafsun      (:)     !< Net carbon uptake at shaded leaves [mg C s-1 m-2_leaf)]
  real, allocatable :: sh_can          (:,:,:) !< Total sensible heat flux per canopy vertical level [W m-3]
  real, allocatable :: le_can          (:,:,:) !< Total latent heat flux per canopy vertical level[W m-3]
  real, allocatable :: Fco2_can        (:,:,:) !< Total CO2 flux per canopy vertical level [mg CO2 m-3 s-1]
  real, allocatable :: S_theta         (:,:,:) !< Temperature source/sink due to canopy[K s-1]
  real, allocatable :: S_qt            (:,:,:) !< Moisture source/sink due to canopy[Kg_w Kg_a-1 s-1]
  real, allocatable :: S_co2           (:,:,:) !< CO2 source/sink due to canopy [ppb s-1]
  real, allocatable :: gcc_leafshad    (:,:,:) !< carbon stomatal conductance at shaded leaves [m s-1]
  real, allocatable :: gcc_leafshad_old(:,:,:) !< carbon stomatal conductance at shaded leaves in previous timestep[m s-1]
  real, allocatable :: gcc_leafsun     (:,:,:) !< carbon stomatal conductance at sunny leaves [m s-1]
  real, allocatable :: gcc_leafsun_old (:,:,:) !< carbon stomatal conductance at sunny leaves in previous timestep[m s-1]
  real, allocatable :: ci_leafshad     (:,:,:) !< internal carbon concentration at shaded leaves [mg m-3]
  real, allocatable :: ci_leafshad_old (:,:,:) !< internal carbon concentration at shaded leaves in previous timestep[mg m-3]
  real, allocatable :: ci_leafsun      (:,:,:) !< internal carbon concentration at sunny leaves  [mg m-3]
  real, allocatable :: ci_leafsun_old  (:,:,:) !< internal carbon concentration at sunny leaves in previous timestep[ [mg m-3]
  real, allocatable :: tskin_can       (:,:) !< tskin at canopy top
  real, allocatable :: albdir_can      (:,:) !< effective direct SW albedo by canopy
  real, allocatable :: albdif_can      (:,:) !< effective diffuse SW albedo by canopy
  real, allocatable :: albsw_can       (:,:) !< effective global SW albedo by canopy
  real, allocatable :: swdir_can       (:) !< direct SW through canopy
  real, allocatable :: swdif_can       (:) !< diffuse SW through canopy
  real, allocatable :: swu_can         (:) !< upwards SW through canopy
  real, allocatable :: lwd_can         (:) !< diffuse SW through canopy
  real, allocatable :: lwu_can         (:) !< diffuse SW through canopy
contains
!-----------------------------------------------------------------------------------------
  SUBROUTINE initcanopy
    use modmpi,      only : myid, mpi_logical, mpi_integer, my_real, comm3d, mpierr
    use modglobal,   only : kmax, ifnamopt, fname_options, ifinput, cexpnr, zh, dzh, dzf,ih,i1,jh,j1,i2,j2
    use modsurfdata, only : nangle_gauss
    implicit none

    integer ierr, k, kp
    character(80) readstring

    namelist/NAMCANOPY/ lcanopy, ncanopy, cd, lai_can, lpaddistr, npaddistr, &
                        wth_total, wqt_total, wsv_total, wth_can, wqt_can, wsv_can, &
                        wth_alph, wqt_alph, wsv_alph, &
                        lcanopyeb,lwidth,llength,transpiretype,leaf_eps,&
                        lrelaxgc_can,kgc_can,lrelaxci_can,kci_can

    if(myid==0) then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMCANOPY,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMCANOPY'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMCANOPY'
      endif
      write(6 ,NAMCANOPY)
      close(ifnamopt)

      ncanopy = min(ncanopy,kmax)
    endif

    call MPI_BCAST(lcanopy      ,   1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(ncanopy      ,   1, mpi_integer , 0, comm3d, mpierr)
    call MPI_BCAST(cd           ,   1, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(lai_can      ,   1, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(lpaddistr    ,   1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(npaddistr    ,   1, mpi_integer , 0, comm3d, mpierr)
    call MPI_BCAST(wth_total    ,   1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(wqt_total    ,   1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(wsv_total    , 100, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(wth_can      ,   1, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(wqt_can      ,   1, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(wsv_can      , 100, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(wth_alph     ,   1, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(wqt_alph     ,   1, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(wsv_alph     , 100, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(lcanopyeb    ,   1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(lwidth       ,   1, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(llength      ,   1, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(transpiretype,   1, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(leaf_eps     ,   1, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(lrelaxgc_can ,   1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(kgc_can      ,   1, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(lrelaxci_can ,   1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(kci_can      ,   1, my_real     , 0, comm3d, mpierr)

    if (.not. (lcanopy)) return

    if (.not. lpaddistr) npaddistr = 11

    allocate(padfactor (npaddistr))
    allocate(ppad      (npaddistr))
    allocate(zpad      (npaddistr))
    allocate(padtemp   (npaddistr))
    allocate(padf      (ncanopy  ))
    allocate(padh      (ncanopy+1))
    allocate(pai       (ncanopy+1))

    ! Determination of padfactor: relative weighing of plant area distribution inside canopy; equidistant from surface to canopy top
    if (lpaddistr) then  !< Profile prescribed by user in the file paddistr.inp.<expnr>
      if (myid==0) then
        open (ifinput,file='paddistr.inp.'//cexpnr)

        do k=1,npaddistr
          read (ifinput,'(a80)') readstring

          do while (readstring(1:1)=='#')     ! Skip the lines that are commented (like headers)
            read (ifinput,'(a80)') readstring
          end do

          read(readstring,*) padfactor(k)
        end do

        close(ifinput)

        !And now, weigh it such that the array of averages of 2 adjacent padfactor values is on average 1 (just in case the user's array does not fulfill that criterion)
        padfactor = padfactor * (npaddistr-1) * 2 / sum(padfactor(1:(npaddistr-1))+padfactor(2:npaddistr))

        write(*,*) 'Prescribed weighing for plant area density from surface to canopy top (equidistant); normalized if necessary'
        do k=1,npaddistr
          write (*,*) padfactor(k)
        end do
      endif

      call MPI_BCAST(padfactor, npaddistr, my_real , 0, comm3d, mpierr)

    else                 !< Standard profile fron Ned Patton
      padfactor = (/ 0.4666666666666667, &
                     0.5307086614173228, &
                     0.6792650918635170, &
                     0.9548556430446193, &
                     1.3154855643044620, &
                     1.5490813648293960, &
                     1.5916010498687660, &
                     1.5275590551181100, &
                     1.2944881889763780, &
                     0.3236220472440945, &
                     0.0000000000000000  /)
    endif
    f_lai_h = lai_can / zh(1+ncanopy) ! LAI of canopy divided by height of the top of the canopy
   
    ppad    = f_lai_h * padfactor ! prescribed PAD-values
    do k=1,npaddistr
      zpad(k) = zh(1+ncanopy) * real(k-1)/real(npaddistr-1)
    end do

    !interpolate the PAD values to the LES grid
    call spline(zpad,ppad,npaddistr,padtemp)
    do k=1,(1+ncanopy)
      call splint(zpad,ppad,padtemp,npaddistr,zh(k),padh(k))
    end do

    ! Interpolate plant area (index) density to full levels
    do k=1,ncanopy
      kp      = k+1
      padf(k) = ( dzh(kp) * padh(k) + dzh(k) * padh(kp) ) / ( dzh(k) + dzh(kp) )
    end do

    ! Vertically integrate the plant area density to arrive at plant area index
    pai = 0.0
    do k=ncanopy,1,-1
      pai(k) = pai(k+1) + dzf(k) * padf(k)
    end do

    if (.not. (lcanopyeb)) return
    
    allocate(cPARleaf_shad(ncanopy))
    allocate(cPARleaf_allsun(ncanopy))
    allocate(cPARleaf_sun(ncanopy,nangle_gauss))
    allocate(cfSL(ncanopy))
    allocate(iLAI_can(ncanopy))
    !allocate(tairk(ncanopy))
    allocate(humidairpa(ncanopy))
    !allocate(PAD(ncanopy))
    allocate(windsp(ncanopy))
 
    allocate(LWin_leafshad(ncanopy))
    allocate(LWout_leafshad(ncanopy))
    allocate(LWnet_leafshad(ncanopy))
    allocate(tleaf_shad(ncanopy))
    allocate(rs_leafshad(ncanopy))
    allocate(rb_leafshad(ncanopy))
    allocate(sh_leafshad(ncanopy))
    allocate(le_leafshad(ncanopy))
    allocate(An_leafshad(ncanopy))
   ! no 3 leaf LEB, otherwise here we shoud add one extra dim for the angles
    allocate(LWin_leafsun(ncanopy))
    allocate(LWout_leafsun(ncanopy))
    allocate(LWnet_leafsun(ncanopy))
    allocate(tleaf_sun(ncanopy))
    allocate(rs_leafsun(ncanopy))
    allocate(rb_leafsun(ncanopy))
    allocate(sh_leafsun(ncanopy))
    allocate(le_leafsun(ncanopy))
    allocate(An_leafsun(ncanopy))
 
    allocate(sh_can(2-ih:i1+ih,2-jh:j1+jh,ncanopy))
    allocate(le_can(2-ih:i1+ih,2-jh:j1+jh,ncanopy))
    allocate(Fco2_can(2-ih:i1+ih,2-jh:j1+jh,ncanopy))
    allocate(S_theta(2-ih:i1+ih,2-jh:j1+jh,ncanopy))
    allocate(S_qt(2-ih:i1+ih,2-jh:j1+jh,ncanopy))
    allocate(S_co2(2-ih:i1+ih,2-jh:j1+jh,ncanopy))
    allocate(gcc_leafshad(2-ih:i1+ih,2-jh:j1+jh,ncanopy))
    allocate(gcc_leafsun(2-ih:i1+ih,2-jh:j1+jh,ncanopy))
    allocate(gcc_leafshad_old(2-ih:i1+ih,2-jh:j1+jh,ncanopy))
    allocate(gcc_leafsun_old(2-ih:i1+ih,2-jh:j1+jh,ncanopy))
    allocate(ci_leafshad(2-ih:i1+ih,2-jh:j1+jh,ncanopy))
    allocate(ci_leafshad_old(2-ih:i1+ih,2-jh:j1+jh,ncanopy))
    allocate(ci_leafsun(2-ih:i1+ih,2-jh:j1+jh,ncanopy))
    allocate(ci_leafsun_old(2-ih:i1+ih,2-jh:j1+jh,ncanopy))
    allocate(tskin_can (i2,j2)) ! same dims as surface variables
    allocate(albdir_can(i2,j2))
    allocate(albdif_can(i2,j2))
    allocate(albsw_can (i2,j2))
    allocate(swdir_can(ncanopy))
    allocate(swdif_can(ncanopy))
    allocate(swu_can(ncanopy))
    allocate(lwd_can(ncanopy))
    allocate(lwu_can(ncanopy))

    return
    
  end subroutine initcanopy

  subroutine canopy
    use modfields,   only : up,vp,wp,e12p,thlp,qtp,svp
    use modsurfdata, only : thlflux, qtflux, svflux,indCO2
    use modglobal,   only : nsv,i2,j2

    implicit none

    integer n
    real :: zeroar(i2,j2)

    zeroar = 0.0

    if (.not. (lcanopy)) return

!!  Momentum affected by trees
    call canopyu(up)
    call canopyv(vp)
    call canopyw(wp)
!!  TKE affected by trees
    call canopye(e12p)
!!  Emissions of heat, moisture and scalars by trees (effect on center of the grid)
    
    if (lcanopyeb) then
    
    !rename the values already calculated in canopyeb when called at modsurface
      
      thlp(:,:,:ncanopy) = thlp(:,:,:ncanopy) + S_theta(:,:,:)    
      qtp(:,:,:ncanopy) = qtp(:,:,:ncanopy) + S_qt(:,:,:)    
      svp(:,:,:ncanopy,indCO2) = svp(:,:,:ncanopy,indCO2) + S_co2(:,:,:)    
    endif

    if (wth_total) then
      call canopyc(thlp,wth_can,thlflux,wth_alph,pai)
    else
      call canopyc(thlp,wth_can, zeroar,wth_alph,pai)
    endif
    if (wqt_total) then
      call canopyc( qtp,wqt_can,qtflux,wqt_alph,pai)
    else
      call canopyc( qtp,wqt_can,zeroar,wqt_alph,pai)
    endif
    do n=1,nsv
      if (wsv_total(n)) then
        call canopyc(svp(:,:,:,n),wsv_can(n),svflux(:,:,n),wsv_alph(n),pai)
      else
        call canopyc(svp(:,:,:,n),wsv_can(n),       zeroar,wsv_alph(n),pai)
      endif
    end do

    return
  end subroutine canopy

  subroutine exitcanopy
    implicit none

    if (.not. (lcanopy)) return

    deallocate(padfactor)
    deallocate(ppad     )
    deallocate(zpad     )
    deallocate(padtemp  )
    deallocate(padf     )
    deallocate(padh     )
    deallocate(pai      )
    if (.not. (lcanopyeb)) return

    deallocate(cPARleaf_shad)
    deallocate(cPARleaf_allsun)
    deallocate(cPARleaf_sun)
    deallocate(cfSL)
    deallocate(iLAI_can)
      !deallocate(tairk)
    deallocate(humidairpa)
      !deallocate(PAD)
    deallocate(windsp)
 
    deallocate(LWin_leafshad)
    deallocate(LWout_leafshad)
    deallocate(LWnet_leafshad)
    deallocate(tleaf_shad)
    deallocate(rs_leafshad)
    deallocate(rb_leafshad)
    deallocate(sh_leafshad)
    deallocate(le_leafshad)
    deallocate(An_leafshad)
    deallocate(LWin_leafsun)
    deallocate(LWout_leafsun)
    deallocate(LWnet_leafsun)
    deallocate(tleaf_sun)
    deallocate(rs_leafsun)
    deallocate(rb_leafsun)
    deallocate(sh_leafsun)
    deallocate(le_leafsun)
    deallocate(An_leafsun)
 
    deallocate(sh_can)
    deallocate(le_can)
    deallocate(Fco2_can)
    deallocate(S_theta)
    deallocate(S_qt)
    deallocate(S_co2)
    deallocate(gcc_leafshad)
    deallocate(gcc_leafsun)
    deallocate(gcc_leafshad_old)
    deallocate(gcc_leafsun_old)
    deallocate(ci_leafshad)
    deallocate(ci_leafshad_old)
    deallocate(ci_leafsun)
    deallocate(ci_leafsun_old)
    deallocate(tskin_can)
    deallocate(albdir_can)
    deallocate(albdif_can)
    deallocate(albsw_can)
    deallocate(swdir_can)
    deallocate(swdif_can)
    deallocate(swu_can)
    deallocate(lwd_can)
    deallocate(lwu_can)
    return
  end subroutine exitcanopy
  subroutine canopyeb(i,j,ps,rk3coef)                                       ! in
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! This is the main routine for the calculation of a canopy energy balance,
   ! taking into account canopy vertical levels and sunlit and shaded leaves.
   ! Steps in this routine are:
   ! 1-Calculates radiation profiles of both SW and LW
   ! 2-Calculates leaf energy balance at each level for sunlit and shaded leaves
   ! independently, as well as the resultant co2 absorption and heat and
   ! moisture fluxes per leaf area
   ! 3-Calculates the resultant heat, co2 and moisture sources at each gridbox
   ! A big part of the subroutines used in canopyeb were originally developed by Ned
   ! Patton and have been adapted, modified or extended where needed.
   !                                                      Xabier Pedruzo
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use modglobal, only  : j1, i1,cp,rlv,rk3step
    use modsurfdata,only : phitot,weight_g,indCO2,albedo_surf,l3leaves,nangle_gauss,MW_CO2,MW_Air,canopyrad,nuco2q
    use modfields, only  : thl0,rhof,qt0,exnf,u0,v0,presf,svm,tmp0
    use modraddata, only : swdir,swdif,swd,swu,lwu,lwd,tskin_rad,albedo_rad,iradiation,irad_par
   
    implicit none
    
    integer, intent(in) :: i,j
    real, intent(in) :: ps,rk3coef
    integer :: k_can
    real    :: LWin
    real    :: pa ! plant area
    real    :: lw_air,lw_leaflayer 

   !                                                                  !
   !######### STEP 1 - Radiation inside the canopy keeping DALES sign convention#################### 
   !                                                                  !
   
   ! rewrite LAI from LAI per vert layer to integrated LAI from top down to that level
   !do k_can = ncanopy,1,(-1)
   !  iLAI_can(k_can) = sum(LAI_can(1:k_can))
   !end do
      iLAI_can(:) = pai(1:ncanopy)
    call canopyrad(ncanopy,lai_can,iLAI_can,swdir(i,j,ncanopy),swdif(i,j,ncanopy),albedo_surf(i,j),& ! in
                   cPARleaf_shad,cPARleaf_sun,cfSL,                                                & ! out for vegetation
                   albdir_can(i,j),albdif_can(i,j),albsw_can(i,j),                                 & ! out for radiation                         
                   swdir_can(:ncanopy),swdif_can(:ncanopy),swu_can(:ncanopy))                        ! out for radiation
    ! tairk(:) = tmp0(i,j,:)
   ! convert humidity into vapor pressure
    do k_can=1,ncanopy
      humidairpa(k_can) =  watervappres(qt0(i,j,k_can),presf(k_can)) ! check that qt0 is in Kg/Kg
    enddo

    ! we do LW as well:
    do k_can =1,ncanopy
      LWin = unexposedleafLWin(tmp0(i,j,k_can), leaf_eps) ! shouldn we use air eps?
      LWin_leafshad(k_can)  = 2. * LWin
      LWin_leafsun(k_can)   = 0.5 * exposedleafLWin(humidairpa(k_can),tmp0(i,j,k_can)) + 1.5 * LWin ! why .5-1.5, and not1-1
    end do

   
    ! assume for LW that leaves/plants form one horizontal layer in each grid and LW only moves upwards and downwards
    do k_can = 1,ncanopy
      !get plant area (m2leaf per m2 ground in a grid)
      pa = pai(k_can)-pai(k_can+1)
      lw_air = unexposedleafLWin(tmp0(i,j,k_can), leaf_eps)
      !area-weighted average between sunlit and shaded leaves
      lw_leaflayer = (LWout_leafsun(k_can)*cfSL(k_can) + LWout_leafshad(k_can)*(1.-cfSL(k_can)))
      if (pa <= 1.0) then ! area-weighted average between air and leaf
        lwu_can(k_can) = lw_air * (1.0-pa) + lw_leaflayer * pa
      else
        lwu_can(k_can) = lw_leaflayer
      endif  

      ! since we assume that both sides of the leaf are at the same temperature
      lwd_can(k_can) =  -lwu_can(k_can)
    end do  
   !                                                                                                ! 
   !######### STEP 2 - Leaf energy balance per level for sunlit and shaded leaves ####################
   !                                                                                                !
    windsp(:) = sqrt(u0(i,j,1:ncanopy)**2+v0(i,j,1:ncanopy)**2)
   
    do k_can =1,ncanopy
     !shaded
      call leafeb_ags(cPARleaf_shad(k_can), LWin_leafshad(k_can), leaf_eps, transpiretype, lwidth, llength, & ! in
                      tmp0(i,j,k_can), humidairpa(k_can),qt0(i,j,k_can), windsp(k_can), presf(k_can),       & ! in
                      rhof(k_can), svm(i,j,k_can,indCO2), phitot(i,j),                                      & ! in
                      gcc_leafshad_old(i,j,k_can),ci_leafshad_old(i,j,k_can),rk3coef,                       & ! in
                      tleaf_shad(k_can), gcc_leafshad(i,j,k_can), rb_leafshad(k_can),ci_leafshad(i,j,k_can),& ! out
                      sh_leafshad(k_can), le_leafshad(k_can), LWout_leafshad(k_can),An_leafshad(k_can))       ! out
      LWnet_leafshad(k_can) = LWin_leafshad(k_can) - LWout_leafshad(k_can)
      if (lrelaxgc_can) then
        if (gccan_old_set .and. rk3step ==3) then
            gcc_leafshad_old(i,j,k_can) = gcc_leafshad(i,j,k_can)
        else if (.not. gccan_old_set) then
          gcc_leafshad_old(i,j,k_can) = gcc_leafshad(i,j,k_can)
        endif
      endif
      if (lrelaxci_can) then
        if (cican_old_set .and. rk3step ==3) then
            ci_leafshad_old(i,j,k_can) = ci_leafshad(i,j,k_can)
        else if (.not. cican_old_set) then
          ci_leafshad_old(i,j,k_can) = ci_leafshad(i,j,k_can)
        endif
      endif

      !sunny
      if (l3leaves) then
        write(*,*)'trying to d 3 sunny leaves, not available'
        stop
      !do angle=1,nangle_gauss
      !  call leafeb_ags(PARleaf_sun(k_can,angle), LWin_leafsun(k_can,angle), leaf_eps, transpiretype, lwidth, llength, &
      !                  tmp0(i,j,k_can), humidairpa(k_can),qt0(i,j,k_can), windsp(k_can), presf(k_can),           & ! in
      !                  rhof(k_can), svm(i,j,k_can,indCO2), phitot(i,j), & 
      !                  gcc_leafsun_old(i,j,k_can,angle),ci_leafsun_old(i,j,k_can,angle),rk3coef, & ! additonal variables needed (Xabi)
      !                  tleaf_sun(k_can,angle), gcc_leafsun(i,j,k_can,angle), rb_leafsun(k_can,angle),ci_leafshad(i,j,k_can,angle),& !
      !                  sh_leafsun(k_can,angle), le_leafsun(k_can,angle), LWout_leafsun(k_can,angle))
      !   LWnet_leafsun(k_can,angle) = LWin_leafsun(k_can,angle) - LWout_leafsun(k_can,angle)
      !  if (lrelaxgc_can) then
      !    if (gccan_old_set .and. rk3step ==3) then
      !      gcc_leafsun_old(i,j,k_can,angle) = gcc_leafshad(i,j,k_can,angle)
      !     else if (.not. gccan_old_set) then
      !      gcc_leafsun_old(i,j,k_can,angle) = gcc_leafshad(i,j,k_can,angle)
      !    endif
      !  endif
      !  if (lrelaxci_can) then
      !    if (cican_old_set .and. rk3step ==3) then
      !      ci_leafsun_old(i,j,k_can,angle) = ci_leafshad(i,j,k_can,angle)
      !     else if (.not. cican_old_set) then
      !      ci_leafsun_old(i,j,k_can,angle) = ci_leafshad(i,j,k_can,angle)
      !    endif
      !  endif
      !
      !end do
      !An_leafsunaver      = sum(weight_g * An_leafsun(i,j,k_can,1:nangle_gauss))
      !gcc_leafsunaver     = sum(weight_g * gcc_leafsun(i,j,k_can,1:(nangle_gauss)))
      !sh_leafsunaver      = sum(weight_g * sh_leafsun(i,j,k_can,1:nangle_gauss))
      !le_leafsunaver      = sum(weight_g * le_leafsun(i,j,k_can,1:(nangle_gauss)))
      else ! angle-dependent PAR are averaged following 3 point gaussian procedure
        cPARleaf_allsun   = sum(weight_g *  cPARleaf_sun(k_can,1:nangle_gauss))
        call leafeb_ags(cPARleaf_allsun(k_can), LWin_leafsun(k_can), leaf_eps, transpiretype, lwidth, llength, & ! in
                        tmp0(i,j,k_can), humidairpa(k_can),qt0(i,j,k_can), windsp(k_can), presf(k_can),        & ! in
                        rhof(k_can), svm(i,j,k_can,indCO2), phitot(i,j),                                       & ! in
                        gcc_leafsun_old(i,j,k_can),ci_leafsun_old(i,j,k_can),rk3coef,                          & ! in 
                        tleaf_sun(k_can), gcc_leafsun(i,j,k_can), rb_leafsun(k_can),ci_leafsun(i,j,k_can),     & ! out
                        sh_leafsun(k_can), le_leafsun(k_can), LWout_leafsun(k_can),An_leafsun(k_can))            ! out
        LWnet_leafsun(k_can) = LWin_leafsun(k_can) - LWout_leafsun(k_can)
        if (lrelaxgc_can) then
          if (gccan_old_set .and. rk3step ==3) then
              gcc_leafsun_old(i,j,k_can) = gcc_leafsun(i,j,k_can)
          else if (.not. gccan_old_set) then
            gcc_leafsun_old(i,j,k_can) = gcc_leafsun(i,j,k_can)
          endif
        endif
        if (lrelaxci_can) then
          if (cican_old_set .and. rk3step ==3) then
              ci_leafsun_old(i,j,k_can) = ci_leafsun(i,j,k_can)
          else if (.not. cican_old_set) then
            ci_leafsun_old(i,j,k_can) = ci_leafsun(i,j,k_can)
          endif
        endif
      endif !l3leaves

    end do !k_can
   
   !                                                                                                ! 
   !######### STEP 3 - determining source terms for LES ####################
   !                                                                                                !
             ! --- calculate sources of heat, water by scaling
             !        from the leaf to the grid volume
    !we clauclate pad:
    !PAD = LAI_can / dz_can  ! 

    call canopysource(sh_leafsun, sh_leafshad,           &
                      le_leafsun, le_leafshad,           &
                      An_leafsun, An_leafshad,           &
                      cfSL,                         &
                      sh_can(i,j,1:ncanopy),le_can(i,j,1:ncanopy),Fco2_can(i,j,1:ncanopy))                                ! recv

    !convert sources to tendencies
    S_theta(i,j,1:ncanopy)   = sh_can(i,j,1:ncanopy)/(rhof(1:ncanopy)*cp)
    S_qt(i,j,1:ncanopy)      = le_can(i,j,1:ncanopy)/(rhof(1:ncanopy)*rlv)
    S_co2(i,j,1:ncanopy)     = Fco2_can(i,j,1:ncanopy)*(MW_Air/MW_CO2) * (1.0/rhof(1:ncanopy))* 1000 !In  ppb/s
    
   !                                                                                                ! 
   !######### STEP 4 - pass on tskin_can, effective canopy albedo and in-canopy radiation profiles####################
   !   
    
    ! tskin_can as weighted average similar to sources in canopysource
    tskin_can(i,j) = tleaf_sun(ncanopy)*cfSL(ncanopy) + tleaf_shad(ncanopy)*(1.-cfSL(ncanopy))
    tskin_rad(i,j)  = tskin_can(i,j)
    albedo_rad(i,j) = albsw_can(i,j)
                      !albdir_can(i,j),albdif_can(i,j)
    
    if (iradiation==irad_par) then ! all terms must be positive
      swdir(i,j,:ncanopy) = -swdir_can(:)                 
      swdif(i,j,:ncanopy) = -swdif_can(:)               
      swd  (i,j,:ncanopy) = -(swdir_can(:) + swdif_can(:))
      swu  (i,j,:ncanopy) =  swu_can(:)               
      lwd  (i,j,:ncanopy) = -lwd_can(:)
      lwu  (i,j,:ncanopy) =  lwu_can(:)
    else
      swdir(i,j,:ncanopy) = swdir_can(:)                 
      swdif(i,j,:ncanopy) = swdif_can(:)               
      swd  (i,j,:ncanopy) = (swdir_can(:) + swdif_can(:))
      swu  (i,j,:ncanopy) = swu_can(:)               
      lwd  (i,j,:ncanopy) = lwd_can(:)
      lwu  (i,j,:ncanopy) = lwu_can(:)
    endif
    
   !                                                                                                ! 
   !######### STEP 5 - check switches about lags in plant response ####################
   !                                                                                                !


   if (lrelaxgc_can .and. (.not. gccan_old_set) ) then
     if (rk3step == 3) then
       gccan_old_set = .true.
     endif
   endif
   if (lrelaxci_can .and. (.not. cican_old_set) ) then
     if (rk3step == 3) then
       cican_old_set = .true.
     endif
   endif  

  end subroutine canopyeb
  subroutine canopyu (putout)
    use modglobal, only  : i1, ih, j1, j2, jh, k1, cu, cv, dzh, imax, jmax
    use modfields, only  : u0, v0, w0
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ucor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: vcor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ftau  (imax,jmax)
    integer             :: k, kp

    ucor = u0 + cu
    vcor = v0 + cv

    do k=1,ncanopy
      kp   = k+1
      ftau = cd * padf(k) * sqrt(ucor(2:i1,2:j1,k)**2 +  &
                ((vcor(1:imax,2:j1,k)+vcor(2:i1,2:j1,k)+vcor(1:imax,3:j2,k)+vcor(2:i1,3:j2,k))/4)**2 + &
                ((dzh(kp)*(w0(1:imax,2:j1,k)+w0(2:i1,2:j1,k))+dzh(k)*(w0(1:imax,2:j1,kp)+w0(2:i1,2:j1,kp)))/(2*(dzh(k)+dzh(kp))))**2 &
                )

      putout(2:i1,2:j1,k) = putout(2:i1,2:j1,k) - ftau * ucor(2:i1,2:j1,k)
    end do

    return
  end subroutine canopyu

  subroutine canopyv (putout)
    use modglobal, only  : i1, i2, ih, j1, jh, k1, cu, cv, dzh, imax, jmax
    use modfields, only  : u0, v0, w0
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ucor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: vcor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ftau  (imax,jmax)
    integer             :: k, kp

    ucor = u0 + cu
    vcor = v0 + cv

    do k=1,ncanopy
      kp   = k+1
      ftau = cd * padf(k) * sqrt(vcor(2:i1,2:j1,k)**2 +  &
                ((ucor(3:i2,2:j1,k)+ucor(2:i1,2:j1,k)+ucor(3:i2,1:jmax,k)+ucor(2:i1,1:jmax,k))/4)**2 + &
                ((dzh(kp)*(w0(2:i1,1:jmax,k)+w0(2:i1,2:j1,k))+dzh(k)*(w0(2:i1,1:jmax,kp)+w0(2:i1,2:j1,kp)))/(2*(dzh(k)+dzh(kp))))**2 &
                )

      putout(2:i1,2:j1,k) = putout(2:i1,2:j1,k) - ftau * vcor(2:i1,2:j1,k)
    end do

    return
  end subroutine canopyv

  subroutine canopyw (putout)
    use modglobal, only  : i1, i2, ih, j1, j2, jh, k1, cu, cv, dzh, dzf, imax, jmax
    use modfields, only  : u0, v0, w0
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ucor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: vcor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ftau  (imax,jmax)
    integer             :: k, km

    ucor = u0 + cu
    vcor = v0 + cv

    do k=2,(ncanopy+1)
      km   = k-1
      ftau = cd * padh(k) * sqrt(w0(2:i1,2:j1,k)**2 +  &
                ((dzf(km)*(ucor(2:i1,2:j1,k)+ucor(3:i2,2:j1,k))+dzf(k)*(ucor(2:i1,2:j1,km)+ucor(3:i2,2:j1,km)))/(4*dzh(k)))**2 + &
                ((dzf(km)*(vcor(2:i1,2:j1,k)+vcor(2:i1,3:j2,k))+dzf(k)*(vcor(2:i1,2:j1,km)+vcor(2:i1,3:j2,km)))/(4*dzh(k)))**2 &
                )

      putout(2:i1,2:j1,k) = putout(2:i1,2:j1,k) - ftau * w0(2:i1,2:j1,k)
    end do

    return
  end subroutine canopyw

  subroutine canopye (putout)
    use modglobal, only  : i1, i2, ih, j1, j2, jh, k1, cu, cv, dzh, imax, jmax
    use modfields, only  : u0, v0, w0, e120
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ucor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: vcor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ftau  (imax,jmax)
    integer             :: k, kp

    ucor = u0 + cu
    vcor = v0 + cv

    do k=1,ncanopy
      kp   = k+1
      ftau = cd * padf(k) * sqrt(((ucor(3:i2,2:j1,k)+ucor(2:i1,2:j1,k))/2)**2 + &
                                 ((vcor(2:i1,3:j2,k)+vcor(2:i1,2:j1,k))/2)**2 + &
                ((dzh(kp)*w0(2:i1,2:j1,k)+dzh(k)*w0(2:i1,2:j1,kp))/(dzh(k)+dzh(kp)))**2 &
                )

      putout(2:i1,2:j1,k) = putout(2:i1,2:j1,k) - e120(2:i1,2:j1,k) * ftau
    end do

    return
  end subroutine canopye

  subroutine canopyc (putout, flux_top, flux_surf, alpha, pai)
    use modglobal, only  : i1, i2, ih, j1, j2, jh, k1, dzf, imax, jmax
    use modfields, only  : rhobh, rhobf
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real, intent(in   ) :: flux_top
    real, intent(in   ) :: flux_surf(i2,j2)
    real, intent(in   ) :: alpha
    real, intent(in   ) :: pai(ncanopy+1)
    real                :: flux_net (i2,j2)
    integer             :: k
    real                :: integratedcontribution(imax,jmax,ncanopy+1), tendency(imax,jmax,ncanopy)

    flux_net                        = flux_top * rhobh(ncanopy+1) - flux_surf * rhobh(1)
    integratedcontribution(:,:,1)   = 0.0

    do k=2,(ncanopy+1)
      integratedcontribution(:,:,k) = flux_net(2:i1,2:j1) * exp(- alpha * pai(k))
    end do
    do k=1,ncanopy
      tendency(:,:,k) = ( integratedcontribution(:,:,(k+1)) - integratedcontribution(:,:,k) ) / ( rhobf(k) * dzf(k) )
    end do

    putout(2:i1,2:j1,1:ncanopy) = putout(2:i1,2:j1,1:ncanopy) + tendency

    return
  end subroutine canopyc
!subroutine canopysource(layers,pad,             & ! in
subroutine canopysource(     sunleafsh, shadeleafsh,             &
                             sunleaflh, shadeleaflh,             &
                             sunleafAn, shadeleafAn,             &
                             sunfrac,                            &
                             sh,lh,Fco2)                                ! out
! ======================================================================
! Original subroutine from Ned Patton ~ adapted and extended where needed
! ======================================================================

     ! ---- subroutine to calculate total sources in this grid volume by
     !        accounting for the percentage of sunlit vs. shaded leaves
     !        and scaling sources from the leaf to the grid volume

     implicit none

     ! --- incoming variables

 !    integer, intent(in)                             :: layers, &
 !                                                       nclos,  &
 !                                                       nrows

     !integer, intent(in)                             :: layers

     !real, intent(in), dimension(layers)             :: pad            ! plant area density [m2 m-3]

 !   real, intent(in), dimension(nclos,nrows,layers) :: sunleafsh,   & ! sensible heat flux (sunlit) [W m-2]
     real, intent(in), dimension(ncanopy)             :: sunleafsh,   & ! sensible heat flux (sunlit) [W m-2]
                                                        shadeleafsh, & ! sensible heat flux (shaded) [W m-2]
                                                        sunleaflh,   & ! latent heat flux (sunlit) [W m-2]
                                                        shadeleaflh, & ! latent heat flux (shaded) [W m-2]
                                                        sunleafAn,   & ! CO2 flux (sunlit) [mg m-2 s-1]
                                                        shadeleafAn, & ! CO2 flux (shaded) [mg m-2 s-1]
                                                        sunfrac      ! fraction of leaves in sun [-]

     ! --- outgoing variables

     !real, intent(out),dimension(nclos,nrows,layers) :: sh,lh          ! heat and moisture source [W m-3]
     real, intent(out),dimension(ncanopy) :: sh,lh        ! heat and moisture source [W m-3]
     real, intent(out),dimension(ncanopy) :: Fco2         ! CO2 source [mg CO2 m-3 s-1]

     ! --- local variables

     !integer :: i,j,ii
     integer :: ii

     ! --- begin calculation

 !    do j = 1,nrows
 !    do i = 1,nclos

         do ii = 1, ncanopy
            ! ---- sunleafsh,sunleaflh,shadeleafsh, and shadeleaflh are in [W m-2].
            !        multiplying by pad [m2 m-3], which represents the frontal
            !        leaf area occupying this grid volume (note that the sensible
            !        heat values have already been multiplied by 2 to account for
            !        two-sided leaves) scales them to the grid volume and
            !        results in a source of heat and moisture in units of [W m-3].

         !  sh(i,j,ii)  = (sunleafsh(i,j,ii)*sunfrac(i,j,ii)  &
         !               + shadeleafsh(i,j,ii)*(1.-sunfrac(i,j,ii))) * pad(ii)

         !  lh(i,j,ii)  = (sunleaflh(i,j,ii)*sunfrac(i,j,ii)  &
         !               + shadeleaflh(i,j,ii)*(1.-sunfrac(i,j,ii))) * pad(ii)
            sh(ii)  = (sunleafsh(ii)*sunfrac(ii) + shadeleafsh(ii)*(1.-sunfrac(ii))) * padf(ii)

            lh(ii)  = (sunleaflh(ii)*sunfrac(ii) + shadeleaflh(ii)*(1.-sunfrac(ii))) * padf(ii)

            Fco2(ii) = (sunleafAn(ii)*sunfrac(ii) + shadeleafAn(ii)*(1.-sunfrac(ii))) * padf(ii)
         enddo

 !    enddo
 !    enddo

      return
 end subroutine canopysource

subroutine leafeb_ags(i_s, i_r, eps, transpiretype, lwidth, llength, & ! incoming
                          tairk, humidairpa,qtair, ws, pres,         & ! in     
                          !, cantype, vegtyp,                    &
                          rho, CO2air, phi_tot,gcc_old,ci_old,rk3coef,      &  ! additonal in variables needed (XPB)
                          tleaf, gccleaf,rb,ci,                         & ! in/out
                          sh, le, LWout,                             & ! out
                          An)                                        ! out extra by XPB
! ======================================================================
! Original subroutine from Ned Patton ~ adapted where needed
! ======================================================================

      ! ---- this routine calculates the leaf energy balance at an individual vertical canopy location
      !      largely following Leuning et al. (1995), Plant Cell Environment, 18, 1183-1200
      !      original code base comes from Alex Guenther and Michael Boy (MEGAN), but
      !      heavily modified by EGP
      !
      !      Key assumption is that we're using the isothermal form of the Penman-Monteith
      !      combination equation (Jones, 1976) to partition the net radiation absorbed by
      !      the leaves into sensible and latent heat.
      !
      !      Note: 'transpiretype' defines whether the leaves are transpiring from all/both sides
      !       of the leaves or only from one side or if there is some transpiration through the cuticle too.
      !
      !      substantially modified by EGP  4/2013
      !         now uses Paw U, 1987, J. Therm. Biol. formulation for iterating

      !modified by Xabi Pedruzo March 2020 to work with Ags and another
      !canopy radiative scheme (Pedruzo-Bagazgoitia et al. 2017)
       use modsurfdata, only : f_Ags,nuco2q
       implicit none

       ! ---- incoming variables

       real, intent(in)    :: & !ppfd,             &      ! absorbed PAR photon flux density [umol m-2 s-1]
                              i_s,              &      ! absorbed incoming solar radiation [W m-2]
                              i_r,              &      ! absorbed ifrared radiation [W m-2]
                              eps,              &      ! leaf IR emissivity
                              transpiretype,    &      ! what type of transpirer?
                              lwidth,           &      ! leaf width/shoot diameter [m]
                              llength,          &      ! leaf/shoot length [m]
                              tairk,            &      ! air temperature [K]
                              humidairpa,       &      ! atmospheric vapor pressure [Pa]
                              qtair,            &      ! atmospheric water mixing ratio [kg kg-1]
                              ws,               &      ! wind speed [m s-1]
                              pres,             &      ! surface pressure [Pa]
                              rho,              &      ! air density [kg m-3]
                              CO2air,           &      ! atmospheric CO2 concentration [ppb]
                              phi_tot,          &      ! Total soil water content [-]
                              gcc_old,          &      ! gc in previous timestep
                              ci_old,           &      ! ci in previous timestep
                              rk3coef                  ! runge-kutta timestep
                              !dt,               &      ! time-step duration [s]
                              !psir,             &      ! root-zone soil water potential [m]
                              !xw2,              &      ! soil moist. parameter controling slope in Ball-Berry model
                              !xwc4                     ! soil moist. control on A_net

       !integer, intent(in) :: cantype                  ! what type of canopy? (multi-layer characteristics)
       !integer, intent(in) :: vegtyp                   ! what type of canopy? (USGS type)

       ! ---- in/out variables

       real, intent(inout) :: tleaf                    ! leaf temperature from last timestep [K]
                                                       !   updated in this routine

       !real, intent(inout) :: rs                       ! stomatal resistance from last timestep [s m-1]
                                                       !   updated in this routine
       real, intent(inout) :: gccleaf                  ! carbon stomatal conductance from last timestep [s m-1]
                                                       !   updated in this routine
       real, intent(inout) :: ci                       ! leaf internal carbon concentration 
                                                       !   updated in this routine
       real, intent(inout) :: rb                       ! boundary layer resistance from last timestep [s m-1]
                                                       !   updated in this routine

       ! ---- outgoing variables

       real, intent(out)   :: sh,               &      ! sensible heat flux [W m-2]
                              le,               &      ! latent heat flux [W m-2]
                              An,               &      ! net CO2 exchange [mg m-2 s-1]  ! XPB
                              LWout                    ! LW from leaf at this temp [W m-2]

       ! ---- local variables

       !integer, parameter :: wet = 0   ! EGP fix for rainy conditions
       integer :: ii
       real    :: gb
       real    :: wind
       real    :: tdelt
       real    :: tleaf_l, tleaf_r
       real    :: rr_l, rr_r
       !DOUBLE PRECISION :: t_n, t_g, t_p  ! if we want a better iteration
       !DOUBLE PRECISION :: f_n, f_g, f_p
       real :: t_n, t_g, t_p
       real :: f_n, f_g, f_p

       real    :: Fleaf
       real    :: fstr,Am,Rdark,alphac,co2abs,CO2comp,Ds,D0,fmin

       real :: humidairkgm3
       integer, parameter :: num_iterations = 50
       real,    parameter :: wind_mn = 1.e-6     ! minimum wind speed [m s-1]
       real,    parameter :: tiny = 1.e-9

       ! --- check for negative windspeed and force minimum windspeed (wind_mn)

       wind = max(abs(ws), wind_mn)

       ! --- calculate current air density [kg m-3]
       !     where, rho(air) = rho(dry_air) + rho(vapor)

       !rho = (pres / (r_d * tairk)) + (humidairpa / (r_v * tairk))

       ! --- make sure gs_old has last timestep's value (mainly for restart)
       !gs_old = 1. /  rs
       !ci_old  = ci

       ! --- convert air vapor pressure from  [Pa] to [kg m-3]

       humidairkgm3 = converthumiditypa2kgm3(humidairpa, tairk)

       ! ---- iterate using a bisection method where we assume that
       !      the function is linear between the current two guesses
       !      'a' and 'b', and choose the next guess to be the value
       !      where that linear function crosses zero, i.e.:
       !
       !       f(x) - f(a)       f(b) - f(a)
       !      -------------  =  -------------
       !          x - a             b - a
       !
       !     solving for the new guess 'x' when f(x) = 0 yields:
       !
       !                        b - a
       !     x = a - f(a) * -------------
       !                     f(b) - f(a)
       !
       !     check convergence against 'tiny'
       !
       !     this iteration method is something like a Newton-Raphson
       !     method, but where we don't presume to know an analytic form
       !     for the first derivative.  Rather, we presume a linear slope

       ! --- since we know that the leaves will be close to the air
       !     temp, bracket the air temp to bound the root.
       tleaf_l = tairk - 10.
       tleaf_r = tairk + 10.

       ! --- evaluate function at the left e``nd point 'a', where
       !     'a' = 'tleaf_l' and 'f(a)' = 'rr_l'

       tdelt = tleaf_l - tairk            ! --- current delta_T = T_leaf - T_air [K]
      !call gem (tleaf_l,tairk,humidairkgm3,wind,ppfd,pres,dt, &   ! send
      !          vegtyp,cantype,lwidth,llength,                &
      !          psir,xw2,xwc4,                                &
      !          gs_old,                                       &   ! send/recv
      !          gb,rs)                                            ! recv
       !xabi edit
       call f_Ags(CO2air,qtair,rho,tairk,pres,tleaf_l,                & ! in
                  phi_tot,i_s,                                        & ! in
                  lrelaxgc_can,gccan_old_set,kgc_can,gcc_old,rk3coef, & ! in
                  lrelaxci_can,cican_old_set,kci_can,ci_old,          & ! in
                  gccleaf,Fleaf,ci,                                   & ! out
                  fstr,Am,Rdark,alphac,co2abs,CO2comp,Ds,D0,fmin)       ! out

       !get the water  stomatal resistancer from stomatal CO2 conductance:
       !rs = 1.0/(nuco2q*gccleaf)
       !get leaf boundary layer
       gb = leafblc(tleaf_l,tairk,wind,lwidth,llength)
       !end xabi edit

       sh      = leafh(tdelt,gb,rho)                                   ! sensible heat flux at this tleaf [W m-2]
       le      = leafle(tleaf_l,humidairkgm3,gb,gccleaf,transpiretype)  ! latent heat flux at this tleaf  [W m-2]
       LWout   = leafLWout(tleaf_l, eps)                               ! outgoing long wave out at this tleaf [W m-2]

       rr_l = i_s + i_r - LWout - sh - le                              ! leaf energy balance at this tleaf [W m-2]

       ! --- evaluate function at the end point 'b', where
       !     'b' = 'tleaf_r' and 'f(b)' = 'rr_r'

       tdelt = tleaf_r - tairk            ! --- current delta_T = T_leaf - T_air [K]
       !call gem (tleaf_r,tairk,humidairkgm3,wind,ppfd,pres,dt, &   ! send
       !          vegtyp,cantype,lwidth,llength,                &
       !          psir,xw2,xwc4,                                &
       !          gs_old,                                       &   ! send/recv
       !          gb,rs)                                            ! recv
       call f_Ags(CO2air,qtair,rho,tairk,pres,tleaf_r,                & ! in
                  phi_tot,i_s,                                        & ! in
                  lrelaxgc_can,gccan_old_set,kgc_can,gcc_old,rk3coef, & ! in
                  lrelaxci_can,cican_old_set,kci_can,ci_old,          & ! in
                  gccleaf,Fleaf,ci,                                   & ! out
                  fstr,Am,Rdark,alphac,co2abs,CO2comp,Ds,D0,fmin)       ! out

       !get the stomatal resistance to water from stomatal CO2 conductance:
       !rs = 1.0/(nuco2q*gccleaf)
       !get leaf boundary layer
       gb = leafblc(tleaf_r,tairk,wind,lwidth,llength)
       
       sh      = leafh(tdelt,gb,rho)                                   ! sensible heat flux at this tleaf [W m-2]
       le      = leafle(tleaf_r,humidairkgm3,gb,gccleaf,transpiretype)  ! latent heat flux at this tleaf  [W m-2]
       LWout   = leafLWout(tleaf_r, eps)                               ! outgoing long wave out at this tleaf [W m-2]

       rr_r = i_s + i_r - LWout - sh - le                              ! leaf energy balance at this tleaf [W m-2]

       ! ---- to simplify the iterations (and eliminate an if-statement)
       !      decide which of the bracketing end points provides
       !      the positive or negative side of the root
       !
       !      t_n = guess at negative side
       !      f_n = the function evaluated at negative side
       !      t_p = guess at positive side
       !      f_p = the function evaluated at positive side
       !
       !      by doing this we basically flip the slope of the function
       !      to always be one direction (i.e. positive or negative)
       if (rr_r > 0.) then
          t_n = tleaf_l
          f_n = rr_l
          t_p = tleaf_r
          f_p = rr_r
       else
          t_n = tleaf_r
          f_n = rr_r
          t_p = tleaf_l
          f_p = rr_l
       endif

       ! ---- now iterate

       do ii = 1,num_iterations

          ! --- find a guess for tleaf (t_g) which is the zero-crossing
          !     of a linear line between the current f_n and f_p

          t_g = t_n - (f_n * (t_p - t_n) / (f_p - f_n))
          ! --- current delta_T = T_leaf - T_air [K]

          tdelt = t_g - tairk

          !call gem (t_g,tairk,humidairkgm3,wind,ppfd,pres,dt, &   ! send
          !          vegtyp,cantype,lwidth,llength,            &
          !          psir,xw2,xwc4,                            &
          !          gs_old,                                   &   ! send/recv
          !          gb,rs)                                        ! recv
          call f_Ags(CO2air,qtair,rho,tairk,pres,t_g,                    & ! in
                     phi_tot,i_s,                                        & ! in
                     lrelaxgc_can,gccan_old_set,kgc_can,gcc_old,rk3coef, & ! in
                     lrelaxci_can,cican_old_set,kci_can,ci_old,          & ! in
                     gccleaf,Fleaf,ci,                                   & ! out
                     fstr,Am,Rdark,alphac,co2abs,CO2comp,Ds,D0,fmin)       ! out

          !get the stomatal resistance to water from stomatal CO2 conductance:
          !rs = 1.0/(nuco2q*gccleaf)
          !get leaf boundary layer
          gb = leafblc(t_g,tairk,wind,lwidth,llength)

          sh      = leafh(tdelt,gb,rho)                               ! sensible heat flux at this t_g [W m-2]
          le      = leafle(t_g,humidairkgm3,gb,gccleaf,transpiretype)  ! latent heat flux at this t_g  [W m-2]
          LWout   = leafLWout(t_g, eps)                               ! outgoing long wave out at this t_g [W m-2]

          f_g = i_s+ i_r - LWout - sh - le                           ! leaf energy balance at this t_g [W m-2]


          ! --- check for convergence
          !if (abs(f_g) < tiny .or. f_g == 0.) goto 100
          if (abs(f_g) < 0.001 .or. f_g == 0.) goto 100   !less strict  criteria proposed XPB 
         ! if (abs(f_g) < 0.001 .or. f_g == 0.) then   
         !  write(*,*)'convergenge reached!!f_g',f_g
         !  write(*,*)'after iteration',ii
         !  goto 100
         !endif

          ! --- if not converged, replace one of the old end points
          !     with the current guess.  if the current guess results
          !     in a value of the function that is > 0, replace
          !     the positive end point with the current guess.  if
          !     the current guess results in a value of the function
          !     that is < 0, replace the negative end point with
          !     the current guess.

          if (f_g > 0.) then
             t_p = t_g
             f_p = f_g
          else
             t_n = t_g
             f_n = f_g
          endif
       enddo

100    continue

       !An    = Fleaf
       An    = leafco2f(co2abs,ci, gb, gccleaf, transpiretype)
       tleaf = t_g    ! return final guess t_g as tleaf [K]
       rb    = 1./gb  ! boundary layer heat resistance [s m-1]
       return
 end subroutine leafeb_ags
  
real function leafh(tdelt,gb,rho)
! ======================================================================
! subroutine from Ned Patton ~ adapted where obvious
! ======================================================================

   ! ----  sensible heat flux term in leaf energy balance [W m-2]
   !        heat flux from both sides of leaf (i.e. times 2)

   ! 2 sides * bl conductance * temperature gradient * rho * cp
   
   use modglobal, only: cp
   implicit none
   real, intent(in) :: tdelt, gb, rho
   leafh = 2. * gb * tdelt * rho * cp

   return
end function leafh
real function leafle(tleaf, ambvap, gb, gcc, transpiretype)
! ======================================================================
! subroutine from Ned Patton ~ adapted where obvious
! ======================================================================

   ! ---- calculate the latent energy term in leaf energy balance
   use modsurfdata,only : nuco2q
   implicit none

   real, intent(in) :: tleaf,         &    ! leaf temperature [K]
                       ambvap,        &    ! atmospheric vapor density [kg m-3]
                       gb,            &    ! leaf boundary layer heat conductance [m s-1]
                       gcc,           &    ! stomatal carbon conductance [s m-1]
                       transpiretype       ! type of transpirer (1=hypostomatous, 2=amphistomatous,
                                           !    1.25=hypostomatous but with some transpiration through cuticle)

   real        :: leafres, vapdeficit, le,rs
   real        :: lathv               ! latent heat of vaporization at this temp [J kg-1]
   
    rs         = 1.0/(nuco2q*gcc)          ! stomatal resistance to water
    leafres    = 1./(1.075*gb) + rs        ! [s m-1] ! 1.075 factor to account for different between heat and water conductance (Goudriaan and van laar 1994)
    vapdeficit = svdtk(tleaf) - ambvap     ! [kg m-3]
    lathv      = lhv(tleaf)                ! [J kg-1]


    ! latent heat of vap [J kg-1] * vap deficit [kg m-3] / leaf resistance [s m-1]
    le = transpiretype * lathv * vapdeficit / leafres

    ! --- for now, don't allow dew on the leaves

    if (le < 0.) then
       leafle = 0.
    else
       leafle = le
    endif

    return
end function leafle
real function leafco2f(co2abs, ci, gb, gcc, transpiretype)
   ! ---- calculate the co2 flux including leaf boundary layer effects

   implicit none

   real, intent(in) :: co2abs,        &    ! atmospheric co2 concentration [mg C m-3]
                       ci,            &    ! leaf internal co2 concentration [mg C m-3]
                       gb,            &    ! leaf boundary layer heat conductance [m s-1]
                       gcc,           &    ! stomatal carbon conductance [s m-1]
                       transpiretype       ! type of transpirer (1=hypostomatous, 2=amphistomatous,
                                           !    1.25=hypostomatous but with some transpiration through cuticle)

   real          :: gbv ! boundary layer water vapor conductance [m s-1]
   real          :: co2f ! net assimilation rate [mg C m-2 s-1]
   gbv =  1.075 * gb ! Goudriaan and van laar 1994
   co2f = -(co2abs-ci) / ( (1./gcc) + 1.4* (1./gbv))! Eq 8.13,"Modelling potential growth processes",Goudriaan and van Laar 1994.
   !co2f = (co2abs-ci) / ( 1.6*rs + 1.4* (1./gbv))! Eq 8.13,"Modelling potential growth processes",Goudriaan and van Laar 1994.
   
   leafco2f = transpiretype * co2f ! if it transpires water vapor, it also absorbs co2
 
end function leafco2f

real function leafblc(tleaf,tairk,wind,lwidth,llength)
! ======================================================================
! subroutine from Ned Patton
! ======================================================================
 
    ! ---- routine to calculate leaf boundary layer conductance [m s-1]
    !
    !      follows: leuning et al, 1995, plant cell environment, 18, 1183-1200
    !      although, uses the maximum of gb_forced and gb_free
    !          like Nikolov et al, Ecol. Modeling, 1995
 
    implicit none
 
    real,    intent(in) :: tleaf,   &       ! leaf temperature [K] (leaf is assumed isothermal)
                           tairk,   &       ! air temperature [K]
                           wind,    &       ! local wind speed [m s-1]
                           lwidth,  &       ! leaf width / shoot diameter [m]
                           llength          ! leaf length / needle length [m]
 
    real :: gb_forced
    real :: gb_free
    real :: tdelta
 
    tdelta = abs(tleaf - tairk)
    gb_forced = 0.003 * sqrt( wind / lwidth )
 
    if (tdelta >= 0) then
       gb_free = 0.5 * 2.06e-5 * sqrt(sqrt(1.6e8 * tdelta * llength**3)) / llength
    else
       gb_free = 0.
    endif
 
 
  ! --- now decide which one to use (for GEM, take the largest)
 
  ! leafblc = gb_forced + gb_free          ! sum of the two per Leuning (1995)
    leafblc = max(gb_free, gb_forced)      ! maximum per Nikolov et al. (1995)
  ! leafblc = 0.5 * (gb_forced + gb_free)  ! take average of the two
 
    return
end function leafblc



! ======================================================================
! subroutine from Ned Patton ~ adapted where obvious
! ======================================================================
      subroutine spline(x,y,n,y2)
      implicit none

      integer n
      real :: yp1 = 2.0e30, ypn = 2.0e30  !< Values needed for the fit
      real x(n), y(n), y2(n)
      integer i, k
      real p, qn, sig, un
      real, allocatable :: u(:)

      allocate(u(max(n-1,1)))

      if(yp1 > .99e30) then
        y2(1) = 0.0
        u(1)  = 0.0
      else
        y2(1) = -0.5
        u(1)  = (3./(x(2) - x(1)))*((y(2) - y(1))/(x(2) - x(1)) - yp1)
      endif
      if(ypn > .99e+30) then
        qn = 0.0
        un = 0.0
      else
        qn = 0.5
        un = (3.0/(x(n) - x(n-1)))* (ypn - (y(n) - y(n-1))/(x(n) - x(n-1)))
      endif

      do i=2,n-1
        sig   = (x(i) - x(i-1))/(x(i+1) - x(i-1))
        p     = sig*y2(i-1) + 2.0
        y2(i) = (sig - 1.0)/p
        u(i)  = (6.0*((y(i+1) - y(i))/(x(i+1) - x(i)) - (y(i) - y(i-1)) &
              / (x(i) - x(i-1)))/(x(i+1) - x(i-1)) - sig*u(i-1))/p
      enddo

      y2(n) = (un - qn*u(n-1))/(qn*y2(n-1) + 1.0)

      do k=n-1,1,-1
        y2(k) = y2(k)*y2(k+1) + u(k)
      enddo

      deallocate(u)

      return
      end subroutine spline

! ======================================================================
! subroutine from Ned Patton ~ adapted where obvious
! ======================================================================

      subroutine splint(xa,ya,y2a,n,x,y)
      implicit none

      integer n
      real x,y, xa(n), y2a(n), ya(n)
      integer k,khi,klo
      real a,b,h
      klo = 1
      khi = n
    1 continue
      if(khi - klo > 1) then
         k = (khi + klo)/2
         if(xa(k) > x) then
           khi = k
         else
           klo = k
         endif
         go to 1
      endif
      h = xa(khi) - xa(klo)
      a = (xa(khi) - x)/h
      b = (x - xa(klo))/h
      y = a*ya(klo) + b*ya(khi) + ((a**3 - a)*y2a(klo) + (b**3 - b)*y2a(khi))*(h**2)/6.0

      return
      end subroutine splint
real function unexposedleafLWin(tk, eps)
! ======================================================================
! subroutine from Ned Patton ~ adapted where obvious
! ======================================================================

   ! ----  calculate ir into leaf that is shaded from the sky

   implicit none
   real, intent(in) :: eps, tk
   real, parameter  :: boltz = 5.67e-8

   unexposedleafLWin = eps * boltz * (tk**4)  !shouldnt we use air_eps?


 end function unexposedleafLWin

 ! ======================================================================

 real function exposedleafLWin(humidpa, tk)
! ======================================================================
! subroutine from Ned Patton ~ adapted where obvious
! ======================================================================

 ! ----  calculate ir into leaf that is exposed to the sky

   implicit none
   real, intent(in) :: tk, humidpa
   real :: emissatm
   real, parameter  :: boltz = 5.67e-8

   ! apparent atmospheric emissivity for clear skies: function of water
   ! vapor pressure [Pa] and ambient air temperature [K] based on
   ! Brutsaert(1975) referenced in Leuning (1995)

   emissatm        = 0.642 * (humidpa / tk)**(1./7.)
   exposedleafLWin = emissatm * boltz * (tk**4)

   return
 end function exposedleafLWin
real function leafLWout(tleaf, eps)
! ======================================================================
! subroutine from Ned Patton - adapted where obvious
! ======================================================================

   ! ---- ir thermal radiation energy emissiom by leaf
   !
   !      note: leaf is presumed to be two-sided (hence the times 2)
   
   use modglobal,only:  boltz
   implicit none
   real, intent(in)  :: tleaf,eps
   !real, parameter  :: boltz = 5.67e-8

   leafLWout = eps * boltz * (2. * (tleaf**4))

   return
end function leafLWout

 real function watervappres(spec_hum, pres)
! ======================================================================
! subroutine from Ned Patton ~ adapted where obvious
! ======================================================================

   ! ---- convert water specific humidity [kg kg-1] to water vapor pressure [Pa]

   implicit none
   real, intent(in) :: spec_hum       ! specific humidity [kg kg-1]
   real, intent(in) :: pres           ! atmospheric pressure [Pa]
   real, parameter :: waterairratio = 18.016/28.97 ! ratio of molecular weight of water to that of air [-]

   watervappres = (spec_hum * pres) / (waterairratio + (1.-waterairratio)*spec_hum)

   return
 end function watervappres
real function converthumiditypa2kgm3(pa, tk)
! ======================================================================
! subroutine from Ned Patton ~ adapted where obvious
! ======================================================================

   ! ---- convert vapor pressure [Pa] into vapor density  [kg m-3]
   
   use modglobal, only : rv
   implicit none
   real, intent(in) :: pa, tk
   !real:: r_v= 461.6

   converthumiditypa2kgm3 = pa / (rv * tk)

   return
end function converthumiditypa2kgm3
real function svdtk(tk)
! ======================================================================
! subroutine from Ned Patton ~ adapted where obvious
! ======================================================================

      ! ---- calculate saturation vapor density [kg m-3]

      implicit none
      real, intent(in) :: tk

      real :: svp!,e_sat,converthumiditypa2kgm3

      ! ---- determine saturation vapor pressure [Pa]

      svp = e_sat(tk)

      ! --- convert vapor pressure [Pa] to vapor density [kg m-3] (using: p = \rho r_v t)

      svdtk = converthumiditypa2kgm3(svp, tk)

      return
end function svdtk
real function lhv(tk)

  ! ---- latent heat of vaporization [J kg-1] from stull p641 (expects tk in [K])

  implicit none
  real, intent(in)  :: tk
  real, parameter  ::  kelvin  = 273.16


  lhv = (2.501 - (2.37e-3 * (tk - kelvin))) * 1.e6
return
end function lhv

real function e_sat(T)
  implicit none
  real, intent(in) :: T ! absolute temperature [K]
  !real             ::  e_sat  ! saturation vapor pressure [Pa]

  e_sat = 0.611e3 * exp(17.2694 * (T - 273.16) / (T - 35.86))
return
end function e_sat

end module modcanopy
