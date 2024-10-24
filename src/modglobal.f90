!> \file modglobal.f90
!!  Declares the global constants

!>
!!  Declares the global constants
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
module modglobal
use iso_c_binding
use modprecision
use netcdf
implicit none
save

      ! Simulation dimensions (parconst.f90)
      integer :: itot = 64
      integer :: jtot = 64
      integer :: imax
      integer :: jmax
      integer :: kmax = 96
      integer ::  i1
      integer ::  j1
      integer ::  k1
      integer ::  k2
      integer ::  i2
      integer ::  j2
      integer ::  nsv = 0       !< Number of additional scalar fields

      integer ::  ih=3
      integer ::  jh=3
      integer ::  kh=1
      integer ::  kcb=0

      character(256) :: fname_options = 'namoptions'
      logical :: lwarmstart = .false.!<   flag for "cold" or "warm" start
      real    :: trestart  = 3600. !<     * each trestart sec. a restart file is written to disk
      integer(kind=longint) :: itrestart !<     * each trestart sec. a restart file is written to disk
      integer(kind=longint)    :: tnextrestart    !<     * each trestart sec. a restart file is written to disk
      character(50) :: startfile    !<    * name of the restart file
      logical :: lstart_netcdf = .false.      !< Start from a NetCDF file

      logical :: llsadv   = .false. !<  switch for large scale forcings
      integer :: ntimedep = 100     !< maximum number of time points for time-dependent forcings

      !<  Global constants modconst.f90
      !< File numbers

      integer, parameter :: ifinput    = 1
      integer, parameter :: ifoutput   = 2
      integer, parameter :: ifnamopt   = 3

      real,parameter :: pi       = 3.141592653589793116
      real(field_r),parameter :: grav     = 9.81             !<    *gravity acceleration.
      real(field_r),parameter :: rd       = 287.04           !<    *gas constant for dry air.
      real(field_r),parameter :: rv       = 461.5            !<    *gas constant for water vapor.
      real(field_r),parameter :: cp       = 1004.            !<    *specific heat at constant pressure (dry air).
      real(field_r),parameter :: rlv     = 2.53e6           !<    *latent heat for vaporisation J/kg
      real(field_r),parameter :: riv     = 2.84e6           !<    *latent heat for sublimation
      real(field_r),parameter :: tup     = 268.             !<    * Temperature range over which mixed phase occurs (high)
      real(field_r),parameter :: tdn     = 253.             !<    * Temperature range over which mixed phase occurs (low)
      real,parameter :: ep       = rd/rv            !<    0.622
      real,parameter :: ep2      = rv/rd - 1.       !<    0.61
      !< real,parameter :: cv       = cp-rd            !<    716.96
      real,parameter :: rcp      = rd/cp            !<    0.286
      real,parameter :: cpr      = cp/rd            !<    3.50
      real,parameter :: rlvocp   = rlv/cp           !<    2.49
      real,parameter :: mair     = 28.967           !< Molar mass of air
      real,parameter :: rhow     = 0.998e3          !<    * Density of water
      real,parameter :: pref0    = 1.e5             !<    *standard pressure used in exner function.
      real,parameter :: tmelt    = 273.16           !<    *temperature of melting of ice.
      real,parameter :: es0      = 610.78           !<    * constants used for computation
      real,parameter :: at       = 17.27            !<    * of saturation mixing ratio
      real,parameter :: bt       = 35.86            !<    * using Tetens Formula.
      real(field_r),parameter :: ekmin    = 1.e-6            !<    *minimum value for k-coefficient.
      real(field_r),parameter :: e12min   = 5.e-5            !<    *minimum value for TKE.
      real,parameter :: fkar     = 0.4              !<    *Von Karman constant
      real(field_r),parameter :: eps1     = 1.e-10           !<    *very small number*
      real,parameter :: epscloud = 1.e-5            !<    *limit for cloud calculation 0.01 g/kg
      real,parameter :: boltz    = 5.67e-8          !<    *Stefan-Boltzmann constant

      ! Land-surface
      real,parameter :: rho_solid_soil = 2700       !< Density of dry solid soil (kg m-3)
      real,parameter :: rho_C_matrix   = 1.6e6      !< Volumetric soil heat capacity [J m-3 K-1]
      real,parameter :: rho_C_water    = 4.18e6     !< Volumetric water heat capacity [J m-3 K-1]
      real,parameter :: gamma_T_matrix = 3.4293695508945325 !< Heat conductivity soil [J s-1 m-1 K-1]
      real,parameter :: gamma_T_water  = 0.57       !< Heat conductivity water [J s-1 m-1 K-1]

      logical :: lcoriol  = .true.  !<  switch for coriolis force
      logical :: lpressgrad = .true.  !<  switch for horizontal pressure gradient force
      integer       :: igrw_damp = 2 !< switch to enable gravity wave damping
      real(field_r) :: geodamptime = 7200. !< time scale for nudging to geowind in sponge layer, prevents oscillations
      real(field_r) :: uvdamprate = 0.  !< rate for damping mean horizontal wind
      real(field_r) :: om22                       !<    *2.*omega_earth*cos(lat)
      real(field_r) :: om23                       !<    *2.*omega_earth*sin(lat)
      real(field_r) :: om22_gs                       !<    *2.*omega_earth*cos(lat)
      real(field_r) :: om23_gs                       !<    *2.*omega_earth*sin(lat)
      real          :: xlat    = 52.              !<    *latitude  in degrees.
      real          :: xlon    = 0.               !<    *longitude in degrees.
      logical       :: lrigidlid = .false. !< switch to enable simulations with a rigid lid
      real(field_r) :: unudge = 1.0   !< Nudging factor if igrw_damp == -1 (nudging mean wind fields to geostrophic values provided by lscale.inp)

      !Base state
      integer :: ibas_prf = 3
      integer, parameter :: ibas_thv    = 1 !< Theta_v constant (Useful in dry cases)
      integer, parameter :: ibas_bou    = 2 !< Boussinesq-like
      integer, parameter :: ibas_st1    = 3 !< Standard atmosphere, with surface temperature correction
      integer, parameter :: ibas_st2    = 4 !< Standard atmosphere, surface temperature 15 Celsius
      integer, parameter :: ibas_usr    = 5 !< User specified

      !Advection scheme
      integer :: iadv_mom = 5, iadv_tke = -1, iadv_thl = -1,iadv_qt = -1,iadv_sv = -1
      integer, parameter :: iadv_null   = 0
      integer, parameter :: iadv_upw    = 1
      integer, parameter :: iadv_cd2    = 2
      integer, parameter :: iadv_5th    = 5
      integer, parameter :: iadv_cd6    = 6
      integer, parameter :: iadv_62     = 62
      integer, parameter :: iadv_52     = 52
      integer, parameter :: iadv_kappa    = 7
      integer, parameter :: iadv_hybrid   = 55
      integer, parameter :: iadv_hybrid_f = 555

      real :: lambda_crit=100. !< maximum value for the smoothness. This controls if WENO or

      ! Tabulated saturation relation
      real, dimension(1:2000) :: ttab
      real, dimension(1:2000) :: esatltab
      real, dimension(1:2000) :: esatitab
      real, dimension(1:2000) :: esatmtab
      real, dimension(-100:4000) :: mygamma251
      real, dimension(-100:4000) :: mygamma21

      !$acc declare copyin(esatltab, esatitab, esatmtab)

      logical :: lmoist   = .true.  !<   switch to calculate moisture fields
      logical :: lnoclouds = .false. !<   switch to enable/disable thl calculations
      logical :: lfast_thermo = .false. !<   switch to enable faster icethermo scheme
      logical :: lsgbucorr= .false.  !<   switch to enable subgrid buoyancy flux
      logical :: lconstexner = .false.  !<  switch to use the initial pressure profile in the exner function

      ! Poisson solver: modpois / modhypre
      integer :: solver_id = 0       ! Identifier for nummerical solver:    0    1   2     3       4
                                     !                                     FFT  SMG PFMG BiCGSTAB GMRES
      integer :: maxiter = 10000     ! Number of iterations                 .    X   X     X       X
      real(real64):: tolerance = 1E-8! Convergence threshold                .    X   X     X       X
      integer :: n_pre = 1           ! Number of pre and post relaxations   .    X   X     X       X
      integer :: n_post =1           ! Number of pre and post relaxations   .    X   X     X       X
      integer :: precond_id = 1      ! Preconditioner ID                    .    .  12   0189     0189
      integer :: maxiter_precond = 1 ! Number of iterations for precondition per iteration
      integer :: hypre_logging = 1   ! HYPRE logging and print level - set higher value for more messages
      type solver_type
         !integer*8 solver,precond
        type(c_ptr) solver, precond
        integer   solver_id, precond_id, maxiter, n_post, n_pre, maxiter_precond
        real      tolerance
      end type
      type(solver_type) :: psolver

      ! Global variables (modvar.f90)
      integer :: xyear  = 0     !<     * year, only for time units in netcdf
      real :: xday      = 1.    !<     * day number
      real :: xtime     = 0.    !<     * GMT time (in hours)
      real(field_r) :: cu        = 0.    !<     * translation velocity in x-direction
      real(field_r) :: cv        = 0.    !<     * translation velocity in y-direction
      real :: runtime   = 300.  !<     * simulation time in secs
      real :: dtmax     = 20.   !<     * maximum time integration interval
      integer(kind=longint) :: idtmax        !<     * maximum time integration interval
      real :: dtav_glob   = 60.
      real :: timeav_glob = 3600.
      real :: tres     = 0.001
      real :: thres     = 5.e-3 !<     * threshold value for inversion height calculations

      integer(kind=longint) :: dt                !<     * time integration interval
      real(field_r) :: rdt                !<     * time integration interval
      integer               :: dt_reason=0  !< indicates which dt limit was the lowest
      integer(kind=longint) :: timee             !<     * elapsed time since the "cold" start
      real :: rtimee             !<     * elapsed time since the "cold" start
      integer(kind=longint) :: btime             !<     * time of (re)start
      integer :: ntrun          !<     * number of timesteps since the start of the run
      integer(kind=longint) :: timeleft
      logical :: ladaptive   = .false.    !<    * adaptive timestepping on or off
      logical :: ltotruntime = .false. !<    * Whether the runtime is counted since the last cold start (if true) or the last warm start (if false, default)

      real    :: courant = -1
      real    :: peclet  = 0.15
      integer(kind=longint) :: dt_lim


      integer :: rk3step = 0

      integer :: iexpnr = 0     !<     * number of the experiment

      character(3) cexpnr

      logical :: loutdirs = .false.       !< if true, create output directories using myidy
      character(20) :: output_prefix = '' !< prefix for output files e.g. for an output directory

      ! Variables for modopenboundary.f90
      logical :: lopenbc = .false., lsynturb = .false.,linithetero = .false.
      type boundary_type
        integer :: nx1,nx2,nx1patch,nx2patch,nx1u,nx2u,nx1v,nx2v,nx1w,nx2w,nx1turb,nx2turb
        real(field_r), allocatable, dimension(:,:,:) :: u,v,w,thl,qt,e12, &
          & u2,v2,w2,uv,uw,vw,thl2,qt2,wthl,wqt,ci,svturb
        real(field_r), allocatable, dimension(:,:,:,:) :: sv
        real(field_r), allocatable, dimension(:,:) :: radcorr,uphase,uphasesingle, &
          radcorrsingle,uturb,vturb,wturb,thlturb,qtturb,e12turb!,randqt,randthl
        real(field_r), allocatable, dimension(:,:,:,:) :: eigvec
        character (len=:), allocatable :: name
      end type
      type(boundary_type), dimension(5) :: boundary
      logical, dimension(5) :: lboundary = .false.
      logical, dimension(5) :: lperiodic =  (/.true., .true., .true., .true., .false./)
      real(field_r) :: dxint=-1.,dyint=-1.,dzint=-1.,tauh=60.,taum=0.,tau=60.,lambda=-1,lambdas=-1.,lambdas_x=-1.,lambdas_y=-1.,lambdas_z=-1.,dxturb=-1.,dyturb=-1.
      integer :: nmodes=100,ntboundary=1,pbc = 3,iturb=0
      real,dimension(:),allocatable :: tboundary

      ! modphsgrd.f90

      real(field_r) :: dx              !<  grid spacing in x-direction
      real(field_r) :: dy              !<  grid spacing in y-direction
      real(field_r) :: dz              !<  grid spacing in z-direction
      real(field_r) :: dxi             !<  1/dx
      real(field_r) :: dyi             !<  1/dy
      real(field_r) :: dzi             !<  1/dz
      real(field_r) :: dxiq            !<  1/(dx*4)
      real(field_r) :: dyiq            !<  1/(dy*4)
      real(field_r) :: dziq            !<  1/(dz*4)
      real(field_r) :: dxi5            !<  1/(dx*2)
      real(field_r) :: dyi5            !<  1/(dy*2)
      real(field_r) :: dzi5            !<  1/(dz*2)
      real(field_r) :: dx2i            !<  (1/dx)**2
      real(field_r) :: dy2i            !<  (1/dy)**2


      real(field_r) :: ijtot
      real(field_r), allocatable :: dzf(:), dzfi(:) !<  thickness of full level, and inverse
      real(field_r), allocatable :: dzh(:), dzhi(:) !<  thickness of half level, and inverse
      real(field_r), allocatable :: zh(:)           !<  height of half level [m]
      real(field_r), allocatable :: zf(:)           !<  height of full level [m]
      real(field_r) :: xsize    = -1                !<  domain size in x-direction
      real(field_r) :: ysize    = -1                !<  domain size in y-direction
      real(field_r), allocatable :: delta(:)        !<  (dx*dy*dz)**(1/3)
      real(field_r), allocatable :: deltai(:)       !<  (dx*dy*dz)**(-1/3)  or dzf**-1 for anisotropic diffusion

      logical :: leq      = .true.  !<  switch for (non)-equidistant mode.
      logical :: lmomsubs = .false.  !<  switch to apply subsidence on the momentum or not
      logical :: ldrydep = .false.         !< On/Off switch dry deposition

      character(80) :: author='', version='DALES 4.4.2'
contains

!> Initialize global settings.
!!
!! Set courant number, calculate the grid sizes (both computational and physical), and set the coriolis parameter
  subroutine initglobal
    use modmpi, only : nprocx, nprocy, myid,comm3d, mpierr, D_MPI_BCAST
    implicit none

    integer :: advarr(4)
    real phi, colat, silat, omega, omega_gs
    real :: ilratio
    integer :: k, n, m, ierr
    integer :: ncid, height_id
    character(80) chmess

    ! phsgrid
    imax = itot/nprocx
    jmax = jtot/nprocy
    i1=imax+1
    j1=jmax+1
    k1=kmax+1
    k2=kmax+2
    i2=imax+2
    j2=jmax+2
    !set the number of ghost cells. NB: This switch has to run in order of required ghost cells
    advarr = (/iadv_mom,iadv_tke,iadv_thl,iadv_qt/)
    if     (any(advarr==iadv_cd6).or.iadv_sv==iadv_cd6) then
      ih = 3
      jh = 3
      kh = 1
    elseif (any(advarr==iadv_62).or.iadv_sv==iadv_62) then
      ih = 3
      jh = 3
      kh = 1
    elseif (any(advarr==iadv_5th).or.iadv_sv==iadv_5th) then
      ih = 3
      jh = 3
      kh = 1
    elseif (any(advarr==iadv_52).or.iadv_sv==iadv_52) then
      ih = 3
      jh = 3
      kh = 1
    elseif (any(advarr==iadv_hybrid).or.iadv_sv==iadv_hybrid) then
      ih = 3
      jh = 3
      kh = 1
    elseif (any(advarr==iadv_hybrid_f).or.iadv_sv==iadv_hybrid_f) then
      ih = 3
      jh = 3
      kh = 1
    elseif (any(advarr==iadv_kappa).or.iadv_sv==iadv_kappa) then
      ih = 2
      jh = 2
      kh = 1
    elseif (any(advarr==iadv_cd2).or.iadv_sv==iadv_cd2) then
      ih = 1
      jh = 1
      kh = 1
    end if

    ! Global constants

    ! esatltab(m) gives the saturation vapor pressure over water at T corresponding to m
    ! esatitab(m) is the same over ice
    ! esatmtab(m) is interpolated between the ice and liquid values with ilratio
    ! http://www.radiativetransfer.org/misc/atmlabdoc/atmlab/h2o/thermodynamics/e_eq_water_mk.html
    ! Murphy and Koop 2005 parameterization formula.
    do m=1,2000
      ttab(m)=150.+0.2*m
      esatltab(m)=exp(54.842763-6763.22/ttab(m)-4.21*log(ttab(m))+0.000367*ttab(m)+&
           tanh(0.0415*(ttab(m)-218.8))*(53.878-1331.22/ttab(m)-9.44523*log(ttab(m))+ 0.014025*ttab(m)))

      esatitab(m)=exp(9.550426-5723.265/ttab(m)+3.53068*log(ttab(m))-0.00728332*ttab(m))
      ilratio = max(0.,min(1.,(ttab(m)-tdn)/(tup-tdn)))
      esatmtab(m) = ilratio*esatltab(m) + (1-ilratio)*esatitab(m)
    end do

    mygamma251(-100)=0.
    mygamma21(-100)=0.
    do m=-99,4000
      mygamma251(m)=max(lacz_gamma(m/100._dp+2.5_dp)/lacz_gamma(m/100._dp+1._dp)*( ((m/100.+3)*(m/100.+2)*(m/100.+1))**(-1./2.) ),0.)
      mygamma21(m)=max(lacz_gamma(m/100._dp+2._dp)/lacz_gamma(m/100._dp+1._dp)*( ((m/100.+3)*(m/100.+2)*(m/100.+1))**(-1./3.) ),0.)
    end do

    ! Select advection scheme for scalars. If not set in the options file, the momentum scheme is used
    if (iadv_tke<0) iadv_tke = iadv_mom
    if (iadv_thl<0) iadv_thl = iadv_mom
    if (iadv_qt<0)  iadv_qt  = iadv_mom

    !CvH remove where
    !where (iadv_sv<0)  iadv_sv  = iadv_mom
    if (iadv_sv < 0) then
      iadv_sv = iadv_mom
      write (*,*) "Warning: iadv_sv not specified, defaulting to iadv_mom = ", iadv_mom
    end if

    phi    = xlat*pi/180.
    colat  = cos(phi)
    silat  = sin(phi)
    if (lcoriol) then
      omega = 7.292e-5
      omega_gs = 7.292e-5
    else
      omega = 0.
      omega_gs = 0.
    end if
    om22   = 2.*omega*colat
    om23   = 2.*omega*silat
    om22_gs   = 2.*omega_gs*colat
    om23_gs   = 2.*omega_gs*silat

    ! Variables
    write(cexpnr,'(i3.3)') iexpnr


    ! Create the physical grid variables
    allocate(dzf(k1), dzfi(k1))
    allocate(dzh(k1), dzhi(k1))
    allocate(zh(k1))
    allocate(zf(k1))
    allocate(delta(k1),deltai(k1))

    ijtot = real(itot*jtot)

    dx = xsize / float(itot)
    dy = ysize / float(jtot)

    ! MPI

    ! Note, that the loop for reading zf and calculating zh
    ! has been split so that reading is only done on PE 1

    if (lstart_netcdf) then
      ierr = nf90_open('init.'//cexpnr//'.nc', NF90_NOWRITE, ncid)
      if (ierr /= nf90_noerr) call abort
      ierr = nf90_inq_varid(ncid, 'zh', height_id)
      if (ierr /= nf90_noerr) call abort
      ierr = nf90_get_var(ncid, height_id, zf, start=(/ 1 /), count=(/ kmax /))
      if (ierr /= nf90_noerr) call abort
      ierr = nf90_close(ncid)
      if (ierr /= nf90_noerr) call abort
    else
      if(myid==0)then
        open (ifinput,file='prof.inp.'//cexpnr,status='old',iostat=ierr)
        if (ierr /= 0) then
           write(6,*) 'Cannot open the file ', 'prof.inp.'//cexpnr
           STOP
        end if
        read(ifinput,'(a72)') chmess
        read(ifinput,'(a72)') chmess

        do k=1,kmax
          read(ifinput,*) zf(k)
        end do
        close(ifinput)

      end if ! end if myid==0
    end if

  ! MPI broadcast kmax elements from zf

    call D_MPI_BCAST(zf, kmax, 0, comm3d, mpierr)

    zh(1) = 0.0
    do k=1,kmax
      zh(k+1) = zh(k) + 2.0*(zf(k)-zh(k))
    end do
    zf(k1)  = zf(kmax)+ 2.0*(zh(k1)-zf(kmax))


    do  k=1,kmax
      dzf(k) = zh(k+1) - zh(k)
    end do
    dzf(k1) = dzf(kmax)

    dzh(1) = 2*zf(1)
    do k=2,k1
      dzh(k) = zf(k) - zf(k-1)
    end do

    dzfi = 1.0 / dzf
    dzhi = 1.0 / dzh

    do k=1,k1

       delta(k) = (dx*dy*dzf(k))**(1./3.)
       deltai(k) = 1./delta(k)     !can be overruled in modsubgrid in case anisotropic diffusion is applied
    end do

  !--------------------------------------------------
  ! *** Check whether the grid is equidistant *****
  !--------------------------------------------------

    leq=.true.
    dz = dzf(1)
    do k=1,k1
      if (abs(dzf(k)-dz)/dz>eps1) then
        leq = .false.
      end if
    end do

  ! MPI

    if(myid==0)then
      if (.not.leq) then
        write(6,*) &
            'WARNING, You are working with a non-equidistant grid!!!!'
      end if
    end if ! end if myid==0

    dxi     = 1./dx
    dyi     = 1./dy
    dzi     = 1./dz

    dxiq    = 0.25*dxi
    dyiq    = 0.25*dyi
    dziq    = 0.25*dzi

    dx2i    = dxi*dxi
    dy2i    = dyi*dyi

    dxi5    = 0.5*dxi
    dyi5    = 0.5*dyi
    dzi5    = 0.5*dzi

    if(myid==0)then
      write (6,*) 'lev    dz     zf      zh       dzh    delta'
      do k=k1,1,-1
        write(6,'(i4,5f10.2)') k,dzf(k),zf(k),zh(k),dzh(k),delta(k)
      end do

      do k=1,k1
        if (dzf(k) <= 0 .or. dzh(k) <= 0) then
          stop "Zero or negative level spacing found in dzf or dzh."
        end if
      end do
    end if
!     tnextrestart = trestart/tres
!     timeleft=ceiling(runtime/tres)

    !$acc enter data copyin(dzf, dzh, dzfi,dzhi, zh, zf, delta, deltai, &
    !$acc&                  esatmtab, esatitab, esatltab, mygamma251, mygamma21)

  end subroutine initglobal
!> Clean up when leaving the run
  subroutine exitglobal
    !$acc exit data delete(dzf, dzh, zh, zf, delta, deltai, &
    !$acc&                 esatmtab, esatitab, esatltab, mygamma251, mygamma21)

    deallocate(dzf,dzh,dzfi,dzhi,zh,zf,delta,deltai)
  end subroutine exitglobal

FUNCTION LACZ_GAMMA(X) RESULT(fn_val)

! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00

!----------------------------------------------------------------------

! This routine calculates the GAMMA function for a real argument X.
!   Computation is based on an algorithm outlined in reference 1.
!   The program uses rational functions that approximate the GAMMA
!   function to at least 20 significant decimal digits.  Coefficients
!   for the approximation over the interval (1,2) are unpublished.
!   Those for the approximation for X .GE. 12 are from reference 2.
!   The accuracy achieved depends on the arithmetic system, the
!   compiler, the intrinsic functions, and proper selection of the
!   machine-dependent constants.

!*******************************************************************

! Explanation of machine-dependent constants.  Let

! beta   - radix for the floating-point representation
! maxexp - the smallest positive power of beta that overflows

! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.

! XBIG   - the largest argument for which GAMMA(X) is representable
!          in the machine, i.e., the solution to the equation
!                  GAMMA(XBIG) = beta**maxexp
! XINF   - the largest machine representable floating-point number;
!          approximately beta**maxexp
! EPS    - the smallest positive floating-point number such that
!          1.0+EPS .GT. 1.0
! XMININ - the smallest positive floating-point number such that
!          1/XMININ is machine representable

!     Approximate values for some important machines are:

!                            beta       maxexp        XBIG

! CRAY-1         (S.P.)        2         8191        966.961
! Cyber 180/855
!   under NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, etc.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, etc.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-Format   (D.P.)        2          127        34.844
! VAX G-Format   (D.P.)        2         1023        171.489

!                            XINF         EPS        XMININ

! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! Cyber 180/855
!   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308

!*******************************************************************

! Error returns

!  The program returns the value XINF for singularities or
!     when overflow would occur.  The computation is believed
!     to be free of underflow and overflow.


!  Intrinsic functions required are:

!     INT, DBLE, EXP, LOG, REAL, SIN


! References: "An Overview of Software Development for Special
!              Functions", W. J. Cody, Lecture Notes in Mathematics,
!              506, Numerical Analysis Dundee, 1975, G. A. Watson
!              (ed.), Springer Verlag, Berlin, 1976.

!              Computer Approximations, Hart, Et. Al., Wiley and
!              sons, New York, 1968.

!  Latest modification: March 12, 1992

!  Authors: W. J. Cody and L. Stoltz
!           Applied Mathematics Division
!           Argonne National Laboratory
!           Argonne, IL 60439

!----------------------------------------------------------------------

IMPLICIT NONE

REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

! Local variables
INTEGER    :: i, n
LOGICAL    :: parity
REAL (dp)  :: fact, sum, xden, xnum, y, y1, ysq, z
!----------------------------------------------------------------------
!  Mathematical constants
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: one = 1.0_dp, half = 0.5_dp, twelve = 12.0_dp,  &
                         two = 2.0_dp, zero = 0.0_dp,  &
                         sqrtpi = 0.9189385332046727417803297_dp,  &
                         pi = 3.1415926535897932384626434_dp
!----------------------------------------------------------------------
!  Machine dependent parameters
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: xbig = 171.624_dp, xminin = 2.23E-308_dp,   &
                         eps = 2.22E-16_dp, xinf = 1.79E308_dp
!----------------------------------------------------------------------
!  Numerator and denominator coefficients for rational minimax
!     approximation over (1,2).
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: P(8) =  &
           (/ -1.71618513886549492533811E+0_dp,  2.47656508055759199108314E+1_dp,  &
              -3.79804256470945635097577E+2_dp,  6.29331155312818442661052E+2_dp,  &
               8.66966202790413211295064E+2_dp, -3.14512729688483675254357E+4_dp,  &
              -3.61444134186911729807069E+4_dp,  6.64561438202405440627855E+4_dp /)
REAL (dp), PARAMETER  :: Q(8) =  &
           (/ -3.08402300119738975254353E+1_dp,  3.15350626979604161529144E+2_dp,  &
              -1.01515636749021914166146E+3_dp, -3.10777167157231109440444E+3_dp,  &
               2.25381184209801510330112E+4_dp,  4.75584627752788110767815E+3_dp,  &
              -1.34659959864969306392456E+5_dp, -1.15132259675553483497211E+5_dp /)
!----------------------------------------------------------------------
!  Coefficients for minimax approximation over (12, INF).
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: c(7) =  &
           (/ -1.910444077728E-03_dp, 8.4171387781295E-04_dp,  &
              -5.952379913043012E-04_dp, 7.93650793500350248E-04_dp,  &
              -2.777777777777681622553E-03_dp, 8.333333333333333331554247E-02_dp,  &
               5.7083835261E-03_dp /)
!----------------------------------------------------------------------

parity = .false.
fact = one
n = 0
y = x
IF (y <= zero) THEN
!----------------------------------------------------------------------
!  Argument is negative
!----------------------------------------------------------------------
  y = -x
  y1 = AINT(y)
  fn_val = y - y1
  IF (fn_val /= zero) THEN
    IF (y1 /= AINT(y1*half)*two) parity = .true.
    fact = -pi / SIN(pi*fn_val)
    y = y + one
  ELSE
    fn_val = xinf
    GO TO 900
  END IF
END IF
!----------------------------------------------------------------------
!  Argument is positive
!----------------------------------------------------------------------
IF (y < eps) THEN
!----------------------------------------------------------------------
!  Argument < EPS
!----------------------------------------------------------------------
  IF (y >= xminin) THEN
    fn_val = one / y
  ELSE
    fn_val = xinf
    GO TO 900
  END IF
ELSE IF (y < twelve) THEN
  y1 = y
  IF (y < one) THEN
!----------------------------------------------------------------------
!  0.0 < argument < 1.0
!----------------------------------------------------------------------
    z = y
    y = y + one
  ELSE
!----------------------------------------------------------------------
!  1.0 < argument < 12.0, reduce argument if necessary
!----------------------------------------------------------------------
    n = INT(y) - 1
    y = y - n
    z = y - one
  END IF
!----------------------------------------------------------------------
!  Evaluate approximation for 1.0 < argument < 2.0
!----------------------------------------------------------------------
  xnum = zero
  xden = one
  DO  i = 1, 8
    xnum = (xnum+p(i)) * z
    xden = xden * z + q(i)
  END DO
  fn_val = xnum / xden + one
  IF (y1 < y) THEN
!----------------------------------------------------------------------
!  Adjust result for case  0.0 < argument < 1.0
!----------------------------------------------------------------------
    fn_val = fn_val / y1
  ELSE IF (y1 > y) THEN
!----------------------------------------------------------------------
!  Adjust result for case  2.0 < argument < 12.0
!----------------------------------------------------------------------
    DO  i = 1, n
      fn_val = fn_val * y
      y = y + one
    END DO
  END IF
ELSE
!----------------------------------------------------------------------
!  Evaluate for argument .GE. 12.0,
!----------------------------------------------------------------------
  IF (y <= xbig) THEN
    ysq = y * y
    sum = c(7)
    DO  i = 1, 6
      sum = sum / ysq + c(i)
    END DO
    sum = sum / y - y + sqrtpi
    sum = sum + (y-half) * LOG(y)
    fn_val = EXP(sum)
  ELSE
    fn_val = xinf
    GO TO 900
  END IF
END IF
!----------------------------------------------------------------------
!  Final adjustments and return
!----------------------------------------------------------------------
IF (parity) fn_val = -fn_val
IF (fact /= one) fn_val = fact / fn_val
900 RETURN
! ---------- Last line of GAMMA ----------
END FUNCTION LACZ_GAMMA

! Subroutine for detecting and reporting namelist errors.
! Prints the last line read before failiure, as debugging help.
subroutine checknamelisterror (ierr, ifnamopt, namelist)
  implicit none
  integer, intent(in) :: ierr, ifnamopt
  character(*), intent(in) :: namelist
  character(len=1000) :: line

  if (ierr > 0) then
     print *, 'Problem in namoptions ', namelist
     print *, 'iostat error: ', ierr
     backspace(ifnamopt)
     read(ifnamopt,fmt='(A)') line
     print *, 'Invalid line: '//trim(line)
     stop 'ERROR: Problem in namelist.'
  endif
end subroutine checknamelisterror

! Subroutine to check for netcdf error
subroutine handle_err(errcode)
use netcdf
implicit none
integer errcode

write(6,*) 'Error: ', nf90_strerror(errcode)
call abort
end subroutine handle_err

end module modglobal
