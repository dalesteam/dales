!> \file modibmdata.f90
!! By Michael Koene (MK), TU Delft, section Atmospheric Physics, 28 January 2019
!! cstep: subroutine airslabsum  moved from modmpi to here to avoid mutual dependencies

module modibmdata

  implicit none
  save

  logical :: lapply_ibm     = .false.        !< Switch to enable immersed boundary method 
  logical :: lreadfile_obstacles  = .false.  !< Switch to read positions and model level height from a file
  logical :: lwallheat      = .false.        !< Switch to apply lateral heat flux from buildings
  logical :: lpoislast      = .true.         !< Switch to use the Poisson solver after the Immersed boundary method
                                             !  .false. will set the order to: ZeroVelocity -> PoissonSolver -> IBM

  real    :: thlwall        = 293.           !< Wall temperature for temperature flux at the sides of the buildings, needed for lateral flux
  real    :: thlroof        = 293.           !< Obstacle roof (top) temperature
  real    :: qtroof         = 0.             !< Obstacle roof specific humidity
  real    :: thlibm         = 293            !< Interior potential temperature of obstacle
  real    :: qtibm          = 0.             !< Wall specific humidity for the latent heat flux at the sides and top of the buildings
                                             !< In modsurface it will be set to the saturation value (but this needs to be adapted)
  real    :: z0m_wall       = 0.03           !< compare with 0.03 m for open flat terrain, grass, few isolated obstacles

  !< Boolean for applying IBM
  logical, allocatable :: libm (:,:,:)       !< Is true inside a building point
  !< Mask field for advection. IBM works with 2nd order advection scheme
  integer, allocatable :: ibm_adv_mask (:,:,:)
  !< Number of grid points in a slab excluding obstacles, and the number of obstacle points
  integer, allocatable :: Nair (:)
  integer :: kibm_max                        !< index of vertical layer that contains the highest obstacle        


end module modibmdata
