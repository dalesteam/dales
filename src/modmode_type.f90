! Definition of the mode type for the M7 scheme
module modmode_type

  use modprecision, only: field_r
  use modtracers,   only: T_tracer

  implicit none

  type mode_t
    !> General properties of the mode
    character(3)  :: name               !< Short name of the mode
    character(50) :: long_name          !< Full name of the mode
    integer       :: nspecies = -1      !< Number of species in the mode
    integer       :: idx_start, idx_end !< Index range in sv0
    integer       :: idx_numb           !< Index of the number concentration in sv0
    logical       :: lactivation        !< Mode has activation
    real(field_r) :: sigma_g            !< Geometric standard deviation lognormal distribution
                                        !! Has to be flexible at some point, probably 

    !> List of species, which are just tracers
    type(T_tracer), allocatable :: species

    !> Properties of the species
    real(field_r), allocatable :: rho_s(:)   !< Density
    real(field_r), allocatable :: kappa_s(:) !< Hygroscopicity
    
    !> Fields
    !! TODO: reduce memory footprint
    real(field_r), allocatable :: aer_conc(:,:,:,:) !< Mass/number concentration
    real(field_r), allocatable :: aer_tend(:,:,:,:) !< Total tendency
    real(field_r), allocatable :: aer_acti(:,:,:,:) !< Activation
    real(field_r), allocatable :: aer_scvc(:,:,:,:) !< Scavenging, cloud
    real(field_r), allocatable :: aer_evpc(:,:,:,:) !< Evaporation, cloud
    real(field_r), allocatable :: aer_scvr(:,:,:,:) !< Scavenging, rain
    real(field_r), allocatable :: aer_evpr(:,:,:,:) !< Evaporation, rain
    real(field_r), allocatable :: aer_auto(:,:,:,:) !< Autoconversion
    real(field_r), allocatable :: aer_accr(:,:,:,:) !< Accretion
    real(field_r), allocatable :: aer_slfc(:,:,:,:) !< Selfcollection, cloud
    real(field_r), allocatable :: aer_slfr(:,:,:,:) !< Selfcollection, rain
    real(field_r), allocatable :: aer_sedr(:,:,:,:) !< Sedimentation, rain
    real(field_r), allocatable :: aer_sedc(:,:,:,:) !< Sedimentation, cloud

  contains
    procedure, pass :: init => mode_init !< For memory allocation

  end type mode_t

contains

  !> Allocates all allocatable arrays of the mode.
  !!
  !! Has to be called after nspecies is set.
  !!
  subroutine mode_init(self)
    use modglobal, only : i1, ih, j1, jh, k1

    class(mode_t), intent(inout) :: self

    integer :: i, j, k, n

    if (self%nspecies < 0) then
      stop "M7 error: invalid number of species set"
    end if

    allocate(self%rho_s(self%nspecies))
    allocate(self%kappa_s(self%nspecies))

    allocate(self%aer_tend(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    allocate(self%aer_acti(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    allocate(self%aer_scvc(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    allocate(self%aer_evpc(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    allocate(self%aer_scvr(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    allocate(self%aer_evpr(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    allocate(self%aer_auto(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    allocate(self%aer_accr(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    allocate(self%aer_slfc(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    allocate(self%aer_slfr(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    allocate(self%aer_sedr(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    allocate(self%aer_sedc(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))

    self%aer_tend(:,:,:,:) = 0
    self%aer_acti(:,:,:,:) = 0
    self%aer_scvc(:,:,:,:) = 0
    self%aer_evpc(:,:,:,:) = 0
    self%aer_scvr(:,:,:,:) = 0
    self%aer_evpr(:,:,:,:) = 0
    self%aer_auto(:,:,:,:) = 0
    self%aer_accr(:,:,:,:) = 0 
    self%aer_slfc(:,:,:,:) = 0   
    self%aer_slfr(:,:,:,:) = 0
    self%aer_sedr(:,:,:,:) = 0
    self%aer_sedc(:,:,:,:) = 0

    !$acc enter data copyin(self%aer_tend(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    !$acc enter data copyin(self%aer_acti(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    !$acc enter data copyin(self%aer_scvc(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    !$acc enter data copyin(self%aer_evpc(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    !$acc enter data copyin(self%aer_scvr(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    !$acc enter data copyin(self%aer_evpr(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    !$acc enter data copyin(self%aer_auto(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    !$acc enter data copyin(self%aer_accr(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    !$acc enter data copyin(self%aer_slfc(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    !$acc enter data copyin(self%aer_slfr(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    !$acc enter data copyin(self%aer_sedr(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))
    !$acc enter data copyin(self%aer_sedc(2-ih:i1+ih,2-jh:j1+jh,k1,self%nspecies))

  end subroutine mode_init

end module modmode_type
