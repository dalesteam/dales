module modaerosol_mode_t

  use modaerosol_common, only: maxspecies, aerosol_densities, aerosol_names
  use modglobal,         only: i1, j1, k1, ih, jh
  use modprecision,      only: field_r
  use modtracers,        only: get_tracer_index
  use modmpi,            only: print_info_stderr

  implicit none

  character(len=*), parameter :: modname = 'modaerosol_mode_t'

  type, public :: mode_t
    logical                :: &
      lactive                   !< Mode is active.
    character(len=3)       :: &
      name                      !< Mode name.
    integer                :: &
      nspecies,               & !< Number of aerosol species in mode.
      itrac_n,                & !< Index of number in tracer container.
      itrac_q(maxspecies),    & !< Indices of mass in tracer container.
      itype(maxspecies)         !< Aerosol type
    real(field_r)          :: &
      sig_g,                  & !< Geometric standard deviation.
      rho(maxspecies)           !< Aerosol densities
    real(field_r), pointer :: &
      n(:,:,:),               & !< Number concentration.
      np(:,:,:),              & !< Tendency of number concentration.
      q(:,:,:,:),             & !< Mass concentrations.
      qp(:,:,:,:)               !< Tendency of mass concentrations.
  contains
    procedure :: construct => mode_construct
    procedure :: add_aerosol => mode_add_aerosol
    procedure :: allocate => mode_allocate
    procedure :: prepare => mode_prepare
    procedure :: finalize => mode_finalize
  end type mode_t

  integer, parameter :: &
    iNUS = 1, &
    iAIS = 2, &
    iACS = 3, &
    iCOS = 4, &
    iAII = 5, &
    iACI = 6, &
    iCOI = 7, &
    iSO4 = 1, &
    iSS = 2,  &
    iPOM = 3, &
    iBC = 4,  &
    iDU = 5

contains

  !> Construct a mode
  !!
  !! @param[in] name Mode short name.
  !! @param[in] sig_g Geometric standard deviation.
  subroutine mode_construct(this, name, sig_g)
    class(mode_t),    intent(inout) :: this
    character(len=3), intent(in)    :: name
    real(field_r),    intent(in)    :: sig_g

    this%name = name
    this%sig_g = sig_g

  end subroutine mode_construct

  !> Add an aerosol species to a mode.
  !!
  !! @param[in] itype Aerosol type.
  subroutine mode_add_aerosol(this, itype)
    class(mode_t), intent(inout) :: this
    integer,       intent(in)    :: itype

    character(len=*), parameter :: routine = modname//'/mode_add_aerosol'

    integer :: isv
    character(len=3) :: specname

    if (.not. this%lactive) then
      this%lactive = .true.
      this%nspecies = 0

      ! Find the tracer for the number concentration
      this%itrac_n = get_tracer_index(this%name//'_n')
      
      if (this%itrac_n < 1) then
        call print_info_stderr(routine, &
          'could not find existing tracer for '//this%name//'_n')
        error stop
      end if
    end if

    this%nspecies = this%nspecies + 1
    this%itype(this%nspecies) = itype
    this%rho(this%nspecies) = aerosol_densities(itype)

    specname = aerosol_names(itype)
    isv = get_tracer_index(trim(specname)//'_'//this%name)

    if (isv < 1) then
      call print_info_stderr(routine, &
        'could not find existing tracer for '//trim(specname)//&
        &' in mode '//this%name)
      error stop
    end if

    this%itrac_q(this%nspecies) = isv

  end subroutine mode_add_aerosol

  !> Allocate memory for mass/number concentrations and tendencies.
  subroutine mode_allocate(this)
    class(mode_t), intent(inout) :: this

    if (.not. this%lactive) return

    allocate(this%n(2:i1,2:j1,1:k1), this%np(2:i1,2:j1,1:k1), &
      this%q(2:i1,2:j1,1:k1,this%nspecies), &
      this%qp(2:i1,2:j1,1:k1,this%nspecies))

    this%n(:,:,:) = 0
    this%np(:,:,:) = 0
    this%q(:,:,:,:) = 0
    this%qp(:,:,:,:) = 0

    !$acc enter data copyin(this%n, this%np, this%q, this%qp)

  end subroutine mode_allocate

  !> Copies aerosol fields to work space.
  !!
  !! @param[in] sv Tracer fields.
  subroutine mode_prepare(this, sv)
    class(mode_t), intent(inout) :: this
    real(field_r), intent(in)    :: sv(2-ih:,2-jh:,1:,1:)

    integer :: i, j, k, s

    !$acc parallel loop gang vector collapse(3) default(present) async
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          this%n(i,j,k) = max(sv(i,j,k,this%itrac_n), 0.0_field_r)
        end do
      end do
    end do

    !$acc parallel loop gang vector collapse(4) default(present) async
    do s = 1, this%nspecies
      do k = 1, k1
        do j = 2, j1
          do i = 2, i1
            this%q(i,j,k,s) = max(sv(i,j,k,this%itrac_q(s)), 0.0_field_r)
          end do
        end do
      end do
    end do

  end subroutine mode_prepare

  !> Copies computed tendencies to svp fields.
  !!
  !! @param[inout] svp Tracer tendency fields.
  subroutine mode_finalize(this, svp)
    class(mode_t), intent(inout) :: this
    real(field_r), intent(inout) :: svp(2-ih:,2-jh:,1:,1:)

    integer :: i, j, k, s

    !$acc parallel loop gang vector collapse(3) default(present) async
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          svp(i,j,k,this%itrac_n) = svp(i,j,k,this%itrac_n) + this%np(i,j,k)
          this%np(i,j,k) = 0
        end do
      end do
    end do

    !$acc parallel loop gang vector collapse(4) default(present) async
    do s = 1, this%nspecies
      do k = 1, k1
        do j = 2, j1
          do i = 2, i1
            svp(i,j,k,this%itrac_q(s)) = svp(i,j,k,this%itrac_q(s)) + &
                                         this%qp(i,j,k,s)
            this%qp(i,j,k,s) = 0
          end do
        end do
      end do
    end do

  end subroutine mode_finalize

end module modaerosol_mode_t