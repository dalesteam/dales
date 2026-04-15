!> Type definitions for M7 modes.
module modaerosol_mode_t

  use modaerosol_common, only: maxspecies, aerosol_densities, aerosol_names, &
                               aerosol_stdnames, mode_names, mode_longnames, &
                               sigma_g_modes
  use modglobal,         only: i1, j1, k1, kmax, ih, jh
  use modprecision,      only: field_r
  use modtracers,        only: get_tracer_index, add_tracer
  use modmpi,            only: print_info_stderr
  use modtimer,          only: timer_tic, timer_toc

  implicit none

  private

  character(len=*), parameter :: modname = 'modaerosol_mode_t'

  public :: mode_container_t
  public :: mode_t
  public :: aerosol_mode_t
  public :: hydrometeor_mode_t
  public :: mode_connection_t
  public :: connect_modes

  !> Mode container type, allowing per-element type variation in the list of
  !! modes.
  type :: mode_container_t
    class(mode_t), pointer :: p
  end type mode_container_t

  type :: mode_connection_t
    integer, allocatable :: cnct(:,:)
  end type mode_connection_t

  !> Base mode type
  type, abstract :: mode_t
    character(len=3) :: &
      name
    character(len=64) :: &
      longname = 'dummy'
    integer ::      &
      nspecies = 0, &
      itype(maxspecies)
    real(field_r) ::  &
      sig_g,          &
      rho(maxspecies)
  contains
    procedure(mode_t_init), deferred :: init
    procedure(mode_t_prepare), deferred :: prepare
    procedure(mode_t_finish), deferred :: finish
  end type mode_t

  interface
    subroutine mode_t_init(this, imode_type, lspecies)
      import :: mode_t
      class(mode_t), intent(inout) :: this
      integer, intent(in) :: imode_type
      logical, intent(in) :: lspecies(5)
    end subroutine mode_t_init
  end interface

  interface
    subroutine mode_t_prepare(this, sv)
      import :: mode_t, field_r
      class(mode_t), intent(inout) :: this
      real(field_r), intent(in) :: sv(:,:,:,:)
    end subroutine mode_t_prepare
  end interface

  interface
    subroutine mode_t_finish(this, svp, svm, delt)
      import :: mode_t, field_r
      class(mode_t), intent(inout) :: this
      real(field_r), intent(inout) :: svp(:,:,:,:)
      real(field_r), intent(in) :: svm(:,:,:,:)
      real(field_r), intent(in) :: delt
    end subroutine mode_t_finish
  end interface

  !> Free aerosol mode.
  type, extends(mode_t) :: aerosol_mode_t
    integer :: &
      itrac_n, &
      itrac_q(maxspecies)
    real(field_r), allocatable :: &
      n(:,:,:),               &
      np(:,:,:),              &
      q(:,:,:,:),             &
      qp(:,:,:,:)
    type(mode_connection_t) :: &
      to_hydro
  contains
    procedure :: init => aerosol_mode_init
    procedure :: prepare => aerosol_mode_prepare
    procedure :: finish => aerosol_mode_finish
  end type aerosol_mode_t

  !> Mode representing in-hydrometeor aerosol.
  !!
  !! Contrary to the regular mode type, this mode does not have a number
  !! number concentration associated with it, since this is handled by the
  !! microphysical scheme.
  type, extends(mode_t) :: hydrometeor_mode_t
    integer :: &
      itrac_q(maxspecies)
    real(field_r), pointer :: &
      q(:,:,:,:),             &
      qp(:,:,:,:) 
  contains
    procedure :: init => hydrometeor_mode_init
    procedure :: prepare => hydrometeor_mode_prepare
    procedure :: finish => hydrometeor_mode_finish
  end type hydrometeor_mode_t

contains

  !> Initialize a free aerosol mode.
  !!
  !! @param[in] imode_type Mode type identifier (see modaerosol_common.f90 for definitions).
  !! @param[in] lspecies List of enabled species for this mode.
  subroutine aerosol_mode_init(this, imode_type, lspecies)

    class(aerosol_mode_t), intent(inout) :: &
      this
    
    integer, intent(in) :: &
      imode_type

    logical, intent(in) :: &
      lspecies(maxspecies)

    integer :: &
      i, s

    this%name = mode_names(imode_type)
    this%longname = mode_longnames(imode_type)
    this%sig_g = sigma_g_modes(imode_type)
    this%nspecies = count(lspecies)

    this%itrac_q(:) = -1

    i = 1
    do s = 1, maxspecies
      if (lspecies(s)) then
        this%itype(i) = s
        i = i + 1
      end if
    end do

    ! Define tracers
    if (this%nspecies > 0) then
      call add_tracer(this%name//'_n', isv=this%itrac_n)

      do s = 1, this%nspecies
        this%rho(s) = aerosol_densities(this%itype(s))
        call add_tracer(trim(aerosol_names(this%itype(s)))//'_'//this%name, &
                        isv=this%itrac_q(s))
      end do

      ! Allocate memory, to be replaced by pointers to sv0 array?
      allocate(this%n(2:i1,2:j1,1:k1), &
               this%np(2:i1,2:j1,1:k1), &
               this%q(1:this%nspecies,2:i1,2:j1,1:k1), &
               this%qp(1:this%nspecies,2:i1,2:j1,1:k1))

    end if

    !$acc enter data copyin(this)
      !$acc enter data create(this%n(2:i1,2:j1,1:k1), &
      !$acc                   this%np(2:i1,2:j1,1:k1), &
      !$acc                   this%q(1:this%nspecies,2:i1,2:j1,1:k1), &
      !$acc                   this%qp(1:this%nspecies,2:i1,2:j1,1:k1))

  end subroutine aerosol_mode_init

  !> Initialize temporary memory before aerosol dynamics.
  !!
  !! @param[in] sv Tracer array.
  subroutine aerosol_mode_prepare(this, sv)

    class(aerosol_mode_t), intent(inout) :: &
      this

    real(field_r), intent(in) :: &
      sv(2-ih:,2-jh:,:,:)

    character(len=*), parameter :: routine = modname//'/aerosol_mode_prepare'

    integer :: &
      i, j, k, s

    call timer_tic(routine, 3)

    if (this%nspecies > 0) then

      !$acc parallel loop collapse(3) default(present) async wait(1)
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1
            this%n(i,j,k) = max(sv(i,j,k,this%itrac_n), 0.0_field_r)
            this%np(i,j,k) = 0
          end do
        end do
      end do

      !$acc parallel loop collapse(4) default(present) async wait(1)
      do s = 1, this%nspecies 
        do k = 1, kmax 
          do j = 2, j1
            do i = 2, i1
              this%q(s,i,j,k) = max(sv(i,j,k,this%itrac_q(s)), 0.0_field_r)
              this%qp(s,i,j,k) = 0
            end do
          end do
        end do
      end do

    end if

    call timer_toc(routine)

  end subroutine aerosol_mode_prepare

  !> Copy out tendencies.
  !!
  !! @param[inout] svp Tracer tendency array.
  subroutine aerosol_mode_finish(this, svp, svm, delt)

    class(aerosol_mode_t), intent(inout) :: &
      this

    real(field_r), intent(inout) :: &
      svp(2-ih:,2-jh:,:,:)

    real(field_r), intent(in) :: &
      svm(2-ih:,2-jh:,:,:)

    real(field_r), intent(in) :: delt

    character(len=*), parameter :: routine = modname//"/aerosol_mode_finish"

    integer :: &
      i, j, k, s ! Loop indices

    real(field_r) :: sv_cor

    call timer_tic(routine, 3)

    if (this%nspecies > 0) then

      !$acc parallel loop collapse(3) default(present) async wait(1)
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1
            sv_cor = min(svp(i,j,k,this%itrac_n) + this%np(i,j,k) &
                         + (svm(i,j,k,this%itrac_n) / delt), &
                         0.0_field_r)
            svp(i,j,k,this%itrac_n) = svp(i,j,k,this%itrac_n) + this%np(i,j,k) - sv_cor
          end do
        end do
      end do

      !$acc parallel loop collapse(4) default(present) async wait(1)
      do s = 1, this%nspecies
        do k = 1, kmax
          do j = 2, j1
            do i = 2, i1
              sv_cor = min(svp(i,j,k,this%itrac_q(s)) + this%qp(s,i,j,k) &
                           + (svm(i,j,k,this%itrac_q(s)) / delt), &
                           0.0_field_r)
              svp(i,j,k,this%itrac_q(s)) = svp(i,j,k,this%itrac_q(s)) &
                                           + this%qp(s,i,j,k) - sv_cor
            end do
          end do
        end do
      end do

    end if

    call timer_toc(routine)

  end subroutine aerosol_mode_finish

  !> Initialize an in-hydrometeor aerosol mode.
  !!
  !! @param[in] imode_type Mode type identifier (see modaerosol_common.f90 for definitions).
  !! @param[in] lspecies List of enabled species for this mode.
  subroutine hydrometeor_mode_init(this, imode_type, lspecies)

    class(hydrometeor_mode_t), intent(inout) :: &
      this

    integer, intent(in) :: &
      imode_type

    logical, intent(in) :: &
      lspecies(maxspecies)

    integer :: &
      i, s

    this%name = mode_names(imode_type)
    this%longname = mode_longnames(imode_type)
    this%sig_g = sigma_g_modes(imode_type)
    this%nspecies = count(lspecies)

    i = 1
    do s = 1, maxspecies    
      if (lspecies(s)) then
        this%itype(i) = s
        i = i + 1
      end if
    end do

    ! Define tracers
    do s = 1, this%nspecies
      this%rho(s) = aerosol_densities(this%itype(s))
      call add_tracer(trim(aerosol_names(this%itype(s)))//'_'//this%name, &
                      isv=this%itrac_q(s))
    end do

    ! Allocate memory, to be replaced by pointers to sv0 array
    allocate(this%q(1:this%nspecies,2:i1,2:j1,1:k1), &
             this%qp(1:this%nspecies,2:i1,2:j1,1:k1))

    this%q(:,:,:,:) = 0
    this%qp(:,:,:,:) = 0

    !$acc enter data copyin(this, this%q(1:this%nspecies,2:i1,2:j1,1:k1), &
    !$acc                   this%qp(1:this%nspecies,2:i1,2:j1,1:k1))

  end subroutine hydrometeor_mode_init

  !> Initialize temporary memory before aerosol dynamics.
  !!
  !! @param[in] sv Tracer array.
  subroutine hydrometeor_mode_prepare(this, sv)

    class(hydrometeor_mode_t), intent(inout) :: &
      this

    real(field_r), intent(in) :: &
      sv(2-ih:,2-jh:,:,:)

    character(len=*), parameter :: routine = &
      modname//'/hydrometeor_mode_prepare'

    integer :: &
      i, j, k, s

    call timer_tic(routine, 3)

    !$acc parallel loop collapse(4) default(present) async wait(1)
    do s = 1, this%nspecies 
      do k = 1, kmax 
        do j = 2, j1
          do i = 2, i1
            this%q(s,i,j,k) = max(sv(i,j,k,this%itrac_q(s)), 0.0_field_r)
            this%qp(s,i,j,k) = 0
          end do
        end do
      end do
    end do

    call timer_toc(routine)

  end subroutine hydrometeor_mode_prepare
  
  !> Copy out tendencies.
  !!
  !! @param[inout] svp Tracer tendency array.
  subroutine hydrometeor_mode_finish(this, svp, svm, delt)

    class(hydrometeor_mode_t), intent(inout) :: &
      this

    real(field_r), intent(inout) :: &
      svp(2-ih:,2-jh:,:,:)
    real(field_r), intent(in) :: &
      svm(2-ih:,2-jh:,:,:)

    real(field_r), intent(in) :: delt

    character(len=*), parameter :: routine = modname//'/hydrometeor_mode_finish'

    integer :: &
      i, j, k, s

    real(field_r) :: sv_cor

    call timer_tic(routine, 3)

    !$acc parallel loop collapse(4) default(present) async wait(1)
    do s = 1, this%nspecies
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1
            sv_cor = min(svp(i,j,k,this%itrac_q(s)) + this%qp(s,i,j,k) &
                         + (svm(i,j,k,this%itrac_q(s)) / delt), &
                         0.0_field_r)
            svp(i,j,k,this%itrac_q(s)) = svp(i,j,k,this%itrac_q(s)) &
                                         + this%qp(s,i,j,k) - sv_cor
          end do
        end do
      end do
    end do

    call timer_toc(routine)

  end subroutine hydrometeor_mode_finish

  !> Setup a mode_connection between two modes.
  !!
  !! @param[in] from Source mode.
  !! @param[in] to Target mode.
  !! @param[out] connection Connection between the modes.
  subroutine connect_modes(from, to, connection)

    class(mode_t), intent(in) :: &
      from, to

    type(mode_connection_t), intent(out) :: &
      connection

    integer :: &
      sf, st ! Source and target indices

    allocate(connection%cnct(2,from%nspecies))

    do sf = 1, from%nspecies
      do st = 1, to%nspecies 
        if (from%itype(sf) == to%itype(st)) then
          connection%cnct(1,sf) = sf
          connection%cnct(2,sf) = st
        end if
      end do
    end do

  end subroutine connect_modes

end module modaerosol_mode_t
