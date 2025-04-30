module modaerosol
  use modglobal,      only: ifnamopt, fname_options, checknamelisterror, &
                            cexpnr, i1, j1, k1, ih, jh, pi, nsv, rhow, kmax
  use modmath,        only: inv_sqrt_two
  use modmicrodata,   only: qcmin
  use modfields,      only: sv0, svp, svm
  use modmicrodata,   only: qcmin, delt
  use modmpi,         only: myid, D_MPI_BCAST, commwrld, mpierr
  use modprecision,   only: field_r
  use modtracers,     only: add_tracer, allocate_tracers, tracer_prop, &
                            get_tracer_index
  use modstat_nc
  use go,             only: goSplitString_s
  use modtimer,       only: timer_tic, timer_toc
  use modlookuptable, only: LT2_t, LT2_create, LT2_set_col, LT2_get_col, &
                            LT2_get_col_inline

  implicit none

  private

  character(len=*), parameter :: modname = 'modaerosol'

  public :: laerosol

  ! Metadata for the modes and aerosols
  integer,           parameter :: &
    maxmodes = 7, &
    maxspecies = 5
  character(len=3),  parameter :: modenames(maxmodes) = &
    ['nus', 'ais', 'acs', 'cos', 'aii', 'aci', 'coi' ]
  character(len=22), parameter :: longnames(maxmodes) = [ &
    'soluble nucleation    ', &
    'soluble Aitken        ', &
    'soluble accumulation  ', &
    'soluble coarse        ', &
    'insoluble Aitken      ', &
    'insoluble accumulation', &
    'insoluble coarse      ' ]
  character(len=3),  parameter :: aerosol_names(maxspecies) = &
    [ 'so4', 'ss ', 'pom', 'bc ', 'du ' ]
  character(len=26), parameter :: aerosol_longnames(maxspecies) = [ &
    'sulfate                   ', &
    'sea salt                  ', &
    'particulate organic matter', &
    'black carbon              ', &
    'dust                      ' ]
  real(field_r), parameter :: eps = 1e-18

  ! Static indices
  integer, parameter :: &
    iNUS = 1, &
    iAIS = 2, &
    iACS = 3, &
    iCOS = 4, &
    iAII = 5, &
    iACI = 6, &
    iCOI = 7
  integer, parameter :: &
    iSO4 = 1, &
    iSS = 2,  &
    iPOM = 3, &
    iBC = 4,  &
    iDU = 5

  type :: mode_connection_t
    logical              :: ldoshift   !< Do shift for this mode
    integer, allocatable :: itarget(:) !< Index of species in target mode
  end type mode_connection_t

  type, public :: mode_t
    logical                 :: lactive     !< Mode is active
    character(len=3)        :: name        !< Short name
    character(len=27)       :: longname    !< Long name
    integer                 :: nspecies    !< Number of aerosol species
    real(field_r)           :: sig_g       !< Geometric standard deviation
    real(field_r), pointer  :: rho(:)      !< Aerosol densities
    real(field_r), pointer  :: itrac(:)    !< Tracer index
    real(field_r), pointer  :: itype(:)    !< Aerosol type
    real(field_r), pointer  :: n(:,:,:)    !< Number concentration
    real(field_r), pointer  :: q(:,:,:,:)  !< Mass concentration
    real(field_r), pointer  :: np(:,:,:)   !< Number concentration tendency
    real(field_r), pointer  :: qp(:,:,:,:) !< Mass concentration tendency
    type(mode_connection_t) :: shift2cloud
    type(mode_connection_t) :: shift2larger
    type(mode_connection_t) :: shift2free
    type(mode_connection_t) :: shift2soluble
  contains
    procedure :: construct => mode_construct
    procedure :: add_aerosol => mode_add_aerosol
    procedure :: allocate => mode_allocate
    procedure :: prepare => mode_prepare
    procedure :: finalize => mode_finalize
  end type mode_t

  ! Public variables
  logical :: &
    laerosol, & !< Switch for enabling aerosol scheme
    lso4,     & !< Switch for enabling sulphur
    lss,      & !< Switch for enabling sea salt
    lpom,     & !< Switch for enabling particulate organic matter
    lbc,      & !< Switch for enabling black carbon
    ldu         !< Switch for enabling dust

  ! Data
  type(mode_t) :: &
    modes(maxmodes)
  ! TODO: make a mode type for these again
  ! and separate soluble from insoluble
  logical :: &
    species_active(maxspecies)
  integer :: &
    n_species_active,     &
    inc_type(maxspecies), &
    inc_idx(maxspecies)
  real(field_r), pointer :: &
    qa_inc(:,:,:,:), &
    qa_inr(:,:,:,:), &
    qap_inc(:,:,:,:), &
    qap_inr(:,:,:,:)

  ! Lookup tables for scavenging
  type(LT2_t) :: &
    inc_tab_m, &
    inc_tab_n, &
    blc_tab_m, &
    blc_tab_n

contains

  function aerosol_get_index_in_mode(itype, mode) result(idx)

    integer,      intent(in) :: itype
    type(mode_t), intent(in) :: mode
    !$acc routine seq

    integer :: idx

    ! Linear search
    do idx = 1, mode%nspecies
      if (itype == mode%itype(idx)) return
    end do

    idx = -1

  end function aerosol_get_index_in_mode

  !> Get array index of aerosol in the in-cloud and in-rain categories
  function aerosol_get_index_in_cloud(itype) result(idx)

    integer, intent(in) :: itype
    !$acc routine seq

    integer :: idx

    idx = inc_idx(itype)

  end function aerosol_get_index_in_cloud

  function aerosol_get_type_in_cloud(idx) result(itype)

    integer, intent(in) :: idx
    !$acc routine seq

    integer :: itype

    itype = inc_type(idx)

  end function aerosol_get_type_in_cloud

  subroutine init_aerosol()

    character(len=*), parameter :: routine = modname//"/initaerosol"

    integer       :: imod, ierr, ncid, nvars, iaer, mode_loc
    character(3)  :: name
    character(64) :: long_name
    character(27) :: modes_str
    real(field_r) :: rho, kappa
    character(3)  :: modes_list(maxmodes)
    integer       :: aero_idx_in_mode
    real(field_r), parameter :: sigma_g(maxmodes) = (/ 1.59, 1.59, 1.59, 2.00, &
      1.59, 1.59, 2.00/)
    real(field_r), parameter :: cldrad(10) = log([5., 10., 15., 20., 25., 30., &
      35., 40., 45., 50.])
    real(field_r), parameter :: rainrate(5) = log([0.01, 0.1, 1., 10., 100.])

    integer, allocatable :: varids(:)

    ! Values for lookup tables
    include "scavenging.inc"

    namelist /NAMAEROSOL/ laerosol, lso4, lss, lpom, lbc, ldu

    ! Read input
    if (myid == 0) then
      ! Namelist
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read(ifnamopt, NAMAEROSOL, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMAEROSOL')
      close(ifnamopt)
    end if

    call d_mpi_bcast(laerosol, 1, 0, commwrld, mpierr)
    call d_mpi_bcast(lso4, 1, 0, commwrld, mpierr)
    call d_mpi_bcast(lss, 1, 0, commwrld, mpierr)
    call d_mpi_bcast(lpom, 1, 0, commwrld, mpierr)
    call d_mpi_bcast(lbc, 1, 0, commwrld, mpierr)
    call d_mpi_bcast(ldu, 1, 0, commwrld, mpierr)

    if (.not. laerosol) return

    ! Setup the modes
    do imod = 1, maxmodes
      call mode_construct(modes(imod), name=modenames(imod), &
        long_name=longnames(imod), sig_g=sigma_g(imod))
    end do

    !                | NUS | AIS | ACS | COS | AII | ACI | COI
    ! SO4            |  x  |  x  |  x  |  x  |     |     |
    ! Sea Salt       |     |     |  x  |  x  |     |     |
    ! Organic Matter |     |  x  |  x  |  x  |  x  |     |
    ! Black Carbon   |     |  x  |  x  |  x  |  x  |     |
    ! Dust           |     |     |  x  |  x  |     |  x  |  x

    ! Sulphuric acid
    if (lso4) then
      block
        integer :: my_modes(4) = [iNUS, iAIS, iACS, iCOS]
        do imod = 1, size(my_modes)
          call mode_add_aerosol(modes(my_modes(imod)), itype=iso4)
        end do
        call add_tracer('so4_c', long_name='so4 in-cloud mass concentration', unit='kg/kg')
        call add_tracer('so4_r', long_name='so4 in-rain mass concentration', unit='kg/kg')
      end block
    end if

    ! Sea salt
    if (lss) then
      block
        integer :: my_modes(2) = [iACS, iCOS]
        do imod = 1, size(my_modes)
          call mode_add_aerosol(modes(my_modes(imod)), itype=iss)
        end do
        call add_tracer('ss_c', long_name='sea salt in-cloud mass concentration', unit='kg/kg')
        call add_tracer('ss_r', long_name='sea salt in-rain mass concentration', unit='kg/kg')
      end block
    end if

    ! Particulate organic matter
    if (lpom) then
      block
        integer :: my_modes(4) = [iAIS, iACS, iCOS, iAII]
        do imod = 1, size(my_modes)
          call mode_add_aerosol(modes(my_modes(imod)), itype=ipom)
        end do
        call add_tracer('pom_c', long_name='organic matter in-cloud mass concentration', unit='kg/kg')
        call add_tracer('pom_r', long_name='organic matter in-rain mass concentration', unit='kg/kg')
      end block
    end if

    ! Black carbon
    if (lbc) then
      block
        integer :: my_modes(4) = [iAIS, iACS, iCOS, iAII]
        do imod = 1, size(my_modes)
          call mode_add_aerosol(modes(my_modes(imod)), itype=ibc)
        end do
        call add_tracer('bc_c', long_name='black carbon in-cloud mass concentration', unit='kg/kg')
        call add_tracer('bc_r', long_name='black carbon in-rain mass concentration', unit='kg/kg')
      end block
    end if

    ! Mineral dust
    if (ldu) then
      block
        integer :: my_modes(4) = [iACS, iCOS, iACI, iCOI]
        do imod = 1, size(my_modes)
          call mode_add_aerosol(modes(my_modes(imod)), itype=idu)
        end do
        call add_tracer('du_c', long_name='mineral dust in-cloud mass concentration', unit='kg/kg')
        call add_tracer('du_r', long_name='mineral dust in-rain mass concentration', unit='kg/kg')
      end block
    end if

    ! Compute indices of aerosols in the array of in-cloud mass
    block

      integer :: i, j
      integer :: index

      species_active(1) = lso4
      species_active(2) = lss
      species_active(3) = lpom
      species_active(4) = lbc
      species_active(5) = ldu

      n_species_active = count(species_active, dim=1)

      do i = 1, maxspecies
        index = 0
        if (.not. species_active(i)) cycle
        do j = 1, i
          index = index + merge(1, 0, species_active(j))
        end do
        inc_idx(i) = index
        if (index > 0) inc_type(index) = i
      end do

    end block

    ! Finally, allocate memory
    do imod = 1, maxmodes
      call mode_allocate(modes(imod))
    end do

    ! "Temporary" arrays for in-cloud and in-rain categories of aerosol
    allocate(qa_inc(2:i1,2:j1,k1,n_species_active), &
             qa_inr(2:i1,2:j1,k1,n_species_active), &
             qap_inc(2:i1,2:j1,k1,n_species_active), &
             qap_inr(2:i1,2:j1,k1,n_species_active))

    qa_inc = 0
    qa_inr = 0
    qap_inc = 0
    qap_inr = 0

    ! Setup the lookup tables for scavenging routines
    inc_tab_m = LT2_create([cldrad(1), aerrad(1)], [cldrad(10), aerrad(60)], &
                           [10, 60], 1)
    call LT2_set_col(inc_tab_m, 1, cldrad, aerrad, scavenging_eff_incloud_m)

    inc_tab_n = LT2_create([cldrad(1), aerrad(1)], [cldrad(10), aerrad(60)], &
                           [10, 60], 1)
    call LT2_set_col(inc_tab_n, 1, cldrad, aerrad, scavenging_eff_incloud_n)

    blc_tab_m = LT2_create([rainrate(1), aerrad_blc(1)], &
                           [rainrate(5), aerrad_blc(100)], &
                           [5, 100], 1)
    call LT2_set_col(blc_tab_m, 1, rainrate, aerrad_blc, &
                     scavenging_eff_belowcloud_m)

    blc_tab_n = LT2_create([rainrate(1), aerrad_blc(1)], &
                           [rainrate(5), aerrad_blc(100)], &
                           [5, 100], 1)
    call LT2_set_col(blc_tab_n, 1, rainrate, aerrad_blc, &
                     scavenging_eff_belowcloud_n)

  end subroutine init_aerosol

  ! Prepares aerosol fields for microphysics calculations
  subroutine aerosol_prepare

    character(len=*), parameter :: routine = modname//'/aerosol_prepare'

    integer :: i, j, k, s, imod
    integer :: itype, idx_c, idx_r

    if (.not. laerosol) return

    call timer_tic(routine, 1)

    ! TODO: Possible optimization: replace this with pointers
    ! need to make sure that the in-cloud species are contiguous in sv array
    ! or: copy them while transposing to (s,k,j,i)

    ! Copy ambient mass and number concentrations to temp fields
    do imod = 1, maxmodes
      call modes(imod)%prepare(sv0)
    end do

    do s = 1, n_species_active
      itype = aerosol_get_type_in_cloud(s)
      idx_c = get_tracer_index(trim(aerosol_names(itype))//"_c")
      idx_r = get_tracer_index(trim(aerosol_names(itype))//"_r")
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1
            ! Copy mass concentrations
            qa_inc(i,j,k,s) = max(sv0(i,j,k,idx_c), 0.0_field_r)
            qa_inr(i,j,k,s) = max(sv0(i,j,k,idx_r), 0.0_field_r)
            ! Reset tendency fields
            qap_inc(i,j,k,s) = 0
            qap_inr(i,j,k,s) = 0
          end do
        end do
      end do
    end do

    call timer_toc(routine)

  end subroutine aerosol_prepare

  subroutine aerosol_finalize

    character(len=*), parameter :: routine = modname//'/aerosol_finalize'

    integer :: i, j, k, s, imod
    integer :: itype, idx_c, idx_r

    if (.not. laerosol) return

    call timer_tic(routine, 1)

    if (laerosol) then
      do imod = 1, maxmodes
        call modes(imod)%finalize(svp, svm, delt)
      end do
    end if

    do s = 1, n_species_active
      itype = aerosol_get_type_in_cloud(s)
      idx_c = get_tracer_index(trim(aerosol_names(itype))//"_c")
      idx_r = get_tracer_index(trim(aerosol_names(itype))//"_r")
      do k = 1, k1
        do j = 2, j1
          do i = 2, i1
            svp(i,j,k,idx_c) = svp(i,j,k,idx_c) + qap_inc(i,j,k,s)
            svp(i,j,k,idx_r) = svp(i,j,k,idx_r) + qap_inr(i,j,k,s)
          end do
        end do
      end do
    end do

    call timer_toc(routine)

  end subroutine aerosol_finalize

  !> Construct a mode
  !!
  !! \param name Short name.
  !! \param long_name Long name.
  !! \param sigma_g Geometric standard deviation.
  subroutine mode_construct(self, name, long_name, sig_g)
    class(mode_t), intent(inout) :: self
    character(3),  intent(in)    :: name
    character(*),  intent(in)    :: long_name
    real(field_r), intent(in)    :: sig_g

    self % name = trim(name)
    self % longname = trim(long_name)
    self % sig_g = sig_g

  end subroutine mode_construct

  subroutine mode_add_aerosol(self, itype)
    class(mode_t), intent(inout) :: self
    integer, intent(in) :: itype

    integer :: isv

    self%lactive = .true.

    self%nspecies = self%nspecies + 1
    self%itype(self%nspecies) = itype

    ! Setup a tracer for the mass concentration
    call add_tracer(trim(aerosol_names(itype))//"_"//self%name, &
      long_name=aerosol_longnames(itype), isv=isv)

    self%itrac(self%nspecies) = isv

  end subroutine mode_add_aerosol

  !> Allocate memory for mass/number concentrations and tendencies.
  subroutine mode_allocate(self)
    class(mode_t), intent(inout) :: self

    ! Make sure static data is available on GPU, even if we don't use this mode
    !$acc enter data copyin(self, self%enabled)

    if (self%nspecies < 1) return

    allocate(self%n(2:i1,2:j1,1:k1), self%np(2:i1,2:j1,1:k1), &
      self%q(2:i1,2:j1,1:k1,self%nspecies), &
      self%qp(2:i1,2:j1,1:k1,self%nspecies))

    self%n(:,:,:) = 0
    self%np(:,:,:) = 0
    self%q(:,:,:,:) = 0
    self%qp(:,:,:,:) = 0

    !$acc enter data copyin(self%n, self%np, self%q, self%qp)

  end subroutine mode_allocate

  !> Copies aerosol fields to work space.
  !!
  !! \param sv Tracer fields.
  subroutine mode_prepare(self, sv)
    class(mode_t), intent(inout) :: self
    real(field_r), intent(in)    :: sv(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv)

    integer :: iaer, sv_idx
    integer :: i, j, k

    if (self % nspecies < 1) return

    sv_idx = get_tracer_index(self%name//'_n')

    !$acc parallel loop gang vector collapse(3) default(present) async
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          self%n(i,j,k) = max(sv(i,j,k,sv_idx), 0.0_field_r)
        end do
      end do
    end do

    ! And the mass concentrations
    do iaer = 1, self%nspecies
      sv_idx = get_tracer_index(aerosol_names(self%itype(iaer))//'_'//self%name)
      !$acc parallel loop gang vector collapse(3) default(present) async
      do k = 1, k1
        do j = 2, j1
          do i = 2, i1
            self%q(i,j,k,iaer) = max(sv(i,j,k,sv_idx), 0.0_field_r)
          end do
        end do
      end do
    end do

    !$acc wait

  end subroutine mode_prepare

  !> Copies computed tendencies to svp fields.
  !!
  !! \param svp Tracer tendency fields.
  subroutine mode_finalize(self, svp, svm, delt)
    class(mode_t), intent(inout) :: self
    real(field_r), intent(inout) :: svp(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv)
    real(field_r), intent(in)    :: svm(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv)
    real(field_r), intent(in)    :: delt

    integer :: iaer, sv_idx
    integer :: i, j, k
    character(3) :: name

    if (self % nspecies < 1) return

    sv_idx = get_tracer_index(self%name//'_n')

    !$acc parallel loop gang vector collapse(3) default(present) async
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          svp(i,j,k,sv_idx) = svp(i,j,k,sv_idx) + &
            max(self%np(i,j,k), -svm(i,j,k,sv_idx) / delt)
          self%np(i,j,k) = 0
        end do
      end do
    end do

    do iaer = 1, self%nspecies
      sv_idx = get_tracer_index(aerosol_names(self%itype(iaer))//'_'//self%name)
      !$acc parallel loop gang vector collapse(3) default(present) async
      do k = 1, k1
        do j = 2, j1
          do i = 2, i1
            svp(i,j,k,sv_idx) = svp(i,j,k,sv_idx) + &
              max(self%qp(i,j,k,iaer), -svm(i,j,k,sv_idx) / delt)
            self%qp(i,j,k,iaer) = 0
          end do
        end do
      end do
    end do

    !$acc wait

  end subroutine mode_finalize

end module modaerosol