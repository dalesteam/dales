module modaerosol

  use, intrinsic :: iso_fortran_env

  use modglobal,      only: ifnamopt, fname_options, checknamelisterror, &
                            cexpnr, i1, j1, k1, ih, jh, pi, nsv, rhow, kmax
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

  public :: aerosol_read_namelist
  public :: init_aerosol
  public :: aerosol_prepare
  public :: aerosol_finalize
  public :: activation
  public :: aerosol_cloud_to_rain
  public :: aerosol_redistribute

  public :: scavenging_cloud

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

  interface erfcinv
    module procedure :: erfcinv_real32
    module procedure :: erfcinv_real64
  end interface

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
    real(field_r)           :: rho(maxspecies)     !< Aerosol densities
    integer                 :: itrac(maxspecies)    !< Tracer index
    integer                 :: itype(maxspecies)    !< Aerosol type
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
  type(mode_t), target :: &
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

  include 'erfcinv.inc'

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

  subroutine aerosol_read_namelist(nml_filename)

    character(len=*), intent(in) :: nml_filename

    integer :: ierr

    namelist /NAMAEROSOL/ laerosol, lso4, lss, lpom, lbc, ldu

    if (myid == 0) then
      ! Namelist
      open(ifnamopt, file=nml_filename, status='old', iostat=ierr)
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

  end subroutine aerosol_read_namelist

  subroutine init_aerosol()

    character(len=*), parameter :: routine = modname//"/initaerosol"

    integer       :: imod, ierr, ncid, nvars, iaer, mode_loc
    character(3)  :: name
    character(64) :: long_name
    character(27) :: modes_str
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
        real(field_r) :: rho = 1841
        do imod = 1, size(my_modes)
          call mode_add_aerosol(modes(my_modes(imod)), iso4, rho)
        end do
        call add_tracer('so4_c', long_name='so4 in-cloud mass concentration', unit='kg/kg')
        call add_tracer('so4_r', long_name='so4 in-rain mass concentration', unit='kg/kg')
      end block
    end if

    ! Sea salt
    if (lss) then
      block
        integer :: my_modes(2) = [iACS, iCOS]
        real(field_r) :: rho = 2165
        do imod = 1, size(my_modes)
          call mode_add_aerosol(modes(my_modes(imod)), iss, rho)
        end do
        call add_tracer('ss_c', long_name='sea salt in-cloud mass concentration', unit='kg/kg')
        call add_tracer('ss_r', long_name='sea salt in-rain mass concentration', unit='kg/kg')
      end block
    end if

    ! Particulate organic matter
    if (lpom) then
      block
        integer :: my_modes(4) = [iAIS, iACS, iCOS, iAII]
        real(field_r) :: rho = 1800
        do imod = 1, size(my_modes)
          call mode_add_aerosol(modes(my_modes(imod)), ipom, rho)
        end do
        call add_tracer('pom_c', long_name='organic matter in-cloud mass concentration', unit='kg/kg')
        call add_tracer('pom_r', long_name='organic matter in-rain mass concentration', unit='kg/kg')
      end block
    end if

    ! Black carbon
    if (lbc) then
      block
        integer :: my_modes(4) = [iAIS, iACS, iCOS, iAII]
        real(field_r) :: rho = 1300
        do imod = 1, size(my_modes)
          call mode_add_aerosol(modes(my_modes(imod)), ibc, rho)
        end do
        call add_tracer('bc_c', long_name='black carbon in-cloud mass concentration', unit='kg/kg')
        call add_tracer('bc_r', long_name='black carbon in-rain mass concentration', unit='kg/kg')
      end block
    end if

    ! Mineral dust
    if (ldu) then
      block
        integer :: my_modes(4) = [iACS, iCOS, iACI, iCOI]
        real(field_r) :: rho = 2560
        do imod = 1, size(my_modes)
          call mode_add_aerosol(modes(my_modes(imod)), idu, rho)
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

  subroutine mode_add_aerosol(self, itype, rho)
    class(mode_t), intent(inout) :: self
    integer, intent(in) :: itype
    real(field_r), intent(in) :: rho

    integer :: isv

    self%lactive = .true.

    self%nspecies = self%nspecies + 1
    self%itype(self%nspecies) = itype

    ! Setup a tracer for the mass concentration
    call add_tracer(trim(aerosol_names(itype))//"_"//self%name, &
      long_name=aerosol_longnames(itype), isv=isv)

    self%itrac(self%nspecies) = isv
    self%rho(self%nspecies) = rho

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
      sv_idx = get_tracer_index(trim(aerosol_names(self%itype(iaer)))//'_'//self%name)
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
      sv_idx = get_tracer_index(trim(aerosol_names(self%itype(iaer)))//'_'//self%name)
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


  !> \brief Aerosol activation based on updraft velocity
  !!
  !! \see https://doi.org/10.5194/acp-15-9217-2015
  !!
  !! \param m_ais Aitken soluble mode.
  !! \param m_acs Accumulation soluble mode.
  !! \param m_cos Coarse soluble mode.
  !! \param m_inc In-cloud mode.
  !! \param w Vertical velocity.
  !! \param delt Time step size.
  subroutine activation(ql, w, nc, delt, ncp)

    real(field_r), intent(in) :: &
      ql(2-ih:,2-jh:,:),         &
      w(2-ih:,2-jh:,:),          &
      nc(2:,2:,:),               &
      delt
    real(field_r), intent(inout) :: &
      ncp(2:,2:,:)

    character(*),  parameter :: routine = modname//"::activation_pn15"

    real(field_r), parameter :: r_crit = 35E-9 ! Critical radius for activation.

    integer :: &
      i, j, k, s, & ! Loop indices.
      my_target     ! Target index in in-cloud mode.

    real(field_r) :: &
      dm,            & ! Median diameter.
      n_act,         & ! Number of activated particles.
      dncdt,         & ! Potential activation tendency.
      fn,            & ! Activated fraction of number concentration.
      fm,            & ! Activated fraciton of mass concentrattion.
      tend_n,        & ! Real tendency of number concentration.
      tend_m,        & ! Real tendency of mass concentration.
      w0               ! Updraft velocity.

    call timer_tic(routine, 1)

    associate(m_cos => modes(iCOS), m_acs => modes(iACS), m_ais => modes(iAIS))

    do k = 1, kmax
      do j = 2, j1
        do i = 2, i1
          if (ql(i,j,k) > qcmin) then
            if (m_ais%lactive) then
              dm = calc_median_diameter(m_ais%n(i,j,k), m_ais%q(i,j,k,:), &
                                        m_ais%rho, m_ais%sig_g)
              if (dm > 0) then
              fn = 1 - 0.5_field_r * erfc(-log(2 * r_crit / &
                   dm + eps) / (sqrt(2.0_field_r) * log(m_ais%sig_g)))
            end if
            end if

            n_act = 1E-6 * (m_acs%n(i,j,k) + m_cos%n(i,j,k) + &
                            fn * m_ais%n(i,j,k))

            w0 = max(0.0_field_r, w(i,j,k))
            dncdt = 1E6 / delt * &
                      (0.1 * (w0 * 100 * n_act / &
                              (w0 * 100 + 0.023_field_r * n_act + eps)))**1.27_field_r &
              - 1E-6 * Nc(i,j,k)
            dncdt = max(dncdt, 0.0_field_r)

            fn = dncdt * delt / (m_cos%n(i,j,k)+ eps)
              fn = max(min(fn, 1.0_field_r), 0.0_field_r)
            if (fn >= 1.0_field_r) then
              fm = 1.0_field_r
            elseif (fn <= 0.0_field_r) then
              fm = 0.0_field_r
            else
              fm = 1 - 0.5_field_r * erfc(erfcinv(2 * fn) &
                - 3 * log(m_cos%sig_g) / sqrt(2.0_field_r))
            end if

            tend_n = fn * m_cos%n(i,j,k) / delt
              tend_n = max(0.0_field_r, tend_n)
              m_cos%np(i,j,k) = m_cos%np(i,j,k) - tend_n
              Ncp(i,j,k) = Ncp(i,j,k) + tend_n

              do s = 1, m_cos % nspecies
              tend_m = max(0.0_field_r, fm * m_cos%q(i,j,k,s) / delt)
                my_target = inc_idx(m_cos%itype(s))
                m_cos%qp(i,j,k,s) = m_cos%qp(i,j,k,s) - tend_m
                qap_inc(i,j,k,my_target) = qap_inc(i,j,k,my_target) + tend_m
              end do

            dncdt = dncdt - tend_n
            dncdt = max(0.0_field_r, dncdt)

            if (m_acs%lactive .and. dncdt > 0) then
              fn = dncdt * delt / (m_acs%n(i,j,k) + eps)
              fn = max(min(fn, 1.0_field_r), 0.0_field_r)
              if (fn >= 1.0_field_r) then
                fm = 1.0_field_r
              elseif (fn <= 0.0_field_r) then
                fm = 0.0_field_r
              else
                fm = 1 - 0.5_field_r * erfc(erfcinv(2 * fn) &
                  - 3 * log(m_acs%sig_g) / sqrt(2.0_field_r))
              end if

              tend_n = fn * m_acs%n(i,j,k) / delt
              tend_n = max(0.0_field_r, tend_n)

              m_acs%np(i,j,k) = m_acs%np(i,j,k) - tend_n
              Ncp(i,j,k) = Ncp(i,j,k) + tend_n

              do s = 1, m_acs % nspecies
                tend_m = fm * m_acs%q(i,j,k,s) / delt
                tend_m = max(0.0_field_r, tend_m)
                my_target = inc_idx(m_acs%itype(s))
                m_acs%qp(i,j,k,s) = m_acs%qp(i,j,k,s) - tend_m
                qap_inc(i,j,k,my_target) = qap_inc(i,j,k,my_target) + tend_m
              end do

              dncdt = dncdt - tend_n
              dncdt = max(0.0_field_r, dncdt)
            end if

            if (m_ais%lactive .and. dncdt > 0) then
              fn = dncdt * delt / (m_ais%n(i,j,k) + eps)
              fn = max(min(fn, 1.0_field_r), 0.0_field_r)
              fm = 1 - 0.5_field_r * &
                erfc(erfcinv(2 * fn) - 3 * log(m_ais%sig_g) / sqrt(2.0_field_r))
              fm = merge(1.0_field_r, fm, fn > 1.0_field_r)

              tend_n = fn * m_ais%n(i,j,k) / delt
              tend_n = max(0.0_field_r, tend_n)

              m_ais%np(i,j,k) = m_ais%np(i,j,k) - tend_n
              Ncp(i,j,k) = Ncp(i,j,k) + tend_n

              do s = 1, m_ais%nspecies
                tend_m = fm * m_ais%q(i,j,k,s) / delt
                tend_m = max(0.0_field_r, tend_m)
                my_target = inc_idx(m_ais%itype(s))
                m_ais%qp(i,j,k,s) = m_ais%qp(i,j,k,s) - tend_m
                qap_inc(i,j,k,my_target) = qap_inc(i,j,k,my_target) + tend_m
              end do
            end if
          end if
        end do
      end do
    end do

    end associate

    call timer_toc(routine)

  end subroutine activation

  !> Move aerosol from in-cloud to in-rain mode based on the tendency of some
  !! rain generation process.
  !!
  !! @param[in] qc Cloud water content.
  !! @param[in] qrp Tendency of rain water content.
  subroutine aerosol_cloud_to_rain(qc, qrp)

    real(field_r), intent(in) :: &
      qc(2-ih:,2-jh:,:),         &
      qrp(2:,2:,:)

    character(len=*), parameter :: routine = modname//'/aerosol_cloud_to_rain'

    integer :: &
      i, j, k, s ! Loop indices

    real(field_r) :: &
      dqadt ! Tendency of in-rain aerosol

    call timer_tic(routine, 1)

    !$acc parallel loop collapse(4) default(present) private(dqadt)
    do s = 1, n_species_active
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1
            if (qrp(i,j,k) > 0) then
            dqadt = qrp(i,j,k) / qc(i,j,k) * qa_inc(i,j,k,s)
            qap_inc(i,j,k,s) = qap_inc(i,j,k,s) - dqadt
            qap_inr(i,j,k,s) = qap_inr(i,j,k,s) + dqadt
            end if
          end do
        end do
      end do
    end do

    call timer_toc(routine)

  end subroutine aerosol_cloud_to_rain

  !> Resuspend aerosols from evaporated rain drops.
  !!
  !! Aerosols are resuspended over the ACS and COS modes. Only rain drops that 
  !! fully evaporate should resuspend an aerosol particle, so the resuspended
  !! aerosol mass is corrected using a correction factor.
  !!
  !! @see https://doi.org/10.1016/j.atmosres.2005.10.012
  !!
  !! @param[in] qr Rain water content.
  !! @param[in] qrp Tendency of rain water content from evaporation.
  !! @param[in] nrp Tendency of rain number concentration from evaporation.
  !! @param[in] delt Time step size.
  subroutine aerosol_resuspend_rain(qr, qrp, nrp, delt)

    real(field_r), intent(in) :: &
      qr(2:,2:,:),               &
      qrp(2:,2:,:),              &
      nrp(2:,2:,:),              &
      delt

    character(len=*), parameter :: routine = modname//'/aero_redistribute'
    real(field_r),    parameter :: Dc = 1E-9 ! Median diameter of resuspended aerosol.

    integer :: &
      i, j, k, s, & ! Loop indices.
      itype,      & ! Aerosol type.
      target_idx    ! Index in target mode.

    real(field_r) :: &
      f_evp,             & ! Fraction of evaporated rain water.
      eps,               & ! Correction factor.
      evapm(maxspecies), & ! Mass of evaporated aerosol.
      evapn,             & ! Number of evaporated aerosol.
      dn,                & ! Number median diameter of evaporated aerosol.
      dm,                & ! Mass median diameter of evaporated aerosol.
      fn,                & ! Number fraction of aerosol resuspended in ACS mode.
      fm                   ! Mass fraction of aerosol resuspended in COS mode.

    call timer_tic(routine, 1)

    associate(m_acs => modes(iACS), m_cos => modes(iCOS))

    do k = 1, kmax
      do j = 2, j1
        do i = 2, i1
          if (qr(i,j,k) > 0) then
            f_evp = (-qrp(i,j,k) * delt) / (qr(i,j,k) + 1E-40)
          f_evp = max(min(f_evp, 1.0_field_r), 0.0_field_r)

            ! Correction factor from Gong et al. (2006).
            ! Evaluates to 1 for f_evp = 1
          eps = (1 - exp(-2 * sqrt(f_evp)) * (1 + 2 * sqrt(f_evp) &
            + 2 * f_evp + (4.0_field_r/3) * f_evp**(3.0_field_r/2))) &
            * (1 - f_evp) + f_evp * f_evp

            evapm(:) = eps * f_evp * qa_inr(i,j,k,:) / delt
            evapn = -1 * nrp(i,j,k)

            ! Compute the median diameter of the resuspended aerosol.
            dn = calc_median_diameter(evapn, evapm(:), rho, 1.5_field_r)
            dm = dn * exp(3 * log(1.5_field_r)**2)

            fn = 0.5_field_r * erfc(-log(dc/(dn + 1E-40)) &
                                    / (log(1.5_field_r) * sqrt(2.0_field_r)))
            fm = 0.5_field_r * erfc(-log(dc/(dm + 1E-40)) &
                                    / (log(1.5_field_r) * sqrt(2.0_field_r)))

            m_acs%np(i,j,k) = m_acs%np(i,j,k) + fn * evapn
            m_cos%np(i,j,k) = m_cos%np(i,j,k) + (1 - fn) * evapn

            do s = 1, n_species_active
              itype = aerosol_get_type_in_cloud(s)
              target_idx = aerosol_get_index_in_mode(itype, m_acs)
              m_acs%qp(i,j,k,target_idx) = m_acs%qp(i,j,k,target_idx) &
                                           + fm * evapm(s)
              target_idx = aerosol_get_index_in_mode(itype, m_cos)
              m_cos%qp(i,j,k,target_idx) = m_cos%qp(i,j,k,target_idx) &
                                           + (1 - fm) * evapm(s)
              qap_inr(i,j,k,s) = qap_inr(i,j,k,s) - evapm(s)
            end do
          end if
        end do
      end do
          end do

    end associate

    call timer_toc(routine)

  end subroutine aerosol_resuspend_rain

          dn = 1E6 * (6 * m_evp / (pi * n_evp * rho_evp))**(1.0_field_r / 3) &
            * exp(-(3.0_field_r / 2) * log(1.5_field_r)**2)
          dm = dn * exp(3 * log(1.5_field_r)**2)

          fn = 0.5_field_r * erfc(-log(dc/dn) / log(1.5_field_r) / sqrt(2.0_field_r))
          fm = 0.5_field_r * erfc(-log(dc/dm) / log(1.5_field_r) / sqrt(2.0_field_r))

          m_acs%np(i,j,k) = m_acs%np(i,j,k) + Fn * n_evp
          m_cos%np(i,j,k) = m_cos%np(i,j,k) + (1 - fn) * n_evp

          do s = 1, n_species_active
            evapt = eps * f_evp * qa_inr(i,j,k,s) / delt
            itype = aerosol_get_type_in_cloud(s)
            target_idx = aerosol_get_index_in_mode(itype, m_acs)
            m_acs%qp(i,j,k,target_idx) = m_acs%qp(i,j,k,target_idx) + fm * evapt
            target_idx = aerosol_get_index_in_mode(itype, m_cos)
            m_cos%qp(i,j,k,target_idx) = m_cos%qp(i,j,k,target_idx) + (1 - fm) * evapt
          end do
        end do
      end do
    end do

    end associate

  end subroutine aerosol_redistribute

  subroutine scavenging_cloud(ql, Nc, rhof)

    real(field_r), intent(in) :: &
      ql(:,:,:) !> Test
    real(field_r), intent(in) :: Nc(:,:,:)
    real(field_r), intent(in) :: rhof(:)
       
    integer               :: &
      i, j, k, m, s,         & !< Loop indices
      target_idx               !< temp
    real(field_r)         :: &
      m_mass,                & !< Mode mean mass.
      m_dens,                & !< Mode mean density.
      dia_a,                 & !< Mean diameter of aerosol.
      dia_c,                 & !< Mean diameter of cloud droplets.
      tend_m,                & !< Tendency of mass concentration.
      tend_n,                & !< Tendency of number concentration.
      fs_m,                  & !< Fraction of mass that is washed out.
      fs_n                     !< Fraction of number that is washed out.
    type(mode_t), pointer :: &
      mode                     !< Pointer to current mode.

    do m = 1, size(modes)
      mode => modes(m)
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1
            if (ql(i,j,k) > qcmin) then
              ! Mode mean properties
              m_mass = 0
              m_dens = 0

              do s = 1, mode%nspecies
                m_mass = m_mass + mode%q(i,j,k,s)
                m_dens = m_dens + (mode%q(i,j,k,s) / mode%rho(s))
              end do

              m_dens = m_mass / (m_dens + eps)

              ! Compute mean cloud droplet diameter and rain rate.
              ! Make sure both stay within the bounds of the lookup table
              dia_c = max( &
                1E6_field_r * (3 * ql(i,j,k) * rhof(k) / &
                  (4 * pi * nc(i,j,k) * rhow + eps)), &
                5.001_field_r &
                )
              dia_c = min(dia_c, 49.999_field_r)

              ! Compute mean aerosol radius in this mode
              dia_a = 0.5 * (6 * m_mass / (pi * mode%n(i,j,k) * m_dens + eps)) &
                      **(1.0_field_r / 3) * exp((-3 * log(mode%sig_g)**2) / 2)
              dia_a = min(100 * dia_a, 8E-3_field_r)
              dia_a = max(dia_a, 1E-8_field_r)

              ! Compute how much aerosol is washed out (number and mass)
              fs_m = LT2_get_col_inline(inc_tab_m, 1, log(dia_c), log(dia_a))
              fs_m = 1E-6 * Nc(i,j,k) * fs_m

              fs_n = LT2_get_col_inline(inc_tab_n, 1, log(dia_c), log(dia_a))
              fs_n = 1E-6 * Nc(i,j,k) * fs_n

              fs_m = merge(1 / delt, fs_m, fs_m * delt > 1 .or. fs_n * delt > 1)
              fs_n = merge(1 / delt, fs_n, fs_m * delt > 1 .or. fs_n * delt > 1)

              ! Remove aerosol from the free modes
              tend_n = fs_n * max(0.0_field_r, mode%n(i,j,k))
              mode%np(i,j,k) = mode%np(i,j,k) - tend_n

              do s = 1, mode%nspecies
                tend_m  = fs_m * max(0.0_field_r, mode%q(i,j,k,s))
                mode%qp(i,j,k,s) = mode%qp(i,j,k,s) - tend_m
                target_idx = inc_idx(mode%itype(s))
                qap_inc(i,j,k,target_idx) = qap_inc(i,j,k,target_idx) + tend_m
              end do
            end if
          end do
        end do
      end do
    end do

  end subroutine scavenging_cloud

  subroutine scavenging_rain(qr, sed_qr)

    real(field_r), intent(in) :: qr(:,:,:)     !< Rain water mixing ratio.
    real(field_r), intent(in) :: sed_qr(:,:,:) !< Sedimentation rate.

    integer :: &
      i, j, k, m, s, & !< Loop indices.
      target_idx 
    real(field_r) :: &
      m_mass, &
      m_dens, &
      rainrate, &
      radi_a, &
      fs_m, &
      fs_n, &
      tend_m, &
      tend_n
    type(mode_t), pointer :: &
      mode

    do m = 1, size(modes)
      mode => modes(m)
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1 
            if (qr(i,j,k) > 0 .and. sed_qr(i,j,k)*3600 > 0.01_field_r) then
              ! Mode mean properties
              m_mass = 0
              m_dens = 0

              do s = 1, mode%nspecies
                m_mass = m_mass + mode%q(i,j,k,s)
                m_dens = m_dens + (mode%q(i,j,k,s) / mode%rho(s))
              end do

              m_dens = m_mass / (m_dens + eps)

              rainrate = sed_qr(i,j,k) * 3600
              rainrate = max(rainrate, 0.01001_field_r)
              rainrate = min(rainrate, 99.999_field_r)

              ! Compute mean aerosol radius in this mode
              radi_a = 0.5 * (6 * m_mass / &
                (pi * mode%n(i,j,k)* m_dens + eps)) &
                **(1.0_field_r / 3) &
                * exp((-3 * log(mode%sig_g)**2) / 2)

              radi_a = (3 * m_mass / (pi * mode%n(i,j,k) * m_dens + eps))**(1.0_field_r / 3) &
                       * exp((-3 * log(mode%sig_g)**2) / 2)

              radi_a = min(0.9999E3_field_r, radi_a * 1E6_field_r)
              radi_a = max(radi_a, 1.001E-3_field_r)

              ! Compute how much aerosol is washed out (number and mass)
              fs_m = LT2_get_col_inline(blc_tab_m, 1, log(rainrate), log(radi_a))
              fs_m = merge(1 / delt, fs_m, fs_m * delt > 1 .or. fs_n * delt > 1)

              fs_n = LT2_get_col_inline(blc_tab_n, 1, log(rainrate), log(radi_a))
              fs_n = merge(1 / delt, fs_n, fs_m * delt > 1 .or. fs_n * delt > 1)

              ! Remove aerosol from the free modes
              tend_n = fs_n * max(0.0_field_r, mode%n(i,j,k))
              mode%np(i,j,k) = mode%np(i,j,k) - tend_n

              do s = 1, mode%nspecies
                tend_m  = fs_m * max(0.0_field_r, mode%q(i,j,k,s))
                target_idx = inc_idx(mode%itype(s))
                mode%qp(i,j,k,s) = mode%qp(i,j,k,s) - tend_m
                qap_inr(i,j,k,target_idx) = qap_inr(i,j,k,target_idx) + tend_m
              end do
            end if
          end do
        end do
      end do
    end do

  end subroutine scavenging_rain

end module modaerosol