!> DALES NAMLCM handling for standalone LCM configuration.
#ifdef USE_LCM
module modlcm_namelist
  use iso_fortran_env, only : real32
  use lcm_host_interface, only : lcm_config_t, lcm_init_config
  use modglobal, only : ifnamopt, checknamelisterror
  use modmpi, only : myid, D_MPI_BCAST, comm3d

  implicit none

  private

  public :: lcm_read_namelist
  public :: lcm_apply_namelist_config

  logical :: lcm_namelist_is_set = .false.

  ! DALES-side storage for NAMLCM. The names intentionally match the
  ! SAM-LCM MICRO_LAGRANGE namelist and are mapped to lcm_config_t.
  integer :: np_xy, nzl, nstep_lcm_init, nstep_cut_supersat
  integer :: n_uvw_iteration, start_sedimentation, start_collection
  integer :: breakup_conserve_mode, n_start_source_aero, n_end_source_aero
  integer :: n_start_source_in, n_end_source_in, start_ice
  integer :: spec_nsave, spec_nsavestart, n_track_particles
  integer :: particle_ntrack, particle_ntrackstart, particle_ntrackstop
  integer :: nz_b_source_aero(4), nz_t_source_aero(4)
  integer :: ny_s_source_aero(4), ny_n_source_aero(4)
  integer :: nx_w_source_aero(4), nx_e_source_aero(4)
  integer :: nz_b_source_in(4), nz_t_source_in(4)
  integer :: ny_s_source_in(4), ny_n_source_in(4)
  integer :: nx_w_source_in(4), nx_e_source_in(4), type_source_in(4)
  logical :: switch_sgs_velocities, switch_sedimentation, switch_collection
  logical :: switch_linear_sampling_collection, switch_breakup
  logical :: switch_analytical_micro, switch_radiation_micro
  logical :: switch_lem_micro, switch_koehler_micro, switch_ffc_micro
  logical :: switch_steady_pichamber, switch_cloud_tracking
  logical :: switch_phase_sources, switch_particle_source
  logical :: switch_walls_pichamber, switch_reflective_walls_pichamber
  logical :: switch_ice, switch_habit_ice, switch_icebreak, switch_rs
  logical :: switch_ds, switch_moving_source, switch_diskice_sizethresh
  logical :: switch_steady_aerosol, switch_spec_aero
  logical :: track_some_particles, horiz_average_spec
  real(real32) :: r_separate_micro, supersat_threshold, alpha_spec
  real(real32) :: r_start_spec, r_end_spec, frac_ice, salinity_source
  real(real32) :: sizethresh_chen, sizethresh_diskice, inject_bottom, inject_top
  real(real32) :: n_aero(4), rm_aero(4), sigma_aero(4), shares_aero(4)
  real(real32) :: n_ice(4), rm_in(4), sigma_in(4), rho_aero(4)
  real(real32) :: molecular_weight_aero(4), vanthoff_aero(4)
  real(real32) :: n_source_aero(4), n_source_in(4)
  real(real32) :: r_source_aero(4), r_source_in(4)
  real(real32) :: sigma_source_aero(4), sigma_source_in(4)
  real(real32) :: p_source_aero(4), p_source_in(4)
  real(real32) :: frac_in_source_particles(4), frac_in_source_weight(4)
  real(real32) :: vx_source_aero(4), vy_source_aero(4), vz_source_aero(4)
  real(real32) :: vx_source_in(4), vy_source_in(4), vz_source_in(4)
  character(len=6) :: init_aero
  character(len=12) :: unit_aero, unit_source, type_source
  character(len=10) :: habitparam, ds_param
  character(len=8) :: icebreak_param

contains

  subroutine lcm_read_namelist(nml_filename)
    use fortran_support, only : nnml_output

    character(len=*), intent(in) :: nml_filename
    integer :: ierr

    namelist /NAMLCM/                                                     &
      np_xy, nzl, nstep_lcm_init, n_uvw_iteration,                        &
      switch_sgs_velocities, r_separate_micro, n_aero, rm_aero,           &
      shares_aero, sigma_aero, n_ice, rm_in, sigma_in,                    &
      switch_sedimentation, spec_nsave, spec_nsavestart,                  &
      switch_lem_micro, switch_collection, start_collection,               &
      start_sedimentation, init_aero, switch_radiation_micro,             &
      switch_analytical_micro, switch_steady_pichamber,                   &
      switch_cloud_tracking, switch_phase_sources, switch_particle_source,&
      nz_b_source_aero, nz_t_source_aero, ny_s_source_aero,               &
      ny_n_source_aero, nx_w_source_aero, nx_e_source_aero,               &
      nz_b_source_in, nz_t_source_in, ny_s_source_in, ny_n_source_in,     &
      nx_w_source_in, nx_e_source_in, n_start_source_aero,                &
      n_end_source_aero, n_source_aero, n_start_source_in,                &
      n_end_source_in, n_source_in, r_source_aero, r_source_in,           &
      alpha_spec, r_start_spec, switch_walls_pichamber,                   &
      switch_reflective_walls_pichamber, p_source_aero, p_source_in,      &
      switch_koehler_micro, switch_linear_sampling_collection,            &
      sigma_source_aero, sigma_source_in, type_source_in,                 &
      switch_ffc_micro, horiz_average_spec, frac_ice,                     &
      frac_in_source_particles, frac_in_source_weight, switch_ice,        &
      switch_habit_ice, habitparam, start_ice, switch_icebreak,           &
      switch_rs, switch_ds, ds_param, icebreak_param, unit_source,        &
      unit_aero, switch_moving_source, vx_source_aero, vy_source_aero,    &
      vz_source_aero, vx_source_in, vy_source_in, vz_source_in,           &
      salinity_source, type_source, r_end_spec, switch_diskice_sizethresh,&
      sizethresh_chen, sizethresh_diskice, switch_breakup,                &
      breakup_conserve_mode, switch_steady_aerosol, inject_bottom,        &
      inject_top, switch_spec_aero, track_some_particles,                 &
      n_track_particles, particle_ntrack, particle_ntrackstart, rho_aero, &
      vanthoff_aero, molecular_weight_aero, particle_ntrackstop,          &
      nstep_cut_supersat, supersat_threshold

    call initialize_lcm_namelist_defaults()

    if (myid == 0) then
      open(ifnamopt, file=nml_filename, status='old', iostat=ierr)
      read(ifnamopt, NAMLCM, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMLCM')
      write(nnml_output, NAMLCM)
      close(ifnamopt)
    end if

    call broadcast_lcm_namelist(ierr)
    lcm_namelist_is_set = .true.
  end subroutine lcm_read_namelist

  subroutine lcm_apply_namelist_config(config)
    type(lcm_config_t), intent(inout) :: config

    if (.not. lcm_namelist_is_set) call initialize_lcm_namelist_defaults()
    call apply_lcm_namelist_to_config(config)
  end subroutine lcm_apply_namelist_config

  subroutine initialize_lcm_namelist_defaults()
    type(lcm_config_t) :: default_config

    call lcm_init_config(default_config)

    np_xy = default_config%np_xy
    nzl = default_config%nzl
    nstep_lcm_init = default_config%nstep_lcm_init
    nstep_cut_supersat = default_config%nstep_cut_supersat
    n_uvw_iteration = default_config%n_uvw_iteration
    start_sedimentation = default_config%start_sedimentation
    start_collection = default_config%start_collection
    breakup_conserve_mode = default_config%breakup_conserve_mode
    n_start_source_aero = default_config%n_start_source_aero
    n_end_source_aero = default_config%n_end_source_aero
    n_start_source_in = default_config%n_start_source_in
    n_end_source_in = default_config%n_end_source_in
    start_ice = default_config%start_ice
    spec_nsave = default_config%spec_nsave
    spec_nsavestart = default_config%spec_nsavestart
    n_track_particles = default_config%n_track_particles
    particle_ntrack = default_config%particle_ntrack
    particle_ntrackstart = default_config%particle_ntrackstart
    particle_ntrackstop = default_config%particle_ntrackstop
    nz_b_source_aero = default_config%nz_b_source_aero
    nz_t_source_aero = default_config%nz_t_source_aero
    ny_s_source_aero = default_config%ny_s_source_aero
    ny_n_source_aero = default_config%ny_n_source_aero
    nx_w_source_aero = default_config%nx_w_source_aero
    nx_e_source_aero = default_config%nx_e_source_aero
    nz_b_source_in = default_config%nz_b_source_in
    nz_t_source_in = default_config%nz_t_source_in
    ny_s_source_in = default_config%ny_s_source_in
    ny_n_source_in = default_config%ny_n_source_in
    nx_w_source_in = default_config%nx_w_source_in
    nx_e_source_in = default_config%nx_e_source_in
    type_source_in = default_config%type_source_in
    switch_sgs_velocities = default_config%switch_sgs_velocities
    switch_sedimentation = default_config%switch_sedimentation
    switch_collection = default_config%switch_collection
    switch_linear_sampling_collection = default_config%switch_linear_sampling_collection
    switch_breakup = default_config%switch_breakup
    switch_analytical_micro = default_config%switch_analytical_micro
    switch_radiation_micro = default_config%switch_radiation_micro
    switch_lem_micro = default_config%switch_lem_micro
    switch_koehler_micro = default_config%switch_koehler_micro
    switch_ffc_micro = default_config%switch_ffc_micro
    switch_steady_pichamber = default_config%switch_steady_pichamber
    switch_cloud_tracking = default_config%switch_cloud_tracking
    switch_phase_sources = default_config%switch_phase_sources
    switch_particle_source = default_config%switch_particle_source
    switch_walls_pichamber = default_config%switch_walls_pichamber
    switch_reflective_walls_pichamber = default_config%switch_reflective_walls_pichamber
    switch_ice = default_config%switch_ice
    switch_habit_ice = default_config%switch_habit_ice
    switch_icebreak = default_config%switch_icebreak
    switch_rs = default_config%switch_rs
    switch_ds = default_config%switch_ds
    switch_moving_source = default_config%switch_moving_source
    switch_diskice_sizethresh = default_config%switch_diskice_sizethresh
    switch_steady_aerosol = default_config%switch_steady_aerosol
    switch_spec_aero = default_config%switch_spec_aero
    track_some_particles = default_config%track_some_particles
    horiz_average_spec = default_config%horiz_average_spec
    r_separate_micro = default_config%r_separate_micro
    supersat_threshold = default_config%supersat_threshold
    alpha_spec = default_config%alpha_spec
    r_start_spec = default_config%r_start_spec
    r_end_spec = default_config%r_end_spec
    frac_ice = default_config%frac_ice
    salinity_source = default_config%salinity_source
    sizethresh_chen = default_config%sizethresh_chen
    sizethresh_diskice = default_config%sizethresh_diskice
    inject_bottom = default_config%inject_bottom
    inject_top = default_config%inject_top
    n_aero = default_config%n_aero
    rm_aero = default_config%rm_aero
    sigma_aero = default_config%sigma_aero
    shares_aero = default_config%shares_aero
    n_ice = default_config%n_ice
    rm_in = default_config%rm_in
    sigma_in = default_config%sigma_in
    rho_aero = default_config%rho_aero
    molecular_weight_aero = default_config%molecular_weight_aero
    vanthoff_aero = default_config%vanthoff_aero
    n_source_aero = default_config%n_source_aero
    n_source_in = default_config%n_source_in
    r_source_aero = default_config%r_source_aero
    r_source_in = default_config%r_source_in
    sigma_source_aero = default_config%sigma_source_aero
    sigma_source_in = default_config%sigma_source_in
    p_source_aero = default_config%p_source_aero
    p_source_in = default_config%p_source_in
    frac_in_source_particles = default_config%frac_in_source_particles
    frac_in_source_weight = default_config%frac_in_source_weight
    vx_source_aero = default_config%vx_source_aero
    vy_source_aero = default_config%vy_source_aero
    vz_source_aero = default_config%vz_source_aero
    vx_source_in = default_config%vx_source_in
    vy_source_in = default_config%vy_source_in
    vz_source_in = default_config%vz_source_in
    init_aero = default_config%init_aero
    unit_aero = default_config%unit_aero
    unit_source = default_config%unit_source
    type_source = default_config%type_source
    habitparam = default_config%habitparam
    ds_param = default_config%ds_param
    icebreak_param = default_config%icebreak_param
  end subroutine initialize_lcm_namelist_defaults

  subroutine apply_lcm_namelist_to_config(config)
    type(lcm_config_t), intent(inout) :: config

    config%np_xy = np_xy
    config%nzl = nzl
    config%nstep_lcm_init = nstep_lcm_init
    config%nstep_cut_supersat = nstep_cut_supersat
    config%n_uvw_iteration = n_uvw_iteration
    config%start_sedimentation = start_sedimentation
    config%start_collection = start_collection
    config%breakup_conserve_mode = breakup_conserve_mode
    config%n_start_source_aero = n_start_source_aero
    config%n_end_source_aero = n_end_source_aero
    config%n_start_source_in = n_start_source_in
    config%n_end_source_in = n_end_source_in
    config%start_ice = start_ice
    config%spec_nsave = spec_nsave
    config%spec_nsavestart = spec_nsavestart
    config%n_track_particles = n_track_particles
    config%particle_ntrack = particle_ntrack
    config%particle_ntrackstart = particle_ntrackstart
    config%particle_ntrackstop = particle_ntrackstop
    config%nz_b_source_aero = nz_b_source_aero
    config%nz_t_source_aero = nz_t_source_aero
    config%ny_s_source_aero = ny_s_source_aero
    config%ny_n_source_aero = ny_n_source_aero
    config%nx_w_source_aero = nx_w_source_aero
    config%nx_e_source_aero = nx_e_source_aero
    config%nz_b_source_in = nz_b_source_in
    config%nz_t_source_in = nz_t_source_in
    config%ny_s_source_in = ny_s_source_in
    config%ny_n_source_in = ny_n_source_in
    config%nx_w_source_in = nx_w_source_in
    config%nx_e_source_in = nx_e_source_in
    config%type_source_in = type_source_in
    config%switch_sgs_velocities = switch_sgs_velocities
    config%switch_sedimentation = switch_sedimentation
    config%switch_collection = switch_collection
    config%switch_linear_sampling_collection = switch_linear_sampling_collection
    config%switch_breakup = switch_breakup
    config%switch_analytical_micro = switch_analytical_micro
    config%switch_radiation_micro = switch_radiation_micro
    config%switch_lem_micro = switch_lem_micro
    config%switch_koehler_micro = switch_koehler_micro
    config%switch_ffc_micro = switch_ffc_micro
    config%switch_steady_pichamber = switch_steady_pichamber
    config%switch_cloud_tracking = switch_cloud_tracking
    config%switch_phase_sources = switch_phase_sources
    config%switch_particle_source = switch_particle_source
    config%switch_walls_pichamber = switch_walls_pichamber
    config%switch_reflective_walls_pichamber = switch_reflective_walls_pichamber
    config%switch_ice = switch_ice
    config%switch_habit_ice = switch_habit_ice
    config%switch_icebreak = switch_icebreak
    config%switch_rs = switch_rs
    config%switch_ds = switch_ds
    config%switch_moving_source = switch_moving_source
    config%switch_diskice_sizethresh = switch_diskice_sizethresh
    config%switch_steady_aerosol = switch_steady_aerosol
    config%switch_spec_aero = switch_spec_aero
    config%track_some_particles = track_some_particles
    config%horiz_average_spec = horiz_average_spec
    config%r_separate_micro = r_separate_micro
    config%supersat_threshold = supersat_threshold
    config%alpha_spec = alpha_spec
    config%r_start_spec = r_start_spec
    config%r_end_spec = r_end_spec
    config%frac_ice = frac_ice
    config%salinity_source = salinity_source
    config%sizethresh_chen = sizethresh_chen
    config%sizethresh_diskice = sizethresh_diskice
    config%inject_bottom = inject_bottom
    config%inject_top = inject_top
    config%n_aero = n_aero
    config%rm_aero = rm_aero
    config%sigma_aero = sigma_aero
    config%shares_aero = shares_aero
    config%n_ice = n_ice
    config%rm_in = rm_in
    config%sigma_in = sigma_in
    config%rho_aero = rho_aero
    config%molecular_weight_aero = molecular_weight_aero
    config%vanthoff_aero = vanthoff_aero
    config%n_source_aero = n_source_aero
    config%n_source_in = n_source_in
    config%r_source_aero = r_source_aero
    config%r_source_in = r_source_in
    config%sigma_source_aero = sigma_source_aero
    config%sigma_source_in = sigma_source_in
    config%p_source_aero = p_source_aero
    config%p_source_in = p_source_in
    config%frac_in_source_particles = frac_in_source_particles
    config%frac_in_source_weight = frac_in_source_weight
    config%vx_source_aero = vx_source_aero
    config%vy_source_aero = vy_source_aero
    config%vz_source_aero = vz_source_aero
    config%vx_source_in = vx_source_in
    config%vy_source_in = vy_source_in
    config%vz_source_in = vz_source_in
    config%init_aero = init_aero
    config%unit_aero = unit_aero
    config%unit_source = unit_source
    config%type_source = type_source
    config%habitparam = habitparam
    config%ds_param = ds_param
    config%icebreak_param = icebreak_param
  end subroutine apply_lcm_namelist_to_config

  subroutine broadcast_lcm_namelist(ierr)
    integer, intent(out) :: ierr

    call D_MPI_BCAST(np_xy, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(nzl, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(nstep_lcm_init, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(nstep_cut_supersat, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(n_uvw_iteration, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(start_sedimentation, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(start_collection, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(breakup_conserve_mode, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(n_start_source_aero, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(n_end_source_aero, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(n_start_source_in, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(n_end_source_in, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(start_ice, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(spec_nsave, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(spec_nsavestart, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(n_track_particles, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(particle_ntrack, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(particle_ntrackstart, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(particle_ntrackstop, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(nz_b_source_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(nz_t_source_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(ny_s_source_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(ny_n_source_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(nx_w_source_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(nx_e_source_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(nz_b_source_in, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(nz_t_source_in, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(ny_s_source_in, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(ny_n_source_in, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(nx_w_source_in, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(nx_e_source_in, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(type_source_in, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_sgs_velocities, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_sedimentation, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_collection, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_linear_sampling_collection, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_breakup, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_analytical_micro, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_radiation_micro, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_lem_micro, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_koehler_micro, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_ffc_micro, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_steady_pichamber, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_cloud_tracking, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_phase_sources, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_particle_source, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_walls_pichamber, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_reflective_walls_pichamber, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_ice, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_habit_ice, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_icebreak, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_rs, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_ds, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_moving_source, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_diskice_sizethresh, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_steady_aerosol, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(switch_spec_aero, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(track_some_particles, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(horiz_average_spec, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(r_separate_micro, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(supersat_threshold, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(alpha_spec, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(r_start_spec, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(r_end_spec, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(frac_ice, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(salinity_source, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(sizethresh_chen, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(sizethresh_diskice, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(inject_bottom, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(inject_top, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(n_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(rm_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(sigma_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(shares_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(n_ice, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(rm_in, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(sigma_in, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(rho_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(molecular_weight_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(vanthoff_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(n_source_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(n_source_in, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(r_source_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(r_source_in, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(sigma_source_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(sigma_source_in, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(p_source_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(p_source_in, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(frac_in_source_particles, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(frac_in_source_weight, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(vx_source_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(vy_source_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(vz_source_aero, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(vx_source_in, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(vy_source_in, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(vz_source_in, 4, 0, comm3d, ierr)
    call D_MPI_BCAST(init_aero, len(init_aero), 0, comm3d, ierr)
    call D_MPI_BCAST(unit_aero, len(unit_aero), 0, comm3d, ierr)
    call D_MPI_BCAST(unit_source, len(unit_source), 0, comm3d, ierr)
    call D_MPI_BCAST(type_source, len(type_source), 0, comm3d, ierr)
    call D_MPI_BCAST(habitparam, len(habitparam), 0, comm3d, ierr)
    call D_MPI_BCAST(ds_param, len(ds_param), 0, comm3d, ierr)
    call D_MPI_BCAST(icebreak_param, len(icebreak_param), 0, comm3d, ierr)
  end subroutine broadcast_lcm_namelist

end module modlcm_namelist
#endif
