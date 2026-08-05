!> Temporary LCM validation prints for DALES-LCM smoke experiments.
#ifdef LCM_VALIDATION_CHECKS
module modlcm_validation_debug
  use iso_fortran_env, only : int64, real64
  use lcm_host_interface, only : lcm_config_t
  use lcm_initialization_state, only : lcm_initial_particle_count_local, &
                                       lcm_initial_particle_count_global
  use lcm_particle_state, only : cell_particles
  use modglobal, only : imax, jmax, kmax, dx, dy, zh
  use modmpi, only : myidx, myidy

  implicit none

  private

  public :: print_lcm_validation_checks

contains

  subroutine print_lcm_validation_checks(lcm_config)
    type(lcm_config_t), intent(in) :: lcm_config

    integer :: initialized_levels
    integer(int64) :: counted_particles
    logical :: cells_are_allocated

    initialized_levels = min(lcm_config%nzl, kmax)
    cells_are_allocated = allocated(cell_particles)

    write(*,*) 'lcm_validation_check particle_cells_allocated=', &
               cells_are_allocated
    write(*,*) 'lcm_validation_check initial_count_nonzero=', &
               lcm_initial_particle_count_local > 0_int64
    write(*,*) 'lcm_validation_check count_local=', &
               lcm_initial_particle_count_local
    write(*,*) 'lcm_validation_check count_global=', &
               lcm_initial_particle_count_global
    write(*,*) 'lcm_validation_check initialized_levels=', &
               initialized_levels

    if (.not. cells_are_allocated) return

    counted_particles = count_lcm_initialized_particles()
    write(*,*) 'lcm_validation_check counted_initialized_particles=', &
               counted_particles
    write(*,*) 'lcm_validation_check particles_only_at_k_le_nzl=', &
               particles_are_only_at_initialized_levels(initialized_levels)

    call print_lcm_validation_cell(representative_cell=1, i=1, j=1, k=1)
    call print_lcm_validation_cell(representative_cell=2, &
                                   i=max(1, imax / 2), &
                                   j=max(1, jmax / 2), &
                                   k=max(1, initialized_levels))
    call print_lcm_validation_cell(representative_cell=3, &
                                   i=imax, j=jmax, &
                                   k=max(1, initialized_levels))
  end subroutine print_lcm_validation_checks

  integer(int64) function count_lcm_initialized_particles() result(counted)
    integer :: i, j, k

    counted = 0_int64
    do k = 1, kmax
      do j = 1, jmax
        do i = 1, imax
          counted = counted + int(cell_particles(i,j,k)%count, int64)
        end do
      end do
    end do
  end function count_lcm_initialized_particles

  logical function particles_are_only_at_initialized_levels(initialized_levels) &
    result(is_valid)
    integer, intent(in) :: initialized_levels

    integer :: i, j, k

    is_valid = .true.
    do k = initialized_levels + 1, kmax
      do j = 1, jmax
        do i = 1, imax
          if (cell_particles(i,j,k)%count /= 0) is_valid = .false.
        end do
      end do
    end do
  end function particles_are_only_at_initialized_levels

  subroutine print_lcm_validation_cell(representative_cell, i, j, k)
    integer, intent(in) :: representative_cell
    integer, intent(in) :: i
    integer, intent(in) :: j
    integer, intent(in) :: k

    integer :: n
    logical :: position_inside_cell
    logical :: physical_weight
    logical :: physical_mass_aero
    logical :: physical_mass_water
    logical :: physical_temperature
    logical :: physical_qv

    write(*,*) 'lcm_validation_check representative_cell=', representative_cell, &
               ' i=', i, ' j=', j, ' k=', k, &
               ' computational_particles=', cell_particles(i,j,k)%count

    if (cell_particles(i,j,k)%count <= 0) return

    n = 1
    position_inside_cell = particle_position_is_inside_cell(i, j, k, n)
    physical_weight = cell_particles(i,j,k)%particles(n)%weight > 0_int64
    physical_mass_aero = cell_particles(i,j,k)%particles(n)%mass_aero > 0.0_real64
    physical_mass_water = cell_particles(i,j,k)%particles(n)%mass_water >= 0.0_real64
    physical_temperature = is_physical_lcm_temperature( &
      cell_particles(i,j,k)%particles(n)%temperature)
    physical_qv = is_physical_lcm_qv(cell_particles(i,j,k)%particles(n)%qv)

    write(*,*) 'lcm_validation_check representative_particle=', &
               representative_cell, &
               ' n=', n, &
               ' position_inside_cell=', position_inside_cell, &
               ' weight_valid=', physical_weight, &
               ' mass_aero_valid=', physical_mass_aero, &
               ' mass_water_valid=', physical_mass_water, &
               ' temperature_valid=', physical_temperature, &
               ' qv_valid=', physical_qv
    write(*,*) 'lcm_validation_check representative_particle_values=', &
               representative_cell, &
               ' x=', cell_particles(i,j,k)%particles(n)%x, &
               ' y=', cell_particles(i,j,k)%particles(n)%y, &
               ' z=', cell_particles(i,j,k)%particles(n)%z, &
               ' weight=', cell_particles(i,j,k)%particles(n)%weight, &
               ' mass_aero=', cell_particles(i,j,k)%particles(n)%mass_aero, &
               ' mass_water=', cell_particles(i,j,k)%particles(n)%mass_water, &
               ' temperature=', cell_particles(i,j,k)%particles(n)%temperature, &
               ' qv=', cell_particles(i,j,k)%particles(n)%qv
  end subroutine print_lcm_validation_cell

  logical function particle_position_is_inside_cell(i, j, k, n) &
    result(is_inside)
    integer, intent(in) :: i
    integer, intent(in) :: j
    integer, intent(in) :: k
    integer, intent(in) :: n

    real(real64) :: x_lower, x_upper
    real(real64) :: y_lower, y_upper
    real(real64) :: z_lower, z_upper

    x_lower = real(myidx * imax + i - 1, real64) * real(dx, real64)
    x_upper = real(myidx * imax + i, real64) * real(dx, real64)
    y_lower = real(myidy * jmax + j - 1, real64) * real(dy, real64)
    y_upper = real(myidy * jmax + j, real64) * real(dy, real64)
    z_lower = real(zh(k), real64)
    z_upper = real(zh(k + 1), real64)

    is_inside = &
      cell_particles(i,j,k)%particles(n)%x >= x_lower .and. &
      cell_particles(i,j,k)%particles(n)%x < x_upper .and. &
      cell_particles(i,j,k)%particles(n)%y >= y_lower .and. &
      cell_particles(i,j,k)%particles(n)%y < y_upper .and. &
      cell_particles(i,j,k)%particles(n)%z >= z_lower .and. &
      cell_particles(i,j,k)%particles(n)%z < z_upper
  end function particle_position_is_inside_cell

  logical function is_physical_lcm_temperature(value) result(is_physical)
    real(real64), intent(in) :: value

    is_physical = value == value .and. value > 150.0_real64 .and. &
                  value < 350.0_real64
  end function is_physical_lcm_temperature

  logical function is_physical_lcm_qv(value) result(is_physical)
    real(real64), intent(in) :: value

    is_physical = value == value .and. value >= 0.0_real64 .and. &
                  value < 0.1_real64
  end function is_physical_lcm_qv

end module modlcm_validation_debug
#endif
