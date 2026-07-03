!> Lateral sponge layer: smooth nudging of scalar fields at the boundaries.
module modlatsponge

  use modfields,       only: svp, sv0
  use modglobal,       only: nsv, itot, jtot, i1, j1, kmax, pi, rdt, ifnamopt, &
                             checknamelisterror
  use modmpi,          only: myid, myidx, myidy, nprocx, nprocy, d_mpi_bcast, &
                             comm3d
  use modprecision,    only: field_r
  use modtimer,        only: timer_tic, timer_toc
  use fortran_support, only: nnml_output

  implicit none

  private

  public :: lateral_sponge_read_namelist
  public :: lateral_sponge

  character(len=*), parameter :: modname = 'modlatsponge'

  logical :: llateral_sponge = .false. !< Whether to apply the lateral sponge layer
  integer :: nudgedepth = 10           !< Number of nudge grid points

  !$acc declare create(nudgedepth)
!$omp declare target (nudgedepth)

contains

  !> Read the namelist for the lateral sponge layer.
  subroutine lateral_sponge_read_namelist(nml_filename)

    character(len=*), intent(in) :: nml_filename

    integer :: ierr

    namelist /lateral_sponge/ llateral_sponge, nudgedepth

    if (myid == 0) then
      open(ifnamopt, file=nml_filename, status='old', iostat=ierr)
      read(ifnamopt, lateral_sponge, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'lateral_sponge')
      write(nnml_output, lateral_sponge)
      close(ifnamopt)
    end if

    call d_mpi_bcast(llateral_sponge, 1, 0, comm3d, ierr)
    call d_mpi_bcast(nudgedepth, 1, 0, comm3d, ierr)

    !$acc update device(nudgedepth)
!!$omp target update to(nudgedepth)

  end subroutine lateral_sponge_read_namelist

  !> Compute the nudging factor.
  pure function nudge_fac(i) result(fac)

    integer, intent(in) :: i !< Index of grid point in the nudging layer (1 at the boundary, increasing towards the interior)

    real(field_r) :: fac !< Nudging factor

    fac = 0.5 + 0.5 * cos((pi / (nudgedepth - 1)) * (i - 1))

  end function nudge_fac

  !> Nudge scalar fields towards zero at the lateral boundaries.
  subroutine lateral_sponge()

    character(len=*), parameter :: routine = modname//'/lateral_sponge'

    integer :: i, j, k, s

    if (.not. llateral_sponge) return

    call timer_tic(routine, 0)

    ! North
    if (myidy == 0) then
      !$acc parallel loop gang vector collapse(4) default(present) async
!!$omp target teams loop collapse(4) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
      do s = 1, nsv
        do k = 1, kmax
          do j = 1, nudgedepth
            do i = 2, i1
              svp(i,j,k,s) = svp(i,j,k,s) &
                             + nudge_fac(j) * (0 - sv0(i,j,k,s)) / rdt
            end do
          end do
        end do
      end do
    end if

    ! South
    if (myidy == nprocy - 1) then
      !$acc parallel loop gang vector collapse(4) default(present) async
!!$omp target teams loop collapse(4) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
      do s = 1, nsv
        do k = 1, kmax
          do j = j1 - nudgedepth + 1, j1
            do i = 2, i1
              svp(i,j,k,s) = svp(i,j,k,s) &
                             + nudge_fac(j1 - j + 1) * (0 - sv0(i,j,k,s)) / rdt
            end do
          end do
        end do
      end do
    end if

    ! East
    if (myidx == nprocx - 1) then
      !$acc parallel loop gang vector collapse(4) default(present) async
!!$omp target teams loop collapse(4) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
      do s = 1, nsv
        do k = 1, kmax
          do j = 2, j1
            do i = i1 - nudgedepth + 1, i1
              svp(i,j,k,s) = svp(i,j,k,s) &
                             + nudge_fac(i1 - i + 1) * (0 - sv0(i,j,k,s)) / rdt
            end do
          end do
        end do
      end do
    end if

    ! West
    if (myidx == 0) then
      !$acc parallel loop gang vector collapse(4) default(present) async
!!$omp target teams loop collapse(4) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
      do s = 1, nsv
        do k = 1, kmax
          do j = 2, j1
            do i = 1, nudgedepth
              svp(i,j,k,s) = svp(i,j,k,s) &
                             + nudge_fac(i) * (0 - sv0(i,j,k,s)) / rdt
            end do
          end do
        end do
      end do
    end if

    !$acc wait

    call timer_toc(routine)

  end subroutine lateral_sponge

end module modlatsponge
