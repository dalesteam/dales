! This file is part of DALES.
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
! Copyright 1993-2025, The DALES Team.
!

!> Routines for computing (masked) slab averaged profiles.
module modslabaverage

  use, intrinsic :: iso_fortran_env, only: real32, real64

  use modglobal, only: imax, jmax, ijtot
  use modmpi,    only: mpi_allreduce, mpi_in_place, mpi_real4, mpi_real8, mpi_integer, &
                       mpi_sum, comm3d, mpierr

  implicit none

  private

  character(len=*), parameter :: modname = 'modslabaverage'

  public :: slabavg

  interface slabavg
    module procedure slabavg_r4
    module procedure slabavg_r8
    module procedure slabavg_r4_masked
    module procedure slabavg_r8_masked
  end interface slabavg

contains

  !> Compute slab averaged profile of a field.
  !!
  !! \param field Variable to compute slab average of.
  !! \param nh Number of ghost cells.
  !!           \note For arrays with no ghost cells, pass nh=0!
  !! \param avg Slab averaged profile.
  !! \param local Only compute average on local MPI domain.
  subroutine slabavg_r4(field, nh, avg, local)

    real(real32), intent(in)  :: field(:,:,:)
    integer,      intent(in)  :: nh
    real(real32), intent(out) :: avg(:)

    logical, optional, intent(in) :: local

    integer      :: i, j, k
    integer      :: is, ie, js, je, ks, ke
    logical      :: do_global
    real(real32) :: fld_sum, norm_fac

    is = lbound(field, dim=1) + nh
    ie = ubound(field, dim=1) - nh
    js = lbound(field, dim=2) + nh
    je = ubound(field, dim=2) - nh
    ks = lbound(field, dim=3)
    ke = ubound(field, dim=3)

    if (present(local)) then
      do_global = .not. local
    else
      do_global = .true.
    end if

    if (do_global) then
      norm_fac = 1.0_real32 / ijtot
    else
      norm_fac = 1.0_real32 / (imax * jmax)
    end if

    !$acc parallel loop gang default(present)
    do k = ks, ke
      fld_sum = 0
      !$acc loop vector collapse(2) reduction(+: fld_sum)
      do j = js, je
        do i = is, ie
          fld_sum = fld_sum + field(i,j,k)
        end do
      end do
      avg(k) = fld_sum * norm_fac
    end do

    ! TODO: experiment with non-blocking allreduce
    if (do_global) then
      !$acc host_data use_device(avg)
      call mpi_allreduce(mpi_in_place, avg, ke, mpi_real4, mpi_sum, &
                         comm3d, mpierr)
      !$acc end host_data
    end if

  end subroutine slabavg_r4

  !> Compute slab averaged profile of a field.
  !!
  !! \param field Variable to compute slab average of.
  !! \param nh Number of ghost cells.
  !!           \note For arrays with no ghost cells, pass nh=0!
  !! \param avg Slab averaged profile.
  !! \param local Only compute average on local MPI domain.
  subroutine slabavg_r8(field, nh, avg, local)

    real(real64), intent(in)  :: field(:,:,:)
    integer,      intent(in)  :: nh
    real(real64), intent(out) :: avg(:)

    logical, optional, intent(in) :: local

    integer      :: i, j, k
    integer      :: is, ie, js, je, ks, ke
    logical      :: do_global
    real(real64) :: fld_sum, norm_fac

    is = lbound(field, dim=1) + nh
    ie = ubound(field, dim=1) - nh
    js = lbound(field, dim=2) + nh
    je = ubound(field, dim=2) - nh
    ks = lbound(field, dim=3)
    ke = ubound(field, dim=3)

    if (present(local)) then
      do_global = .not. local
    else
      do_global = .true.
    end if

    if (do_global) then
      norm_fac = 1.0_real64 / ijtot
    else
      norm_fac = 1.0_real64 / (imax * jmax)
    end if

    !$acc parallel loop gang default(present)
    do k = ks, ke
      fld_sum = 0
      !$acc loop vector collapse(2) reduction(+: fld_sum)
      do j = js, je
        do i = is, ie
          fld_sum = fld_sum + field(i,j,k)
        end do
      end do
      avg(k) = fld_sum * norm_fac
    end do

    if (do_global) then
      !$acc host_data use_device(avg)
      call mpi_allreduce(mpi_in_place, avg, ke, mpi_real8, mpi_sum, &
                         comm3d, mpierr)
      !$acc end host_data
    end if

  end subroutine slabavg_r8

  !> Compute masked slab averaged profile of a field.
  !!
  !! \param field Variable to compute slab average of.
  !! \param mask Mask to apply. Values marked as .TRUE. will be included in the average.
  !! \param nh Number of ghost cells.
  !!           \note For arrays with no ghost cells, pass nh=0!
  !! \param avg Slab averaged profile.
  !! \param local Only compute average on local MPI domain.
  !! \param fillvalue Value to use for fully masked layers. Default is 0.
  subroutine slabavg_r4_masked(field, mask, nh, avg, local, fillvalue)

    real(real32), intent(in)  :: field(:,:,:)
    logical,      intent(in)  :: mask(:,:,:)
    integer,      intent(in)  :: nh
    real(real32), intent(out) :: avg(:)

    logical,      optional, intent(in) :: local
    real(real32), optional, intent(in) :: fillvalue

    integer      :: i, j, k
    integer      :: is, ie, js, je, ks, ke
    integer      :: test, n_cells
    logical      :: do_global
    real(real32) :: fld_sum, fillvalue_

    is = lbound(field, dim=1) + nh
    ie = ubound(field, dim=1) - nh
    js = lbound(field, dim=2) + nh
    je = ubound(field, dim=2) - nh
    ks = lbound(field, dim=3)
    ke = ubound(field, dim=3)

    if (present(local)) then
      do_global = .not. local
    else
      do_global = .true.
    end if

    fillvalue_ = 0
    if (present(fillvalue)) fillvalue_ = fillvalue

    ! Use a block to place n_cells_tot on the stack
    block
      integer :: n_cells_tot(ks:ke)

      !$acc data create(n_cells_tot)

      !$acc parallel loop gang default(present)
      do k = ks, ke
        fld_sum = 0
        n_cells = 0
        !$acc loop vector collapse(2) reduction(+: fld_sum, n_cells)
        do j = js, je
          do i = is, ie
            test = merge(1, 0, mask(i,j,k))
            fld_sum = fld_sum + test * field(i,j,k)
            n_cells = n_cells + test
          end do
        end do
        avg(k) = fld_sum
        n_cells_tot(k) = n_cells
      end do

      if (do_global) then
        !$acc host_data use_device(avg, n_cells_tot)
        call mpi_allreduce(mpi_in_place, avg, ke, mpi_real4, mpi_sum, comm3d, &
                           mpierr)
        call mpi_allreduce(mpi_in_place, n_cells_tot, ke, mpi_integer, &
                           mpi_sum, comm3d, mpierr)
        !$acc end host_data
      end if

      !$acc parallel loop gang default(present)
      do k = ks, ke
        avg(k) = merge(avg(k) / n_cells_tot(k), fillvalue_, n_cells_tot(k) > 0)
      end do

      !$acc end data

    end block

  end subroutine slabavg_r4_masked

  !> Compute masked slab averaged profile of a field.
  !!
  !! \param field Variable to compute slab average of.
  !! \param mask Mask to apply. Values marked as .TRUE. will be included in the average.
  !! \param nh Number of ghost cells.
  !!           \note For arrays with no ghost cells, pass nh=0!
  !! \param avg Slab averaged profile.
  !! \param local Only compute average on local MPI domain.
  !! \param fillvalue Value to use for fully masked layers. Default is 0.
  subroutine slabavg_r8_masked(field, mask, nh, avg, local, fillvalue)

    real(real64), intent(in)  :: field(:,:,:)
    logical,      intent(in)  :: mask(:,:,:)
    integer,      intent(in)  :: nh
    real(real64), intent(out) :: avg(:)

    logical,      optional, intent(in) :: local
    real(real64), optional, intent(in) :: fillvalue

    integer      :: i, j, k
    integer      :: is, ie, js, je, ks, ke
    integer      :: test, n_cells
    logical      :: do_global
    real(real64) :: fld_sum, fillvalue_

    is = lbound(field, dim=1) + nh
    ie = ubound(field, dim=1) - nh
    js = lbound(field, dim=2) + nh
    je = ubound(field, dim=2) - nh
    ks = lbound(field, dim=3)
    ke = ubound(field, dim=3)

    if (present(local)) then
      do_global = .not. local
    else
      do_global = .true.
    end if

    fillvalue_ = 0
    if (present(fillvalue)) fillvalue_ = fillvalue

    ! Use a block to place n_cells_tot on the stack
    block
      integer :: n_cells_tot(ks:ke)

      !$acc data create(n_cells_tot)

      !$acc parallel loop gang default(present)
      do k = ks, ke
        fld_sum = 0
        n_cells = 0
        !$acc loop vector collapse(2) reduction(+: fld_sum, n_cells)
        do j = js, je
          do i = is, ie
            test = merge(1, 0, mask(i,j,k))
            fld_sum = fld_sum + test * field(i,j,k)
            n_cells = n_cells + test
          end do
        end do
        avg(k) = fld_sum
        n_cells_tot(k) = n_cells
      end do

      if (do_global) then
        !$acc host_data use_device(avg, n_cells_tot)
        call mpi_allreduce(mpi_in_place, avg, ke, mpi_real8, mpi_sum, comm3d, &
                           mpierr)
        call mpi_allreduce(mpi_in_place, n_cells_tot, ke, mpi_integer, &
                           mpi_sum, comm3d, mpierr)
        !$acc end host_data
      end if

      !$acc parallel loop gang default(present)
      do k = ks, ke
        avg(k) = merge(avg(k) / n_cells_tot(k), fillvalue_, n_cells_tot(k) > 0)
      end do

      !$acc end data

    end block

  end subroutine slabavg_r8_masked

end module modslabaverage
