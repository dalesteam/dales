!> \file modhypre.f90
!!  Layer to deal with the HYPRE library
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
!  Copyright 2019 Netherlands eScience Center
!
!  Based on hypre/src/test/f77_struct.f

module modhypre

implicit none
private
public :: inithypre
save

  integer   mpi_comm_hypre

  integer   ierr
  integer*8 grid
  integer*8 stencil
  integer*8 matrixA

contains
  subroutine inithypre
    use mpi
    use modmpi, only : myidx, myidy, nprocx, nprocy

    use modglobal, only : imax, jmax, kmax, dzf, dx, dy, &
                          itot, jtot


    real      values(7*imax*jmax)
    integer   ilower(3), iupper(3), periodic(3)
    integer   offsets(3,7), stencil_indices(7)
    integer   i, j, k, zero

    data zero / 0 /

    write (*,*) 'inithypre'
    ! Have hypre reuse the comm world
    mpi_comm_hypre = MPI_COMM_WORLD

    !-----------------------------------------------------------------------
    !     1. Set up the grid.  Here we use only one part.  Each processor
    !     describes the piece of the grid that it owns.
    !-----------------------------------------------------------------------

    ! Set up the grid
    call HYPRE_StructGridCreate(mpi_comm_hypre, 3, grid, ierr)

    ilower(1) = myidx * nprocx
    ilower(2) = myidy * nprocy
    ilower(3) = 0

    iupper(1) = ilower(1) + imax - 1
    iupper(2) = ilower(2) + jmax - 1
    iupper(3) = kmax - 1
    call HYPRE_StructGridSetExtents(grid, ilower, iupper, ierr)

    periodic(1) = itot
    periodic(2) = jtot
    periodic(3) = 0
    call HYPRE_StructGridSetPeriodic(grid, periodic, ierr)

    ! This is a collective call finalizing the grid assembly
    call HYPRE_StructGridAssemble(grid, ierr)

    !-----------------------------------------------------------------------
    !     2. Define the discretization stencil
    !
    !                     5| /3
    !                      |/
    !                 0----6----1
    !                     /|
    !                   2/ |4
    !-----------------------------------------------------------------------

    ! Create an empty 2D, 7-pt stencil object
    ! TODO: have a look at 17 or 27 point stencils?

    offsets(1,1) = -1
    offsets(2,1) = 0
    offsets(3,1) = 0

    offsets(1,2) = 1
    offsets(2,2) = 0
    offsets(3,2) = 0

    offsets(1,3) = 0
    offsets(2,3) = -1
    offsets(3,3) = 0

    offsets(1,4) = 0
    offsets(2,4) = 1
    offsets(3,4) = 0

    offsets(1,5) = 0
    offsets(2,5) = 0
    offsets(3,5) = -1

    offsets(1,6) = 0
    offsets(2,6) = 0
    offsets(3,6) = 1

    offsets(1,7) = 0
    offsets(2,7) = 0
    offsets(3,7) = 0

    call HYPRE_StructStencilCreate(3, 7, stencil, ierr)
    do i = 1, 7
      call HYPRE_StructStencilSetElement(stencil, i-1, offsets(1,i), ierr) 
    enddo

    !-----------------------------------------------------------------------
    !     3. Setting up the Struct Matrix
    !-----------------------------------------------------------------------

    call HYPRE_StructMatrixCreate(mpi_comm_hypre, grid, stencil, matrixA, ierr)
    call HYPRE_StructMatrixInitialize(matrixA, ierr)

    !-----------------------------------------------------------------------
    !     4. Set the coefficients for the grid
    !-----------------------------------------------------------------------

    stencil_indices(1) = 0
    stencil_indices(2) = 1
    stencil_indices(3) = 2
    stencil_indices(4) = 3
    stencil_indices(5) = 4
    stencil_indices(6) = 5
    stencil_indices(7) = 6

    ! enter stencil values per horizontal slab

    ! start with the non-boundary points
    do k = 2,kmax-1
      ilower(3) = k - 1 ! Fortran to C indexing
      iupper(3) = k - 1 ! Fortran to C indexing
      do i = 1, imax*jmax*7, 7
        values(i+0) = -1.0 / (dx * dx) ! west
        values(i+1) = -1.0 / (dx * dx) ! east

        values(i+2) = -1.0 / (dy * dy) ! south
        values(i+3) = -1.0 / (dy * dy) ! north 

        ! TODO: grid spacing in the vertical. should this be:
        ! a) the full level thickness at level k
        ! b) the distance to the next center of the cell dzf
        ! let's start with (a) 
        values(i+4) = -1.0  /(dzf(k) * dzf(k)) ! below
        values(i+5) = -1.0  /(dzf(k) * dzf(k)) ! above

        ! NOTE: becuase we sum the values(), we have an extra minus sign
        values(i+6) = -2.0 * (values(i+0) + values(i+1) + values(i+2)) ! center
      enddo

      call HYPRE_StructMatrixSetBoxValues(matrixA, &
           ilower(1), iupper(1), 7, stencil_indices, values, ierr)
    enddo

    ! wipe out stencil values outside of domain
    ! bottom
    ilower(3) = 0
    iupper(3) = 0
    do i = 1, imax*jmax*7, 7
      values(i+0) = -1.0 / (dx * dx) ! west
      values(i+1) = -1.0 / (dx * dx) ! east
      values(i+2) = -1.0 / (dy * dy) ! south
      values(i+3) = -1.0 / (dy * dy) ! north 
      values(i+4) = 0 ! below
      values(i+5) = -1.0  /(dzf(1) * dzf(1)) ! above
      values(i+6) = -4.0 * (values(i+0) + values(i+1) + values(i+2)) ! center
    enddo

    call HYPRE_StructMatrixSetBoxValues(matrixA, &
         ilower(1), iupper(1), 7, stencil_indices, values, ierr)

    ! top
    ilower(3) = kmax - 1
    iupper(3) = kmax - 1
    do i = 1, imax*jmax*7, 7
      values(i+0) = -1.0 / (dx * dx) ! west
      values(i+1) = -1.0 / (dx * dx) ! east
      values(i+2) = -1.0 / (dy * dy) ! south
      values(i+3) = -1.0 / (dy * dy) ! north 
      values(i+4) = -1.0  /(dzf(kmax) * dzf(kmax)) ! below
      values(i+5) = 0 ! above
      values(i+6) = -4.0 * (values(i+0) + values(i+1) + values(i+2)) ! center
    enddo

    call HYPRE_StructMatrixSetBoxValues(matrixA, &
         ilower(1), iupper(1), 7, stencil_indices, values, ierr)


    call HYPRE_StructMatrixAssemble(matrixA, ierr)
    call HYPRE_StructMatrixPrint(matrixA, zero, ierr)
    write (*,*) 'inithypre - done'

  end subroutine

end module
