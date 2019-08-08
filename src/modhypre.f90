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
public :: inithypre, solve_hypre, exithypre
save

  integer   mpi_comm_hypre

  integer   ierr
  integer*8 grid, stencil, solver
  integer*8 matrixA, vectorX, vectorB ! Solve Ax = b
  integer   ilower(3), iupper(3), periodic(3)
  integer   maxiter, n_pre, n_post
  real      tol

  integer   solver_id, zero, one

  data solver_id  / 0 /
  data zero / 0 /
  data one  / 1 /

contains
  subroutine inithypre
    use mpi
    use modmpi, only : myid, myidx, myidy, nprocx, nprocy

    use modglobal, only : imax, jmax, kmax, dzf, dzh, dx, dy, &
                          itot, jtot

    use modfields, only : rhobf, rhobh

    implicit none

    ! NOTE: we assume kmax is larger than 7,
    ! when we make the stencil entries
    real      values(kmax*imax*jmax)
    real      cx, cy, cz_down, cz_up, cc

    integer   offsets(3,7), stencil_indices(7)
    integer   i, j, k

    ! Have hypre reuse the comm world
    mpi_comm_hypre = MPI_COMM_WORLD

    !-----------------------------------------------------------------------
    !     1. Set up the grid.
    !        Each processor describes the piece of the grid that it owns.
    !-----------------------------------------------------------------------

    ! Set up the grid
    call HYPRE_StructGridCreate(mpi_comm_hypre, 3, grid, ierr)

    ilower(1) = myidx * imax
    ilower(2) = myidy * jmax
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

    write (*,*) 'HYPRE Setup grid: myid, ilower, iupper', myid, ilower, iupper

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
    !        Each processor describes the piece of the grid that it owns.
    !-----------------------------------------------------------------------

    stencil_indices(1) = 0
    stencil_indices(2) = 1
    stencil_indices(3) = 2
    stencil_indices(4) = 3
    stencil_indices(5) = 4
    stencil_indices(6) = 5
    stencil_indices(7) = 6

    ! enter stencil values per horizontal slab
    do k = 1,kmax
      ilower(3) = k - 1
      iupper(3) = k - 1

      cx = rhobf(k)/(dx*dx)
      cy = rhobf(k)/(dy*dy)

      ! variable:   corresponds to:
      !  cz_down       a(k)
      !  cz_up         c(k)
      !  cc            b(k)
      if (k == 1) then
        cz_down = 0.0 ! a(1) = 0
        cz_up =   rhobh(k+1)/(dzf(k) * dzh(k+1))
        cc = -(cz_up + 2.0 * (cx + cy))
      else if (k == kmax) then
        cz_down = rhobh(k  )/(dzf(k) * dzh(k  ))
        cz_up =   .0  ! c(kmax) = 0
      else
        cz_down = rhobh(k  )/(dzf(k) * dzh(k  ))
        cz_up =   rhobh(k+1)/(dzf(k) * dzh(k+1))
      endif
      cc = -(cz_down + cz_up + 2.0 * (cx + cy))

      do i = 1, imax*jmax*7, 7
        values(i+0) = cx ! west
        values(i+1) = cx ! east
        values(i+2) = cy ! south
        values(i+3) = cy ! north
        values(i+4) = cz_down ! below
        values(i+5) = cz_up ! above
        values(i+6) = cc ! center
      enddo

      call HYPRE_StructMatrixSetBoxValues(matrixA, &
           ilower, iupper, 7, stencil_indices, values, ierr)
    enddo

    call HYPRE_StructMatrixAssemble(matrixA, ierr)
    call HYPRE_StructMatrixPrint(matrixA, zero, ierr)

    !-----------------------------------------------------------------------
    !     5. Set up the rhs and initial guess
    !        Each processor describes the piece of the grid that it owns.
    !-----------------------------------------------------------------------

    call HYPRE_StructVectorCreate(mpi_comm_hypre, grid, vectorB, ierr)
    call HYPRE_StructVectorCreate(mpi_comm_hypre, grid, vectorX, ierr)

    call HYPRE_StructVectorInitialize(vectorB, ierr)
    call HYPRE_StructVectorInitialize(vectorX, ierr)

    ! initialize some values as starting point for the iterative solver
    do i=1,imax*jmax
        values(i) = 1e-5
    enddo
    do k=1,kmax
      ilower(3) = k - 1
      iupper(3) = k - 1
      call HYPRE_StructVectorSetBoxValues(vectorX, ilower, iupper, values, ierr)
    enddo
    call HYPRE_StructVectorAssemble(vectorX, ierr)

    !-----------------------------------------------------------------------
    !     5. Choose a solver and initialize it
    !-----------------------------------------------------------------------

    if (solver_id .eq. 0) then
      ! Solve the system using SMG
      maxiter = 50
      tol = 1e-9
      n_pre = 1
      n_post = 1

      call HYPRE_StructSMGCreate(mpi_comm_hypre, solver, ierr)
      call HYPRE_StructSMGSetMemoryUse(solver, zero, ierr)
      call HYPRE_StructSMGSetMaxIter(solver, maxiter, ierr)
      call HYPRE_StructSMGSetTol(solver, tol, ierr)
      call HYPRE_StructSMGSetRelChange(solver, zero, ierr)
      call HYPRE_StructSMGSetNumPreRelax(solver, n_pre, ierr)
      call HYPRE_StructSMGSetNumPostRelax(solver, n_post, ierr)
      call HYPRE_StructSMGSetPrintLevel(solver, one, ierr)
      call HYPRE_StructSMGSetLogging(solver, one, ierr)
      call HYPRE_StructSMGSetup(solver, matrixA, vectorB, vectorX, ierr)
      write (*,*) 'Selected solver 1 (SMG) with parameters:', maxiter, tol, n_pre, n_post
    else
      write (*,*) 'Invalid solver in inithypre', solver
      call exit(-1)
    endif

    write (*,*) 'inithypre - done'
  end subroutine

  subroutine solve_hypre
    use modmpi, only : myid
    use modglobal, only : i1, j1, ih, jh, imax, jmax, kmax
    use modpois, only : p

    implicit none

    ! real, intent(inout) :: p(2-ih:i1+ih,2-jh:j1+jh,kmax)
    real values(imax,jmax)

    real final_res_norm
    integer i,j,k, num_iterations
    real totalsum

    !-----------------------------------------------------------------------
    !     1. Set up the rhs
    !-----------------------------------------------------------------------
    do k=1,kmax
      do j=1,jmax
        do i=1,imax
          values(i,j) = p(i+1,j+1,k)
        enddo
      enddo

      ilower(3) = k - 1
      iupper(3) = k - 1
      call HYPRE_StructVectorSetBoxValues(vectorB, ilower, iupper, values, ierr)
    enddo
    call HYPRE_StructVectorAssemble(vectorB, ierr)

    ! use current values (ie. the solution to the previous call) as starting point
    !do j=1,jmax
    !  do i=1,imax
    !    values(i,j) = 1e-5
    !  enddo
    !enddo
    !do k=1,kmax
    !  ilower(3) = k - 1
    !  iupper(3) = k - 1
    !  call HYPRE_StructVectorSetBoxValues(vectorX, ilower, iupper, values, ierr)
    !enddo
    !call HYPRE_StructVectorAssemble(vectorX, ierr)

    !-----------------------------------------------------------------------
    !     2. Call a solver
    !-----------------------------------------------------------------------

    if (solver_id .eq. 0) then
      ! Solve the system using SMG
      call HYPRE_StructSMGSolve(solver, matrixA, vectorB, vectorX, ierr)
      if (myid == 0) then
        write (*,*) 'HYPRE solver status (ierr)', ierr
      endif

      call HYPRE_StructSMGGetNumIterations(solver, num_iterations, ierr)
      call HYPRE_StructSMGGetFinalRelative(solver, final_res_norm, ierr)

      if (myid == 0) then
        write (*,*) 'HYPRE Number of iterations / max iteration', num_iterations, maxiter
        write (*,*) 'HYPRE Final residual norm / target', final_res_norm, tol
      endif
    endif

    !-----------------------------------------------------------------------
    !     3. Copy solution
    !-----------------------------------------------------------------------

    do k=1,kmax
      ilower(3) = k - 1
      iupper(3) = k - 1
      call HYPRE_StructVectorGetBoxValues(vectorX, ilower, iupper, values, ierr)

      do j=1,jmax
        do i=1,imax
          p(i+1,j+1,k) = values(i,j)
        enddo
      enddo
    enddo
  end subroutine

  subroutine exithypre
    call HYPRE_StructSMGDestroy(solver, ierr)
    call HYPRE_StructGridDestroy(grid, ierr)
    call HYPRE_StructStencilDestroy(stencil, ierr)
    call HYPRE_StructMatrixDestroy(matrixA, ierr)
    call HYPRE_StructVectorDestroy(vectorB, ierr)
    call HYPRE_StructVectorDestroy(vectorX, ierr)
  end subroutine

end module
