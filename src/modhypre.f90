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
!  2019 Jisk Attema, Netherlands eScience Center
!
!  Based on hypre/src/test/f77_struct.f

module modhypre
use iso_c_binding
implicit none
private
public :: inithypre, solve_hypre, exithypre, set_initial_guess
save

  integer(c_int) mpi_comm_hypre

  integer(c_int) ierr
  type(c_ptr) grid, stencil, solver, precond_solver
  type(c_ptr) matrixA, vectorX, vectorB ! Solve Ax = b
  integer   ilower(3), iupper(3), periodic(3)

  integer   zero, one
  data zero / 0 /
  data one  / 1 /

! We need to append the _ and add bind(c, name) to call the c interface defined
! for fortran
interface
  subroutine HYPRE_StructBiCGSTABCreate ( comm, solver, ierr) &
    bind(c, name="hypre_structbicgstabcreate_")
    use iso_c_binding
    implicit none
    integer(c_int) comm
    type(c_ptr) :: solver
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructSMGCreate ( comm, solver, ierr) &
    bind(c, name="hypre_structsmgcreate_")
    use iso_c_binding
    implicit none
    integer(c_int) comm
    type(c_ptr) :: solver
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructPFMGCreate ( comm, solver, ierr) &
    bind(c, name="hypre_structpfmgcreate_")
    use iso_c_binding
    implicit none
    integer(c_int) comm
    type(c_ptr) :: solver
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructJacobiCreate ( comm, solver, ierr) &
    bind(c, name="hypre_structjacobicreate_")
    use iso_c_binding
    implicit none
    integer(c_int) comm
    type(c_ptr) :: solver
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructPCGCreate ( comm, solver, ierr) &
    bind(c, name="hypre_structpcgcreate_")
    use iso_c_binding
    implicit none
    integer(c_int) comm
    type(c_ptr) :: solver
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructGMRESCreate ( comm, solver, ierr) &
    bind(c, name="hypre_structgmrescreate_")
    use iso_c_binding
    implicit none
    integer(c_int) comm
    type(c_ptr) :: solver
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructLGMRESCreate ( comm, solver, ierr) &
    bind(c, name="hypre_structlgmrescreate_")
    use iso_c_binding
    implicit none
    integer(c_int) comm
    type(c_ptr) :: solver
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructSMGSetRelChange ( solver, rel_change, ierr) &
    bind(c, name="hypre_structsmgsetrelchange_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) rel_change, ierr
  end subroutine
  subroutine HYPRE_StructPFMGSetRelChange ( solver, rel_change, ierr) &
    bind(c, name="hypre_structpfmgsetrelchange_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) rel_change, ierr
  end subroutine
  subroutine HYPRE_StructBiCGSTABSetTol ( solver, tol, ierr) &
    bind(c, name="hypre_structbicgstabsettol_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr
    real(c_double) tol
  end subroutine
  subroutine HYPRE_StructSMGSetTol ( solver, tol, ierr) &
    bind(c, name="hypre_structsmgsettol_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr
    real(c_double) tol
  end subroutine
  subroutine HYPRE_StructPFMGSetTol ( solver, tol, ierr) &
    bind(c, name="hypre_structpfmgsettol_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr
    real(c_double) tol
  end subroutine
  subroutine HYPRE_StructGMRESSetTol ( solver, tol, ierr) &
    bind(c, name="hypre_structgmressettol_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr
    real(c_double) tol
  end subroutine
  subroutine HYPRE_StructLGMRESSetTol ( solver, tol, ierr) &
    bind(c, name="hypre_structlgmressettol_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr
    real(c_double) tol
  end subroutine
  subroutine HYPRE_StructPCGSetTol ( solver, tol, ierr) &
    bind(c, name="hypre_structpcgsettol_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr
    real(c_double) tol
  end subroutine
  subroutine HYPRE_StructBiCGSTABSetMaxIter ( solver, max_iter, ierr) &
    bind(c, name="hypre_structbicgstabsetmaxiter_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr, max_iter
  end subroutine
  subroutine HYPRE_StructSMGSetMaxIter ( solver, max_iter, ierr) &
    bind(c, name="hypre_structsmgsetmaxiter_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr, max_iter
  end subroutine
  subroutine HYPRE_StructPFMGSetMaxIter ( solver, max_iter, ierr) &
    bind(c, name="hypre_structpfmgsetmaxiter_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr, max_iter
  end subroutine
  subroutine HYPRE_StructJacobiSetMaxIter ( solver, max_iter, ierr) &
    bind(c, name="hypre_structjacobisetmaxiter_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr, max_iter
  end subroutine
  subroutine HYPRE_StructPCGSetMaxIter ( solver, max_iter, ierr) &
    bind(c, name="hypre_structpcgsetmaxiter_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr, max_iter
  end subroutine
  subroutine HYPRE_StructGMRESSetMaxIter ( solver, max_iter, ierr) &
    bind(c, name="hypre_structgmressetmaxiter_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr, max_iter
  end subroutine
  subroutine HYPRE_StructLGMRESSetMaxIter ( solver, max_iter, ierr) &
    bind(c, name="hypre_structlgmressetmaxiter_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr, max_iter
  end subroutine
  subroutine HYPRE_StructPCGSetPrecond ( solver, precond_id, precond_solver, ierr) &
    bind(c, name="hypre_structpcgsetprecond_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver, precond_solver
    integer(c_int) ierr, precond_id
  end subroutine
  subroutine HYPRE_StructGMRESSetPrecond ( solver, precond_id, precond_solver, ierr) &
    bind(c, name="hypre_structgmressetprecond_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver, precond_solver
    integer(c_int) ierr, precond_id
  end subroutine
  subroutine HYPRE_StructBiCGSTABSetPrecond ( solver, precond_id, precond_solver, ierr) &
    bind(c, name="hypre_structbicgstabsetprecond_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver, precond_solver
    integer(c_int) ierr, precond_id
  end subroutine
  subroutine HYPRE_StructLGMRESSetPrecond ( solver, precond_id, precond_solver, ierr) &
    bind(c, name="hypre_structlgmressetprecond_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver, precond_solver
    integer(c_int) ierr, precond_id
  end subroutine
  subroutine HYPRE_StructSMGSetLogging ( solver, logging, ierr) &
    bind(c, name="hypre_structsmgsetlogging_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr, logging
  end subroutine
  subroutine HYPRE_StructPFMGSetLogging ( solver, logging, ierr) &
    bind(c, name="hypre_structpfmgsetlogging_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr, logging
  end subroutine
  subroutine HYPRE_StructPCGSetLogging ( solver, logging, ierr) &
    bind(c, name="hypre_structpcgsetlogging_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr, logging
  end subroutine
  subroutine HYPRE_StructGMRESSetLogging ( solver, logging, ierr) &
    bind(c, name="hypre_structgmressetlogging_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr, logging
  end subroutine
  subroutine HYPRE_StructBiCGSTABSetLogging ( solver, logging, ierr) &
    bind(c, name="hypre_structbicgstabsetlogging_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr, logging
  end subroutine
  subroutine HYPRE_StructLGMRESSetLogging ( solver, logging, ierr) &
    bind(c, name="hypre_structlgmressetlogging_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr, logging
  end subroutine
  subroutine HYPRE_StructPCGSetPrintLevel ( solver, level, ierr) &
    bind(c, name="hypre_structpcgsetprintlevel_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr, level
  end subroutine
  subroutine HYPRE_StructLGMRESSetPrintLevel ( solver, level, ierr) &
    bind(c, name="hypre_structlgmressetprintlevel_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) ierr, level
  end subroutine
  subroutine HYPRE_StructSMGSetup ( solver, matrixA, vectorB, vectorX, ierr) &
    bind(c, name="hypre_structsmgsetup_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver, matrixA, vectorB, vectorX
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructPFMGSetup ( solver, matrixA, vectorB, vectorX, ierr) &
    bind(c, name="hypre_structpfmgsetup_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver, matrixA, vectorB, vectorX
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructPCGSetup ( solver, matrixA, vectorB, vectorX, ierr) &
    bind(c, name="hypre_structpcgsetup_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver, matrixA, vectorB, vectorX
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructGMRESSetup ( solver, matrixA, vectorB, vectorX, ierr) &
    bind(c, name="hypre_structgmressetup_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver, matrixA, vectorB, vectorX
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructBiCGSTABSetup ( solver, matrixA, vectorB, vectorX, ierr) &
    bind(c, name="hypre_structbicgstabsetup_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver, matrixA, vectorB, vectorX
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructLGMRESSetup ( solver, matrixA, vectorB, vectorX, ierr) &
    bind(c, name="hypre_structlgmressetup_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver, matrixA, vectorB, vectorX
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructPCGSetTwoNorm ( solver, two_norm, ierr ) &
    bind(c, name="hypre_structpcgsettwonorm_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) two_norm, ierr
  end subroutine
  subroutine HYPRE_StructSMGSetNumPostRelax ( solver, num_post_relax, ierr) &
    bind(c, name="hypre_structsmgsetnumpostrelax_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) num_post_relax, ierr
  end subroutine
  subroutine HYPRE_StructPFMGSetNumPostRelax ( solver, num_post_relax, ierr) &
    bind(c, name="hypre_structpfmgsetnumpostrelax_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) num_post_relax, ierr
  end subroutine
  subroutine HYPRE_StructSMGSetNumPreRelax ( solver, num_pre_relax, ierr) &
    bind(c, name="hypre_structsmgsetnumprerelax_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) num_pre_relax, ierr
  end subroutine
  subroutine HYPRE_StructPFMGSetNumPreRelax ( solver, num_pre_relax, ierr) &
    bind(c, name="hypre_structpfmgsetnumprerelax_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) num_pre_relax, ierr
  end subroutine
  subroutine HYPRE_StructPFMGSetRelaxType ( solver, relax_type, ierr) &
    bind(c, name="hypre_structpfmgsetrelaxtype_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) relax_type, ierr
  end subroutine
  subroutine HYPRE_StructSMGSetMemoryUse ( solver, memory_use, ierr) &
    bind(c, name="hypre_structsmgsetmemoryuse_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) memory_use, ierr
  end subroutine
  subroutine HYPRE_StructLGMRESGetFinalRel ( solver, norm, ierr) &
    bind(c, name="hypre_structlgmresgetfinalrel_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    real(c_double) :: norm
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructBiCGSTABGetNumItera ( solver, num_iterations, ierr) &
    bind(c, name="hypre_structbicgstabgetnumitera_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) num_iterations, ierr
  end subroutine
  subroutine HYPRE_StructGMRESGetNumIteratio ( solver, num_iterations, ierr) &
    bind(c, name="hypre_structgmresgetnumiteratio_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) num_iterations, ierr
  end subroutine
  subroutine HYPRE_StructPCGGetNumIterations ( solver, num_iterations, ierr) &
    bind(c, name="hypre_structpcggetnumiterations_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) num_iterations, ierr
  end subroutine
  subroutine HYPRE_StructSMGGetNumIterations ( solver, num_iterations, ierr) &
    bind(c, name="hypre_structsmggetnumiterations_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) num_iterations, ierr
  end subroutine
  subroutine HYPRE_StructPFMGGetNumIteration ( solver, num_iterations, ierr) &
    bind(c, name="hypre_structpfmggetnumiteration_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) num_iterations, ierr
  end subroutine
  subroutine HYPRE_StructLGMRESGetNumIter ( solver, num_iterations, ierr) &
    bind(c, name="hypre_structlgmresgetnumiter_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    integer(c_int) num_iterations, ierr
  end subroutine
  subroutine HYPRE_StructPCGSolve ( solver, matrixA, vectorB, vectorX, ierr) &
    bind(c, name="hypre_structpcgsolve_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver, matrixA, vectorB, vectorX
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructBiCGSTABSolve ( solver, matrixA, vectorB, vectorX, ierr) &
    bind(c, name="hypre_structbicgstabsolve_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver, matrixA, vectorB, vectorX
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructSMGSolve ( solver, matrixA, vectorB, vectorX, ierr) &
    bind(c, name="hypre_structsmgsolve_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver, matrixA, vectorB, vectorX
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructPFMGSolve ( solver, matrixA, vectorB, vectorX, ierr) &
    bind(c, name="hypre_structpfmgsolve_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver, matrixA, vectorB, vectorX
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructGMRESSolve ( solver, matrixA, vectorB, vectorX, ierr) &
    bind(c, name="hypre_structgmressolve_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver, matrixA, vectorB, vectorX
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructLGMRESSolve ( solver, matrixA, vectorB, vectorX, ierr) &
    bind(c, name="hypre_structlgmressolve_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver, matrixA, vectorB, vectorX
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructVectorAssemble ( vector, ierr ) &
    bind(c, name="hypre_structvectorassemble_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: vector
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructVectorSetBoxValues ( vector, ilower, iupper, values, ierr ) &
    bind(c, name="hypre_structvectorsetboxvalues_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: vector
    integer(c_int) :: ilower(*), iupper(*)
    real(c_double) :: values(*)
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructVectorGetBoxValues ( vector, ilower, iupper, values, ierr ) &
    bind(c, name="hypre_structvectorgetboxvalues_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: vector
    integer(c_int) :: ilower(*), iupper(*)
    real(c_double) :: values(*)
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructVectorInitialize ( vector, ierr ) &
    bind(c, name="hypre_structvectorinitialize_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: vector
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructVectorCreate ( comm, grid, vector, ierr ) &
    bind(c, name="hypre_structvectorcreate_")
    use iso_c_binding
    implicit none
    integer(c_int) :: comm
    type(c_ptr) :: grid, vector
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructMatrixAssemble ( matrix, ierr ) &
    bind(c, name="hypre_structmatrixassemble_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: matrix
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructMatrixSetBoxValues ( matrix, ilower, iupper, &
      num_stencil_indices, stencil_indices, values, ierr )            &
    bind(c, name="hypre_structmatrixsetboxvalues_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: matrix
    integer(c_int) :: ilower(*), iupper(*), stencil_indices(*)
    real(c_double) :: values(*)
    integer(c_int) num_stencil_indices, ierr
  end subroutine
  subroutine HYPRE_StructMatrixInitialize ( matrix, ierr ) &
    bind(c, name="hypre_structmatrixinitialize_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: matrix
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructMatrixCreate ( comm, grid, stencil, matrix, ierr ) &
    bind(c, name="hypre_structmatrixcreate_")
    use iso_c_binding
    implicit none
    integer(c_int) :: comm
    type(c_ptr) :: grid, stencil, matrix
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructStencilSetElement ( stencil, element_index, offset, ierr ) &
    bind(c, name="hypre_structstencilsetelement_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: stencil
    integer(c_int) :: offset(*)
    integer(c_int) element_index, ierr
  end subroutine
  subroutine HYPRE_StructStencilCreate ( dim, size, stencil, ierr ) &
    bind(c, name="hypre_structstencilcreate_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: stencil
    integer(c_int) dim, size, ierr
  end subroutine
  subroutine HYPRE_StructGridAssemble ( grid, ierr ) &
    bind(c, name="hypre_structgridassemble_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: grid
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructGridSetPeriodic ( grid, periodic, ierr ) &
    bind(c, name="hypre_structgridsetperiodic_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: grid
    integer(c_int) :: periodic(*)
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructGridSetExtents ( grid, ilower, iupper, ierr ) &
    bind(c, name="hypre_structgridsetextents_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: grid
    integer(c_int) :: ilower(*), iupper(*)
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructGridCreate ( comm, dim, grid, ierr ) &
    bind(c, name="hypre_structgridcreate_")
    use iso_c_binding
    implicit none
    integer(c_int) :: comm, dim
    type(c_ptr) :: grid
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructGMRESGetFinalRelati ( solver, norm, ierr ) &
    bind(c, name="hypre_structgmresgetfinalrelati_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    real(c_double) :: norm
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructBiCGSTABGetFinalRel ( solver, norm, ierr ) &
    bind(c, name="hypre_structbicgstabgetfinalrel_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    real(c_double) :: norm
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructPFMGGetFinalRelativ ( solver, norm, ierr ) &
    bind(c, name="hypre_structpfmggetfinalrelativ_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    real(c_double) :: norm
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructPCGGetFinalRelative ( solver, norm, ierr ) &
    bind(c, name="hypre_structpcggetfinalrelative_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    real(c_double) :: norm
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructSMGGetFinalRelative ( solver, norm, ierr ) &
    bind(c, name="hypre_structsmggetfinalrelative_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: solver
    real(c_double) :: norm
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructVectorDestroy ( obj , ierr ) &
    bind(c, name="hypre_structvectordestroy_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: obj
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructMatrixDestroy ( obj , ierr ) &
    bind(c, name="hypre_structmatrixdestroy_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: obj
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructStencilDestroy ( obj , ierr ) &
    bind(c, name="hypre_structstencildestroy_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: obj
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructGridDestroy ( obj , ierr ) &
    bind(c, name="hypre_structgriddestroy_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: obj
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructLGMRESDestroy ( obj , ierr ) &
    bind(c, name="hypre_structlgmresdestroy_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: obj
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructGMRESDestroy ( obj , ierr ) &
    bind(c, name="hypre_structgmresdestroy_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: obj
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructJacobiDestroy ( obj , ierr ) &
    bind(c, name="hypre_structjacobidestroy_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: obj
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructBiCGSTABDestroy ( obj , ierr ) &
    bind(c, name="hypre_structbicgstabdestroy_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: obj
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructPFMGDestroy ( obj , ierr ) &
    bind(c, name="hypre_structpfmgdestroy_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: obj
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructSMGDestroy ( obj , ierr ) &
    bind(c, name="hypre_structsmgdestroy_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: obj
    integer(c_int) ierr
  end subroutine
  subroutine HYPRE_StructPCGDestroy ( obj , ierr ) &
    bind(c, name="hypre_structpcgdestroy_")
    use iso_c_binding
    implicit none
    type(c_ptr) :: obj
    integer(c_int) ierr
  end subroutine

end interface

contains
  subroutine initprecond
    ! Setup a preconditioner for the BiCGSTAB or GMRES solvers
    use modglobal, only : solver_id, precond, maxiter, n_pre, n_post

    implicit none

    ! The precond_id flags mean :
    ! 0 - setup a smg preconditioner
    ! 1 - setup a pfmg preconditioner
    ! 7 - setup a jacobi preconditioner
    ! 8 - setup a ds preconditioner
    ! 9 - dont setup a preconditioner
    if (precond == 0) then
      call HYPRE_StructSMGCreate(mpi_comm_hypre, precond_solver, ierr)
      call HYPRE_StructSMGSetMemoryUse(precond_solver, zero, ierr)
      call HYPRE_StructSMGSetMaxIter(precond_solver, maxiter, ierr)
      call HYPRE_StructSMGSetNumPreRelax(precond_solver, n_pre, ierr)
      call HYPRE_StructSMGSetNumPostRelax(precond_solver, n_post, ierr)
      call HYPRE_StructSMGSetTol(precond_solver, 0.0d0, ierr)
    else if (precond == 1) then
      call HYPRE_StructPFMGCreate(mpi_comm_hypre, precond_solver, ierr)
      call HYPRE_StructPFMGSetMaxIter(precond_solver, maxiter, ierr)
      ! weighted Jacobi = 1; red-black GS = 2
      call HYPRE_StructPFMGSetRelaxType(precond_solver, 2, ierr)
      call HYPRE_StructPFMGSetNumPreRelax(precond_solver, n_pre, ierr)
      call HYPRE_StructPFMGSetNumPostRelax(precond_solver, n_post, ierr)
      call HYPRE_StructPFMGSetTol(precond_solver, 0.0d0, ierr)
      ! call HYPRE_StructPFMGSetDxyz(precond_solver, dxyz, ierr)
      ! call HYPRE_StructPFMGSetSkipRelax(precond, 1, ierr)
      ! call HYPRE_StructPFMGSetRAPType(precond, 1, ierr)
    else if (precond == 7 .and. solver_id == 5) then
      call HYPRE_StructJacobiCreate(mpi_comm_hypre, precond_solver, ierr)
      call HYPRE_StructJacobiSetMaxIter(precond_solver, maxiter, ierr)
    else if (precond == 8) then
      precond_solver = c_null_ptr
    else if (precond == 9) then
      precond_solver = c_null_ptr
    else
      write (*,*) 'Invalid preconditioner in inithypre', precond
      write (*,*) 'Possbile values are (0) SMG (1) PFMG (8) DS (9) None'
      call exit(-1)
    endif
  end subroutine

  subroutine exitprecond
    use modglobal, only : precond

    implicit none

    if (precond == 0) then
      call HYPRE_StructSMGDestroy(precond_solver, ierr)
    else if (precond == 1) then
      call HYPRE_StructPFMGDestroy(precond_solver, ierr)
    else if (precond == 7) then
      call HYPRE_StructJacobiDestroy(precond_solver, ierr)
    endif
  end subroutine

  subroutine inithypre
    use mpi
    use modmpi, only : myid, myidx, myidy, nprocx, nprocy
    use modglobal, only : imax, jmax, kmax, dzf, dzh, dx, dy, itot, jtot, &
      solver_id, maxiter, n_pre, n_post, tolerance, precond

    use modfields, only : rhobf, rhobh

    implicit none

    ! NOTE: we assume kmax is larger than 7,
    ! when we make the stencil entries
    real, allocatable ::  values(:)
    real      cx, cy, cz_down, cz_up, cc

    ! integer   num_ghost(6)
    integer   offsets(3,7), stencil_indices(7)
    integer   i, j, k

    real wtime

    ! data num_ghost / 3, 3, 3, 3, 0, 0 /

    allocate(values(imax*jmax*kmax))

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

    WRITE (*,*) 'HYPRE Setup grid: myid, ilower, iupper', myid, ilower, iupper

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
    !call HYPRE_StructMatrixSetNumGhost(matrixA, num_ghost, ierr)
    call HYPRE_StructMatrixInitialize(matrixA, ierr)

    !-----------------------------------------------------------------------
    !     4. Set the coefficients for the grid
    !        Each processor describes the piece of the grid that it owns.
    !-----------------------------------------------------------------------

    wtime = MPI_Wtime()

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
    ! call HYPRE_StructMatrixPrint(matrixA, zero, ierr)

    wtime = MPI_Wtime() - wtime
    if (myid == 0) then
      write (*,*) 'HYPRE assembling struct matrix took', wtime
    endif

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

    if (solver_id == 1) then
      ! Solve the system using SMG
      if (myid == 0) then
        write (*,*) 'Selected solver 1 (SMG) with parameters:', maxiter, tolerance, n_pre, n_post
      endif
      call HYPRE_StructSMGCreate(mpi_comm_hypre, solver, ierr)
      call HYPRE_StructSMGSetMemoryUse(solver, zero, ierr)
      call HYPRE_StructSMGSetMaxIter(solver, maxiter, ierr)
      call HYPRE_StructSMGSetTol(solver, tolerance, ierr)
      call HYPRE_StructSMGSetRelChange(solver, zero, ierr)
      call HYPRE_StructSMGSetNumPreRelax(solver, n_pre, ierr)
      call HYPRE_StructSMGSetNumPostRelax(solver, n_post, ierr)
      call HYPRE_StructSMGSetLogging(solver, 1, ierr)
      call HYPRE_StructSMGSetup(solver, matrixA, vectorB, vectorX, ierr)
    else if (solver_id == 2) then
      ! Solve the system using PFMG
      if (myid == 0) then
        write (*,*) 'Selected solver 2 (PFMG) with parameters:', maxiter, tolerance, n_pre, n_post, precond
      endif
      call HYPRE_StructPFMGCreate(mpi_comm_hypre, solver, ierr)
      call HYPRE_StructPFMGSetMaxIter(solver, maxiter, ierr)
      call HYPRE_StructPFMGSetTol(solver, tolerance, ierr)
      call HYPRE_StructPFMGSetRelChange(solver, zero, ierr)

      ! weighted Jacobi = 1; red-black GS = 2
      if (precond == 1 .or. precond == 2) then
        call HYPRE_StructPFMGSetRelaxType(solver, precond, ierr)
      else
        write (*,*) 'Invalid preconditioner in inithypre', precond
        write (*,*) 'Possbile values are (1) weighted jacobi (2) red-black GS'
        call exit(-1)
      endif

      call HYPRE_StructPFMGSetNumPreRelax(solver, n_pre, ierr)
      call HYPRE_StructPFMGSetNumPostRelax(solver, n_post, ierr)
      ! call HYPRE_StructPFMGSetSkipRelax(solver, one)
      ! call HYPRE_StructPFMGSetDxyz(solver, dxyz, ierr)

      call HYPRE_StructPFMGSetLogging(solver, 2, ierr)
      call HYPRE_StructPFMGSetup(solver, matrixA, vectorB, vectorX, ierr)

    else if (solver_id == 3) then
      ! Solve the system using BiCGSTAB
      if (myid == 0) then
        write (*,*) 'Selected solver 3 (BiCGSTAB) with parameters:', maxiter, tolerance, n_pre, n_post, precond
      endif
      call HYPRE_StructBiCGSTABCreate(mpi_comm_hypre, solver, ierr)
      call HYPRE_StructBiCGSTABSetMaxIter(solver, maxiter, ierr)
      call HYPRE_StructBiCGSTABSetTol(solver, tolerance, ierr)
      call initprecond
      call HYPRE_StructBiCGSTABSetPrecond(solver, precond, precond_solver, ierr)

      call HYPRE_StructBiCGSTABSetLogging(solver, 2, ierr)
      call HYPRE_StructBiCGSTABSetup(solver, matrixA, vectorB, vectorX, ierr)

    else if (solver_id == 4) then
      ! Solve the system using GMRES
      if (myid == 0) then
        write (*,*) 'Selected solver 4 (GMRES) with parameters:', maxiter, tolerance, n_pre, n_post, precond
      endif
      call HYPRE_StructGMRESCreate(mpi_comm_hypre, solver, ierr)
      call HYPRE_StructGMRESSetTol(solver, tolerance, ierr)
      call HYPRE_StructGMRESSetMaxIter(solver, maxiter, ierr)
      call initprecond
      call HYPRE_StructGMRESSetPrecond(solver, precond, precond_solver, ierr)

      call HYPRE_StructGMRESSetLogging(solver, 2, ierr)
      call HYPRE_StructGMRESSetup(solver, matrixA, vectorB, vectorX, ierr)

    else if (solver_id == 5) then
      ! Solve the system using PCG
      if (myid == 0) then
        write (*,*) 'Selected solver 5 (PCG) with parameters:', maxiter, tolerance, n_pre, n_post, precond
      endif
      call HYPRE_StructPCGCreate(mpi_comm_hypre, solver, ierr)
      call HYPRE_StructPCGSetTol(solver, tolerance, ierr)
      call HYPRE_StructPCGSetTwoNorm(solver, 1, ierr)
      call HYPRE_StructPCGSetMaxIter(solver, maxiter, ierr)
      call initprecond
      call HYPRE_StructPCGSetPrecond(solver, precond, precond_solver, ierr)

      call HYPRE_StructPCGSetLogging(solver, 2, ierr)
      call HYPRE_StructPCGSetPrintLevel(solver, 0, ierr)
      call HYPRE_StructPCGSetup(solver, matrixA, vectorB, vectorX, ierr)

    else if (solver_id == 6) then
      ! Solve the system using LGMRES
      if (myid == 0) then
        write (*,*) 'Selected solver 6 (LGMRES) with parameters:', maxiter, tolerance, n_pre, n_post, precond
      endif
      call HYPRE_StructLGMRESCreate(mpi_comm_hypre, solver, ierr)
      call HYPRE_StructLGMRESSetTol(solver, tolerance, ierr)
      call HYPRE_StructLGMRESSetMaxIter(solver, maxiter, ierr)
      call initprecond
      call HYPRE_StructLGMRESSetPrecond(solver, precond, precond_solver, ierr)

      call HYPRE_StructLGMRESSetLogging(solver, 2, ierr)
      call HYPRE_StructLGMRESSetPrintLevel(solver, 0, ierr)
      call HYPRE_StructLGMRESSetup(solver, matrixA, vectorB, vectorX, ierr)

    else
      if (myid == 0) then
        write (*,*) 'Invalid solver in inithypre', solver
      endif
      call exit(-1)
    endif

    deallocate(values)

  end subroutine

  subroutine set_initial_guess(p)
    use modglobal, only : i1, j1, ih, jh, imax, jmax, kmax

    implicit none

    real, intent(inout) :: p(2-ih:i1+ih,2-jh:j1+jh,kmax)
    real values(imax,jmax)

    integer i,j,k
    do k=1,kmax
      do j=1,jmax
        do i=1,imax
          values(i,j) = p(i,j,k)
        enddo
      enddo
      ilower(3) = k - 1
      iupper(3) = k - 1
      call HYPRE_StructVectorSetBoxValues(vectorX, ilower, iupper, values, ierr)
    enddo
    call HYPRE_StructVectorAssemble(vectorX, ierr)
  end subroutine

  subroutine solve_hypre(p, converged)
    use modmpi, only : myid
    use modglobal, only : i1, j1, ih, jh, imax, jmax, kmax, solver_id, maxiter

    implicit none

    real, intent(inout) :: p(2-ih:i1+ih,2-jh:j1+jh,kmax)
    logical, intent(out) :: converged
    real values(imax,jmax)

    real final_res_norm
    integer i,j,k, num_iterations, stat
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

    !-----------------------------------------------------------------------
    !     2. Call a solver
    !-----------------------------------------------------------------------

    if (solver_id == 1) then
      ! Solve the system using SMG
      call HYPRE_StructSMGSolve(solver, matrixA, vectorB, vectorX, stat)
      call HYPRE_StructSMGGetNumIterations(solver, num_iterations, ierr)
      call HYPRE_StructSMGGetFinalRelative(solver, final_res_norm, ierr)

    else if (solver_id == 2) then
      ! Solve the system using PFMG
      call HYPRE_StructPFMGSolve(solver, matrixA, vectorB, vectorX, stat)
      call HYPRE_StructPFMGGetNumIteration(solver, num_iterations, ierr)
      call HYPRE_StructPFMGGetFinalRelativ(solver, final_res_norm, ierr)

    else if (solver_id == 3) then
      ! Solve the system using BiCGSTAB
      call HYPRE_StructBiCGSTABSolve(solver, matrixA, vectorB, vectorX, stat)
      call HYPRE_StructBiCGSTABGetNumItera(solver, num_iterations, ierr)
      call HYPRE_StructBiCGSTABGetFinalRel(solver, final_res_norm, ierr)

    else if (solver_id == 4) then
      ! Solve the system using GMRES
      call HYPRE_StructGMRESSolve(solver, matrixA, vectorB, vectorX, stat)
      call HYPRE_StructGMRESGetNumIteratio(solver, num_iterations, ierr)
      call HYPRE_StructGMRESGetFinalRelati(solver, final_res_norm, ierr)

    else if (solver_id == 5) then
      ! Solve the system using PCG
      call HYPRE_StructPCGSolve(solver, matrixA, vectorB, vectorX, stat)
      call HYPRE_StructPCGGetNumIterations(solver, num_iterations, ierr)
      call HYPRE_StructPCGGetFinalRelative(solver, final_res_norm, ierr)

    else if (solver_id == 6) then
      ! Solve the system using LGMRES
      call HYPRE_StructLGMRESSolve(solver, matrixA, vectorB, vectorX, stat)
      call HYPRE_StructLGMRESGetNumIter(solver, num_iterations, ierr)
      call HYPRE_StructLGMRESGetFinalRel(solver, final_res_norm, ierr)

    else

      write (*,*) 'Invalid solver in solve_hypre', solver_id
      call exit(-1)
    endif

    if (myid == 0) then
      write (*,*) 'HYPRE Iterations, residual norm', num_iterations, final_res_norm
    endif
    if (stat /= 0) then
      if (stat == 1) then
        write (*,*) 'HYPRE solver status (ierr)', stat, 'generic error'
      else if (stat == 2) then
        write (*,*) 'HYPRE solver status (ierr)', stat, 'unable to allocate memory'
      else if (stat == 4) then
        write (*,*) 'HYPRE solver status (ierr)', stat, 'argument error'
      else if (stat == 256) then
        write (*,*) 'HYPRE solver status (ierr)', stat, 'method did not converge as expected'
        converged = .false.
        return
      endif
      call exit(-1)
    endif

    if (num_iterations >= maxiter) then
      converged = .false.
      return
    else
      converged = .true.
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
    use modglobal, only : solver_id, precond

    implicit none

    if (solver_id == 1) then
      call HYPRE_StructSMGDestroy(solver, ierr)
    else if (solver_id == 2) then
      call HYPRE_StructPFMGDestroy(solver, ierr)
    else if (solver_id == 3) then
      if (precond == 0) then
        call HYPRE_StructSMGDestroy(precond_solver, ierr)
      else if (precond == 1) then
        call HYPRE_StructPFMGDestroy(precond_solver, ierr)
      endif
      call HYPRE_StructBiCGSTABDestroy(solver, ierr)
    else if (solver_id == 4) then
      call exitprecond
      call HYPRE_StructGMRESDestroy(solver, ierr)
    else if (solver_id == 5) then
      call exitprecond
      call HYPRE_StructPCGDestroy(solver, ierr)
    else if (solver_id == 6) then
      call exitprecond
      call HYPRE_StructLGMRESDestroy(solver, ierr)
    else
      write (*,*) 'Invalid solver in exit_hypre', solver_id
      call exit(-1)
    endif

    call HYPRE_StructGridDestroy(grid, ierr)
    call HYPRE_StructStencilDestroy(stencil, ierr)
    call HYPRE_StructMatrixDestroy(matrixA, ierr)
    call HYPRE_StructVectorDestroy(vectorB, ierr)
    call HYPRE_StructVectorDestroy(vectorX, ierr)
  end subroutine

end module
