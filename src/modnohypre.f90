!> \file modhypre.f90
!! Dummy iterative solver
!
!  Dummy functions for compiling DALES without the HYPRE library.
!  All functions will print an error message and abort.
!
!  2019 Jisk Attema, Netherlands eScience Center
!

module modhypre
use modglobal, only : solver_type
implicit none
private
public :: inithypre_grid, inithypre_solver, solve_hypre, exithypre_grid, exithypre_solver, set_initial_guess, set_zero_guess

contains

  subroutine inithypre_grid
    implicit none

    call error_and_exit()
  end subroutine

  subroutine inithypre_solver(solver,solver_id,maxiter,tolerance,precond_id,n_pre,n_post)
    implicit none
    type(solver_type), intent(inout) :: solver
    real, intent(in) :: tolerance
    integer, intent(in) :: solver_id, maxiter, precond_id,n_pre,n_post

    call error_and_exit()
  end subroutine

  subroutine exithypre_grid
    implicit none

    call error_and_exit()
  end subroutine

  subroutine exithypre_solver(solver)
    implicit none
    type(solver_type), intent(inout) :: solver

    call error_and_exit()
  end subroutine

  subroutine set_initial_guess(p)
    use modglobal, only : i1, j1, ih, jh, kmax
    implicit none

    real, intent(inout) :: p(2-ih:i1+ih,2-jh:j1+jh,kmax)

    call error_and_exit()
  end subroutine

  subroutine set_zero_guess()
    implicit none

    call error_and_exit()
  end subroutine

  subroutine solve_hypre(solver, p, converged)
    use modglobal, only : i1, j1, ih, jh, kmax
    implicit none

    real, intent(inout) :: p(2-ih:i1+ih,2-jh:j1+jh,kmax)
    logical, intent(out) :: converged
    type(solver_type), intent(inout) :: solver

    call error_and_exit()
    converged = .false. ! suppress warnings about intent(out) variable
                        ! not being assigned
  end subroutine

  subroutine error_and_exit
    implicit none

    write (*,*) 'DALES was compiled without HYPRE.'
    write (*,*) 'Use the poisson solver (solver_id=0),'
    write (*,*) 'or recompile DALES with the option USE_HYPRE.'

    call exit(-1)
  end subroutine

end module
