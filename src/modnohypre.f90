!> \file modhypre.f90
!! Dummy iterative solver
!
!  Dummy functions for compiling DALES without the HYPRE library.
!  All functions will print an error message and abort.
!
!  2019 Jisk Attema, Netherlands eScience Center
!

module modhypre

implicit none
private
public :: inithypre, solve_hypre, exithypre, set_initial_guess

contains

  subroutine inithypre
    implicit none

    call error_and_exit()
  end subroutine

  subroutine exithypre
    implicit none

    call error_and_exit()
  end subroutine

  subroutine set_initial_guess(p)
    use modglobal, only : i1, j1, ih, jh, kmax

    implicit none

    real, intent(inout) :: p(2-ih:i1+ih,2-jh:j1+jh,kmax)

    call error_and_exit()
  end subroutine

  subroutine solve_hypre(p)
    use modglobal, only : i1, j1, ih, jh, kmax

    implicit none

    real, intent(inout) :: p(2-ih:i1+ih,2-jh:j1+jh,kmax)

    call error_and_exit()
  end subroutine

  subroutine error_and_exit
    implicit none

    write (*,*) 'DALES was compiled without HYPRE.'
    write (*,*) 'Use the poisson solver (solver_id=0),'
    write (*,*) 'or recompile DALES with the option USE_HYPRE.'

    call exit(-1)
  end subroutine

end module
