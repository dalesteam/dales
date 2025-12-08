!> Common declarations for FFT-based and HYPRE Poisson solvers.
module modpois_data
  
  use, intrinsic :: iso_c_binding

  use modprecision, only: real64, pois_r

  implicit none

  integer :: solver_id = 0

  real(pois_r), pointer     :: p(:,:,:)  !< Pressure fluctuations in real space.
  real(pois_r), pointer     :: Fp(:,:,:) !< Pressure fluctuations in Fourier space.
  real(pois_r), allocatable :: d(:,:,:)  !< Work array for tridiagonal solver.
  real(pois_r), allocatable :: xyrt(:,:) !< Eigenvalues.

  real(pois_r), allocatable :: pup(:,:,:) !< Work array for rhs.
  real(pois_r), allocatable :: pvp(:,:,:) !< Work array for rhs.
  real(pois_r), allocatable :: pwp(:,:,:) !< Work array for rhs.
  real(pois_r), allocatable :: a(:)       !< Work array for solver.
  real(pois_r), allocatable :: b(:)       !< Work array for solver.
  real(pois_r), allocatable :: c(:)       !< Work array for solver.

  integer :: ps, pe, qs, qe !< Start and end indices of Fourier space matrices.

  ! HYPRE-specific variables
  integer      :: maxiter = 10000     !< Number of iterations.
  real(real64) :: tolerance = 1E-8    !< Convergence threshold.
  integer      :: n_pre = 1           !< Number of pre relaxations.
  integer      :: n_post = 1          !< Number of post relaxations.
  integer      :: precond_id = 1      !< Preconditioner ID.
  integer      :: maxiter_precond = 1 !< Number of iterations for preconditioner.
  integer      :: hypre_logging = 1   !< HYPRE logging and print level.

  type solver_type
    type(c_ptr)  :: solver, precond
    integer      :: solver_id, precond_id, maxiter, n_post, n_pre, &
                    maxiter_precond
    real(real64) :: tolerance
  end type

  type(solver_type) :: psolver

end module modpois_data