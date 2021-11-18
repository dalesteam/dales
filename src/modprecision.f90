!> \file modprecision.f90

!>
!! Declare all precisions used throughout DALES
!>
!! We will start with the precision on a per-file basis 
!!  \author Victor Azizi, Escience Center
!!  \par Revision list

!! Notes
!! Real(4) does not (necessarily) means a 4-byte single precision floating point
!! this is up to the compiler to decide. -freal-4-real8 breaks this for example

! Preprocessor macro to concatenate 2 strings
! Gfortran and Nag use traditional c-preprocessor, therefore the macro is
! different for them
#if defined(__GFORTRAN__) || defined(NAGFOR)
#define PASTE(a) a
#define CAT(a,b) PASTE(a)b
#else
#define PASTE(a,b) a ## b
#define CAT(a,b) PASTE(a,b)
#endif

module modprecision
use iso_fortran_env, only : real32, real64, int32, int64

integer, parameter :: longint = int64

!< double precision kind for lacz_gamma
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)

!< Parameter kinds, for rrtmg radiation scheme
integer, parameter :: kind_rb  = real64     !selected_real_kind(12) ! 8 byte real
integer, parameter :: kind_im  = int32      !selected_int_kind(6)   ! 4 byte integer
integer, parameter :: SHR_KIND_R4 = real32  !selected_real_kind( 6) ! 4 byte real
integer, parameter :: SHR_KIND_IN = kind(1) ! native integer

!! Module specific precisions
  ! Precision for the most common fields u, v, w
  ! And all other fields that do not have
  ! their own kind and need to interoperate
integer, parameter :: field_r = CAT(real,FIELD_PRECISION)
                                      
  ! Precision for the poisson solver,
  ! HYPRE is double-only, so real32 doesn't
  ! work with -DUSE_HYPRE=T 
  ! Precision for p and Fp fields
  ! Precision for all internal fields in
  ! the old fft and fftw module
integer, parameter :: pois_r  = CAT(real,POIS_PRECISION)
                                      
end module
