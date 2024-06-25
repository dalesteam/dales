#!/usr/bin/env python
"""This is a simple example for tuning Fortran OpenACC code with the kernel tuner"""

from kernel_tuner import tune_kernel
from kernel_tuner.utils.directives import (
    Code,
    OpenACC,
    Fortran,
    extract_directive_signature,
    extract_directive_code,
    extract_preprocessor,
    generate_directive_function,
    extract_directive_data,
    allocate_signature_memory,
    extract_initialization_code
)

code = """
#define jmax 128
#define imax 128
#define kmax 80
#define vert_size 3180
#define horz_size 6400

subroutine advecu_2nd
  !$tuner initialize
  real(8) :: dxiq, dyiq, dziq
  real(8), dimension(imax+2,jmax+2,kmax+1) :: a_in, a_out, u0, v0, w0
  real(8), dimension(kmax+1) :: rhobf
  integer :: i,j,k,ip,im,jp,jm,kp,km,i1,j1
  i1 = imax + 1
  j1 = jmax + 1
  a_in = 0.
  a_out = 0.
  u0 = 1.
  v0 = 1.5
  w0 = 0.5
  rhobf = 1.15
  dxiq = 1 / (horz_size / imax)**4
  dyiq = 1 / (horz_size / jmax)**4
  dziq = 1 / (vert_size / kmax)**4
  !$acc enter data copyin(a_in(:imax+2,:jmax+2,:kmax+1))
  !$acc enter data copyin(a_out(:imax+2,:jmax+2,:kmax+1))
  !$acc enter data copyin(u0(:imax+2,:jmax+2,:kmax+1))
  !$acc enter data copyin(v0(:imax+2,:jmax+2,:kmax+1))
  !$acc enter data copyin(w0(:imax+2,:jmax+2,:kmax+1))
  !$acc enter data copyin(rhobf(:kmax+1))
  !$tuner stop

  !$tuner start advecu_2nd 
  !$acc parallel loop collapse(ncoll) default(present) vector_length(nthreads)
  do k = 1, kmax
    do j = 2, jmax + 1
      do i = 2, imax + 1
        a_out(i,j,k)  = a_out(i,j,k)- ( &
                ( &
                (a_in(i,j,k)+a_in(i+1,j,k))*(u0(i,j,k)+u0(i+1,j,k)) &
                -(a_in(i,j,k)+a_in(i-1,j,k))*(u0(i,j,k)+u0(i-1,j,k)) &
                )*dxiq &
                +(  &
                (a_in(i,j,k)+a_in(i,j+1,k))*(v0(i,j+1,k)+v0(i-1,j+1 ,k)) &
                -(a_in(i,j,k)+a_in(i,j-1,k))*(v0(i,j  ,k)+v0(i-1,j  ,k)) &
                )*dyiq )
      end do
    end do
  end do
  !$tuner stop
  !$tuner deinitialize
  !$acc exit data delete(a_in(:imax+2,:jmax+2,:kmax+1))
  !$acc exit data delete(a_out(:imax+2,:jmax+2,:kmax+1))
  !$acc exit data delete(u0(:imax+2,:jmax+2,:kmax+1))
  !$acc exit data delete(v0(:imax+2,:jmax+2,:kmax+1))
  !$acc exit data delete(w0(:imax+2,:jmax+2,:kmax+1))
  !$acc exit data delete(rhobf(:kmax+1))
  !$tuner stop
end subroutine advecu_2nd
"""

# Extract tunable directive
app = Code(OpenACC(), Fortran())
preprocessor = extract_preprocessor(code)
signature = extract_directive_signature(code, app)
body = extract_directive_code(code, app)
# Allocate memory on the host
data = extract_directive_data(code, app)
# Generate kernel string
init = extract_initialization_code(code, app)

args = allocate_signature_memory(data["advecu_2nd"], preprocessor)
kernel_string = generate_directive_function(
    preprocessor, signature["advecu_2nd"], body["advecu_2nd"], app, data=data["advecu_2nd"], initialization=init
)
tune_params = dict()
tune_params["nthreads"] = [32 * i for i in range(1, 33)]
tune_params["ncoll"] = [1, 2]
metrics = dict()
#metrics["GB/s"] = lambda x: ((2 * 4 * len(args[0])) + (4 * len(args[0]))) / (x["time"] / 10**3) / 10**9

#answer = [None, None, args[0] + args[1], None, None]

tune_kernel(
    "advecu_2nd",
    kernel_string,
    0,
    args,
    tune_params,
    metrics=metrics,
#    answer=answer,
    compiler_options=["-fast", "-acc=gpu", "-gpu=cc80"],
    compiler="nvfortran",
)
