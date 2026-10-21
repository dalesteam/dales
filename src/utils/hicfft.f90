!> Interface to either cuFFT or hipFFT. Prevents a preprocessor mess in modgpufft.
module hicfft
  use iso_c_binding
  use modprecision,    only: pois_r
  use fortran_support, only: finish
#if USE_CUDA
  use cufft
#elif USE_HIP
  use hipfort_hipfft
#endif 

  implicit none

  private

  public :: plan_t
  public :: HICFFT_FWD_TYPE
  public :: HICFFT_BWD_TYPE
  public :: hicfftDestroy
  public :: hicfftExecForward
  public :: hicfftExecBackward
  public :: hicfftGetSize
  public :: hicfftPlanMany
  public :: hicfftSetAutoAllocation
  public :: hicfftSetWorkArea
  
  character(len=*), parameter :: modname = "hicfft"

  type plan_t
#if USE_CUDA
    integer :: handle
#else
    type(c_ptr) :: handle
#endif
  end type plan_t

#if USE_CUDA
#if POIS_PRECISION==32
  integer, parameter :: HICFFT_FWD_TYPE = CUFFT_R2C
  integer, parameter :: HICFFT_BWD_TYPE = CUFFT_C2R
#else
  integer, parameter :: HICFFT_FWD_TYPE = CUFFT_D2Z
  integer, parameter :: HICFFT_BWD_TYPE = CUFFT_Z2D
#endif
#elif USE_HIP
#if POIS_PRECISION==32
  integer, parameter :: HICFFT_FWD_TYPE = HIPFFT_R2C
  integer, parameter :: HICFFT_BWD_TYPE = HIPFFT_C2R
#else
  integer, parameter :: HICFFT_FWD_TYPE = HIPFFT_D2Z
  integer, parameter :: HICFFT_BWD_TYPE = HIPFFT_Z2D
#endif
#else
  integer, parameter :: HICFFT_FWD_TYPE = 0
  integer, parameter :: HICFFT_BWD_TYPE = 0
#endif

contains

  subroutine hicfftDestroy(plan)
    type(plan_t), intent(inout) :: plan
    character(len=*), parameter :: routine = modname//"hicfftDestroy"
#if USE_CUDA
    call check(cufftDestroy(plan%handle), routine)
#elif USE_HIP
    call check(hipfftDestroy(plan%handle), routine)
#endif
  end subroutine hicfftDestroy

  subroutine hicfftExecForward(plan, idata, odata)
    type(plan_t), intent(in) :: plan
    real(pois_r), pointer, intent(in) :: idata(:,:,:)
    real(pois_r), pointer, intent(out) :: odata(:,:,:)
    character(len=*), parameter :: routine = modname//"hicfftExecForward"
    !$acc host_data use_device(idata, odata)
    !$omp target data use_device_addr(idata, odata)
#if POIS_PRECISION==32
#if USE_CUDA
    call check(cufftExecR2C(plan%handle, idata, odata), routine)
#elif USE_HIP
    call check(hipfftExecR2C(plan%handle, c_loc(idata), c_loc(odata)), routine)
#endif
#else
#if USE_CUDA
    call check(cufftExecD2Z(plan%handle, idata, odata), routine)
#elif USE_HIP
    call check(hipfftExecD2Z(plan%handle, c_loc(idata), c_loc(odata)), routine)
#endif
#endif
    !$acc end host_data
    !$omp end target data
  end subroutine hicfftExecForward

  subroutine hicfftExecBackward(plan, idata, odata)
    type(plan_t), intent(in) :: plan
    real(pois_r), pointer, intent(in) :: idata(:,:,:)
    real(pois_r), pointer, intent(out) :: odata(:,:,:)
    character(len=*), parameter :: routine = modname//"hicfftExecBackward"
    !$acc host_data use_device(idata, odata)
    !$omp target data use_device_addr(idata, odata)
#if POIS_PRECISION==32
#if USE_CUDA
    call check(cufftExecC2R(plan%handle, idata, odata), routine)
#elif USE_HIP
    call check(hipfftExecC2R(plan%handle, c_loc(idata), c_loc(odata)), routine)
#endif
#else
#if USE_CUDA
    call check(cufftExecZ2D(plan%handle, idata, odata), routine)
#elif USE_HIP
    call check(hipfftExecZ2D(plan%handle, c_loc(idata), c_loc(odata)), routine)
#endif
#endif
    !$acc end host_data
    !$omp end target data
  end subroutine hicfftExecBackward

  subroutine hicfftGetSize(plan, workSize)
    type(plan_t), intent(in) :: plan
    integer(c_intptr_t), intent(out) :: workSize
    character(len=*), parameter :: routine = modname//"hicfftGetSize"
#if USE_CUDA
    call check(cufftGetSize(plan%handle, workSize), routine)
#elif USE_HIP
    call check(hipfftGetSize(plan%handle, workSize), routine)
#endif
  end subroutine hicfftGetSize

  subroutine hicfftPlanMany(plan, rank, n, inembed, istride, idist, onembed, &
                          ostride, odist, myType, batch)
    type(plan_t), intent(inout) :: plan
    integer, intent(in) :: rank, n, inembed, istride, idist, onembed, ostride, &
                           odist, myType, batch
    character(len=*), parameter :: routine = modname//"hicfftPlanMany"
#if USE_CUDA
    call check(cufftPlanMany(plan%handle, rank, n, inembed, istride, idist, &
                             onembed, ostride, odist, myType, batch), routine)
#elif USE_HIP
    call check(hipfftPlanMany(plan%handle, rank, n, inembed, istride, idist, &
                              onembed, ostride, odist, myType, batch), routine)
#endif
  end subroutine hicfftPlanMany

  subroutine hicfftSetAutoAllocation(plan, autoAllocate)
    type(plan_t), intent(inout) :: plan
    integer, intent(in) :: autoAllocate
    character(len=*), parameter :: routine = modname//"hicfftSetAutoAllocation"
#if USE_CUDA
    call check(cufftSetAutoAllocation(plan%handle, autoAllocate), routine)
#elif USE_HIP
    call check(hipfftSetAutoAllocation(plan%handle, autoAllocate), routine)
#endif
  end subroutine hicfftSetAutoAllocation

  subroutine hicfftSetWorkArea(plan, workspace)
    type(plan_t), intent(inout) :: plan
    real(pois_r), allocatable, target, intent(in) :: workspace(:)
    character(len=*), parameter :: routine = modname//"hicfftSetWorkArea"
    !$acc host_data use_device(workspace)
    !$omp target data use_device_addr(workspace)
#if USE_CUDA
    call check(cufftSetWorkArea(plan%handle, workspace), routine)
#elif USE_HIP
    call check(hipfftSetWorkArea(plan%handle, c_loc(workspace)), routine)
#endif
    !$acc end host_data
    !$omp end target data
  end subroutine hicfftSetWorkArea

  subroutine check(istat, routine)

    integer,          intent(in) :: istat
    character(len=*), intent(in) :: routine

    character(len=64) :: msg

    if (istat /= 0) then
#if USE_CUDA
      write(msg, *) "cuFFT returned non-zero exit code", istat
#elif USE_HIP
      write(msg, *) "hipFFT returned non-zero exit code", istat
#endif
      call finish(routine, msg)
    end if

  end subroutine check

end module hicfft
