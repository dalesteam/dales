!> Interface to either cuFFT or hipFFT. Prevents a preprocessor mess in modgpufft.
module hicfft
  use iso_c_binding
  use modprecision, only: pois_r
#if HAVE_CUDA
  use cufft
#elif HAVE_HIP
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

  type plan_t
#if HAVE_HIP
    type(c_ptr) :: handle
#else
    integer :: handle
#endif
  end type plan_t

#if HAVE_CUDA
#if POIS_PRECISION==32
  integer, parameter :: HICFFT_FWD_TYPE = CUFFT_R2C
  integer, parameter :: HICFFT_BWD_TYPE = CUFFT_C2R
#else
  integer, parameter :: HICFFT_FWD_TYPE = CUFFT_D2Z
  integer, parameter :: HICFFT_BWD_TYPE = CUFFT_Z2D
#endif
#elif HAVE_HIP
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

  function hicfftDestroy(plan) result(istat)
    type(plan_t), intent(inout) :: plan
    integer :: istat
#if HAVE_CUDA
    istat = cufftDestroy(plan%handle)
#elif HAVE_HIP
    istat = hipfftDestroy(plan%handle)
#endif
  end function hicfftDestroy

  function hicfftExecForward(plan, idata, odata) result(istat)
    type(plan_t), intent(in) :: plan
    real(pois_r), intent(in) :: idata(:,:,:)
    real(pois_r), intent(out) :: odata(:,:,:)
    integer :: istat
#if POIS_PRECISION==32
#if HAVE_CUDA
    istat = cufftExecR2C(plan%handle, idata, odata)
#elif HAVE_HIP
    istat = hipfftExecR2C(plan%handle, c_loc(idata), c_loc(odata))
#endif
#else
#if HAVE_CUDA
    istat = cufftExecD2Z(plan%handle, idata, odata)
#elif HAVE_HIP
    istat = hipfftExecD2Z(plan%handle, c_loc(idata), c_loc(odata))
#endif
#endif
  end function hicfftExecForward

  function hicfftExecBackward(plan, idata, odata) result(istat)
    type(plan_t), intent(in) :: plan
    real(pois_r), intent(in) :: idata(:,:,:)
    real(pois_r), intent(out) :: odata(:,:,:)
    integer :: istat
#if POIS_PRECISION==32
#if HAVE_CUDA
    istat = cufftExecC2R(plan%handle, idata, odata)
#elif HAVE_HIP
    istat = hipfftExecC2R(plan%handle, c_loc(idata), c_loc(odata))
#endif
#else
#if HAVE_CUDA
    istat = cufftExecZ2D(plan%handle, idata, odata)
#elif HAVE_HIP
    istat = hipfftExecZ2D(plan%handle, c_loc(idata), c_loc(odata))
#endif
#endif
  end function hicfftExecBackward

  function hicfftGetSize(plan, workSize) result(istat)
    type(plan_t), intent(in) :: plan
    integer(8), intent(out) :: workSize
    integer :: istat
#if HAVE_CUDA
    istat = cufftGetSize(plan%handle, workSize)
#elif HAVE_HIP
    istat = hipfftGetSize(plan%handle, c_loc(workSize))
#endif
  end function hicfftGetSize

 function hicfftPlanMany(plan, rank, n, inembed, istride, idist, onembed, &
                          ostride, odist, myType, batch) result(istat)
    type(plan_t), intent(inout) :: plan
    integer, intent(in) :: rank, n, inembed, istride, idist, onembed, ostride, &
                           odist, myType, batch
    integer :: istat 
#if HAVE_CUDA
    istat = cufftPlanMany(plan%handle, rank, n, inembed, istride, idist, &
                          onembed, ostride, odist, myType, batch)
#elif HAVE_HIP
    istat = hipfftPlanMany(plan%handle, rank, n, inembed, istride, idist, &
                           onembed, ostride, odist, myType, batch)
#endif
  end function hicfftPlanMany

  function hicfftSetAutoAllocation(plan, autoAllocate) result(istat)
    type(plan_t), intent(inout) :: plan
    integer, intent(in) :: autoAllocate
    integer :: istat
#if HAVE_CUDA
    istat = cufftSetAutoAllocation(plan%handle, autoAllocate)
#elif HAVE_HIP
    istat = hipfftSetAutoAllocation(plan%handle, autoAllocate)
#endif
  end function hicfftSetAutoAllocation

  function hicfftSetWorkArea(plan, workspace) result(istat)
    type(plan_t), intent(inout) :: plan
    real(pois_r), intent(in) :: workspace(:)
    integer :: istat
    !$acc host_data use_device(workspace)
#if HAVE_CUDA
    istat = cufftSetWorkArea(plan%handle, workspace)
#elif HAVE_HIP
    istat = hipfftSetWorkArea(plan%handle, c_loc(workspace))
#endif
    !$acc end host_data
  end function hicfftSetWorkArea

end module hicfft
