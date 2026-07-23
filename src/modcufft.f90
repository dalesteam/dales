module modcufft
  use, intrinsic :: iso_c_binding 

  use modtimer
  use modmpi
  use modglobal, only: itot, jtot, kmax, i1, j1, &
                       imax, jmax, ih, jh, dxi, dyi, pi, ijtot
  use modprecision, only: pois_r
  use modtranspose, only: t_transposer
  use modgpu,       only: workspace_0
  use hicfft,       only: plan_t, HICFFT_FWD_TYPE, HICFFT_BWD_TYPE, hicfftDestroy, &
                          hicfftExecForward, hicfftExecBackward, hicfftGetSize, &
                          hicfftPlanMany, hicfftSetAutoAllocation, &
                          hicfftsetWorkArea

  implicit none

  character(len=*), parameter :: modname = 'modcufft'

#if defined(DALES_GPU)

  save
    real :: norm_fac !< Normalization factor

    real(pois_r), allocatable, target :: p_halo(:) !< Pressure with halos
    real(pois_r), allocatable, target :: p_nohalo(:)
    real(pois_r), pointer :: px(:,:,:), py(:,:,:)

    integer :: nphix, nphiy
    integer :: konx, kony, iony, jonx

    type(plan_t) :: planx, planxi, plany, planyi !< Plan handles

    integer(c_intptr_t) :: worksize, max_worksize !< Size of the required workspace

    type(t_transposer) :: transposer

  contains
    !< Setup plans, workspace, etc
    subroutine cufftinit(p, Fp, d, xyrt, ps, pe, qs, qe)
      use modgpu, only: workspace_0, allocate_workspace

      implicit none

      real(pois_r), pointer :: p(:,:,:) !< Pressure, spatial domain, with halos
      real(pois_r), pointer :: Fp(:,:,:) !< Pressure, spectral domain, with halos
      real(pois_r), allocatable :: xyrt(:,:) !< Array of eigenvalues
      real(pois_r), allocatable :: d(:,:,:)
      integer, intent(out) :: ps, pe, qs, qe

      integer(kind=8) :: sz
      integer :: fftsize, inembed, onembed, idist, odist, istride, ostride

      ! Dimensions of the transposes
      ! For explanation of the variables, see modfftw.f90/fftwinit

      transposer = t_transposer()
      sz = transposer%get_buffer_size()

      !$acc enter data copyin(transposer)
      !$omp target enter data map(to:transposer)

      konx = transposer%konx
      iony = transposer%iony
      jonx = transposer%jonx

      ! Number of complex coefficients
      nphix = itot/2 + 1
      nphiy = jtot/2 + 1
      
      sz = max(kmax * imax * jmax, &
               konx * (2 * nphix) * jmax, &
               konx * imax * (2 * nphiy))

      ! Allocate memory for the pressure
      allocate(p_halo(1:(imax+2*ih)*(jmax+2*jh)*kmax))
      allocate(p_nohalo(sz))

      !$acc enter data create(p_halo, p_nohalo)
      !$omp target enter data map(alloc:p_halo,p_nohalo)

      p(2-ih:i1+ih,2-jh:j1+jh,1:kmax) => p_halo(1:(imax+2*ih)*(jmax+2*jh)*kmax) ! z-aligned
      px(1:nphix*2,1:jmax,1:konx) => p_nohalo(1:konx*jmax*(nphix*2)) ! x-aligned
      py(1:nphiy*2,1:konx,1:iony) => p_nohalo(1:konx*iony*(nphiy*2)) ! y-aligned

      if (nprocs == 1) then
        Fp(2-ih:i1+ih,2-jh:j1+jh,1:kmax) => p_halo(1:(imax+2*ih)*(jmax+2*jh)*kmax) ! z-aligned
      else
        Fp(1:iony,1:jonx,1:kmax) => p_halo(1:iony*jonx*kmax)
      end if

      ! x-direction
      fftsize = itot
      inembed = itot
      onembed = nphiy
      idist = nphix*2
      odist = nphix
      istride = 1
      ostride = 1

      call hicfftPlanMany( &
        planx, &
        1, &
        fftsize, &
        inembed, &
        istride, &
        idist, &
        onembed, &
        ostride, &
        odist, &
        HICFFT_FWD_TYPE, &
        jmax*konx &
      )
      call hicfftSetAutoAllocation(planx, 0)

      call hicfftPlanMany( &
        planxi, &
        1, &
        fftsize, &
        onembed, &
        ostride, &
        odist, &
        inembed, &
        istride, &
        idist, &
        HICFFT_BWD_TYPE, &
        jmax*konx &
      )
      call hicfftSetAutoAllocation(planxi, 0)

      ! y-direction

      fftsize = jtot
      inembed = jtot
      onembed = nphiy
      idist = nphiy*2
      odist = nphiy
      istride = 1
      ostride = 1

      call hicfftPlanMany( &
        plany, &
        1, &
        fftsize, &
        inembed, &
        istride, &
        idist, &
        onembed, &
        ostride, &
        odist, &
        HICFFT_FWD_TYPE, &
        iony*konx&
      )
      call hicfftSetAutoAllocation(plany, 0)

      call hicfftPlanMany( &
        planyi, &
        1, &
        fftsize, &
        onembed, &
        ostride, &
        odist, &
        inembed, &
        istride, &
        idist, &
        HICFFT_BWD_TYPE, &
        iony*konx&
      )
      call hicfftSetAutoAllocation(planyi, 0)

      ! Determine the workspace needed for FFTs and transposes
      max_worksize = -1
      
      call hicfftGetSize(planx, worksize)
      max_worksize = max(max_worksize, worksize)
      call hicfftGetSize(planxi, worksize)
      max_worksize = max(max_worksize, worksize)
      call hicfftGetSize(plany, worksize)
      max_worksize = max(max_worksize, worksize)
      call hicfftGetSize(planyi, worksize)
      max_worksize = max(max_worksize, worksize)
      
      ! max_worksize is in bytes, so convert it to number of elements by dividing by the size of a real number
      worksize = max_worksize / (storage_size(1._pois_r) / 8)

      worksize = max(worksize, sz)

      call allocate_workspace(int(worksize))

      call hicfftSetWorkArea(planx, workspace_0)
      call hicfftSetWorkArea(planxi, workspace_0)
      call hicfftSetWorkArea(plany, workspace_0)
      call hicfftSetWorkArea(planyi, workspace_0)

      if (nprocs == 1) then
        allocate(xyrt(2-ih:i1+ih,2-jh:j1+jh))
        allocate(d(2-ih:i1+ih,2-jh:j1+jh,kmax))
        ps = 2
        pe = i1
        qs = 2
        qe = j1
      else
        allocate(xyrt(iony,jonx))
        allocate(d(iony,jonx,kmax))
        ps = 1
        pe = iony
        qs = 1
        qe = jonx
      end if

      call init_factors(xyrt)

      norm_fac = 1 / real((itot*jtot))

      !$acc enter data copyin(xyrt, d)
      !$omp target enter data map(to:xyrt,d)

    end subroutine cufftinit

    !< Exit routine
    subroutine cufftexit(p, Fp, d, xyrt)

      implicit none

      real(pois_r), pointer :: p(:,:,:), Fp(:,:,:)
      real(pois_r), allocatable :: d(:,:,:), xyrt(:,:)


      deallocate(d, xyrt, p_halo, p_nohalo)

      nullify(p, Fp)

      call hicfftDestroy(planx)
      call hicfftDestroy(planxi)
      call hicfftDestroy(plany)
      call hicfftDestroy(planyi)
      
    end subroutine cufftexit

    subroutine init_factors(xyrt)
      implicit none

      real(pois_r), allocatable :: xyrt(:,:)
      real(pois_r) :: xrt(itot), yrt(jtot)
      integer :: iswap(itot), jswap(jtot)
      integer i,j,nh
      
      ! cuFFT orders the Fourier coefficients like this:
      !   
      !   r[0],i[0],r[1],i[1],r[2],i[2],...,r[n/2],i[n/2],r[n/2+1],i[n/2+1]
      ! 
      ! i[0] and i[n/2+1] are 0, so data is reordered like this:
      ! 
      !   r[0],r[n/2+1],r[1],i[1],...,r[n/2],i[n/2]
      !
      ! TODO: this needs to work for uneven number of grid points too

      do i=1,itot
        xrt(i) = -4.*dxi*dxi*(sin(float(i-1)*pi/itot))**2
      end do

      ! Swap order
      nh = (itot+1) / 2
      iswap(1) = 1
      iswap(2) = nh + (1-mod(itot,2))
      do i=2,itot-1
        if (i <= nh) then
          iswap(2*i-1) = i
        else
          iswap(itot-2*(i-(nh+1))-mod(itot,2)) = i+1
        end if
      end do
      xrt(:) = xrt(iswap(:))

      do j=1,jtot
        yrt(j) = -4.*dxi*dxi*(sin(float(j-1)*pi/jtot))**2
      end do

      ! Swap order
      nh = (jtot+1) / 2
      jswap(1) = 1
      jswap(2) = nh + (1-mod(jtot,2))
      do j=2,jtot-1
        if (j <= nh) then
          jswap(2*j-1) = j
        else
          jswap(jtot-2*(j-(nh+1))-mod(jtot,2)) = j+1
        end if
      end do
      yrt(:) = yrt(jswap(:))

      xyrt = 0

      if (nprocs == 1) then
        do j=2,j1
          do i=2,i1
            xyrt(i,j) = (xrt(i-1) + yrt(j-1))
          end do
        end do
      else
        do j = 1, jonx
          do i = 1, iony
            xyrt(i,j) = xrt(myidy*iony+i) + yrt(myidx*jonx+j)
          end do
        end do
      end if

    end subroutine init_factors

    !< Forward transforms 
    subroutine cufftf(p, Fp)

      real(pois_r), pointer :: p(:,:,:), Fp(:,:,:)
      integer :: i, j, k, ii

      call timer_tic('modcufft/cufftf', 1)

      call transposer%z_to_x(p, px, workspace_0)

      call hicfftExecForward(planx, px, px)

      call postprocess_f_fft(px, (/2*nphix, jmax, konx/), itot)
      call transposer%x_to_y(px, py, workspace_0)

      call hicfftExecForward(plany, py, py)

      call postprocess_f_fft(py, (/2*nphiy, konx, iony/), jtot)

      call transposer%y_to_z(py, Fp, workspace_0)

      call timer_toc('modcufft/cufftf')

    end subroutine cufftf

    !< Backward transforms
    subroutine cufftb(p, Fp)
      
      real(pois_r), pointer :: p(:,:,:), Fp(:,:,:)
      integer :: i, j, k, ii

      call timer_tic('modcufft/cufftb', 1)

      call transposer%z_to_y(Fp, py, workspace_0)
      call preprocess_b_fft(py, (/2*nphiy, konx, iony/), jtot)

      call hicfftExecBackward(planyi, py, py)

      call transposer%y_to_x(py, px, workspace_0)
      call preprocess_b_fft(px, (/2*nphix, jmax, konx/), itot)

      call hicfftExecBackward(planxi, px, px)

      call transposer%x_to_z(px, p, workspace_0)
      
      !$acc parallel loop collapse(3) default(present)
      !$omp target teams loop collapse(3) defaultmap(present:aggregate)&
      !$omp defaultmap(present:allocatable)
      do k=1,kmax
        do j=2,j1
          do i=2,i1
            p(i,j,k) = p(i,j,k) * norm_fac
          end do
        end do
      end do

      call timer_toc('modcufft/cufftb')

    end subroutine cufftb

    !> Postprocess signal after forward FFT
    subroutine postprocess_f_fft(arr, dim, len)
      implicit none

      real(pois_r), pointer, intent(inout) :: arr(:,:,:)
      integer, intent(in) :: dim(:)
      integer, intent(in) :: len

      integer :: j, k, sz_2, sz_3

      sz_2 = dim(2)
      sz_3 = dim(3)

      !$acc parallel loop collapse(2) default(present)
      !$omp target teams loop collapse(2) defaultmap(present:aggregate)&
      !$omp defaultmap(present:allocatable)
      do k = 1, sz_3
        do j = 1, sz_2
          arr(2,j,k) = arr(len+1,j,k)
        end do
      end do
    
    end subroutine postprocess_f_fft

    !< Preprocess signal before inverse FFT
    subroutine preprocess_b_fft(arr, dim, len)
      implicit none

      real(pois_r), pointer, intent(inout) :: arr(:,:,:)
      integer, intent(in) :: dim(:)
      integer, intent(in) :: len

      integer :: i, j, k, sz_2, sz_3

      sz_2 = dim(2)
      sz_3 = dim(3)

      !$acc parallel loop collapse(2) default(present)
      !$omp target teams loop collapse(2) defaultmap(present:aggregate)&
      !$omp defaultmap(present:allocatable)
      do k = 1, sz_3
        do j = 1, sz_2
          arr(len+1,j,k) = arr(2,j,k)
          arr(2,j,k) = 0.
        end do
      end do

    end subroutine preprocess_b_fft

#else
  contains

    subroutine cufftinit(p, Fp, d, xyrt, ps, pe, qs, qe)
      real(pois_r), pointer :: p(:,:,:)
      real(pois_r), pointer :: Fp(:,:,:)
      real(pois_r), allocatable :: d(:,:,:)
      real(pois_r), allocatable :: xyrt(:,:)
      integer, intent(out) :: ps, pe, qs, qe
      call error_and_exit()
      ps=0
      pe=0
      qs=0
      qe=0
      p => NULL()
      Fp => NULL()
      xyrt = 0 ! to suppress warning about unused argument.
      d = 0    ! These are not allocated, so would be bad to actually do
    end subroutine cufftinit

    subroutine cufftexit(p, Fp, d, xyrt)
      real(pois_r), pointer :: p(:,:,:)
      real(pois_r), pointer :: Fp(:,:,:)
      real(pois_r), allocatable :: d(:,:,:)
      real(pois_r), allocatable :: xyrt(:,:)
      call error_and_exit()
      p => NULL()
      Fp => NULL()
      xyrt = 0 ! to suppress warning about unused argument.
      d = 0    ! These are not allocated, so would be bad to actually do
    end subroutine cufftexit

    subroutine cufftf(p, Fp)
      real(pois_r), pointer :: p(:,:,:)
      real(pois_r), pointer :: Fp(:,:,:)
      call error_and_exit()
      p = 0 ! to suppress warning about unused argument.
      Fp = 0
    end subroutine cufftf

    subroutine cufftb(p, Fp)
      real(pois_r), pointer :: p(:,:,:)
      real(pois_r), pointer :: Fp(:,:,:)
      call error_and_exit()
      p = 0 ! to suppress warning about unused argument.
      Fp = 0
    end subroutine cufftb

    subroutine error_and_exit
      write(*,*) "DALES was compiled without GPU support, but cuFFT solver was selected"
      write(*,*) "Use another solver (solver_id), or compile with SYST=NV-OpenACC"
      call exit(-1)
    end subroutine

#endif

end module modcufft
