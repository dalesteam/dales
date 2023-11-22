module modcufft
  use, intrinsic :: iso_c_binding 

  use modmpi
  use modglobal, only: itot, jtot, kmax, i1, j1, &
                       imax, jmax, ih, jh, dxi, dyi, pi, ijtot
  use modprecision, only: pois_r

  implicit none

#if defined(_OPENACC)

  save
    real :: norm_fac !< Normalization factor
    integer :: istat !< cuFFT return status 

    real(pois_r), allocatable, target :: p_halo(:) !< Pressure with halos
    real(pois_r), allocatable, target :: p_nohalo(:)
    real(pois_r), pointer :: px(:,:,:), py(:,:,:)

    integer :: nphix, nphiy
    integer :: konx, kony, iony, jonx

    integer :: planx, planxi, plany, planyi !< Plan handles

    integer(int_ptr_kind()) :: worksize, max_worksize !< Size of the required workspace

  contains
    !< Setup plans, workspace, etc
    subroutine cufftinit(p, Fp, d, xyrt, ps, pe, qs, qe)
      use cufft
      use modgpu, only: workspace_0, allocate_workspace

      implicit none

      real(pois_r), pointer :: p(:,:,:) !< Pressure, spatial domain, with halos
      real(pois_r), pointer :: Fp(:,:,:) !< Pressure, spectral domain, with halos
      real(pois_r), allocatable :: xyrt(:,:) !< Array of eigenvalues
      real(pois_r), allocatable :: d(:,:,:)
      integer, intent(out) :: ps, pe, qs, qe

      integer(kind=8) :: sz
      integer :: fftsize, inembed, onembed, idist, odist, istride, ostride
      integer :: CUFFT_FWD_TYPE, CUFFT_BWD_TYPE

      ! Dimensions of the transposes
      ! For explanation of the variables, see modfftw.f90/fftwinit

      konx = kmax / nprocx
      if (mod(kmax, nprocx) > 0) then
        konx = konx + 1
      end if

      kony = kmax / nprocy
      if (mod(kmax, nprocy) > 0) then
        kony = kony + 1
      end if

      iony = itot / nprocy
      if (mod(itot, nprocy) > 0) then
        iony = iony + 1
      end if

      jonx = jtot / nprocx
      if (mod(jtot, nprocx) > 0) then
        jonx = jonx + 1
      end if

      ! Number of complex coefficients
      nphix = itot/2 + 1
      nphiy = jtot/2 + 1
      
      sz = max(imax * jmax * konx * nprocx, & ! z-aligned
               iony * jmax * konx * nprocy, & ! x-aligned
               iony * jonx * konx * nprocx)   ! y-aligned 

      ! Allocate memory for the pressure
      allocate(p_halo(1:(imax+2*ih)*(jmax+2*jh)*kmax))
      allocate(p_nohalo(kmax*(nphix*2)*(nphiy*2)))

      !$acc enter data create(p_halo, p_nohalo)

      p(2-ih:i1+ih,2-jh:j1+jh,1:kmax) => p_halo(1:(imax+2*ih)*(jmax+2*jh)*kmax) ! z-aligned
      px(1:nphix*2,1:jmax,1:konx) => p_nohalo(1:konx*jmax*(nphix*2)) ! x-aligned
      py(1:nphiy*2,1:konx,1:iony) => p_nohalo(1:konx*iony*(nphiy*2)) ! y-aligned

      if (nprocs == 1) then
        Fp(2-ih:i1+ih,2-jh:j1+jh,1:kmax) => p_halo(1:(imax+2*ih)*(jmax+2*jh)*kmax) ! z-aligned
      else
        Fp(1:iony,1:jonx,1:kmax) => p_halo(1:iony*jonx*kmax)
      end if
        

      ! Precision
#if POIS_PRECISION==32
      CUFFT_FWD_TYPE = CUFFT_R2C
      CUFFT_BWD_TYPE = CUFFT_C2R
#else
      CUFFT_FWD_TYPE = CUFFT_D2Z
      CUFFT_BWD_TYPE = CUFFT_Z2D
#endif

      ! x-direction
      fftsize = itot
      inembed = itot
      onembed = nphiy
      idist = nphix*2
      odist = nphix
      istride = 1
      ostride = 1

      istat = cufftSetAutoAllocation(planx, 0)
      istat = cufftPlanMany( &
        planx, &
        1, &
        fftsize, &
        inembed, &
        istride, &
        idist, &
        onembed, &
        ostride, &
        odist, &
        CUFFT_FWD_TYPE, &
        jmax*konx &
      )
      call check_exitcode(istat)

      istat = cufftSetAutoAllocation(planxi, 0)
      istat = cufftPlanMany( &
        planxi, &
        1, &
        fftsize, &
        onembed, &
        ostride, &
        odist, &
        inembed, &
        istride, &
        idist, &
        CUFFT_BWD_TYPE, &
        jmax*konx &
      )

      call check_exitcode(istat)

      ! y-direction

      fftsize = jtot
      inembed = jtot
      onembed = nphiy
      idist = nphiy*2
      odist = nphiy
      istride = 1
      ostride = 1

      istat = cufftSetAutoAllocation(plany, 0)
      istat = cufftPlanMany( &
        plany, &
        1, &
        fftsize, &
        inembed, &
        istride, &
        idist, &
        onembed, &
        ostride, &
        odist, &
        CUFFT_FWD_TYPE, &
        iony*konx&
      )

      call check_exitcode(istat)

      istat = cufftSetAutoAllocation(planyi, 0)
      istat = cufftPlanMany( &
        planyi, &
        1, &
        fftsize, &
        onembed, &
        ostride, &
        odist, &
        inembed, &
        istride, &
        idist, &
        CUFFT_BWD_TYPE, &
        iony*konx&
      )

      call check_exitcode(istat)

      ! Determine the workspace needed for FFTs and transposes
      max_worksize = -1
      
      istat = cufftGetSize(planx, worksize)
      max_worksize = max(max_worksize, worksize)
      istat = cufftGetSize(planxi, worksize)
      max_worksize = max(max_worksize, worksize)
      istat = cufftGetSize(plany, worksize)
      max_worksize = max(max_worksize, worksize)
      istat = cufftGetSize(planyi, worksize)
      max_worksize = max(max_worksize, worksize)
      
      ! max_worksize is in bytes, so convert it to number of elements by dividing by the size of a real number
      worksize = max_worksize / (storage_size(1._pois_r) / 8)

      worksize = max(worksize, ((nphix*2)*(nphiy*2)*kmax))

      call allocate_workspace(int(worksize))

      !$acc host_data use_device(workspace_0)
      istat = cufftSetWorkArea(planx, workspace_0)
      istat = cufftSetWorkArea(planxi, workspace_0)
      istat = cufftSetWorkArea(plany, workspace_0)
      istat = cufftSetWorkArea(planyi, workspace_0)
      !$acc end host_data

      call check_exitcode(istat)

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

    end subroutine cufftinit

    !< Exit routine
    subroutine cufftexit(p, Fp, d, xyrt)
      use cufft

      implicit none

      real(pois_r), pointer :: p(:,:,:), Fp(:,:,:)
      real(pois_r), allocatable :: d(:,:,:), xyrt(:,:,:)


      deallocate(d, xyrt, p_halo, p_nohalo)

      nullify(p, Fp)

      istat = cufftDestroy(planx)
      istat = cufftDestroy(planxi)
      istat = cufftDestroy(plany)
      istat = cufftDestroy(planyi)
      
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
      use cufft

      implicit none

      real(pois_r), pointer :: p(:,:,:), Fp(:,:,:)
      integer :: i, j, k, ii
      
      call transpose_a1(p, px)

      !$acc host_data use_device(px)
#if POIS_PRECISION==32
      istat = cufftExecR2C(planx, px, px)
#else
      istat = cufftExecD2Z(planx, px, px)
#endif
      !$acc end host_data
      
      call postprocess_f_fft(px, (/2*nphix, jmax, konx/), itot)
      call transpose_a2(px, py)
      
      !$acc host_data use_device(py)
#if POIS_PRECISION==32
      istat = cufftExecR2C(plany, py, py)
#else
      istat = cufftExecD2Z(plany, py, py)
#endif
      !$acc end host_data
      call postprocess_f_fft(py, (/2*nphiy, konx, iony/), jtot)

      call transpose_a3(py, Fp)

    end subroutine cufftf

    !< Backward transforms
    subroutine cufftb(p, Fp)
      use cufft

      implicit none
      
      real(pois_r), pointer :: p(:,:,:), Fp(:,:,:)
      integer :: i, j, k, ii

      call transpose_a3inv(py, Fp)
      call preprocess_b_fft(py, (/2*nphiy, konx, iony/), jtot)

      !$acc host_data use_device(py)
#if POIS_PRECISION==32
      istat = cufftExecC2R(planyi, py, py)
#else
      istat = cufftExecZ2D(planyi, py,  py)
#endif
      !$acc end host_data

      call check_exitcode(istat)
      call transpose_a2inv(px, py)
      call preprocess_b_fft(px, (/2*nphix, jmax, konx/), itot)

      !$acc host_data use_device(px)
#if POIS_PRECISION==32
      istat = cufftExecC2R(planxi, px, px)
#else
      istat = cufftExecZ2D(planxi, px, px)
#endif
      !$acc end host_data

      call check_exitcode(istat)
      call transpose_a1inv(p, px)
      
      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=2,j1
          do i=2,i1
            p(i,j,k) = p(i,j,k) * norm_fac
          end do
        end do
      end do

    end subroutine cufftb

    subroutine transpose_a1(p, px)
      use modgpu, only: workspace_0, workspace_1
      implicit none

      real(pois_r), pointer, intent(in) :: p(:,:,:)
      real(pois_r), pointer, intent(out) :: px(:,:,:)

      integer :: i, j, k, n, ii

      if (nprocs == 1) then
        !$acc parallel loop collapse(3) default(present)
        do k=1,kmax
          do j=1,jtot
            do i=1,itot
              px(i,j,k) = p(i+1,j+1,k)
            end do
          end do
        end do
      else
        ! TODO: clean this one up
        !$acc parallel loop collapse(4) default(present) private(ii)
        do n = 0, nprocx-1
          do k = 1, konx
            do j = 2, j1
              do i = 2, i1
                ii = (i-1) + (j-2)*imax + (k-1)*imax*jmax + n*imax*jmax*konx
                if (k+n*konx <= kmax) workspace_0(ii) = p(i,j,k+n*konx) 
              end do
            end do
          end do
        end do

        !$acc host_data use_device(workspace_0, workspace_1)
        call D_MPI_ALLTOALL(workspace_0, imax*jmax*konx, &
                            workspace_1, imax*jmax*konx, &
                            commrow, mpierr)
        !$acc end host_data

        !$acc parallel loop collapse(4) default(present) private(ii)
        do n = 0, nprocx-1
          do k = 1, konx
            do j = 1, jmax
              do i = 1, imax
                ii = i + (j-1)*imax + (k-1)*imax*jmax + n*imax*jmax*konx
                px(i+n*imax,j,k) = workspace_1(ii)
              end do
            end do
          end do
        end do
      end if

    end subroutine transpose_a1

    subroutine transpose_a1inv(p, px)
      use modgpu, only: workspace_0, workspace_1
      implicit none

      real(pois_r), pointer, intent(in) :: px(:,:,:)
      real(pois_r), pointer, intent(out) :: p(:,:,:)

      integer :: i, j, k, n, ii

      if (nprocs == 1) then
        !$acc parallel loop collapse(3) default(present)
        do k = 1, kmax
          do j = 1, jtot
            do i = 1, itot
              p(i+1,j+1,k) = px(i,j,k)
            end do
          end do
        end do
      else
        !$acc parallel loop collapse(4) default(present) private(ii)
        do n = 0, nprocx-1
          do k = 1, konx
            do j = 1, jmax
              do i = 1, imax
                ii = i + (j-1)*imax + (k-1)*imax*jmax + n*imax*jmax*konx
                workspace_0(ii) = px(i+n*imax,j,k)
              end do
            end do
          end do
        end do

        !$acc host_data use_device(workspace_0, workspace_1)
        call D_MPI_ALLTOALL(workspace_0, imax*jmax*konx, &
                            workspace_1, imax*jmax*konx, &
                            commrow, mpierr)
        !$acc end host_data

        !$acc parallel loop collapse(4) default(present) private(ii)
        do n = 0, nprocx-1
          do k = 1, konx
            do j = 2, j1
              do i = 2, i1
                ii = (i-1) + (j-2)*imax + (k-1)*imax*jmax + n*imax*jmax*konx
                if (k+n*konx <= kmax) p(i,j,k+n*konx) = workspace_1(ii)
              end do
            end do
          end do
        end do
      end if

    end subroutine transpose_a1inv

    subroutine transpose_a2(px, py)
      use modgpu, only: workspace_0, workspace_1
      implicit none

      real(pois_r), pointer, intent(in) :: px(:,:,:)
      real(pois_r), pointer, intent(out) :: py(:,:,:)

      integer :: i, j, k, n, ii

      if (nprocs == 1) then
        !$acc parallel loop collapse(3) default(present) private(ii)
        do k = 1, kmax
          do j = 1, jtot
            do i = 1, itot
              ii = i + (j-1)*itot + (k-1)*itot*jtot
              workspace_0(ii) = px(i,j,k)
            end do
          end do
        end do

        !$acc parallel loop collapse(3) default(present) private(ii)
        do k = 1, kmax
          do j = 1, jtot
           do i = 1, itot
              ii = i + (j-1)*itot + (k-1)*itot*jtot
              py(j,k,i) = workspace_0(ii)
            end do
          end do
        end do
      else
        !$acc parallel loop collapse(4) default(present) private(ii)
        do n = 0, nprocy-1
          do k = 1, konx
            do j = 1, jmax
              do i = 1, iony
                ii = i + (j-1)*iony + (k-1)*iony*jmax + n*iony*jmax*konx
                if (i <= itot) workspace_0(ii) = px(i+n*iony,j,k)
              end do
            end do
          end do
        end do

        !$acc host_data use_device(workspace_0, workspace_1)
        call D_MPI_ALLTOALL(workspace_0, iony*jmax*konx, &
                            workspace_1, iony*jmax*konx, &
                            commcol, mpierr)
        !$acc end host_data

        !$acc parallel loop collapse(4) default(present) private(ii)
        do n = 0, nprocy-1
          do k = 1, konx
            do i = 1, iony
              do j = 1, jmax
                ii = i + (j-1)*iony + (k-1)*iony*jmax + n*iony*jmax*konx
                py(j+n*jmax,k,i) = workspace_1(ii)
              end do
            end do
          end do
        end do
                
      end if

    end subroutine transpose_a2

    subroutine transpose_a2inv(px, py)
      use modgpu, only: workspace_0, workspace_1
      implicit none

      real(pois_r), pointer, intent(in) :: px(:,:,:)
      real(pois_r), pointer, intent(out) :: py(:,:,:)

      integer :: i, j, k, n, ii

      if (nprocs == 1) then
        !$acc parallel loop collapse(3) default(present) private(ii)
        do k = 1, kmax
          do j = 1, jtot
            do i = 1, itot
              ii = j + (i-1)*jtot + (k-1)*itot*jtot
              workspace_0(ii) = py(j,k,i)
            end do
          end do
        end do

        !$acc parallel loop collapse(3) default(present) private(ii)
        do k = 1, kmax
          do j = 1, jtot
            do i = 1, itot
              ii = j + (i-1)*jtot + (k-1)*itot*jtot
              px(i,j,k) = workspace_0(ii)
            end do
          end do
        end do
      else
        !$acc parallel loop collapse(4) default(present) private(ii)
        do n = 0, nprocy-1
          do k = 1, konx
            do i = 1, iony
              do j = 1, jmax
                ii = i + (j-1)*iony + (k-1)*iony*jmax + n*iony*jmax*konx
                workspace_0(ii) = py(j+n*jmax,k,i)
              end do
            end do
          end do
        end do

        !$acc host_data use_device(workspace_0, workspace_1)
        call D_MPI_ALLTOALL(workspace_0, iony*jmax*konx, &
                            workspace_1, iony*jmax*konx, &
                            commcol, mpierr)
        !$acc end host_data

        !$acc parallel loop collapse(4) default(present) private(ii)
        do n = 0, nprocy-1
          do k = 1, konx
            do j = 1, jmax
              do i = 1, iony
                ii = i + (j-1)*iony + (k-1)*iony*jmax + n*iony*jmax*konx
                if (i+n*iony <= itot) px(i+n*iony,j,k) = workspace_1(ii)
              end do
            end do
          end do
        end do
      end if

    end subroutine transpose_a2inv

    subroutine transpose_a3(py, Fp)
      use modgpu, only: workspace_0, workspace_1
      implicit none

      real(pois_r), pointer, intent(in) :: py(:,:,:)
      real(pois_r), pointer, intent(out) :: Fp(:,:,:)

      integer :: i, j, k, n, ii

      if (nprocs == 1) then
        !$acc parallel loop collapse(3) default(present)
        do k = 1, kmax
          do j = 1, jtot
            do i = 1, itot
              Fp(i+1,j+1,k) = py(j,k,i)
            end do
          end do
        end do
      else
        !$acc parallel loop collapse(4) default(present) private(ii)
        do n = 0, nprocx-1
          do k = 1, konx
            do i = 1, iony
              do j = 1, jonx
                ii = j + (i-1)*jonx + (k-1)*iony*jonx + n*iony*jonx*konx
                if (j+n*jonx <= jtot) workspace_0(ii) = py(j+n*jonx,k,i)
              end do
            end do
          end do
        end do

        !$acc host_data use_device(workspace_0, workspace_1)
        call D_MPI_ALLTOALL(workspace_0, iony*jonx*konx, &
                            workspace_1, iony*jonx*konx, &
                            commrow, mpierr)
        !$acc end host_data

        !$acc parallel loop collapse(4) default(present) private(ii)
        do n = 0, nprocx-1
          do k = 1, konx
            do j = 1, jonx
              do i = 1, iony
                ii = j + (i-1)*jonx + (k-1)*iony*jonx + n*iony*jonx*konx
                if (k+n*konx <= kmax) Fp(i,j,k+n*konx) = workspace_1(ii)
              end do
            end do
          end do
        end do

      end if

    end subroutine transpose_a3

    subroutine transpose_a3inv(py, Fp)
      use modgpu, only: workspace_0, workspace_1
      implicit none

      real(pois_r), pointer, intent(in) :: Fp(:,:,:)
      real(pois_r), pointer, intent(out) :: py(:,:,:)

      integer :: i, j, k, n, ii

      if (nprocs == 1) then
        !$acc parallel loop collapse(3) default(present)
        do k=1,kmax
          do j=1,jtot
            do i=1,itot
              py(j,k,i) = Fp(i+1,j+1,k)
            end do
          end do
        end do
      else
        !$acc parallel loop collapse(4) default(present) private(ii)
        do n = 0, nprocx-1
          do k = 1, konx
            do j = 1, jonx
              do i = 1, iony
                ii = j + (i-1)*jonx + (k-1)*iony*jonx + n*iony*jonx*konx
                if (k+n*konx <= kmax) workspace_0(ii) = Fp(i,j,k+n*konx)
              end do
            end do
          end do
        end do

        !$acc host_data use_device(workspace_0, workspace_1)
        call D_MPI_ALLTOALL(workspace_0, iony*jonx*konx, &
                            workspace_1, iony*jonx*konx, &
                            commrow, mpierr)
        !$acc end host_data

        !$acc parallel loop collapse(4) default(present) private(ii)
        do n = 0, nprocx-1
          do k = 1, konx
            do i = 1, iony
              do j = 1, jonx
                ii = j + (i-1)*jonx + (k-1)*iony*jonx + n*iony*jonx*konx
                if (j+n*jonx <= jtot) py(j+n*jonx,k,i) = workspace_1(ii)
              end do
            end do
          end do
        end do
      end if

    end subroutine transpose_a3inv

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
      do k = 1, sz_3
        do j = 1, sz_2
          arr(len+1,j,k) = arr(2,j,k)
          arr(2,j,k) = 0.
        end do
      end do

    end subroutine preprocess_b_fft

    !< Checks the exitcode of cuFFT calls
    subroutine check_exitcode(istat)
      implicit none
      integer, intent(in) :: istat
      
      if ( istat /= 0 ) then
        write(*,*) "cuFFT returned nonzero exitcode: ", istat
        stop
      end if

    end subroutine check_exitcode

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
    end subroutine cufftinit

    subroutine cufftexit(p, Fp, d, xyrt)
      real(pois_r), pointer :: p(:,:,:)
      real(pois_r), pointer :: Fp(:,:,:)
      real(pois_r), allocatable :: d(:,:,:)
      real(pois_r), allocatable :: xyrt(:,:)
      call error_and_exit()
    end subroutine cufftexit

    subroutine cufftf(p, Fp)
      real(pois_r), pointer :: p(:,:,:)
      real(pois_r), pointer :: Fp(:,:,:)
      call error_and_exit()
    end subroutine cufftf

    subroutine cufftb(p, Fp)
      real(pois_r), pointer :: p(:,:,:)
      real(pois_r), pointer :: Fp(:,:,:)
      call error_and_exit()
    end subroutine cufftb

    subroutine error_and_exit
      write(*,*) "DALES was compiled without GPU support, but cuFFT solver was selected"
      write(*,*) "Use another solver (solver_id), or compile with SYST=NV-OpenACC"
      call exit(-1)
    end subroutine

#endif

end module modcufft
