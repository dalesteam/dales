module modcufft
  use, intrinsic :: iso_c_binding 

  use modmpi, only: nprocx, nprocy 
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

    integer :: planx, planxi, plany, planyi !< Plan handles

    integer(int_ptr_kind()) :: worksize, max_worksize !< Size of the required workspace

  contains
    !< Setup plans, workspace, etc
    subroutine cufftinit(p, Fp, d, xyrt, ps, pe, qs, qe)
      use cufft
      use modgpu, only: workspace, allocate_workspace

      implicit none

      real(pois_r), pointer :: p(:,:,:) !< Pressure, spatial domain, with halos
      real(pois_r), pointer :: Fp(:,:,:) !< Pressure, spectral domain, with halos
      real(pois_r), allocatable :: xyrt(:,:) !< Array of eigenvalues
      real(pois_r), allocatable :: d(:,:,:)
      integer, intent(out) :: ps, pe, qs, qe

      integer :: itot12, jtot12
      integer :: fftsize, inembed, onembed, idist, odist, istride, ostride
      integer :: CUFFT_FWD_TYPE, CUFFT_BWD_TYPE
      
      ! Setup pressure array
      allocate(p_halo(1:(imax+2*ih)*(jmax+2*jh)*kmax))
      !$acc enter data create(p_halo)
      p(2-ih:i1+ih,2-jh:j1+jh,1:kmax) => p_halo(1:(imax+2*ih)*(jmax+2*jh)*kmax)
      Fp(2-ih:i1+ih,2-jh:j1+jh,1:kmax) => p_halo(1:(imax+2*ih)*(jmax+2*jh)*kmax)

      ! Number of complex coefficients
      itot12 = itot/2 + 1
      jtot12 = jtot/2 + 1

      ! Allocate workspace
      ! Keep in mind that cuFFT does real to complex transforms, where one complex number
      ! consists of two real numbers
      allocate(p_nohalo(kmax*(itot12*2)*(jtot12*2)))
      !$acc enter data create(p_nohalo)

      px(1:itot12*2,1:jtot,1:kmax) => p_nohalo(1:kmax*jtot*(itot12*2))
      py(1:jtot12*2,1:itot,1:kmax) => p_nohalo(1:kmax*itot*(jtot12*2))

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
      inembed = itot12
      onembed = itot12
      idist = itot12*2
      odist = itot12
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
        jtot*kmax &
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
        jtot*kmax &
      )

      call check_exitcode(istat)

      ! y-direction

      fftsize = jtot
      inembed = jtot
      onembed = jtot12
      idist = jtot12*2
      odist = jtot12
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
        itot*kmax &
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
        itot*kmax &
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

      worksize = max(worksize, ((itot12*2)*(jtot12*2)*kmax))

      call allocate_workspace(int(worksize))

      !$acc host_data use_device(workspace)
      istat = cufftSetWorkArea(planx, workspace)
      istat = cufftSetWorkArea(planxi, workspace)
      istat = cufftSetWorkArea(plany, workspace)
      istat = cufftSetWorkArea(planyi, workspace)
      !$acc end host_data

      call check_exitcode(istat)

      allocate(xyrt(2-ih:i1+ih,2-jh:j1+jh))
      allocate(d(2-ih:i1+ih,2-jh:j1+jh,kmax))

      call init_factors(xyrt)

      ps = 2
      pe = i1
      qs = 2
      qe = j1

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

      do j=2,j1
        do i=2,i1
          xyrt(i,j) = (xrt(i-1) + yrt(j-1))
        end do
      end do

    end subroutine init_factors

    !< Forward transforms 
    subroutine cufftf(p, Fp)
      use cufft
      use modgpu, only: workspace

      implicit none

      real(pois_r), pointer :: p(:,:,:), Fp(:,:,:)
      integer :: i, j, k, ii
      
      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            px(i,j,k) = p(i+1,j+1,k)
          end do
        end do
      end do

      !$acc host_data use_device(px)
#if POIS_PRECISION==32
      istat = cufftExecR2C(planx, px, px)
#else
      istat = cufftExecD2Z(planx, px, px)
#endif
      !$acc end host_data
      
      ! Reorder
      !$acc parallel loop collapse(2) default(present)
      do k=1,kmax
        do j=1,jtot
          px(2,j,k) = px(itot+1,j,k)
        end do
      end do

      ! Transpose to y-contiguous
      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            ii = i + ((j-1)*2*(itot/2+1)) + ((k-1)*2*(itot/2+1)*2*(jtot/2+1))
            workspace(ii) = px(i,j,k)
          end do
        end do
      end do

      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
         do i=1,itot
            ii = i + ((j-1)*2*(itot/2+1)) + ((k-1)*2*(itot/2+1)*2*(jtot/2+1))
            py(j,i,k) = workspace(ii)
          end do
        end do
      end do
      
      !$acc host_data use_device(py)
#if POIS_PRECISION==32
      istat = cufftExecR2C(plany, py, py)
#else
      istat = cufftExecD2Z(plany, py, py)
#endif
      !$acc end host_data

      ! Reorder
      !$acc parallel loop collapse(2) default(present)
      do k=1,kmax
        do i=1,itot
          py(2,i,k) = py(jtot+1,i,k)
        end do
      end do

      ! Transpose
      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            ii = j + ((i-1)*2*(jtot/2+1)) + ((k-1)*2*(itot/2+1)*2*(jtot/2+1))
            workspace(ii) = py(j,i,k)
          end do
        end do
      end do

      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            ii = j + ((i-1)*2*(jtot/2+1)) + ((k-1)*2*(itot/2+1)*2*(jtot/2+1))
            px(i,j,k) = workspace(ii)
          end do
        end do
      end do

      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            p(i+1,j+1,k) = px(i,j,k)
          end do
        end do
      end do

    end subroutine cufftf

    !< Backward transforms
    subroutine cufftb(p, Fp)
      use cufft
      use modgpu, only: workspace

      implicit none
      
      real(pois_r), pointer :: p(:,:,:), Fp(:,:,:)
      integer :: i, j, k, ii
      
      ! 1. Fill in workspace
      ! 2. Reorder data to original cuFFT format
      ! 3. Backwards transform
      ! 4. Transpose
      ! 5. Reorder data
      ! 6. Backwards transform
      ! 7. Transpose
      ! 8. Fill in pressure array
      
      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            px(i,j,k) = Fp(i+1,j+1,k)
          end do
        end do
      end do

      !$acc parallel loop collapse(2) default(present)
      do k=1,kmax
        do j=1,jtot
          px(itot+1,j,k) = px(2,j,k)
          px(2,j,k) = 0.
        end do
      end do

      !$acc host_data use_device(px)
#if POIS_PRECISION==32
      istat = cufftExecC2R(planxi, px, px)
#else
      istat = cufftExecZ2D(planxi, px, px)
#endif
      !$acc end host_data

      call check_exitcode(istat)

      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            ii = i + ((j-1)*2*(itot/2+1)) + ((k-1)*2*(itot/2+1)*2*(jtot/2+1))
            workspace(ii) = px(i,j,k)
          end do
        end do
      end do

      ! TODO: Hier klopt natuurlijk helemaal niets van
      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            ii = i + ((j-1)*2*(itot/2+1)) + ((k-1)*2*(itot/2+1)*2*(jtot/2+1))
            py(j,i,k) = workspace(ii)
          end do
        end do
      end do

      !$acc parallel loop collapse(2) default(present)
      do k=1,kmax
        do i=1,itot
          py(jtot+1,i,k) = py(2,i,k)
          py(2,i,k) = 0
        end do
      end do

      !$acc host_data use_device(py)
#if POIS_PRECISION==32
      istat = cufftExecC2R(planyi, py, py)
#else
      istat = cufftExecZ2D(planyi, py,  py)
#endif
      !$acc end host_data

      call check_exitcode(istat)

      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
         do i=1,itot
            ii = j + ((i-1)*2*(jtot/2+1)) + ((k-1)*2*(itot/2+1)*2*(jtot/2+1))
            workspace(ii) = py(j,i,k)
          end do
        end do
      end do
      
      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            ii = j + ((i-1)*2*(jtot/2+1)) + ((k-1)*2*(itot/2+1)*2*(jtot/2+1))
            px(i,j,k) = workspace(ii)
          end do
        end do
      end do

      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            px(i,j,k) = px(i,j,k) * norm_fac
          end do
        end do
      end do

      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            p(i+1,j+1,k) = px(i,j,k)
          end do
        end do
      end do
      
    end subroutine cufftb

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
