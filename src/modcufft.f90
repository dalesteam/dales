module modcufft
  use, intrinsic :: iso_c_binding 

  use modmpi, only: nprocx, nprocy 
  use modglobal, only: itot, jtot, ktot, kmax, &
                       ih, jh, dxi, dyi, pi, ijtot
  use modprecision, only: pois_r

  implicit none

  save
    integer :: method

    real :: norm_fac !< Normalization factor

    integer :: istat !< cuFFT return status 
    
    ! 2D transforms
    integer :: planxy, planxyi !< Plan handles
    integer(int_ptr_kind()) :: worksize !< Pointer to the size of the required workspace
    real(pois_r), allocatable, target :: p_in_h(:,:,:) !< Pressure, with halo
    real(pois_r), allocatable :: p_in_nh(:,:,:) !< Pressure, no halo
    real(pois_r), allocatable :: p_out(:,:,:) !< Fourier coefficients

  contains
    !< Setup plans, workspace, etc
    subroutine cufftinit(p, Fp, d, xyrt, ps, pe, qs, qe)
      use cufft

      implicit none

      real(pois_r), pointer :: p(:,:,:) !< Pressure, spatial domain
      real(pois_r), pointer :: Fp(:,:,:) !< Pressure, spectral domain
      real(pois_r), allocatable :: xyrt(:,:) !< Array of eigenvalues
      real(pois_r), allocatable :: d(:,:,:)
      integer, intent(out) :: ps, pe, qs, qe

      integer :: jtot12
      integer :: idist, odist, istride, ostride
      integer :: inembed(2), onembed(2)

      if (nprocx == 1 .and. nprocy == 1) then
        method = 2 ! Single GPU, can do 2D transforms
      else
        method = 1 ! Multi GPU, do 1D transforms + transposes
      end if
      
      if (method == 2) then
        call setup_plan_2D(p, Fp, d, xyrt, ps, pe, qs, qe)
      else
        stop
      end if

      call init_factors(xyrt)

      norm_fac = 1/sqrt(ijtot)

    end subroutine cufftinit

    !< Exit routine
    subroutine cufftexit(p, Fp, d, xyrt)
      use cufft

      implicit none

      real(pois_r), pointer :: p(:,:,:), Fp(:,:,:)
      real(pois_r), allocatable :: d(:,:,:), xyrt(:,:,:)

      deallocate(d, xyrt)

      nullify(p, Fp)

      if (method == 2) then
        istat = cufftDestroy(planxy)
        istat = cufftDestroy(planxyi)
      end if
      
    end subroutine cufftexit

    subroutine init_factors(xyrt)
      implicit none

      real(pois_r), allocatable :: xyrt(:,:)
      real(pois_r) :: xrt(itot), yrt(jtot)

      integer i,j,ii,jj
      
      ! cuFFT orders the Fourier coefficients like this:
      !   
      !   r[0],r[1],i[1],r[2],i[2],...,r[n/2],i[n/2],r[n/2+1]
      ! 
      ! So we order the eigenvalues in the same way here
      !
      ! TODO: this needs to work for uneven number of grid points too

      ! x direction
      xrt(1) = 0
      do i=2,(itot/2)
        ii = (i-2)*2+2
        xrt(ii) = -4.*dxi*dxi*(sin((i-1)*pi/itot))**2 
        xrt(ii+1) = xrt(ii)
      end do
      if (mod(itot,2) == 0) then
        xrt(itot) = -4.*dxi*dxi
      end if

      ! y direction
      yrt(1) = 0
      do j=2,(jtot/2)
        jj = (j-2)*2+2
        yrt(jj) = -4.*dyi*dyi*(sin((j-1)*pi/jtot))**2
        yrt(jj+1) = yrt(jj) 
      end do
      if (mod(jtot,2) == 0) then
        yrt(jtot) = -4.*dyi*dyi
      end if

      if (method == 2) then
        do j=2,jtot
          do i=2,itot
            xyrt(i,j) = (xrt(i-1)+yrt(j-1))
          end do
        end do
      else
        stop "init_factors: illegal method"
      end if

    end subroutine init_factors

    !< Forward transforms 
    subroutine cufftf(p, Fp)
      use cufft

      implicit none

      real(pois_r), pointer :: p(:,:,:), Fp(:,:,:)
      integer :: i, j, k

      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            p_in_nh(i,j,k) = p_in_h(i+1,j+1,k)
          end do
        end do
      end do


      if (method == 2) then
        !$acc host_data use_device(p_in_nh, p_out)
        istat = cufftExecD2Z(planxy, p_in_nh, p_out)
        !$acc end host_data
        call check_exitcode(istat)
      else
        stop
      end if

      ! Normalization
      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            Fp(i,j,k) = Fp(i,j,k) * norm_fac
          end do
        end do
      end do

    end subroutine cufftf

    !< Backward transforms
    subroutine cufftb(p, Fp)
      use cufft

      implicit none
      
      real(pois_r), pointer :: p(:,:,:), Fp(:,:,:)
      integer :: i, j, k

      ! Normalization
      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            Fp(i,j,k) = Fp(i,j,k) * norm_fac
          end do
        end do
      end do

      if (method == 2) then
        !$acc host_data use_device(p_in_nh, p_out)
        istat = cufftExecZ2D(planxyi, p_out, p_in_nh)
        !$acc end host_data
        call check_exitcode(istat)
      else
        stop
      end if

      !$acc parallel loop collapse(3) default(present)
      do k=1,kmax
        do j=1,jtot
          do i=1,itot
            p_in_h(i+1,j+1,k) = p_in_nh(i,j,k)
          end do
        end do
      end do


    end subroutine cufftb

    subroutine setup_plan_2D(p, Fp, d, xyrt, ps, pe, qs, qe)
      use cufft

      implicit none

      real(pois_r), pointer :: p(:,:,:)
      !double complex, pointer :: Fp(:,:,:)
      real(pois_r), pointer :: Fp(:,:,:)
      real(pois_r), allocatable :: xyrt(:,:), d(:,:,:)
      integer, intent(out) :: ps, pe, qs, qe

      integer :: fftsize(2)
      integer :: jtot12
      integer :: inembed(2), onembed(2), idist, odist, istride, ostride

      ! Setup arrays for pressure fluctuations
      allocate(p_in_h(2-ih:itot+ih,2-jh:jtot+jh,1:kmax))
      allocate(p_in_nh(1:itot,1:jtot,1:kmax))

      p => p_in_h

      ! Allocate complex array for the Fourier coefficients
      ! From cuFFT documentation: the output of a 2D FFT of N1*N2 values has
      ! N1*(floor(N2/2)+1) elements 
      !jtot12 = jtot/2 + 1
      allocate(p_out(itot,jtot,kmax))
        
      ! Data layout
      fftsize = (/ itot, jtot /)
      inembed = (/ itot, jtot /)
      onembed = (/ itot, jtot /)
      idist = itot*jtot
      odist = itot*jtot
      istride = 1
      ostride = 1

      ! Forward transforms 
      istat = cufftCreate(planxy)
      istat = cufftPlanMany(&
        planxy, &    ! Plan handle
        2, &         ! Dimensionality
        fftsize, &   ! Sizes
        inembed, &   ! inembed
        istride, &   ! istride
        idist, &     ! idist
        onembed, &   ! onembed
        ostride, &   ! ostride
        odist, &     ! odist
        CUFFT_D2Z, & ! Double real to double complex transform
        kmax &       ! Batch size
      )
      
      call check_exitcode(istat)

      ! Backward transforms
      istat = cufftCreate(planxyi)
      istat = cufftPlanMany(&
        planxyi, &    ! Plan handle
        2, &         ! Dimensionality
        fftsize, &   ! Sizes
        onembed, &   ! inembed
        ostride, &   ! istride
        odist, &     ! idist
        inembed, &   ! onembed
        istride, &   ! ostride
        odist, &     ! odist
        CUFFT_Z2D, & ! Double real to double complex transform
        kmax &       ! Batch size
      )

      call check_exitcode(istat)

      allocate(xyrt(2-ih:itot+ih,2-jh:jtot+jh))
      allocate(d(2-ih:itot+ih,2-jh:jtot+jh,1:kmax))

      ps = 2
      pe = itot
      qs = 2
      qe = jtot

    end subroutine setup_plan_2D

    !< Checks the exitcode of cuFFT calls
    subroutine check_exitcode(istat)
      implicit none
      integer, intent(in) :: istat
      
      if ( istat /= 0 ) then
        write(*,*) "cuFFT returned nonzero exitcode: ", istat
        stop
      end if

    end subroutine check_exitcode

end module modcufft
