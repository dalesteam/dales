module modcufft
  use, intrinsic :: iso_c_binding 

  use modmpi, only: nprocx, nprocy 
  use modglobal, only: itot, jtot, ktot, i1, j1, kmax, &
                       ih, jh, dxi, dyi, pi
  use modprecision, only: pois_r

  implicit none

  save
    integer :: method

    integer, allocatable, dimension(:) :: fftsize

    integer :: istat !< cuFFT return status 
    
    ! 2D transforms
    integer :: planxy, planxyi !< Plan handles
    integer(int_ptr_kind()) :: worksize !< Pointer to the size of the required workspace

  contains
    !< Setup plans, workspace, etc
    subroutine cufftinit
      use cufft

      implicit none

      real(pois_r), allocatable :: xyrt(:,:) !< Array of eigenvalues

      if (nprocx == 1 .and. nprocy == 1) then
        method = 2 ! Single GPU, can do 2D transforms
      else
        method = 1 ! Multi GPU, do 1D transforms + transposes
      end if
      
      if (method == 2) then
        allocate(fftsize(2))
        
        fftsize(1) = itot ! Size of the x-dimension
        fftsize(2) = jtot ! Size of the y-dimension

        ! Forward transforms 

        istat = cufftCreate(planxy)
        
        ! TODO: disable this for cuDecomp
        istat = cufftSetAutoAllocation(planxy, 1)

        istat = cufftMakePlanMany(&
          planxy, &    ! Plan handle
          2, &         ! Dimensionality
          fftsize, &   ! Sizes
          null(), &    ! inembed
          1, &         ! istride
          1, &         ! idist
          null(), &    ! onembed
          1, &         ! ostride
          1, &         ! odist
          CUFFT_D2Z, & ! Double real to double complex transform
          ktot, &      ! Batch size
          worksize &   ! Size of temporary workspace
        )
      
        call check_exitcode(istat)

        ! Backward transforms
        
        istat = cufftCreate(planxyi)
        
        ! TODO: disable this for cuDecomp
        istat = cufftSetAutoAllocation(planxyi, 1)

        istat = cufftMakePlanMany(&
          planxyi, &   ! Plan handle
          2, &         ! Dimensionality
          fftsize, &   ! Sizes
          null(), &    ! inembed
          1, &         ! istride
          1, &         ! idist
          null(), &    ! onembed
          1, &         ! ostride
          1, &         ! odist
          CUFFT_Z2D, & ! Double complex to double real transform
          ktot, &      ! Batch size
          worksize &   ! Size of temporary work space 
        )

        call check_exitcode(istat)

        allocate(xyrt(2-ih:i1+ih,2-jh:j1+jh))

      end if      


      call init_factors(xyrt)

    end subroutine cufftinit

    !< Exit routine
    subroutine cufftexit
      use cufft

      implicit none

      deallocate(fftsize)

      if (method == 2) then
        istat = cufftDestroy(planxy)
        istat = cufftDestroy(planxyi)
      end if
      
    end subroutine cufftexit

    subroutine init_factors(xyrt)
      implicit none

      real(pois_r), allocatable :: xyrt(:,:)
      real(pois_r) :: xrt(itot), yrt(jtot)

      integer i,j
      
      ! x direction, use FFTW order as in modfftw.f90
      xrt(1) = 0
      do i=2,(itot/2)
        xrt(i) = -4.*dxi*dxi*(sin((i-1)*pi/itot))**2 
        xrt(itot-i+2) = xrt(i)
      end do
      if (mod(itot,2) == 0) then
        xrt(1+itot/2) = -4.*dxi*dxi
      end if

      ! y direction
      yrt(1) = 0
      do j=2,(jtot/2)
        yrt(j) = -4.*dyi*dyi*(sin((j-1)*pi/jtot))**2
        yrt(jtot-j+2) = yrt(j) 
      end do
      if (mod(jtot,2) == 0) then
        yrt(1+jtot/2) = -4*dyi*dyi
      end if

      if (method == 2) then
        do j=2,j1
          do i=2,i1
            xyrt(i,j) = (xrt(i-1)+yrt(j-1))
          end do
        end do
      else
        stop "init_factors: illegal method"
      end if

    end subroutine init_factors

    !< Forward transforms 
    subroutine cufftf
    implicit none
    end subroutine cufftf

    !< Backward transforms
    subroutine cufftb
    implicit none
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

end module modcufft
