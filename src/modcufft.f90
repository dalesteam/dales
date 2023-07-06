module modcufft
  use, intrinsic :: iso_c_binding 

  use modmpi, only: nprocx, nprocy 
  use modglobal, only: itot, jtot, ktot

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

      end if      

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
