module modcufft

  implicit none

  contains
    !< Setup plan, workspace, etc
    subroutine cufftinit
    implicit none
    end subroutine cufftinit

    !< Exit routine
    subroutine cufftexit
    implicit none
    end subroutine cufftexit

    !< Forward transforms 
    subroutine cufftf
    implicit none
    end subroutine cufftf

    !< Backward transforms
    subroutine cufftb
    implicit none
    end subroutine cufftb

end module modcufft
