!###############################################################################
!
! WARNING: do not confuse this file with go.F90 from LOTOS-EUROS!
! This is an interface between le_drydepos_gas_depac.f90 and the rest of DALES.
! The DEPAC routines use 'gol', 'goErr' and 'goPr' that were defined in the GO
! routines for LOTOS-EUROS. This file contains an implementation of these 
! functions. This way, le_drydepos_depac.f90 did not have to be adapted to DALES
! and it is directly transferrable between LE and DALES.
!
! L. Geers, June 2022
!###############################################################################

module GO
    implicit none

    public :: gol, goPr, goErr

    ! buffer size:
    integer, parameter :: lengol = 1024

    ! buffer for standard output
    character(len=lengol) :: gol

contains
    !> Print ordinary message to stdout
    !!
    !! Subroutine prints `gol` buffer to stdout.
    subroutine goPr
        ! write buffer to standard output (stdout)
        write(*, '(a)') trim(gol)
        ! clear buffer
        gol = ''
    end subroutine goPr

    !> Print error message. 
    !!
    !! Currently, errors are also sent to stdout. Need to adapt this to stderr.
    subroutine goErr
        call goPr
    end subroutine goErr

    !> Transform lowercase letters to uppercase ones in a string.
    !!
    !! Leave capitals, numbers and punctuation untouched
    !!
    !! @param[in] strIn String to convert to uppercase
    !! @returns Uppercase string
    !!
    !! @see https://stackoverflow.com/questions/10759375/how-can-i-write-a-to-upper-or-to-lower-function-in-f90
    function to_upper(strIn) result(strOut)
        ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
        ! Original author: Clive Page
        
        implicit none

        character(len=*), intent(in) :: strIn
        character(len=len(strIn)) :: strOut
        integer :: i,j

        do i = 1, len(strIn)
            j = iachar(strIn(i:i))
            if (j>= iachar("a") .and. j<=iachar("z") ) then
                strOut(i:i) = achar(iachar(strIn(i:i))-32)
            else
                strOut(i:i) = strIn(i:i)
            end if
        end do
        
    end function to_upper
end module ! GO