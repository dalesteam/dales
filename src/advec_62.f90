!> \file advec_62.f90
!!  Does advection with a 6th order horizontal central differencing scheme.
!!  Vertical advection is done with 2nd order scheme.
!! \par Revision list
!! \par Authors
!! \see Wicker and Skamarock 2002
!!
!! A higher-order accuracy in the calculation of the advection is reached with a
!! sixth order central differencing scheme.
!! \latexonly
!!!! \endlatexonly
!!
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!

module advec_62
contains
!> Advection at cell center
subroutine advecc_62(putin, putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi,dyi,dzf,dzh,leq
  use modfields, only : u0, v0, w0,rhobf

  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the cell centered field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency
  real                                       :: inv2dzfk, rhobf_p, rhobf_m

  integer :: i,j,k

  if (leq) then

  k = 1
  inv2dzfk = 1./(2. * dzf(k))
  rhobf_p = rhobf(k+1)/rhobf(k)
  do j=2,j1
     do i=2,i1
        putout(i,j,k)  = putout(i,j,k)- ( &
             ( &
             u0(i+1,j,k)/60. &
             *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
             -u0(i,j,k)/60. &
             *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
             )*dxi&
             +(&
             v0(i,j+1,k)/60. &
             *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
             -v0(i,j,k)/60. &
             *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
             )* dyi &
             + ( &
             w0(i,j,k+1) * (rhobf_p * putin(i,j,k+1) + putin(i,j,k)) &
             ) * inv2dzfk &
             )
     end do
  end do


  do k=2,kmax
     inv2dzfk = 1./(2. * dzf(k))
     rhobf_p = rhobf(k+1)/rhobf(k)
     rhobf_m = rhobf(k-1)/rhobf(k)

    do j=2,j1
      do i=2,i1


              putout(i,j,k)  = putout(i,j,k)- (  &
                  ( &
                      u0(i+1,j,k)/60. &
                      *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                      -u0(i,j,k)/60. &
                      *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                  )*dxi&
                +(&
                      v0(i,j+1,k)/60. &
                      *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                      -v0(i,j,k)/60. &
                      *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                  )* dyi &
                + ( &
                  w0(i,j,k+1) * (rhobf_p * putin(i,j,k+1) + putin(i,j,k)) &
                  -w0(i,j,k)  * (rhobf_m * putin(i,j,k-1) + putin(i,j,k)) &
                  ) * inv2dzfk &
                  )

      end do
    end do
 end do
else ! non-equidistant grid
 k = 1
  inv2dzfk = 1./(2. * dzf(k))
  rhobf_p = rhobf(k+1)/rhobf(k)
  do j=2,j1
     do i=2,i1
        putout(i,j,k)  = putout(i,j,k)- ( &
             ( &
             u0(i+1,j,k)/60. &
             *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
             -u0(i,j,k)/60. &
             *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
             )*dxi&
             +(&
             v0(i,j+1,k)/60. &
             *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
             -v0(i,j,k)/60. &
             *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
             )* dyi &
             + ( &
             w0(i,j,k+1) * (rhobf_p * putin(i,j,k+1) * dzf(k) +  putin(i,j,k) * dzf(k+1) ) / dzh(k+1) &
             ) * inv2dzfk &
             )
     end do
  end do


  do k=2,kmax
     inv2dzfk = 1./(2. * dzf(k))
     rhobf_p = rhobf(k+1)/rhobf(k)
     rhobf_m = rhobf(k-1)/rhobf(k)

    do j=2,j1
      do i=2,i1


              putout(i,j,k)  = putout(i,j,k)- (  &
                  ( &
                      u0(i+1,j,k)/60. &
                      *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                      -u0(i,j,k)/60. &
                      *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                  )*dxi&
                +(&
                      v0(i,j+1,k)/60. &
                      *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                      -v0(i,j,k)/60. &
                      *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                  )* dyi &
                  + ( &
                   w0(i,j,k+1) * (rhobf_p * putin(i,j,k+1) * dzf(k) +  putin(i,j,k) * dzf(k+1) ) / dzh(k+1) &
                  -w0(i,j,k )  * (rhobf_m * putin(i,j,k-1) * dzf(k) +  putin(i,j,k) * dzf(k-1) ) / dzh(k) &
                  ) * inv2dzfk &
                  )

      end do
    end do
 end do
 end if
end subroutine advecc_62


!> Advection at the u point.
subroutine advecu_62(putin,putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzf,dzh,leq
  use modfields, only : u0, v0, w0,rhobf

  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the u field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency
  real                                       :: inv4dzfk, rhobf_p, rhobf_m

  integer :: i,j,k

  if (leq) then
  
  k = 1
  inv4dzfk = 1./(4. * dzf(k))
  rhobf_p = rhobf(k+1)/rhobf(k)
  do j=2,j1
     do i=2,i1

        putout(i,j,k)  = putout(i,j,k)- ( &
             ( &
             (u0(i+1,j,k)+u0(i,j,k))/60. &
             *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
             -(u0(i,j,k)+u0(i-1,j,k))/60. &
             *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
             )*dxi5&
             +(&
             (v0(i,j+1,k)+v0(i-1,j+1,k))/60. &
             *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
             -(v0(i,j,k)+v0(i-1,j,k))/60. &
             *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
             )* dyi5 &
             + ( &
             (rhobf_p * putin(i,j,k+1) + putin(i,j,k)) *(w0(i,j,k+1)+ w0(i-1,j,k+1)) &
             ) * inv4dzfk &
             )
     end do
  end do


    do k=2,kmax
       inv4dzfk = 1./(4. * dzf(k))
       rhobf_p = rhobf(k+1)/rhobf(k)
       rhobf_m = rhobf(k-1)/rhobf(k)

      do j=2,j1
        do i=2,i1

          putout(i,j,k)  = putout(i,j,k)- ( &
                (&
                    (u0(i+1,j,k)+u0(i,j,k))/60. &
                    *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i-1,j,k))/60. &
                    *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                )*dxi5&
              +(&
                    (v0(i,j+1,k)+v0(i-1,j+1,k))/60. &
                    *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i-1,j,k))/60. &
                    *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                )* dyi5 &
              + ( &
                     (putin(i,j,k) + rhobf_p * putin(i,j,k+1) )*(w0(i,j,k+1)+w0(i-1,j,k+1)) &
                    -(putin(i,j,k) + rhobf_m * putin(i,j,k-1) )*(w0(i,j,k  )+w0(i-1,j,k  )) &
                ) * inv4dzfk &
                )

        end do
      end do
    end do
 else ! non-equidistant grid
    k = 1
  inv4dzfk = 1./(4. * dzf(k))
  rhobf_p = rhobf(k+1)/rhobf(k)
  do j=2,j1
     do i=2,i1

        putout(i,j,k)  = putout(i,j,k)- ( &
             ( &
             (u0(i+1,j,k)+u0(i,j,k))/60. &
             *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
             -(u0(i,j,k)+u0(i-1,j,k))/60. &
             *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
             )*dxi5&
             +(&
             (v0(i,j+1,k)+v0(i-1,j+1,k))/60. &
             *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
             -(v0(i,j,k)+v0(i-1,j,k))/60. &
             *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
             )* dyi5 &
             + ( &
             (putin(i,j,k) * dzf(k+1) + rhobf_p * putin(i,j,k+1)*dzf(k)) * (w0(i,j,k+1)+w0(i-1,j,k+1)) / dzh(k+1) &
             ) * inv4dzfk &
             )
     end do
  end do


    do k=2,kmax
       inv4dzfk = 1./(4. * dzf(k))
       rhobf_p = rhobf(k+1)/rhobf(k)
       rhobf_m = rhobf(k-1)/rhobf(k)

      do j=2,j1
        do i=2,i1

          putout(i,j,k)  = putout(i,j,k)- ( &
                (&
                    (u0(i+1,j,k)+u0(i,j,k))/60. &
                    *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i-1,j,k))/60. &
                    *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                )*dxi5&
              +(&
                    (v0(i,j+1,k)+v0(i-1,j+1,k))/60. &
                    *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i-1,j,k))/60. &
                    *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                )* dyi5 &
              + ( &
                     (putin(i,j,k) * dzf(k+1) + rhobf_p * putin(i,j,k+1) *dzf(k)) * (w0(i,j,k+1)+w0(i-1,j,k+1)) / dzh(k+1) &
                    -(putin(i,j,k) * dzf(k-1) + rhobf_m * putin(i,j,k-1) *dzf(k)) * (w0(i,j,k  )+w0(i-1,j,k  )) / dzh(k)   &                     
                ) * inv4dzfk &
                )

        end do
      end do
    end do
    
  end if

end subroutine advecu_62



!> Advection at the v point.
subroutine advecv_62(putin, putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzf,dzh,leq
  use modfields, only : u0, v0, w0,rhobf
  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the v field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency
  real                                       :: inv4dzfk, rhobf_p, rhobf_m

  integer :: i,j,k

  if (leq) then

  k = 1
  inv4dzfk = 1./(4. * dzf(k))
  rhobf_p = rhobf(k+1)/rhobf(k)

    do j=2,j1
     do i=2,i1

        putout(i,j,k)  = putout(i,j,k)- ( &
             ( &
             (u0(i+1,j,k)+u0(i+1,j-1,k))/60. &
             *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
             -(u0(i,j,k)+u0(i,j-1,k))/60. &
             *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
             )*dxi5&
             +(&
             (v0(i,j+1,k)+v0(i,j,k))/60. &
             *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
             -(v0(i,j,k)+v0(i,j-1,k))/60. &
             *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
             )* dyi5 &
             +( &
             (w0(i,j,k+1)+w0(i,j-1,k+1)) *(rhobf_p * putin(i,j,k+1)+putin(i,j,k)) &                          
             ) * inv4dzfk  &
             )
     end do
  end do
        
    do k=2,kmax
       inv4dzfk = 1./(4. * dzf(k))
       rhobf_p = rhobf(k+1)/rhobf(k)
       rhobf_m = rhobf(k-1)/rhobf(k)

      do j=2,j1
        do i=2,i1
        
          putout(i,j,k)  = putout(i,j,k)- ( &
                ( &
                    (u0(i+1,j,k)+u0(i+1,j-1,k))/60. &
                    *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i,j-1,k))/60. &
                    *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                 )*dxi5&
                +(&
                    (v0(i,j+1,k)+v0(i,j,k))/60. &
                    *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i,j-1,k))/60. &
                    *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                  )* dyi5 &
                + ( &
                (w0(i,j,k+1)+w0(i,j-1,k+1))*(rhobf_p * putin(i,j,k+1) + putin(i,j,k))&
                -(w0(i,j,k) +w0(i,j-1,k))  *(rhobf_m * putin(i,j,k-1) + putin(i,j,k))&
                  ) * inv4dzfk  &
                  )
        end do
      end do
    end do
 else ! non-equidistant grid
    k = 1
  inv4dzfk = 1./(4. * dzf(k))
  rhobf_p = rhobf(k+1)/rhobf(k)

    do j=2,j1
     do i=2,i1

        putout(i,j,k)  = putout(i,j,k)- ( &
             ( &
             (u0(i+1,j,k)+u0(i+1,j-1,k))/60. &
             *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
             -(u0(i,j,k)+u0(i,j-1,k))/60. &
             *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
             )*dxi5&
             +(&
             (v0(i,j+1,k)+v0(i,j,k))/60. &
             *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
             -(v0(i,j,k)+v0(i,j-1,k))/60. &
             *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
             )* dyi5 &
             +( &
             (w0(i,j,k+1)+w0(i,j-1,k+1)) * (rhobf_p * putin(i,j,k+1)*dzf(k) + putin(i,j,k)*dzf(k+1)) / dzh(k+1)&
             ) * inv4dzfk  &
             )
     end do
  end do
        
    do k=2,kmax
       inv4dzfk = 1./(4. * dzf(k))
       rhobf_p = rhobf(k+1)/rhobf(k)
       rhobf_m = rhobf(k-1)/rhobf(k)

      do j=2,j1
        do i=2,i1
        
          putout(i,j,k)  = putout(i,j,k)- ( &
                ( &
                    (u0(i+1,j,k)+u0(i+1,j-1,k))/60. &
                    *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i,j-1,k))/60. &
                    *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                 )*dxi5&
                +(&
                    (v0(i,j+1,k)+v0(i,j,k))/60. &
                    *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i,j-1,k))/60. &
                    *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                  )* dyi5 &
                  + ( &
                  (w0(i,j,k+1)+w0(i,j-1,k+1)) * (rhobf_p * putin(i,j,k+1) * dzf(k) + putin(i,j,k) * dzf(k+1)) / dzh(k+1)&
                 -(w0(i,j,k)  +w0(i,j-1,k))   * (rhobf_m * putin(i,j,k-1) * dzf(k) + putin(i,j,k) * dzf(k-1)) / dzh(k)  &
                  ) * inv4dzfk  &
                  )
        end do
      end do
   end do
end if

end subroutine advecv_62



!> Advection at the w point.
subroutine advecw_62(putin, putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzh
  use modfields, only : u0, v0, w0,rhobh
  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the w field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency
  real                                       :: inv4dzhk, rhobh_p, rhobh_m

  integer :: i,j,k

  !if (leq) then

    do k=2,kmax
       inv4dzhk = 1./(4. * dzh(k))
       rhobh_p = rhobh(k+1)/rhobh(k)
       rhobh_m = rhobh(k-1)/rhobh(k)

      do j=2,j1
        do i=2,i1

           putout(i,j,k)  = putout(i,j,k)- ( &
                 (&
                     (u0(i+1,j,k)+u0(i+1,j,k-1))/60. &
                     *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                     -(u0(i,j,k)+u0(i,j,k-1))/60. &
                     *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                 )*dxi5&
                +(&
                     (v0(i,j+1,k)+v0(i,j+1,k-1))/60. &
                     *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                     -(v0(i,j,k)+v0(i,j,k-1))/60. &
                     *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                  )* dyi5 &
                + ( &
                      (putin(i,j,k)+rhobh_p * putin(i,j,k+1) )*(w0(i,j,k) + w0(i,j,k+1)) &
                     -(putin(i,j,k)+rhobh_m * putin(i,j,k-1) )*(w0(i,j,k) + w0(i,j,k-1)) &
                  )*inv4dzhk &
                  )
       end do
      end do
     end do

end subroutine advecw_62
end module advec_62
