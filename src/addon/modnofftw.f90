!> \file modfft2d.f90
!>
!!  Solve the poisson equation using a FFT + tridiagonal solver. Uses FFTW.
!>
!!  \author Jisk Attema, Netherlands eScience Center
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
!  Copyright 2014 Netherlands eScience Center
!
module modfftw
use modprecision, only : pois_r

implicit none

contains

  subroutine fftwinit(p, Fp, d, xyrt, ps,pe,qs,qe)
    real(pois_r), pointer      :: p(:,:,:)
    real(pois_r), pointer      :: Fp(:,:,:)
    real(pois_r), allocatable  :: d(:,:,:)
    real(pois_r), allocatable  :: xyrt(:,:)
    integer,intent(out)        :: ps,pe,qs,qe
    call error_and_exit()
    ps=0 ! suppress warnings about intent(out) variables not being assigned
    pe=0
    qs=0
    qe=0
 end subroutine

 subroutine fftwexit(p,Fp,d,xyrt)
    real(pois_r), pointer     :: p(:,:,:)
    real(pois_r), pointer     :: Fp(:,:,:)
    real(pois_r), allocatable :: d(:,:,:)
    real(pois_r), allocatable :: xyrt(:,:)
    call error_and_exit()
 end subroutine

  subroutine fftwf(p, Fp)
    real(pois_r), pointer :: p(:,:,:)
    real(pois_r), pointer :: Fp(:,:,:)
    call error_and_exit()
  end subroutine

  subroutine fftwb(p, Fp)
    real(pois_r), pointer :: p(:,:,:)
    real(pois_r), pointer :: Fp(:,:,:)
    call error_and_exit()
  end subroutine

  subroutine error_and_exit
    write (*,*) 'DALES was compiled without FFTW.'
    write (*,*) 'Use the default poisson solver (solver_id=0),'
    write (*,*) 'or recompile DALES with the option USE_FFTW.'

    call exit(-1)
  end subroutine

end module modfftw
