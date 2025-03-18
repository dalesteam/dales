!> \file modscalarpulse.f90
!!  Adds a pulse of qt to a single, horizontal, spatial wavenumber, 
!!  over a predefined height and at a preset time

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

module modscalarpulse

use modglobal, only : longint
implicit none
private
public :: initscalarpulse, scalarpulse, lscalarpulse

save
! pulse governing variables (from namelist)
  logical :: lscalarpulse     = .false.
  logical :: lcpmip           = .false.
  integer(kind=longint) :: timepulse
  real    :: amppulse, kpulse, zminpulse, zmaxpulse, radius
  real, allocatable :: qtav0(:),qtav1(:)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initscalarpulse
    use mpi
    use modmpi,    only :myid,D_MPI_BCAST, &
                         mpierr,commwrld
    use modglobal, only :cexpnr,runtime,ifnamopt,fname_options, &
                         checknamelisterror,tres,k1

    implicit none

    integer ierr
    namelist/NAMSCALARPULSE/ &
    lscalarpulse, lcpmip, timepulse, amppulse, kpulse, zminpulse, zmaxpulse, radius

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMSCALARPULSE,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMSCALARPULSE')
      write(6 ,NAMSCALARPULSE)
      close(ifnamopt)
    end if

    timepulse = timepulse/tres

    call D_MPI_BCAST(lscalarpulse       ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(lcpmip             ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(timepulse          ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(amppulse           ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(kpulse             ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(zminpulse          ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(zmaxpulse          ,1,0,commwrld,mpierr)
    call D_MPI_BCAST(radius             ,1,0,commwrld,mpierr)

    if (lcpmip) allocate(qtav0(k1), qtav1(k1))

  end subroutine initscalarpulse

  subroutine scalarpulse

    use modglobal, only : rk3step,timee,dt_lim
    use modmpi,    only : myid
    implicit none
    if (.not. lscalarpulse) return
    if (rk3step/=3) return

    if(timee<timepulse) then
      dt_lim = minval((/dt_lim,timepulse-timee/))
      return
    end if
    if (timee>=timepulse) then
      if (myid==0) then
        print *, 'Performing scalar pulse...'
      end if
      call do_scalarpulse
      lscalarpulse = .false.
    end if

  end subroutine scalarpulse

  subroutine do_scalarpulse

    use modfields     , only : qtm
    use modglobal     , only : i1,j1,k1,imax,jmax, &
                               itot,jtot,dx,dy,zf,pi, &
                               ih,jh,ijtot
    use modmpi        , only : myidx,myidy,myid,slabsum

    logical kstartflag
    integer kstart, kend, i, j, k
    real    xf, yf, facx, facy, qtpulse, center_x, center_y

    kstartflag = .true.
    xf         = myidx*imax*dx
    yf         = myidy*jmax*dy
    facx       = 2*pi*kpulse/(itot*dx)
    facy       = 2*pi*kpulse/(jtot*dy)
    qtpulse    = 0.
    qtav0      = 0.
    qtav1      = 0.

    ! Calculate the levels to apply the perturbation at
    do k=1,k1
      if (zf(k) >= zminpulse .and. kstartflag) then
        kstart = k
        kstartflag = .false.
      end if
      if (zf(k) > zmaxpulse) then
        kend = k
        exit
      end if
    end do
    if (myid == 0) then
      print *, 'kstart, kend', kstart, kend 
    end if

    ! Apply the perturbation
    if (lcpmip) then
      ! Protocol for cpmip

      ! Calculate domain-mean profile
      call slabsum(qtav0 ,1,k1,qtm ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      qtav0 = qtav0 / ijtot

      ! Find the domain center
      center_x = dx*itot*0.5
      center_y = dy*jtot*0.5

      ! Add the pulse
      do i=2,i1
        xf = myidx*imax*dx + dx*(i - 1.5)
        do j=2,j1
          yf = myidy*jmax*dy + dy*(j - 1.5)
          if ((xf-center_x)**2.0D0+(yf-center_y)**2.0D0 .lt. radius**2.0D0) then
            qtpulse = amppulse * cos( pi/2.0D0 * (sqrt( (xf-center_x)**2.0D0 + (yf-center_y)**2.0D0 ) / radius) )**2.0D0
          else
            qtpulse = 0.
          end if
          if (myid == 0) then
            print *, 'x, y, qtpulse', xf, yf, qtpulse
          end if
          do k=kstart,kend 
            qtm(i,j,k) = qtm(i,j,k) + qtpulse
          end do
        end do
      end do

      ! Calculate domain-mean profile again
      call slabsum(qtav1 ,1,k1,qtm ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      qtav1 = qtav1 / ijtot

      ! And subtract the difference everywhere
      do i=2,i1
        do j=2,j1
          do k=kstart,kend 
            qtm(i,j,k) = qtm(i,j,k) - (qtav1(k) - qtav0(k))
          end do
        end do
      end do

    else
      do i=2,i1
        xf = myidx*imax*dx + dx*(i - 1.5)
        do j=2,j1
          yf = myidy*jmax*dy + dy*(j - 1.5)
          qtpulse = amppulse * cos(facx*xf + pi) &
                             * cos(facy*yf + pi)
          if (myid == 0) then
            print *, 'x, y, qtpulse', xf, yf, qtpulse
          end if
          do k=kstart,kend 
            qtm(i,j,k) = qtm(i,j,k) + qtpulse
          end do
        end do
      end do
    end if

  end subroutine do_scalarpulse

end module modscalarpulse