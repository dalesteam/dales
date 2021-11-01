!> \file modchecksim.f90
!!  Monitors Courant and Peclet numbers, and divergence.

!>
!!  Monitors Courant and Peclet numbers, and divergence.
!>
!!  These numbers are put out to screen either every tcheck seconds, or every time step (if tcheck=0).
!!  \author Thijs Heus,MPI-M
!!  \author Hans Cuijpers, KNMI
!!  \par Revision list
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
!
module modchecksim
  use modglobal, only : longint

  implicit none
  private
  public initchecksim,checksim

  real    :: tcheck = 0.
  integer(kind=longint) :: tnext = 3600.,itcheck
  real    :: dtmn =0.,ndt =0.

  ! explanations for dt_limit, determined in tstep_update()
  character (len=15) :: dt_reasons(0:5) = [character(len=15):: "initial step", "timee", "dt_lim" , "idtmax", "velocity", "diffusion"]
  
  save
contains
!> Initializing Checksim. Read out the namelist, initializing the variables
  subroutine initchecksim
    use modglobal, only : ifnamopt, fname_options,dtmax,ladaptive,btime,tres,checknamelisterror
    use modmpi,    only : myid,comm3d,mpierr, D_MPI_BCAST
    implicit none
    integer :: ierr
    namelist/NAMCHECKSIM/ &
    tcheck

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMCHECKSIM,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMCHECKSIM')
      write(6 ,NAMCHECKSIM)
      close(ifnamopt)

      if (.not. ladaptive .and. tcheck < dtmax) then
        tcheck = dtmax
      end if
    end if

    call D_MPI_BCAST(tcheck     ,1,0,comm3d,mpierr)
    itcheck = floor(tcheck/tres)
    tnext = itcheck+btime


  end subroutine initchecksim
!>Run checksim. Timekeeping, and output
  subroutine checksim
    use modglobal, only : timee,rtimee, rk3step, rdt
    use modmpi,    only : myid
    implicit none
    character(20) :: timeday
    if (timee ==0) return
    if (rk3step/=3) return
    dtmn = dtmn +rdt; ndt =ndt+1.
    if(timee<tnext) return
    tnext = tnext+itcheck
    dtmn  = dtmn / ndt
    if (myid==0) then
      call date_and_time(time=timeday)
      write (*,*) '================================================================='
      write (*,'(3A,F9.2,A,F12.9)') 'Time of Day: ', timeday(1:10),'    Time of Simulation: ', rtimee, '    dt: ',dtmn
    end if
    call calccourant
    call calcpeclet
    call chkdiv
    dtmn  = 0.
    ndt   = 0.

  end subroutine checksim
!>      Calculates the courant number as in max(w)*deltat/deltaz
  subroutine calccourant
    use modglobal, only : i1,j1,kmax,dx,dy,dzh
    use modfields, only : u0,v0,w0
    use modmpi,    only : myid,comm3d,mpierr,mpi_max, D_MPI_ALLREDUCE
    implicit none

    real, allocatable, dimension (:) :: courxl,courx,couryl,coury,courzl,courz,courtotl,courtot
    integer       :: k

    allocate(courxl(kmax),courx(kmax),couryl(kmax),coury(kmax),courzl(kmax),courz(kmax),courtotl(kmax),courtot(kmax))
    courxl = 0.0
    courx  = 0.0
    couryl = 0.0
    coury  = 0.0
    courzl = 0.0
    courz  = 0.0
    courtot  = 0.0
    do k=1,kmax
      courxl(k)=maxval(abs(u0(2:i1,2:j1,k)))*dtmn/dx
      couryl(k)=maxval(abs(v0(2:i1,2:j1,k)))*dtmn/dy
      courzl(k)=maxval(abs(w0(2:i1,2:j1,k)))*dtmn/dzh(k)
      courtotl(k)=maxval(u0(2:i1,2:j1,k)*u0(2:i1,2:j1,k)/(dx*dx)+v0(2:i1,2:j1,k)*v0(2:i1,2:j1,k)/(dy*dy)+&
      w0(2:i1,2:j1,k)*w0(2:i1,2:j1,k)/(dzh(k)*dzh(k)))*dtmn*dtmn
    end do
    call D_MPI_ALLREDUCE(courxl  ,courx  ,kmax,MPI_MAX,comm3d,mpierr)
    call D_MPI_ALLREDUCE(couryl  ,coury  ,kmax,MPI_MAX,comm3d,mpierr)
    call D_MPI_ALLREDUCE(courzl  ,courz  ,kmax,MPI_MAX,comm3d,mpierr)
    call D_MPI_ALLREDUCE(courtotl,courtot,kmax,MPI_MAX,comm3d,mpierr)
    if (myid==0) then
      write(*,'(A,3ES10.2,I5,ES10.2,I5)') 'Courant numbers (x,y,z,tot):',&
      maxval(courx(1:kmax)),maxval(coury(1:kmax)),maxval(courz(1:kmax)),maxloc(courz(1:kmax)),sqrt(maxval(courtot(1:kmax))),maxloc(courtot(1:kmax))
    end if

    deallocate(courxl,courx,couryl,coury,courzl,courz,courtotl,courtot)

    return
  end subroutine calccourant
!> Calculates the cell peclet number as max(ekm) *deltat/deltax**2
  subroutine calcpeclet

    use modglobal, only : i1,j1,k1,kmax,dx,dy,dzh
    use modsubgrid,only : ekm
    use modmpi,    only : myid,comm3d,mpierr,mpi_max, D_MPI_ALLREDUCE
    implicit none


    real, allocatable, dimension (:) :: peclettotl,peclettot
    integer       :: k

    allocate(peclettotl(k1),peclettot(k1))
    peclettotl = 0.
    peclettot  = 0.
    do k=1,kmax
      peclettotl(k)=maxval(ekm(2:i1,2:j1,k))*dtmn/minval((/dzh(k),dx,dy/))**2
    end do

    call D_MPI_ALLREDUCE(peclettotl,peclettot,k1,MPI_MAX,comm3d,mpierr)
    if (myid==0) then
      write(6,'(A,ES10.2,I5)') 'Cell Peclet number:',maxval(peclettot(1:kmax)),maxloc(peclettot(1:kmax))
    end if

    deallocate(peclettotl,peclettot)

    return
  end subroutine calcpeclet
!> Checks local and total divergence
  subroutine chkdiv

    use modglobal, only : i1,j1,kmax,dx,dy,dzf,dt_reason
    use modfields, only : u0,v0,w0,rhobf,rhobh
    use modmpi,    only : myid,comm3d,mpi_sum,mpi_max,mpierr, D_MPI_ALLREDUCE
    implicit none



    real div, divmax, divtot
    real divmaxl, divtotl
    integer i, j, k

    divmax = 0.
    divtot = 0.
    divmaxl= 0.
    divtotl= 0.

    do k=1,kmax
    do j=2,j1
    do i=2,i1
       div = &
                rhobf(k) * (u0(i+1,j,k) - u0(i,j,k) )/dx + &
                rhobf(k) * (v0(i,j+1,k) - v0(i,j,k) )/dy + &
                (rhobh(k+1)*w0(i,j,k+1) - rhobh(k)*w0(i,j,k) )/dzf(k)
      divmaxl = max(divmaxl,abs(div))
      divtotl = divtotl + div*dx*dy*dzf(k)
    end do
    end do
    end do

    call D_MPI_ALLREDUCE(divtotl, divtot, 1,     &
                          MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(divmaxl, divmax, 1,     &
                          MPI_MAX, comm3d,mpierr)

    if(myid==0)then
      write(6 ,'(A,2ES11.2,A,A)')'divmax, divtot = ', divmax, divtot,  '       dt limited by ', dt_reasons(dt_reason)
   end if

   return    
    
  end subroutine chkdiv

end module modchecksim

