!> \file tstep.f90
!!  Performs the time integration

!>
!!  Performs the time integration
!>
!! Tstep uses adaptive timestepping and 3rd order Runge Kutta time integration.
!! The adaptive timestepping chooses it's delta_t according to the courant number
!! and the cell peclet number, depending on the advection scheme in use.
!!
!!
!!  \author Chiel van Heerwaarden, Wageningen University
!!  \author Thijs Heus,MPI-M
!! \see Wicker and Skamarock 2002
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

!> Determine time step size dt in initialization and update time variables
!!
!! The size of the timestep Delta t is determined adaptively, and is limited by both the Courant-Friedrichs-Lewy criterion CFL
!! \latexonly
!! \begin{equation}
!! \CFL = \mr{max}\left(\left|\frac{u_i \Delta t}{\Delta x_i}\right|\right),
!! \end{equation}
!! and the diffusion number $d$. The timestep is further limited by the needs of other modules, e.g. the statistics.
!! \endlatexonly
module tstep
use modtimer
implicit none

  ! Arrays for keeping track of maximum values
!  real, allocatable, dimension(:) :: ummax
!  real, allocatable, dimension(:) :: vmmax
!  real, allocatable, dimension(:) :: wmmax
!  real, allocatable, dimension(:) :: ekhmax

contains
!> Allocate arrays
subroutine inittstep
  use modglobal, only : kmax
  implicit none
end subroutine inittstep

!> Deallocate arrays
subroutine exittstep
  implicit none
end subroutine exittstep

subroutine tstep_update
  use modglobal, only : i1,j1,k1,rk3step,timee,rtimee,dtmax,dt,ntrun,courant,peclet,dt_reason,nsv, &
                        kmax,dx,dy,dzh,dt_lim,ladaptive,timeleft,idtmax,rdt,tres,longint ,lwarmstart
  use modfields, only : um,vm,wm,up,vp,wp,thlp,svp,qtp,e12p
  use modsubgrid,only : ekm,ekh
  use modmpi,    only : comm3d,mpierr,mpi_max,D_MPI_ALLREDUCE
  use modtimer,       only : timer_tic, timer_toc
  implicit none

  integer       :: i, j, k, n
  real,save     :: courtotmax=-1,peclettot=-1
  real          :: courold, cfl_sq_l, cfl_sq, peclettotl, pecletold, pe_ekm, pe_ekh, min_size_sq
  logical,save  :: spinup=.true.

  call timer_tic('tstep/tstep_update', 0)

  if(lwarmstart) spinup = .false.

  rk3step = mod(rk3step,3) + 1
  if(rk3step == 1) then

    ! Initialization
    if (spinup) then
      if (ladaptive) then
        courold = courtotmax
        pecletold = peclettot
        peclettotl = 0.0
        cfl_sq_l = -1.0
        !$acc parallel loop collapse(3) default(present) reduction(max:cfl_sq_l, peclettotl)
        do k = 1, kmax
          do j = 2, j1
            do i = 2, i1
              cfl_sq_l =  (um(i,j,k)*rdt/dx) * (um(i,j,k)*rdt/dx) &
                        + (vm(i,j,k)*rdt/dy) * (vm(i,j,k)*rdt/dy) &
                        + (wm(i,j,k)*rdt/dzh(k)) * (wm(i,j,k)*rdt/dzh(k))
              min_size_sq = min(dzh(k),min(dx,dy))**2
              pe_ekm = ekm(i,j,k)*rdt/min_size_sq
              pe_ekh = ekh(i,j,k)*rdt/min_size_sq
              peclettotl = max(pe_ekm, pe_ekh)
            enddo
          enddo
        enddo
        call D_MPI_ALLREDUCE(cfl_sq_l,cfl_sq,1,MPI_MAX,comm3d,mpierr)
        call D_MPI_ALLREDUCE(peclettotl,peclettot,1,MPI_MAX,comm3d,mpierr)
        courtotmax = sqrt(cfl_sq)

        if ( pecletold>0) then
          dt = min(timee,dt_lim,idtmax,floor(rdt/tres*courant/courtotmax,longint),floor(rdt/tres*peclet/peclettot,longint))
          dt_reason = minloc((/timee,dt_lim,idtmax,floor(rdt/tres*courant/courtotmax,longint),floor(rdt/tres*peclet/peclettot,longint)/),1)
          if (abs(courtotmax-courold)/courold<0.1 .and. (abs(peclettot-pecletold)/pecletold<0.1)) then
            spinup = .false.
          end if
        end if
        rdt = dble(dt)*tres
        dt_lim = timeleft
        timee   = timee  + dt
        rtimee  = dble(timee)*tres
        timeleft=timeleft-dt
        ntrun   = ntrun  + 1
      else
        dt = 2 * dt
        if (dt >= idtmax) then
          dt = idtmax
          spinup = .false.
        end if
        rdt = dble(dt)*tres
      end if

    ! Normal time loop
    else
      if (ladaptive) then
        peclettotl = 1e-5
        cfl_sq_l = -1.0
        !$acc parallel loop collapse(3) default(present) reduction(max:cfl_sq_l, peclettotl)
        do k = 1, kmax
          do j = 2, j1
            do i = 2, i1
              cfl_sq_l =  (um(i,j,k)*rdt/dx) * (um(i,j,k)*rdt/dx) &
                        + (vm(i,j,k)*rdt/dy) * (vm(i,j,k)*rdt/dy) &
                        + (wm(i,j,k)*rdt/dzh(k)) * (wm(i,j,k)*rdt/dzh(k))
              min_size_sq = min(dzh(k),min(dx,dy))**2
              pe_ekm = ekm(i,j,k)*rdt/min_size_sq
              pe_ekh = ekh(i,j,k)*rdt/min_size_sq
              peclettotl = max(pe_ekm, pe_ekh)
            enddo
          enddo
        enddo
        call D_MPI_ALLREDUCE(cfl_sq_l,cfl_sq,1,MPI_MAX,comm3d,mpierr)
        call D_MPI_ALLREDUCE(peclettotl,peclettot,1,MPI_MAX,comm3d,mpierr)
        courtotmax = sqrt(cfl_sq)

        dt = min(timee,dt_lim,idtmax,floor(rdt/tres*courant/courtotmax,longint),floor(rdt/tres*peclet/peclettot,longint))
        dt_reason = minloc((/timee,dt_lim,idtmax,floor(rdt/tres*courant/courtotmax,longint),floor(rdt/tres*peclet/peclettot,longint)/),1)
        rdt = dble(dt)*tres
        timeleft=timeleft-dt
        dt_lim = timeleft
        timee   = timee  + dt
        rtimee  = dble(timee)*tres
        ntrun   = ntrun  + 1
      else
        dt = idtmax
        rdt = dtmax
        ntrun   = ntrun  + 1
        timee   = timee  + dt !ntimee*dtmax
        rtimee  = dble(timee)*tres
        timeleft=timeleft-dt
      end if
    end if
  end if

  ! set all tendencies to zero
  !$acc parallel loop collapse(3) default(present)
  do k = 1, k1
    do j = 2, j1
      do i = 1, i1
        up(i,j,k)=0.
        vp(i,j,k)=0.
        wp(i,j,k)=0.
        thlp(i,j,k)=0.
        e12p(i,j,k)=0.
        qtp(i,j,k)=0.
      enddo
    enddo
  enddo

  ! Scalars
  if (nsv > 0) then
    !$acc parallel loop collapse(4) default(present)
    do n = 1, nsv
      do k = 1, k1
        do j = 2, j1
          do i = 1, i1
            svp(i,j,k,n)=0.
          enddo
        enddo
      enddo
    enddo
  endif

  call timer_toc('tstep/tstep_update')
end subroutine tstep_update


!> Time integration is done by a third order Runge-Kutta scheme.
!!
!! \latexonly
!! With $f^n(\phi^n)$ the right-hand side of the appropriate equation for variable
!! $\phi=\{\fav{u},\fav{v},\fav{w},e^{\smfrac{1}{2}},\fav{\varphi}\}$, $\phi^{n+1}$
!! at $t+\Delta t$ is calculated in three steps:
!! \begin{eqnarray}
!! \phi^{*} &=&\phi^n + \frac{\Delta t}{3}f^n(\phi^n)\nonumber\\\\
!! \phi^{**} &=&\phi^{n} + \frac{\Delta t}{2}f^{*}(\phi^{*})\nonumber\\\\
!! \phi^{n+1} &=&\phi^{n} + \Delta t f^{**}(\phi^{**}),
!! \end{eqnarray}
!! with the asterisks denoting intermediate time steps.
!! \endlatexonly
!! \see Wicker and Skamarock, 2002
subroutine tstep_integrate


  use modglobal, only : rdt,rk3step,e12min,i1,j1,kmax,nsv
  use modfields, only : u0,um,up,v0,vm,vp,w0,wm,wp,&
                        thl0,thlm,thlp,qt0,qtm,qtp,&
                        e120,e12m,e12p,sv0,svm,svp
  implicit none

  integer :: i,j,k,n

  real rk3coef

  call timer_tic('tstep/tstep_integrate', 0)

  rk3coef = rdt / (4. - dble(rk3step))

  if(rk3step /= 3) then
    !$acc parallel loop collapse(3) default(present)
    do k = 1, kmax
      do j = 2, j1
        do i = 1, i1
          u0(i,j,k)   = um(i,j,k)   + rk3coef * up(i,j,k)
          v0(i,j,k)   = vm(i,j,k)   + rk3coef * vp(i,j,k)
          w0(i,j,k)   = wm(i,j,k)   + rk3coef * wp(i,j,k)
          thl0(i,j,k) = thlm(i,j,k) + rk3coef * thlp(i,j,k)
          qt0(i,j,k)  = qtm(i,j,k)  + rk3coef * qtp(i,j,k)
          e120(i,j,k) = max(e12min, e12m(i,j,k) + rk3coef * e12p(i,j,k))
        end do
      end do
    end do

    ! Scalars
    if (nsv > 0) then
      !$acc parallel loop collapse(4) default(present)
      do n = 1, nsv
        do k = 1, kmax
          do j = 2, j1
            do i = 1, i1
              sv0(i,j,k,n) = svm(i,j,k,n) + rk3coef * svp(i,j,k,n)
            end do
          end do
        end do
      end do
    endif
  else ! step 3 - store result in both ..0 and ..m
    !$acc parallel loop collapse(3) default(present)
    do k = 1, kmax
      do j = 2, j1
        do i = 1, i1
          um(i,j,k)   = um(i,j,k)   + rk3coef * up(i,j,k)
          u0(i,j,k)   = um(i,j,k)
          vm(i,j,k)   = vm(i,j,k)   + rk3coef * vp(i,j,k)
          v0(i,j,k)   = vm(i,j,k)
          wm(i,j,k)   = wm(i,j,k)   + rk3coef * wp(i,j,k)
          w0(i,j,k)   = wm(i,j,k)
          thlm(i,j,k) = thlm(i,j,k) + rk3coef * thlp(i,j,k)
          thl0(i,j,k) = thlm(i,j,k)
          qtm(i,j,k)  = qtm(i,j,k)  + rk3coef * qtp(i,j,k)
          qt0(i,j,k)  = qtm(i,j,k)
          e12m(i,j,k) = max(e12min, e12m(i,j,k) + rk3coef * e12p(i,j,k))
          e120(i,j,k) = e12m(i,j,k)
        end do
      end do
    end do

    ! Scalars
    if (nsv > 0) then
      !$acc parallel loop collapse(4) default(present)
      do n = 1, nsv
        do k = 1, kmax
          do j = 2, j1
            do i = 1, i1
              svm(i,j,k,n) = svm(i,j,k,n) + rk3coef * svp(i,j,k,n)
              sv0(i,j,k,n) = svm(i,j,k,n)
            end do
          end do
        end do
      end do
    endif
  end if
  call timer_toc('tstep/tstep_integrate')
end subroutine tstep_integrate

end module tstep
