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

module modtstep

implicit none
save

real,allocatable,dimension(:) :: rkCoef  ! Runge Kutta time integration coefficients (Wicker and Skamarock)
real,allocatable,dimension(:) :: rka,rkb ! Runge Kutta weighing coefficients (e.g. Eq. 9 Carpenter and Kennedy, 1994)

contains

  subroutine inittstep
    use modglobal, only : iTimeInt,rkStep,rkMaxStep
    use modmpi,    only : mpierr,myid
    
    if (iTimeInt==1) then ! Wicker and Skamarock RK3 scheme
      rkMaxStep = 3
      allocate(rkCoef(rkMaxStep))
      rkCoef = (/ 1./3.,1./2.,1. /)

    elseif (iTimeInt==2) then ! Williamson Low Storage RK3 scheme
      rkMaxStep = 3
      allocate(rka(rkMaxStep),rkb(rkMaxStep))
      rka = (/ 0.   ,-5./9. ,-153./128. /)
      rkb = (/ 1./3.,15./16.,8./15.     /)
      ! The following lines are necessary for now, because other modules still use the old Runge Kutta method
      allocate(rkCoef(rkMaxStep))
 !     rkCoef = (/ 1./3.,1./2.,1. /)
      rkCoef = (/ 1./3.,15./16.,8./15.     /)
      
    else
      if(myid==0)then
        write(6,*)'ERROR: Nonvalid timeintegration scheme selection'
      end if
      call MPI_FINALIZE(mpierr)
      stop
    end if 

  end subroutine inittstep


  subroutine tstep_update
  
  
    use modglobal, only : i1,j1,rkStep,rkMaxStep,timee,rtimee,runtime,btime,dtmax,dt,ntimee,ntrun,courant,peclet,&
                          kmax,k1,dx,dy,dzh,dt_lim,ladaptive,timeleft,idtmax,rdt,tres,longint ,lwarmstart
    use modfields, only : u0,v0,w0,rhobf
    use modsubgrid,only : ekm
    use modmpi,    only : myid,comm3d,mpierr,mpi_max,my_real
    implicit none
  
    real, allocatable, dimension (:) :: courtotl,courtot
  
    integer       :: k
    real,save     :: courtotmax=-1,peclettot=-1
    real          :: courold,peclettotl,pecletold
    logical,save  :: spinup=.true.
  
    allocate(courtotl(k1),courtot(k1))
  
    if(lwarmstart) spinup = .false.
  
    rkStep = mod(rkStep,rkMaxStep)+1
    
    if(rkStep == 1) then
      ! Initialization
      if (spinup) then
        if (ladaptive) then
          courold = courtotmax
          pecletold = peclettot
          peclettotl=0.0
          do k=1,kmax
            courtotl(k)=maxval(u0(2:i1,2:j1,k)*u0(2:i1,2:j1,k)/(dx*dx)+v0(2:i1,2:j1,k)*v0(2:i1,2:j1,k)/(dy*dy)+&
            w0(2:i1,2:j1,k)*w0(2:i1,2:j1,k)/(dzh(k)*dzh(k)))*rdt*rdt
          end do
          call MPI_ALLREDUCE(courtotl,courtot,k1,MY_REAL,MPI_MAX,comm3d,mpierr)
          courtotmax=0.0
          do k=1,kmax
            courtotmax=max(courtotmax,courtot(k))
          enddo
          courtotmax=sqrt(courtotmax)
          do k=1,kmax
            peclettotl=max(peclettotl,maxval(ekm(2:i1,2:j1,k)/rhobf(k))*rdt/minval((/dzh(k),dx,dy/))**2)
          end do
          call MPI_ALLREDUCE(peclettotl,peclettot,1,MY_REAL,MPI_MAX,comm3d,mpierr)
          if ( pecletold>0) then
            dt = min(timee,dt_lim,idtmax,floor(rdt/tres*courant/courtotmax,longint),floor(rdt/tres*peclet/peclettot,longint))
            if (abs(courtotmax-courold)/courold<0.1 .and. (abs(peclettot-pecletold)/pecletold<0.1)) then
              spinup = .false.
            end if
          end if
          rdt = dble(dt)*tres
          dt_lim = timeleft
          timee   = timee  + dt
          rtimee  = dble(timee)*tres
          timeleft=timeleft-dt
          ntimee  = ntimee + 1
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
          do k=1,kmax
            courtotl(k)=maxval((u0(2:i1,2:j1,k)*rdt/dx)*(u0(2:i1,2:j1,k)*rdt/dx)+(v0(2:i1,2:j1,k)*rdt/dy)*&
            (v0(2:i1,2:j1,k)*rdt/dy)+(w0(2:i1,2:j1,k)*rdt/dzh(k))*(w0(2:i1,2:j1,k)*rdt/dzh(k)))
          end do      
          call MPI_ALLREDUCE(courtotl,courtot,k1,MY_REAL,MPI_MAX,comm3d,mpierr)
          courtotmax=0.0
          do k=1,kmax
              courtotmax=max(courtotmax,sqrt(courtot(k)))
          enddo
          do k=1,kmax
            peclettotl=max(peclettotl,maxval(ekm(2:i1,2:j1,k)/rhobf(k))*rdt/minval((/dzh(k),dx,dy/))**2)
          end do
          call MPI_ALLREDUCE(peclettotl,peclettot,1,MY_REAL,MPI_MAX,comm3d,mpierr)
          dt = min(timee,dt_lim,idtmax,floor(rdt/tres*courant/courtotmax,longint),floor(rdt/tres*peclet/peclettot,longint))
          rdt = dble(dt)*tres
          timeleft=timeleft-dt
          dt_lim = timeleft
          timee   = timee  + dt
          rtimee  = dble(timee)*tres
          ntimee  = ntimee + 1
          ntrun   = ntrun  + 1
        else
          dt = idtmax
          rdt = dtmax
          ntimee  = ntimee + 1
          ntrun   = ntrun  + 1
          timee   = timee  + dt !ntimee*dtmax
          rtimee  = dble(timee)*tres
          timeleft=timeleft-dt
        end if
      end if
    end if
  
    deallocate(courtotl,courtot)
  
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
  
    use modglobal, only : i1,j1,kmax,nsv,rdt,rkStep,rkMaxStep,e12min,lmoist,k1,ih,jh,rslabs,kcb,iTimeInt
    use modfields, only : u0,um,up,v0,vm,vp,w0,wm,wp,wp_store,&
                          thl0,thlm,thlp,qt0,qtm,qtp,&
                          e120,e12m,e12p,sv0,svm,svp
    implicit none
  
    integer i,j,k,n
    

  
    if (iTimeInt==1) then    ! the Wicker and Skamarock RK3 scheme
     
      do k=1,kmax
        do j=2,j1
          do i=2,i1
    
            u0(i,j,k)   = um(i,j,k)   + rkCoef(rkStep) * rdt * up(i,j,k)
            v0(i,j,k)   = vm(i,j,k)   + rkCoef(rkStep) * rdt * vp(i,j,k)
            w0(i,j,k)   = wm(i,j,k)   + rkCoef(rkStep) * rdt * wp(i,j,k)
            thl0(i,j,k) = thlm(i,j,k) + rkCoef(rkStep) * rdt * thlp(i,j,k)
            qt0(i,j,k)  = qtm(i,j,k)  + rkCoef(rkStep) * rdt * qtp(i,j,k)
            e120(i,j,k) = e12m(i,j,k) + rkCoef(rkStep) * rdt * e12p(i,j,k)
    
            e120(i,j,k) = max(e12min,e120(i,j,k))
            e12m(i,j,k) = max(e12min,e12m(i,j,k))
    
            do n=1,nsv
              sv0(i,j,k,n) = svm(i,j,k,n) + rkCoef(rkStep) * rdt * svp(i,j,k,n)
            end do
    
          end do
        end do
      end do
    
      if(rkStep == rkMaxStep) then
        um = u0
        vm = v0
        wm = w0
        thlm = thl0
        qtm  = qt0
        e12m = e120
        svm = sv0
      end if

      ! Store the vertical velocity tendency for later use in the sampling module
      wp_store = wp
      ! Reset all tendency variables
      up=0.
      vp=0.
      wp=0.
      thlp=0.
      qtp=0.
      svp=0.
      e12p=0.

    elseif (iTimeInt==2) then ! The RK3 scheme proposed by Williamson (1980)
     
      do k=1,kmax
        do j=2,j1
          do i=2,i1

            u0(i,j,k)   = u0(i,j,k)   + rkb(rkStep)*rdt*up(i,j,k)
            v0(i,j,k)   = v0(i,j,k)   + rkb(rkStep)*rdt*vp(i,j,k)
            w0(i,j,k)   = w0(i,j,k)   + rkb(rkStep)*rdt*wp(i,j,k)
            thl0(i,j,k) = thl0(i,j,k) + rkb(rkStep)*rdt*thlp(i,j,k)
            qt0(i,j,k)  = qt0(i,j,k)  + rkb(rkStep)*rdt*qtp(i,j,k)
            e120(i,j,k) = e120(i,j,k) + rkb(rkStep)*rdt*e12p(i,j,k)

            e120(i,j,k) = max(e12min,e120(i,j,k))

            do n=1,nsv
              sv0(i,j,k,n) = sv0(i,j,k,n) + rkb(rkStep)*rdt*svp(i,j,k,n)
            end do

          end do
        end do
      end do

      if (rkStep==rkMaxStep) then
        wp_store = wp
      end if
      ! Set the appropriate offsets for the tendencies
      up  = rka(rkStep)*rdt*up
      vp  = rka(rkStep)*rdt*vp
      wp  = rka(rkStep)*rdt*wp
      thlp= rka(rkStep)*rdt*thlp
      qtp = rka(rkStep)*rdt*qtp
      e12p= rka(rkStep)*rdt*e12p
      svp = rka(rkStep)*rdt*svp

    end if
     
  end subroutine tstep_integrate

end module modtstep
