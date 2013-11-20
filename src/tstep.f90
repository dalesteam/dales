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
!!  \author Johan van der Dussen
!! \see Wicker and Skamarock 2002
!! \see Williamson 1980
!! \see Gottlieb and Shu 1998
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
use modglobal, only : rkStep,rkMaxStep,iTimeInt,iTimeWicker,iTimeLowStor,iTimeTVD

implicit none
save

real,allocatable,dimension(:) :: rka,rkb,rkc ! Runge Kutta weighting coefficients (e.g. Eq. 9 Carpenter and Kennedy, 1994)

contains

  ! This subroutine initializes the coefficients that are used for the
  ! Runge-Kutta time integration schemes.
  subroutine inittstep
    use modmpi,    only : mpierr,myid
    
    select case (iTimeInt)
      case(iTimeWicker) ! Wicker and Skamarock RK3 scheme
        rkMaxStep = 3
        allocate(rka(rkMaxStep),rkb(rkMaxStep),rkc(rkMaxStep))
        rka(:) = 0.
        rkb    = (/ 1./3.,1./2.,1. /)
        rkc    = (/ 1./3.,1./2.,1. /)

      case(iTimeLowStor) ! Williamson Low Storage RK3 scheme
        rkMaxStep = 3
        allocate(rka(rkMaxStep),rkb(rkMaxStep),rkc(rkMaxStep))
        rka = (/-5./9.,-153./128., 0.     /)
        rkb = (/ 1./3.,15./16.   , 8./15. /)
        rkc = (/ 1./3., 3./4.    , 1.     /)
        
      case(iTimeTVD) ! Third order total variation diminishing Runge-Kutta scheme
        rkMaxStep = 3
        allocate(rka(rkMaxStep),rkb(rkMaxStep),rkc(rkMaxStep))
        rka = (/ 1.   , 1./4.,   0.  /)
        rkb = (/ 1.   , 1./4., 2./3. /)
        rkc = (/ 1./4., 1./2., 1.    /)
        
      case default
        if(myid==0)then
          write(6,*)'ERROR: Nonvalid timeintegration scheme selection'
        end if
        call MPI_FINALIZE(mpierr)
        stop
    end select

    ! Set the rkStep to the appropriate value (see tstep_update)
    rkStep = rkMaxStep

  end subroutine inittstep


  subroutine tstep_update
    use modglobal, only : timee,rtimee,runtime,btime,dtmax,dt,ntimee,ntrun,courant,peclet,&
                          dt_lim,ladaptive,timeleft,idtmax,rdt,subDt,tres,longint,lwarmstart
    implicit none
  
    real          :: courantold=-1.,pecletold=-1., &
                     courantmax,pecletmax
    logical,save  :: spinup=.true.
    integer       :: k
  
    if (lwarmstart) spinup = .false.

    ! Determine the new timestep at the end of the previous full RK time step
    if (rkStep==rkMaxStep) then
      if (spinup) then
      !== Initialization =========
        if (ladaptive) then
          courantmax= getCourant()
          pecletmax = getPeclet()
          
          if (pecletold>0) then
            dt = min( timee,  &
                      dt_lim, &
                      idtmax, &
                      floor(rdt/tres*courant/courantmax,longint), &
                      floor(rdt/tres*peclet /pecletmax ,longint)  )
            if (abs(courantmax-courantold)/courantold<0.1 .and. (abs(pecletmax-pecletold)/pecletold<0.1)) then
              spinup = .false.
            end if
          end if
          ! Store the previous Courant and Peclet numbers
          courantold = courantmax
          pecletold  = pecletmax
          ! Set all relevant time variables
          rdt     = dble(dt)*tres
          dt_lim  = timeleft
          timee   = timee  + dt
          rtimee  = dble(timee)*tres
          timeleft= timeleft-dt
          ntimee  = ntimee + 1
          ntrun   = ntrun  + 1
        else ! ladaptive=.false.
          dt = 2*dt
          if (dt >= idtmax) then
            dt = idtmax
            spinup = .false.
          end if
          rdt = dble(dt)*tres
        end if
      !== Normal time loop
      else
        if (ladaptive) then
          courantMax = getCourant()
          pecletMax  = getPeclet()
          ! Determine the maximum stable time step (in millisec, type longint)
          dt      = min( timee,  &
                         dt_lim, &
                         idtmax, &
                         floor(rdt/tres*courant/courantmax,longint), &
                         floor(rdt/tres*peclet /pecletmax ,longint)  )
          rdt     = dble(dt)*tres ! Timestep in seconds
          timeleft= timeleft-dt
          dt_lim  = timeleft
          timee   = timee  + dt
          rtimee  = dble(timee)*tres
          ntimee  = ntimee + 1
          ntrun   = ntrun  + 1
        else ! ladaptive = .false.
          dt      = idtmax
          rdt     = dtmax
          ntimee  = ntimee + 1
          ntrun   = ntrun  + 1
          timee   = timee  + dt !ntimee*dtmax
          rtimee  = dble(timee)*tres
          timeleft=timeleft-dt
        end if ! ladaptive
      end if ! spinup
    end if ! rkStep==rkMaxStep

    ! Update the index of the RK substep and the magnitude of the subtimestep
    rkStep = mod(rkStep,rkMaxStep)+1
    subDt = rkb(rkStep)*rdt
  
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
  
    use modglobal, only : i1,j1,kmax,nsv,rdt,subDt,e12min,ih,jh,rslabs,kcb
    use modfields, only : u0,um,up,v0,vm,vp,w0,wm,wp,&
                          thl0,thlm,thlp,qt0,qtm,qtp,&
                          e120,e12m,e12p,sv0,svm,svp
    implicit none
 
    select case (iTimeInt)
 
      case(iTimeWicker,iTimeTVD)
        u0  (2:i1,2:j1,1:kmax) = um  (2:i1,2:j1,1:kmax) + subDt*up  (2:i1,2:j1,1:kmax)
        v0  (2:i1,2:j1,1:kmax) = vm  (2:i1,2:j1,1:kmax) + subDt*vp  (2:i1,2:j1,1:kmax)
        w0  (2:i1,2:j1,1:kmax) = wm  (2:i1,2:j1,1:kmax) + subDt*wp  (2:i1,2:j1,1:kmax)
        thl0(2:i1,2:j1,1:kmax) = thlm(2:i1,2:j1,1:kmax) + subDt*thlp(2:i1,2:j1,1:kmax)
        qt0 (2:i1,2:j1,1:kmax) = qtm (2:i1,2:j1,1:kmax) + subDt*qtp (2:i1,2:j1,1:kmax)
        e120(2:i1,2:j1,1:kmax) = e12m(2:i1,2:j1,1:kmax) + subDt*e12p(2:i1,2:j1,1:kmax)
        if (nsv>0) sv0(2:i1,2:j1,1:kmax,1:nsv) = svm(2:i1,2:j1,1:kmax,1:nsv) + subDt*svp(2:i1,2:j1,1:kmax,1:nsv)
    
        if (rkStep == rkMaxStep) then
          um = u0
          vm = v0
          wm = w0
          thlm = thl0
          qtm  = qt0
          e12m = e120
          svm = sv0
        end if
  
      case (iTimeLowStor)
        u0  (2:i1,2:j1,1:kmax) = u0  (2:i1,2:j1,1:kmax) + subDt*up  (2:i1,2:j1,1:kmax)
        v0  (2:i1,2:j1,1:kmax) = v0  (2:i1,2:j1,1:kmax) + subDt*vp  (2:i1,2:j1,1:kmax)
        w0  (2:i1,2:j1,1:kmax) = w0  (2:i1,2:j1,1:kmax) + subDt*wp  (2:i1,2:j1,1:kmax)
        thl0(2:i1,2:j1,1:kmax) = thl0(2:i1,2:j1,1:kmax) + subDt*thlp(2:i1,2:j1,1:kmax)
        qt0 (2:i1,2:j1,1:kmax) = qt0 (2:i1,2:j1,1:kmax) + subDt*qtp (2:i1,2:j1,1:kmax)
        e120(2:i1,2:j1,1:kmax) = e120(2:i1,2:j1,1:kmax) + subDt*e12p(2:i1,2:j1,1:kmax)
        if (nsv>0) sv0(2:i1,2:j1,1:kmax,1:nsv) = sv0(2:i1,2:j1,1:kmax,1:nsv) + subDt*svp(2:i1,2:j1,1:kmax,1:nsv)
  
      case default
    end select

    ! Set a lower limit value for e120 (not sure why this is required...)
    where (e120 < e12min)
      e120 = e12min
    end where

    ! Set the appropriate offsets for the tendencies (for Wicker and Skamarock, rka(:)=0 )
    up   = rka(rkStep)*up
    vp   = rka(rkStep)*vp
    wp   = rka(rkStep)*wp
    thlp = rka(rkStep)*thlp
    qtp  = rka(rkStep)*qtp
    e12p = rka(rkStep)*e12p
    if (nsv>0) svp = rka(rkStep)*svp
     
  end subroutine tstep_integrate

  !======== Determine the maximum Courant number
  real function getCourant()
    use modfields, only : u0,v0,w0
    use modglobal, only : i1,j1,kmax,dx,dy,dzh,rdt
    use modmpi,    only : my_real,mpi_max,comm3d,mpierr
    implicit none
    real    :: courantl
    integer :: k

    ! Determine the maximum Courant number
    courantl = 0.
    do k=1,kmax
      courantl = max( courantl, &
                      maxval( u0(2:i1,2:j1,k)*u0(2:i1,2:j1,k)/(dx*dx) + &
                              v0(2:i1,2:j1,k)*v0(2:i1,2:j1,k)/(dy*dy) + &
                              w0(2:i1,2:j1,k)*w0(2:i1,2:j1,k)/(dzh(k)*dzh(k)) )*rdt*rdt )
    end do
    courantl = sqrt(courantl)
    call MPI_ALLREDUCE(courantl,getCourant,1,MY_REAL,MPI_MAX,comm3d,mpierr)

  end function getCourant

  !======== Determine the maximum Peclet number
  real function getPeclet()
    use modfields, only : ekm,ekh,rhobf
    use modglobal, only : i1,j1,kmax,dx,dx2i,dy,dy2i,dzh,rdt
    use modmpi,    only : my_real,mpi_max,comm3d,mpierr
    implicit none
    integer :: k
    real    :: pecletl

    pecletl = 1.e-5
    ! The original version of the Peclet number:
    do k=1,kmax
      pecletl = max( pecletl, &
                     maxval(ekm(2:i1,2:j1,k)/rhobf(k))*rdt / minval((/dzh(k),dx,dy/))**2 )
    end do
    ! New version of the Peclet number. Note that now this version uses ekh,
    ! which is approximately ekm*3. The Peclet number should therefore be 3
    ! times as large.
    !do k=1,kmax
    !  pecletl = max( pecletl, &
    !                 maxval(ekh(2:i1,2:j1,k)/rhobf(k)*rdt*dx2i), &
    !                 maxval(ekh(2:i1,2:j1,k)/rhobf(k)*rdt*dy2i), &
    !                 maxval(ekh(2:i1,2:j1,k)/rhobf(k)*rdt/(dzh(k)*dzh(k))))
    !end do

    call MPI_ALLREDUCE(pecletl,getPeclet,1,MY_REAL,MPI_MAX,comm3d,mpierr)

  end function getPeclet

end module modtstep
