!> \file daleslib.f90
!! Library API for Dales.
!  This file is part of DALESLIB.
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

! Initializes the libary's dales instance. This routine requires a path to the namoptions file to be passed 
! as an argument, and optionally the mpi communicator can be set too (this is particularly useful when 
! nesting into other parallelized codes). If the latter argument is not provided, the communicator will be 
! assumed to be MPI_COMM_WORLD. 

module daleslib

    implicit none

    integer, parameter :: FIELDID_U=1
    integer, parameter :: FIELDID_V=2
    integer, parameter :: FIELDID_W=3
    integer, parameter :: FIELDID_THL=4
    integer, parameter :: FIELDID_QT=5
    
    integer :: my_task,master_task

    real, allocatable :: u_tend(:)
    real, allocatable :: v_tend(:)
    real, allocatable :: thl_tend(:)
    real, allocatable :: qt_tend(:)
    real :: ps_tend
    
    contains

        subroutine initialize(path,mpi_comm)

            !!----------------------------------------------------------------
            !!     0.0    USE STATEMENTS FOR CORE MODULES
            !!----------------------------------------------------------------
            use modmpi,             only : initmpicomm,myid
            use modstartup,         only : startup

            !----------------------------------------------------------------
            !     0.1     USE STATEMENTS FOR ADDONS STATISTICAL ROUTINES
            !----------------------------------------------------------------
            use modcape,            only : initcape
            use modchecksim,        only : initchecksim
            use modstat_nc,         only : initstat_nc
            !use modspectra2,       only : initspectra2
            use modtimestat,        only : inittimestat
            use modgenstat,         only : initgenstat
            use modradstat,         only : initradstat
            use modlsmstat,         only : initlsmstat
            use modsampling,        only : initsampling
            use modquadrant,        only : initquadrant
            use modcrosssection,    only : initcrosssection 
            use modAGScross,        only : initAGScross
            use modlsmcrosssection, only : initlsmcrosssection
            use modcloudfield,      only : initcloudfield
            use modfielddump,       only : initfielddump
            use modsamptend,        only : initsamptend

            use modbulkmicrostat,   only : initbulkmicrostat
            use modbudget,          only : initbudget
            use modheterostats,     only : initheterostats

            ! modules below are disabled by default to improve compilation time
            !use modstress,         only : initstressbudget

            !use modtilt,           only : inittilt
            !use modparticles,      only : initparticles
            use modnudge,           only : initnudge
            !use modprojection,     only : initprojection
            use modchem,            only : initchem
            use modcanopy,          only : initcanopy
            implicit none

            character(len=256), intent(in)  :: path
            integer, intent(in), optional   :: mpi_comm

            !----------------------------------------------------------------
            !     0      INITIALIZE MPI COMMUNICATOR
            !----------------------------------------------------------------

            call initmpicomm(mpi_comm)

            !----------------------------------------------------------------
            !     1      READ NAMELISTS,INITIALISE GRID, CONSTANTS AND FIELDS
            !----------------------------------------------------------------

            call startup(path)

            !---------------------------------------------------------
            !      2     INITIALIZE STATISTICAL ROUTINES AND ADD-ONS
            !---------------------------------------------------------
            call initchecksim
            call initstat_nc   ! Should be called before stat-routines that might do netCDF
            call inittimestat  ! Timestat must preceed all other timeseries that could write in the same netCDF file (unless stated otherwise
            call initgenstat   ! Genstat must preceed all other statistics that could write in the same netCDF file (unless stated otherwise
            !call inittilt
            call initsampling
            call initquadrant
            call initcrosssection
            call initAGScross
            call initlsmcrosssection
            !call initprojection
            call initcloudfield
            call initfielddump
            call initsamptend
            call initradstat
            call initlsmstat
            !call initparticles
            call initnudge
            call initbulkmicrostat
            call initbudget
            !call initstressbudget
            call initchem
            call initheterostats
            call initcanopy

            !call initspectra2
            call initcape

            !Set additional library information
            my_task=myid
            master_task=0

            call initdaleslib
            
        end subroutine initialize


        ! allocate arrays for tendencies set through the daleslib interface        
        subroutine initdaleslib
          use modglobal, only: kmax
          ! use modforces, only: lforce_user
          implicit none

          allocate(u_tend(1:kmax))
          allocate(v_tend(1:kmax))
          allocate(thl_tend(1:kmax))
          allocate(qt_tend(1:kmax))
          u_tend = 0
          v_tend = 0
          thl_tend = 0
          qt_tend = 0
          ! lforce_user = .true.
          ps_tend = 0
        end subroutine initdaleslib

                ! deallocate arrays for tendencies 
        subroutine exitdaleslib
          implicit none

          deallocate(u_tend)
          deallocate(v_tend)
          deallocate(thl_tend)
          deallocate(qt_tend)
          
        end subroutine exitdaleslib

        subroutine force_tendencies
          use modglobal,   only : i1,j1,kmax
          use modfields,   only : up,vp,thlp,qtp

          implicit none
          integer k

          do k=1,kmax
             up  (2:i1,2:j1,k) = up  (2:i1,2:j1,k) + u_tend(k) 
             vp  (2:i1,2:j1,k) = vp  (2:i1,2:j1,k) + v_tend(k) 
             thlp(2:i1,2:j1,k) = thlp(2:i1,2:j1,k) + thl_tend(k)
             qtp (2:i1,2:j1,k) = qtp (2:i1,2:j1,k) + qt_tend(k)            
          enddo
        end subroutine force_tendencies



        
        ! Performs a single time step. This is exactly the code inside the time loop of the main program. 
        ! If the adaptive time stepping is enabled, the new time stamp is not guaranteed to be a dt later.

        subroutine step
            !!----------------------------------------------------------------
            !!     0.0    USE STATEMENTS FOR CORE MODULES
            !!----------------------------------------------------------------
            use modstartup,         only : writerestartfiles
            use modtimedep,         only : timedep
            use modboundary,        only : boundary, grwdamp! JvdD ,tqaver
            use modthermodynamics,  only : thermodynamics
            use modmicrophysics,    only : microsources
            use modsurface,         only : surface
            use modsubgrid,         only : subgrid
            use modforces,          only : forces, coriolis, lstend
            use modradiation,       only : radiation
            use modpois,            only : poisson
            !use modedgecold,       only : coldedge

            !----------------------------------------------------------------
            !     0.1     USE STATEMENTS FOR ADDONS STATISTICAL ROUTINES
            !----------------------------------------------------------------
            use modcape,            only : docape
            use modchecksim,        only : checksim
            !use modspectra2,       only : dospecs, tanhfilter
            use modtimestat,        only : timestat
            use modgenstat,         only : genstat
            use modradstat,         only : radstat
            use modlsmstat,         only : lsmstat
            use modsampling,        only : sampling
            use modquadrant,        only : quadrant
            use modcrosssection,    only : crosssection
            use modAGScross,        only : AGScross
            use modlsmcrosssection, only : lsmcrosssection
            use modcloudfield,      only : cloudfield
            use modfielddump,       only : fielddump
            use modsamptend,        only : samptend,tend_start,tend_adv,tend_subg,tend_force,&
                tend_rad,tend_ls,tend_micro,tend_topbound,tend_pois,tend_addon,tend_coriolis,leibniztend

            use modbulkmicrostat,   only : bulkmicrostat
            use modbudget,          only : budgetstat
            use modheterostats,     only : heterostats

            ! modules below are disabled by default to improve compilation time
            !use modstress,         only : stressbudgetstat

            !use modtilt,           only : tiltedgravity, tiltedboundary
            !use modparticles,      only : particles
            use modnudge,           only : nudge
            !use modprojection,     only : projection
            use modchem,            only : twostep
            use modcanopy,          only : canopy

            implicit none


                      
            call tstep_update                           ! Calculate new timestep
            call timedep
            call update_ps
            
            call samptend(tend_start,firstterm=.true.)

            !-----------------------------------------------------
            !   3.1   RADIATION         
            !-----------------------------------------------------
            call radiation !radiation scheme
            call samptend(tend_rad)

            !-----------------------------------------------------
            !   3.2   THE SURFACE LAYER
            !-----------------------------------------------------
            call surface

            !-----------------------------------------------------
            !   3.3   ADVECTION AND DIFFUSION
            !-----------------------------------------------------
            call advection
            call samptend(tend_adv)
            call subgrid
            call canopy
            call samptend(tend_subg)

            !-----------------------------------------------------
            !   3.4   REMAINING TERMS
            !-----------------------------------------------------
            call coriolis !remaining terms of ns equation
            call samptend(tend_coriolis)
            call forces !remaining terms of ns equation
            call force_tendencies         ! NOTE - not standard DALES, these are our own tendencies
            call samptend(tend_force)

            call lstend !large scale forcings
            call samptend(tend_ls)
            call microsources !Drizzle etc.
            call samptend(tend_micro)

            !------------------------------------------------------
            !   3.4   EXECUTE ADD ONS
            !------------------------------------------------------
            call nudge
            !    call dospecs
            !    call tiltedgravity

            call samptend(tend_addon)

            !-----------------------------------------------------------------------
            !   3.5  PRESSURE FLUCTUATIONS, TIME INTEGRATION AND BOUNDARY CONDITIONS
            !-----------------------------------------------------------------------
            call grwdamp !damping at top of the model
            !JvdD    call tqaver !set thl, qt and sv(n) equal to slab average at level kmax
            call samptend(tend_topbound)
            call poisson
            call samptend(tend_pois,lastterm=.true.)

            call tstep_integrate                        ! Apply tendencies to all variables
            call boundary
            !call tiltedboundary
            !-----------------------------------------------------
            !   3.6   LIQUID WATER CONTENT AND DIAGNOSTIC FIELDS
            !-----------------------------------------------------
            call thermodynamics
            call leibniztend
            !-----------------------------------------------------
            !   3.7  WRITE RESTARTFILES AND DO STATISTICS
            !------------------------------------------------------
            call twostep
            !call coldedge
            call checksim
            call timestat  !Timestat must preceed all other timeseries that could write in the same netCDF file (unless stated otherwise
            call genstat  !Genstat must preceed all other statistics that could write in the same netCDF file (unless stated otherwise
            call radstat
            call lsmstat
            call sampling
            call quadrant
            call crosssection
            call AGScross
            call lsmcrosssection
            !call tanhfilter
            call docape
            !call projection
            call cloudfield
            call fielddump
            !call particles

            call bulkmicrostat
            call budgetstat
            !call stressbudgetstat
            call heterostats

            call writerestartfiles


        end subroutine step

        ! update the surface pressure according to ps_tend
        ! perform the update every full timestep to avoid introducing ps0, psm
        subroutine update_ps
          use modglobal,          only: rdt, rk3step, lmoist
          use modsurfdata,        only : ps, qts
          use modsurface,         only : qtsurf


          if(rk3step == 1) then
             ps = ps + ps_tend * rdt
             
             ! modtimedep calls qtsurf to update surface values after changing ps, so we'll do the same
             if (lmoist .and. ps_tend /= 0) then
                call qtsurf
             else
                qts = 0.
             endif
          endif
        end subroutine update_ps
        




        
        ! Performs repeatedly time stepping until the tstop (measured in seconds after the cold start) is reached. 
        ! Note that the resulting model time may be a bit later than tstop.

        subroutine run_to(tstop)

            use modglobal,          only: timee,longint,rk3step

            implicit none

            integer(kind=longint),intent(in) :: tstop

            do while (timee < tstop .or. rk3step < 3)
                call step
            end do

        end subroutine run_to

        ! Returns the model time, measured in s after the cold start.

        function get_model_time() result(ret)

            use modglobal,          only: timee,longint

            implicit none

            integer(kind=longint) :: ret

            ret=timee

        end function get_model_time

        ! Cleans up; closes file handles ans deallocates arrays. This is a copy of the code 
        ! in the main program after the time loop. 

        subroutine finalize

            !!----------------------------------------------------------------
            !!     0.0    USE STATEMENTS FOR CORE MODULES
            !!----------------------------------------------------------------
            use modstartup,         only : exitmodules

            !----------------------------------------------------------------
            !     0.1     USE STATEMENTS FOR ADDONS STATISTICAL ROUTINES
            !----------------------------------------------------------------
            use modcape,            only : exitcape
            use modgenstat,         only : exitgenstat
            use modradstat,         only : exitradstat
            use modlsmstat,         only : exitlsmstat
            use modsampling,        only : exitsampling
            use modquadrant,        only : exitquadrant
            use modcrosssection,    only : exitcrosssection  
            use modAGScross,        only : exitAGScross
            use modlsmcrosssection, only : exitlsmcrosssection
            use modcloudfield,      only : cloudfield
            use modfielddump,       only : exitfielddump
            use modsamptend,        only : exitsamptend

            use modbulkmicrostat,   only : exitbulkmicrostat
            use modbudget,          only : exitbudget
            use modheterostats,     only : exitheterostats

            ! modules below are disabled by default to improve compilation time
            !use modstress,         only : exitstressbudget

            !use modtilt,           only : exittilt
            !use modparticles,      only : exitparticles
            use modnudge,           only : exitnudge
            use modcanopy,          only : exitcanopy

            implicit none

            !--------------------------------------------------------
            !    4    FINALIZE ADD ONS AND THE MAIN PROGRAM
            !-------------------------------------------------------
            call exitgenstat
            call exitradstat
            call exitlsmstat
            !call exitparticles
            call exitnudge
            call exitsampling
            call exitquadrant
            call exitsamptend
            call exitbulkmicrostat
            call exitbudget
            !call exitstressbudget
            call exitcrosssection
            call exitAGScross
            call exitlsmcrosssection
            call exitcape
            call exitfielddump
            call exitheterostats
            call exitcanopy
            call exitmodules
            call exitdaleslib

        end subroutine finalize

        function grid_shape() result(s)

            use modglobal, only: itot,jtot,kmax

            implicit none

            integer:: s(3)

            s=(/itot,jtot,kmax/)

        end function

        function allocate_z_axis(a) result(ierr)

            use modglobal, only: kmax

            implicit none

            real, allocatable :: a(:)
            integer           :: ierr

            ierr=0
            
            if(allocated(a)) then
                deallocate(a,stat=ierr)
            endif
            
            if(ierr/=0) then
                return
            endif

            allocate(a(kmax),stat=ierr)

        end function

        function allocate_2d(a) result(ierr)

            use modglobal, only: itot,jtot

            implicit none

            real, allocatable :: a(:,:)
            integer           :: ierr

            ierr=0
            
            if(allocated(a)) then
                deallocate(a,stat=ierr)
            endif
            
            if(ierr/=0) then
                return
            endif

            allocate(a(itot,jtot),stat=ierr)

        end function
        
        function allocate_3d(a) result(ierr)

            use modglobal, only: itot,jtot,kmax

            implicit none

            real, allocatable :: a(:,:,:)
            integer           :: ierr

            ierr=0
            
            if(allocated(a)) then
                deallocate(a,stat=ierr)
            endif
            
            if(ierr/=0) then
                return
            endif

            allocate(a(itot,jtot,kmax),stat=ierr)

        end function


    ! Retrieves the z-layer average of the given array
    ! Al : input field, 3d.  
    ! Ag : output vector, 1d, in Z direction. Assumed to be as high as Al
    ! NOTE averages the FULL input array - if less is wanted, pass a slice !
    function gatherlayeravg(Al,Ag) result(ret)
      use mpi
      use modmpi, only: comm3d, my_real, mpierr, myid, nprocs
      real, intent(in)      :: Al(:,:,:)
      real, intent(out)     :: Ag(:)
      integer               :: k, nk, ret

      nk = size(Al, 3)
      Ag = (/ (sum(Al(:,:,k)), k=1,nk) /)     ! sum layers of Al
      
      !in-place reduction
      if (myid == 0) then
         CALL mpi_reduce(MPI_IN_PLACE, Ag, nk, MY_REAL, MPI_SUM, 0, comm3d, ret)
      else
         CALL mpi_reduce(          Ag, Ag, nk, MY_REAL, MPI_SUM, 0, comm3d, ret)
      endif
      
      if (myid == 0) then
         Ag = Ag / (size(Al,1) * size(Al,2) * nprocs)
      endif
    end function gatherlayeravg

    ! map global indices to local
    ! gi = 0...itot-1  <-- zero-based global indices, good for numpy
    ! gj = 0...jtot-1
    ! i = 2...i1  <- two-based local indices, good for accessing the fields
    ! j = 2...j1
    function localindex(gi, gj, gk, i, j, k) result(ret)
      use modglobal, only: imax, jmax, kmax, i1, j1
      use modmpi, only: myidx, myidy
      integer, intent(in) :: gi, gj, gk
      integer, intent(out) :: i, j, k
      integer :: ret

      i = gi - myidx * imax + 2
      j = gj - myidy * jmax + 2 
      k = gk + 1
      
      if (i >= 2 .and. i <= i1 .and. j >= 2 .and. j <= j1 .and. k >= 1 .and. k <= kmax) then
         ret = 1  ! this point belongs to us
      else
         ret = 0
      endif
      
    end function localindex
    
    
    ! getter function for 3D field data
    ! g_i, g_j, g_k are index arrays
    ! data from the specified field is extracted and returned in a
    ! a(r) = field(g_(r), g_j(r), g_k(r))
    !
    ! Assumptions: the indexing arrays are 0-based, as generated by numpy
    !              the data is stored in field(1:imax,1:jmax,1:kmax)
    ! the array should be sliced to contain ONLY the physical cells before calling,
    ! so the data arrays are 1-based
    function gathervol(g_i,g_j,g_k,a,n,field) result(ret)
      use modglobal, only: imax, jmax, kmax, i1, j1
      use mpi
      use modmpi, only: myidx, myidy, comm3d, my_real, myid
      
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_i,g_j,g_k
      real,    dimension(n), intent(out)  :: a
      real,    intent(in)                 :: field (:,:,:)
      integer                             :: r, i, j, k, is, js, ks, ret

      is = size(field, 1)
      js = size(field, 2)
      ks = size(field, 3)

      ! store the data that belongs to "my" chunk or 0
      ! this is reduced with SUM, yielding the full data
      do r = 1,n
         i = g_i(r) - myidx * imax + 1
         j = g_j(r) - myidy * jmax + 1 
         k = g_k(r)                + 1
         if (i >= 1 .and. i <= is .and. j >= 1 .and. j <= js .and. k >= 1 .and. k <= ks) then
            a(r) = field(i,j,k)
            !print *, '(', g_i(r), g_j(r), g_k(r), ') ', i, j, k, '(', myidx, myidy, ')'
         else
            a(r) = 0
         endif
      enddo
      
      !in-place reduction
      if (myid == 0) then
         CALL mpi_reduce(MPI_IN_PLACE, a, n, MY_REAL, MPI_SUM, 0, comm3d, ret)
      else
         CALL mpi_reduce(           a, a, n, MY_REAL, MPI_SUM, 0, comm3d, ret)
      endif
      
    end function gathervol

    ! gather the Liquid Water Path
    ! vertical sum of QL * rho * dz = integral dz * rho * q_l
    ! unit: kg / m^2
    !
    ! g_i, g_j are index arrays
    ! data from the specified field is extracted and returned in a
    ! a(r) = field(g_(r), g_j(r))
    
    ! Assumptions: the indexing arrays are 0-based, as generated by numpy
    !              the data is stored in field(1:imax,1:jmax,1:kmax)
    ! the array should be sliced to contain ONLY the physical cells before calling,
    ! so the data arrays are 1-based
    function gatherLWP(g_i,g_j,a,n,field) result(ret)
      use modglobal, only: imax, jmax, kmax, i1, j1, dzf
      use modfields, only: rhobf
      use mpi
      use modmpi, only: myidx, myidy, comm3d, my_real, myid
      
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_i,g_j
      real,    dimension(n), intent(out)  :: a
      real,    intent(in)                 :: field (:,:,:)
      integer                             :: r, i, j, k, is, js, ks, ret

      is = size(field, 1)
      js = size(field, 2)
      ks = size(field, 3)

      ! store the data that belongs to "my" chunk or 0
      ! this is reduced with SUM, yielding the full data
      do r = 1,n
         i = g_i(r) - myidx * imax + 1
         j = g_j(r) - myidy * jmax + 1 
         if (i >= 1 .and. i <= is .and. j >= 1 .and. j <= js ) then
            a(r) = sum(field(i,j,:) * dzf * rhobf) ! vertical sum of the field, weighted by cell height and density
         else
            a(r) = 0
         endif
      enddo
      
      !in-place reduction
      if (myid == 0) then
         CALL mpi_reduce(MPI_IN_PLACE, a, n, MY_REAL, MPI_SUM, 0, comm3d, ret)
      else
         CALL mpi_reduce(           a, a, n, MY_REAL, MPI_SUM, 0, comm3d, ret)
      endif
      
    end function gatherLWP
    
end module daleslib
