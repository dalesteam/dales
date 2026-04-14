!> \file modibm.f90
!! Grid-conforming Immersed Boundary Method (IBM)

!>
!!  \author Michael Koene, Delft University of Technology, 2018-2019
!!  \author Stephan de Roode, Delft University of Technology, 2018-2024
!!  \author Steven van der Linden, Delft University of Technology, 2025-
!!  \author André van Ginkel, Delft University of Technology, 2025-
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
! Copyright 2025-2026 Delft University of Technology
!

module modibm

  use iso_fortran_env, only : real64
  use modglobal,       only : rd, rv, grav, ijtot, iinput, input_netcdf, ifnamopt, checknamelisterror, fname_options
  use modprecision,    only : field_r
  use modmpi,          only : myid, comm3d, mpierr, myidx, myidy, d_mpi_bcast, excjs, d_mpi_allreduce, &
                            mpi_max, mpi_sum
  use modsurface,      only : psim, psih, calc_obl_iter
  use modibmdata,      only : lapply_ibm,lpoislast, lwallheat, &
                            thlwall, qtwall, thlroof, qtroof, thlibm, qtibm, &
                            z0m_wall, z0h_wall, fluid_mask
  use modtimer
  use modlogging,      only : finish, warning, message
  
  implicit none
  save
  private
  character(len=*), parameter :: modname = 'modibm'
  
  public :: ixw_p, ixw_m, iyw_p, iyw_m, izw_p, iobst
  
  integer :: Nobst, Nobst_wide                                                                    !< Number of obstacles on pure domain and incl. halo cells 
  integer :: Nxwalls_plus, Nywalls_plus, Nzwalls_plus, Nxwalls_min, Nywalls_min, Nzwalls_min      !< Number of walls oriented in positve/negative x,y,z directions
  integer, allocatable :: ixw_p(:,:), ixw_m(:,:), iyw_p(:,:), iyw_m(:,:), izw_p(:,:)!, izw_m(:,:) !< Indices of walls oriented in positve/negative x,y,z directions
  integer, allocatable :: iobst(:,:)                                                              !< Indices of obstacles
  
  real(field_r) :: dx_half, dy_half, Cm_xwall, Cm_ywall, Cd_xwall, Cd_ywall, Cm_zwall, Cd_zwall, z_MO !< Additional variables/parameters

  public :: ibm_read_namelist, initibm, exitibm, applyibm, zerowallvelocity

contains

  !> Read ibm namelist entries and broadcast settings.
  subroutine ibm_read_namelist(nml_filename)
    use fortran_support, only: nnml_output
    
    character(len=*), intent(in) :: nml_filename
    
    character(len=*), parameter :: routine = modname//'/ibm_read_namelist'

    integer :: ierr

    ! Read in NAMOPTIONS parameters related to IBM
    namelist/IBM/ &
      lapply_ibm, lwallheat, thlwall, thlibm, thlroof, qtibm, lpoislast, z0m_wall, z0h_wall

    if (myid == 0) then

      open (ifnamopt, file=fname_options, status='old', iostat=ierr)
      read (ifnamopt, IBM, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'IBM')
      write (nnml_output, IBM)
      close (ifnamopt)

    end if

    call d_mpi_bcast(lapply_ibm         ,    1, 0, comm3d, mpierr)
    call d_mpi_bcast(lwallheat          ,    1, 0, comm3d, mpierr)
    call d_mpi_bcast(thlwall            ,    1, 0, comm3d, mpierr)
    call d_mpi_bcast(thlibm             ,    1, 0, comm3d, mpierr)
    call d_mpi_bcast(qtibm              ,    1, 0, comm3d, mpierr)
    call d_mpi_bcast(thlroof            ,    1, 0, comm3d, mpierr)
    call d_mpi_bcast(qtroof             ,    1, 0, comm3d, mpierr)
    call d_mpi_bcast(lpoislast          ,    1, 0, comm3d, mpierr)
    call d_mpi_bcast(z0m_wall           ,    1, 0, comm3d, mpierr)
    call d_mpi_bcast(z0h_wall           ,    1, 0, comm3d, mpierr)

  end subroutine

  !> Initializes the immersed boundary method by reading obstacles and placing them on the grid.
  subroutine initibm

    use modglobal,        only : zh, zf, itot, jtot, ih, i1, i2, jh, j1, j2, k1, imax, jmax, kmax, cexpnr, ifinput, &
                                nsv, cu, cv, ijtot, &
                                iadv_mom,iadv_tke,iadv_thl,iadv_qt,iadv_sv,iadv_cd2, &
                                ibas_prf, &
                                dx,dy,fkar
    use modsurface,       only : lmostlocal
    use modsubgriddata,   only : lanisotrop, lsmagorinsky
    use fortran_support,  only : nnml_output

    implicit none

    integer         :: i, j, k, ierr, no = 0
    integer         :: advarr(5)
    character(100)  :: readstring

    character(len=*), parameter :: routine = modname//'/initibm'

    ! Temporary fields and profiles for the processing of the IBM input
    real(field_r), allocatable :: bc_height(:,:)                                                                            !< Height of immersed boundary at cell center i,j
    integer,       allocatable :: tiobst(:,:), tixw_p(:,:), tixw_m(:,:), tiyw_p(:,:), tiyw_m(:,:), tizw_p(:,:)!, izw_m(:,:) !< Indices of walls oriented p/n x,y,z

    call timer_tic('modibm/initibm',0)

    if( myid==0 ) then

      ! Do some checks for conflicting settings and warn/stop further execution
      if ( lapply_ibm ) then
        if (abs(cu) > 0) call finish(routine, 'Domain translation not allowed with IBM, set cu to zero')  
        if (abs(cv) > 0) call finish(routine, 'Domain translation not allowed with IBM, set cv to zero')  

        if (ibas_prf .ne. 2) then
          ibas_prf = 2
          call warning(routine, &
              'ibas_pr is overwritten to 2 (Boussinesq approximation with constant density) ' // &
              'height dependent density gives probles with correction of vertical advective '  // &
              'tendencies at the top of obstacles'  &
          )  
        end if

        ! TODO: implement kappa advection for tracers
        advarr = (/iadv_mom,iadv_tke,iadv_thl,iadv_qt,iadv_sv/)
        if (any(advarr/=iadv_cd2)) then
          call finish(routine, &
              'Current IBM implementation only works with 2nd order advection. '  // &
              'Proper check for kappa advection of scalars to be implemented'&
          )  
        end if

        if ( lanisotrop ) then
          call warning(routine, &
              'WARNING: you are using IBM with anisotropic grids in x,y-direction '  // & 
              'while possible, this may cause unexpected results (blending effects)' &
          )  
        end if

        if ( .not.(lsmagorinsky) ) then
          call warning(routine, &
              'WARNING: subgrid TKE (e120) is not (yet) explicitly corrected for walls. '  // &
              'This includes production, destruction and dissipation terms in sgs-tke budget' &
          )  
        end if

        !! Check doesn't work currently, related to order in startup routine
        ! if ( .not. lmostlocal ) then
        !   stop 'ERROR: you must use local monin-obukhov with IBM'
        ! end if

      end if

    end if

    ! Step out of further subroutine when IBM is switched off
    if (.not. (lapply_ibm)) then
       call timer_toc('modibm/initibm')
       return
    endif

    ! Calculate law-of-wall coefficients for vertical walls (constant; no stability correction on vertical walls)
    dx_half = 0.5 * dx
    dy_half = 0.5 * dy
    Cm_xwall = (fkar/(log(dx_half/z0m_wall)))**2
    Cm_ywall = (fkar/(log(dy_half/z0m_wall)))**2

    ! TODO extend later as well with more options: constant flux, reactive wall temperature
    ! Set thermal boundary conditions for vertical walls to either Dirichlet (lwallheat=.true.) or zero flux
    if (lwallheat) then
      Cd_xwall = fkar**2 / log(dx_half/z0m_wall) / log(dx_half/z0h_wall)
      Cd_ywall = fkar**2 / log(dy_half/z0m_wall) / log(dy_half/z0h_wall)
    else
      Cd_xwall = 0.
      Cd_ywall = 0.
    end if

    allocate(bc_height(itot+1,jtot+1))             ! use itot+1, jtot+1 and start writing at index = 2 to conform to field indices
    allocate(fluid_mask(2-ih:i1+ih,2-jh:j1+jh,k1))

    bc_height(:,:)  = 0.
    fluid_mask (:,:,:)    = .true.

    ! Definition of obstacles
    if (myid==0) then
      if (iinput == input_netcdf) then

        call init_ibm_from_nc(bc_height)

      else

       call message(routine, 'Reading inputfile ibm.inp.', cexpnr) 

        open (ifinput,file='ibm.inp.'//cexpnr)
          do k=1,7
            read (ifinput,'(a100)') readstring
            write (6,*) readstring
          end do

          do j=jtot+1,2,-1
            do i=2,itot+1
              read(ifinput,'(F6.1)') bc_height(i,j)
            end do
          end do

        close(ifinput)

        call message(routine, 'Succesfully read inputfile in modibm') 

      end if
    end if 

    call d_mpi_bcast(bc_height, (itot+1)*(jtot+1), 0, comm3d, mpierr)

    ! Determine obstacle cells. Checks if obstacle height is above midpoint of vertical cell (= full levels). Corresponds to >50% of cell being filled.
    Nobst = 0
    do i=2,i1 
      do j=2,j1
         do k=1,kmax
            if (zf(k) <= bc_height(i+myidx*imax,j+myidy*jmax)) then
              fluid_mask(i,j,k) = .false.                             ! Set grid cell to obstacle
              Nobst             = Nobst + 1                           ! Increase counter of obstacles by one
           end if
        end do
      end do
    end do

    call excjs(fluid_mask  , 2,i1,2,j1,1,k1,ih,jh)

    call message(routine, 'Start determination of wall positions')

    !> \brief Identify sidewalls based on fluid_mask
    !!  u-positions (with index i) are to the left of grid center (with index i)
    !!  walls with normal in minus x-direction will therefore be at same index, walls with normal in plus direction at i+1
    !!  see visual example below
    !!
    !! building positions   :   0     X     X     X     0
    !! grid centered index  :  i-2   i-1    i    i+1   i+2
    !! fluid_mask           :   F     T     T     T     F
    !! xwall min yes/no     :   F     T     F     F     F
    !! xwall plus yes/no    :   F     F     F     F     T

    ! Preset (temporary) arrays for internal buildings and wall indices [Nobst+2*ih*k1+2*jh*k1+4*ih*jh] provides upper bound)
    allocate(tiobst(Nobst+2*ih*j1*k1+2*i1*jh*k1+4*ih*jh,3))                                                                     ! for internal building points
    allocate(tixw_p(Nobst+i2*k1,3), tixw_m(Nobst+i2*k1,3), tiyw_p(Nobst+j2*k1,3), tiyw_m(Nobst+j2*k1,3), tizw_p(Nobst+j2*k1,3)) ! for x- and y-walls in positive and negative directions

    Nxwalls_plus = 0; Nxwalls_min = 0; Nywalls_plus = 0; Nywalls_min = 0; Nzwalls_plus = 0

    Nobst_wide = 0;
    do k=1,kmax
      do j=1,j2
        do i=1,i2

          if ( .not.(fluid_mask(i,j,k)) ) then                    ! check if internal obstable point
            Nobst_wide           = Nobst_wide + 1                 ! local counter for obstacle points
            tiobst(Nobst_wide,1) = i
            tiobst(Nobst_wide,2) = j
            tiobst(Nobst_wide,3) = k
          else                                                    ! else check for obstacles next to fluid

            if ( .not.(fluid_mask(i-1,j,k)) ) then                ! obstacle left of fluid -> wall with normal in positive x-direction
              Nxwalls_plus           = Nxwalls_plus + 1
              tixw_p(Nxwalls_plus,1) = i
              tixw_p(Nxwalls_plus,2) = j
              tixw_p(Nxwalls_plus,3) = k
            end if

            if ( .not.(fluid_mask(i+1,j,k)) ) then                ! obstacle right of fluid -> wall with normal in negative x-direction
              Nxwalls_min           = Nxwalls_min + 1
              tixw_m(Nxwalls_min,1) = i
              tixw_m(Nxwalls_min,2) = j
              tixw_m(Nxwalls_min,3) = k
            end if

            if ( .not.(fluid_mask(i,j-1,k)) ) then                ! obstacle in front of fluid -> wall with normal in positive y-direction
              Nywalls_plus           = Nywalls_plus + 1
              tiyw_p(Nywalls_plus,1) = i
              tiyw_p(Nywalls_plus,2) = j
              tiyw_p(Nywalls_plus,3) = k
            end if

            if ( .not.(fluid_mask(i,j+1,k)) ) then                ! obstacle at back of fluid -> wall with normal in negative y-direction
              Nywalls_min           = Nywalls_min + 1
              tiyw_m(Nywalls_min,1) = i
              tiyw_m(Nywalls_min,2) = j
              tiyw_m(Nywalls_min,3) = k
            end if

            if (k == 1)  cycle                                    ! if near surface, go to next loop iteration (as fluid_mask(i,j,0) doesn't exist)

            if ( .not.(fluid_mask(i,j,k-1)) ) then                ! obstacle below fluid -> wall with normal in positive z-direction
              Nzwalls_plus            = Nzwalls_plus  + 1
              tizw_p(Nzwalls_plus, 1) = i
              tizw_p(Nzwalls_plus, 2) = j
              tizw_p(Nzwalls_plus, 3) = k
            end if

            !! Current impementation does not allow overhanging obstacles (i.e. walls in negative z)
            ! if (k == kmax)  cycle                                 ! if near top, go to next loop iteration (as fluid_mask(i,j,kmax+1) doesn't exist)

            ! if ( .not.(fluid_mask(i,j,k+1)) ) then                ! obstacle above fluid -> wall with normal in negative z-direction
            !   Nzwalls_min            = Nzwalls_min  + 1
            !   tizw_m(Nzwalls_min, 1) = i
            !   tizw_m(Nzwalls_min, 2) = j
            !   tizw_m(Nzwalls_min, 3) = k
            ! end if

          end if 

        end do
      end do
    end do

    allocate(iobst(Nobst_wide  ,3))
    allocate(ixw_p(Nxwalls_plus,3))
    allocate(ixw_m(Nxwalls_min, 3))
    allocate(iyw_p(Nywalls_plus,3))
    allocate(iyw_m(Nywalls_min ,3))
    allocate(izw_p(Nzwalls_plus,3))

    ! Copy temporary arrays with indices into final ones
    iobst(1:Nobst_wide,  :) = tiobst(1:Nobst_wide,  :)
    ixw_p(1:Nxwalls_plus,:) = tixw_p(1:Nxwalls_plus,:)
    ixw_m(1:Nxwalls_min ,:) = tixw_m(1:Nxwalls_min ,:)
    iyw_p(1:Nywalls_plus,:) = tiyw_p(1:Nywalls_plus,:)
    iyw_m(1:Nywalls_min ,:) = tiyw_m(1:Nywalls_min ,:)
    izw_p(1:Nzwalls_plus,:) = tizw_p(1:Nzwalls_plus,:)

    deallocate(bc_height)
    deallocate(tiobst, tixw_p, tixw_m, tiyw_p, tiyw_m, tizw_p)

    !$acc enter data copyin(fluid_mask, iobst, ixw_p, ixw_m, iyw_p, iyw_m, izw_p)
!!$omp target enter data map(to:fluid_mask,iobst,ixw_p,ixw_m,iyw_p,&
!!$omp iyw_m,izw_p)

    call timer_toc('modibm/initibm')

    return
  end subroutine initibm

  !> Clears the memory for the immersed boundary method.
  subroutine exitibm

    if (.not. (lapply_ibm)) return

    !$acc exit data delete(fluid_mask, iobst, ixw_p, ixw_m, iyw_p, iyw_m, izw_p)
!!$omp target exit data map(delete:fluid_mask,iobst,ixw_p,ixw_m,iyw_p,&
!!$omp iyw_m,izw_p)

    deallocate(iobst)
    deallocate(ixw_p)
    deallocate(ixw_m)
    deallocate(iyw_p)
    deallocate(iyw_m)
    deallocate(izw_p)

    deallocate(fluid_mask)

    return
  end subroutine exitibm

  !> Helper function to read obstacles from NetCDF.
  subroutine init_ibm_from_nc(bc_height)

    use netcdf
    use modnetcdf,   only : check
    use modglobal,   only : iexpnr
    use modglobal,   only : itot, jtot

    implicit none

    character(len=*), parameter :: routine = modname//'/init_ibm_from_nc'

    real(field_r),  intent(out) :: bc_height(:,:)

    character(32) :: input_file = 'ibm.inp_xxx.nc'
    integer       :: ncid, varid, len_x, len_y

    write(input_file(9:11), '(i3.3)') iexpnr

    call message(routine, "Reading IBM input: ", input_file) 
    call message(routine, "Expecting dimensions x,y and variable bc_height(:,:)")

    call check( nf90_open(input_file, nf90_nowrite, ncid), input_file, __LINE__)

    ! Check if dimensions of ibm.inp_xxx.nc agree with the DALES domain
    call check( nf90_inq_dimid(ncid, 'x', varid), input_file, __LINE__ )
    call check( nf90_inquire_dimension(ncid, varid, len=len_x), input_file, __LINE__ )
    if (len_x /= itot) then
      call finish(routine, "STOPPED. x-dimension of ibm.inp differs from DALES domain: ", len_x, " /=", itot)
    end if

    call check( nf90_inq_dimid(ncid, 'y', varid), input_file, __LINE__ )
    call check( nf90_inquire_dimension(ncid, varid, len=len_y), input_file, __LINE__ )
    if (len_y /= jtot) then
      call finish(routine, "STOPPED. y-dimension of ibm.inp differs from DALES domain: ", len_y, " /=", jtot)
    end if

    ! Get variable bc height from nc file
    call check( nf90_inq_varid( ncid, 'bc_height', varid), input_file, __LINE__ )
    call check( nf90_get_var(ncid, varid, bc_height(2:itot+1,2:jtot+1) , &
                              count = (/itot, jtot/) ), input_file, __LINE__ )
    call check( nf90_close(ncid), input_file, __LINE__ )

    call message(routine, 'Succesfully read netCDF inputfile in modibm')

  end subroutine init_ibm_from_nc
  
  !> Applies the immersed boundary through forcings.
  subroutine applyibm

    use modfields,      only : u0, v0, w0, thl0, qt0, e120, sv0, &
                               up, vp, wp, thlp, qtp, e12p, svp, &
                               thl0av, qt0av, rhobf, rhobh
    use modglobal,      only : rk3step, kmax, i1, j1, k1, ih, jh, rdt, timee, dx, dy, dx2i, dy2i, dzh, dzhi, dzf, dzfi, zf, zh, nsv, e12min, fkar
    use modsurface,     only : lneutral, lmostlocal
    use modsubgriddata, only : ekm, ekh
    use modmpi,         only : excjs

    implicit none

    character(len=*), parameter :: routine = modname//'/applyibm'

    integer           :: i, j, k, nn, nc, retval
    real(field_r)     :: rk3coef, rk3coefi
    real(field_r)     :: emmo, empo, emom, emop, eomm
    real(field_r)     :: u_at_v_min, u_at_v_plus, v_at_u_min, v_at_u_plus
    real(field_r)     :: w_at_v_min, w_at_v_plus, w_at_u_min, w_at_u_plus
    real(field_r)     :: u_at_w_min, u_at_w_plus, v_at_w_min, v_at_w_plus
    real(field_r)     :: uspeed, ucc, vcc, z_MO
    real(field_r)     :: tau_vu_plus, tau_vu_min, tau_wu_min, tau_wu_plus, tau_uv_min, tau_uv_plus, tau_wv_min, tau_wv_plus
    real              :: Lob
    
    if (.not. lapply_ibm) return

    call timer_tic('modibm/applyibm',0)

    rk3coef = rdt / (4. - dble(rk3step))
    rk3coefi = 1. / rk3coef

    ! Set tendencies inside obstacles (i.e., correct for any drift from previous integration step)
    !$acc parallel loop gang vector default(present)
!!$omp target teams loop defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
    do nn = 1,Nobst_wide  !1!< Svdldit werkt niet..
      i = iobst(nn,1)
      j = iobst(nn,2)
      k = iobst(nn,3)

      ! Correction of velocities also corrects walls w. normals in negative x,y-directions (due to staggered grid arrangement)
      up(i,j,k)   = -u0(i,j,k)*rk3coefi
      vp(i,j,k)   = -v0(i,j,k)*rk3coefi
      wp(i,j,k)   = -w0(i,j,k)*rk3coefi
      thlp(i,j,k) = (thlibm - thl0(i,j,k))*rk3coefi
      qtp (i,j,k) = (qtibm  - qt0(i,j,k) )*rk3coefi
      e12p(i,j,k) = (e12min - e120(i,j,k))*rk3coefi   ! Maintain e12min to prevent ekm/ekh going to NaN
                                                      ! Maybe better to explicitly set eddy diffusivities to tiny value
      do nc=1,nsv
          svp(i,j,k,nc) = - sv0(i,j,k,nc)*rk3coefi
      end do

    end do

    ! All corrections consist of 2 parts: cancel the "wrong" diffusion imposed by modsubgrid at fluid points at the wall, and add wall friction/flux.
    ! In this framework, we assume no correction for advection is needed, as velocities at walls should already be close to zero 
    ! TODO: 1. do explicitly remove advective tendencies over walls (mainly for improvement of conservations of scalars)
    !       2. apply IBM conditions explicitly to e12 as well

    ! Correct tendencies for walls in positive z-direction (only works when k>1, which should be the case for vertical walls [see initibm])
    !$acc parallel loop gang vector default(present)
!!$omp target teams loop defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
    do nn = 1,Nzwalls_plus
      i = izw_p(nn,1)
      j = izw_p(nn,2)
      k = izw_p(nn,3)

      ! First set the normal velocity at the wall (correct for any drift)
      wp(i,j,k)   = -w0(i,j,k)*rk3coefi

      ! Remove "wrong" diffusive tendencies
      ! for u(i,j,k):
      emom = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i+1,j,k  ) )  + &
            dzf(k  ) * ( ekm(i,j,k-1) + ekm(i+1,j,k-1) ) ) * &
            ( 0.25_field_r * dzhi(k) )
      up(i,j,k) = up(i,j,k) + 0.5_field_r * rhobh(k)/rhobf(k) * emom * ( (u0(i,j,k) - u0(i,j,k-1)) *dzhi(k)) *dzfi(k)

      ! for u(i+1,j,k):
      emom = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i-1,j,k  ) )  + &
            dzf(k  ) * ( ekm(i,j,k-1) + ekm(i-1,j,k-1) ) ) * &
            ( 0.25_field_r * dzhi(k) )
      up(i+1,j,k) = up(i+1,j,k) + 0.5_field_r * rhobh(k)/rhobf(k) * emom * ( (u0(i+1,j,k) - u0(i+1,j,k-1)) *dzhi(k)) *dzfi(k)

      ! for v(i,j,k):
      eomm = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i,j-1,k ) ) + &
            dzf(k  ) * ( ekm(i,j,k-1) + ekm(i,j-1,k-1) ) ) * &
            ( .25_field_r * dzhi(k) )
      vp(i,j,k) = vp(i,j,k) + 0.5_field_r * rhobh(k)/rhobf(k) * eomm * ( (v0(i,j,k) - v0(i,j,k-1)) * dzhi(k)) * dzfi(k)

      ! for v(i,j+1,k):
      eomm = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i,j+1,k ) ) + &
            dzf(k  ) * ( ekm(i,j,k-1) + ekm(i,j+1,k-1) ) ) * &
            ( .25_field_r * dzhi(k) )
      vp(i,j+1,k) = vp(i,j+1,k) + 0.5_field_r * rhobh(k)/rhobf(k) * eomm * ( (v0(i,j+1,k)-v0(i,j+1,k-1)) * dzhi(k)) * dzfi(k)

      ! Include correct Monin-Obukhov wall drag (with stability correction)
      ucc    = 0.5_real64 * (u0(i,j,k) + u0(i+1,j,k))
      vcc    = 0.5_real64 * (v0(i,j,k) + v0(i,j+1,k))
      uspeed = 0.5_field_r * sqrt( ( u0(i,j,k) + u0(i+1,j,k) )**2 + ( v0(i,j,k) + v0(i,j+1,k) )**2 )
      z_MO   = zf(k) - zh(k)
      Lob    = -1.e10

      if (lneutral) then
        Cm_zwall = fkar**2 / (log(z_MO / z0m_wall))** 2
        Cd_zwall = fkar**2 / (log(z_MO / z0m_wall)) / (log(z_MO / z0h_wall))
      else
        retval = calc_obl_iter(thl0(i,j,k), qt0(i,j,k), real(thlroof), &
                               real(qtroof), z_MO, real(z0m_wall), real(z0h_wall), & 
                               ucc, vcc, Lob)

        Cm_zwall = fkar**2 / (log(z_MO / z0m_wall) - psim(z_MO / Lob) + psim(z0m_wall / Lob))** 2
        Cd_zwall = fkar**2 / (log(z_MO / z0m_wall) - psim(z_MO / Lob) + psim(z0m_wall / Lob)) / (log(z_MO / z0h_wall) - psih(z_MO / Lob) + psih(z0h_wall / Lob))
      end if

      up(i  ,j,k) = up(i  ,j,k) - 0.25_field_r * rhobh(k)/rhobf(k) * Cm_zwall * ( u0(i,j,k) + u0(i+1,j,k) ) * uspeed * dzfi(k)
      up(i+1,j,k) = up(i+1,j,k) - 0.25_field_r * rhobh(k)/rhobf(k) * Cm_zwall * ( u0(i,j,k) + u0(i+1,j,k) ) * uspeed * dzfi(k)
      vp(i,j  ,k) = vp(i,j  ,k) - 0.25_field_r * rhobh(k)/rhobf(k) * Cm_zwall * ( v0(i,j,k) + v0(i,j+1,k) ) * uspeed * dzfi(k)
      vp(i,j+1,k) = vp(i,j+1,k) - 0.25_field_r * rhobh(k)/rhobf(k) * Cm_zwall * ( v0(i,j,k) + v0(i,j+1,k) ) * uspeed * dzfi(k)

      ! tentative fix for vertical diffusion of temperature over z-walls. Also here, 0.5 comes from interpolation of ekm.
      thlp(i,j,k) = thlp(i,j,k) + 0.5_field_r * rhobh(k)/rhobf(k) * ( ( ( dzf(k-1) * ekm(i,j,k) ) + ( dzf(k ) * ekm(i,j,k-1)) )* dzhi(k)  ) * ( thl0(i,j,k) - thl0(i,j,k-1) ) * dzhi(k) * dzfi(k)
      qtp (i,j,k) = qtp (i,j,k) + 0.5_field_r * rhobh(k)/rhobf(k) * ( ( ( dzf(k-1) * ekm(i,j,k) ) + ( dzf(k ) * ekm(i,j,k-1)) )* dzhi(k)  ) * ( qt0 (i,j,k) - qt0 (i,j,k-1) ) * dzhi(k) * dzfi(k)

      ! zero/constant flux of heat/moisture currently not implemented. Implies roof temperatures become diagnostic (requiring energy balance)
      ! doubting whether rhobh(k)/rhobf(k) are truly correct here, i.e., at right positions..
      thlp(i,j,k) = thlp(i,j,k) - rhobh(k)/rhobf(k) * Cd_zwall * ( thl0(i,j,k) - thlroof ) * uspeed * dzfi(k)
      qtp (i,j,k) = qtp (i,j,k) - rhobh(k)/rhobf(k) * Cd_zwall * (  qt0(i,j,k) - qtroof  ) * uspeed * dzfi(k)
      ! no flux of e12 from wall into flow

      ! currently only for zero flux of scalars from walls (i.e., no emission at walls). To be extended later
      do nc=1,nsv
        svp(i,j,k,nc) = svp(i,j,k,nc) - 0 !rhobh(k)/rhobf(k) * Cd_zwall * ( sv0(i,j,k,nc) -svroof(nc) ) * uspeed * dzfi(k) !< currently set to zero
      end do

    end do

    ! Correct tendencies for walls in positive x-direction
    !$acc parallel loop gang vector default(present)
!!$omp target teams loop defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
    do nn = 1,Nxwalls_plus
      i = ixw_p(nn,1)
      j = ixw_p(nn,2)
      k = ixw_p(nn,3)

      ! First set the normal velocity at the wall (correct for any drift)
      up(i,j,k)   = -u0(i,j,k)*rk3coefi

      ! Remove "wrong" diffusive tendencies and replace with wall drag:
      ! essentially every wall facet is split up in two parts, where half the correction goes into v(i,j,k) and half into v(i,j+1,k)
      ! if the x-plane next to the current one is not a wall, no partial correction to v(i,j+1,k) occurs (but a partial contribution of diffusion does)
      ! if this x-plane is also a wall, the total correction on a momentum component will therefore be split over two entries in the wall list. 

      ! for v(i,j,k):
      emmo = 0.25_field_r * ( ekm(i,j,k) + ekm(i,j-1,k) + ekm(i-1,j,k) + ekm(i-1,j-1,k) )
      w_at_v_plus = 0.25_field_r * ( w0(i,j,k) + w0(i,j,k+1) + w0(i,j-1,k) + w0(i,j-1,k+1) )  !at v(i,j,k)
      tau_vu_plus = log_wallaw(v0(i,j,k), w_at_v_plus, Cm_xwall)

      vp(i,j  ,k) = vp(i,j  ,k) + 0.5_field_r * emmo * ( (v0(i,j,k) - v0(i-1,j,k) ) / dx) / dx - 0.5_field_r * tau_vu_plus / dx ! factor 0.5 originates to avoid double correction

      ! for v(i,j+1,k):
      empo = 0.25_field_r * ( ekm(i,j,k) + ekm(i,j+1,k) + ekm(i-1,j,k) + ekm(i-1,j+1,k) )
      w_at_v_plus = 0.25_field_r * (w0(i,j+1,k) + w0(i,j+1,k+1) + w0(i,j,k) + w0(i,j,k+1) )
      tau_vu_plus = log_wallaw(v0(i,j+1,k), w_at_v_plus, Cm_xwall)

      vp(i,j+1,k) = vp(i,j+1,k) + 0.5_field_r * empo * ( (v0(i,j+1,k) - v0(i-1,j+1,k) ) / dx) / dx - 0.5_field_r * tau_vu_plus / dx

      ! for w(i,j,k):
      if (k /= 1) then ! not correctable when at surface (k = 1)
        emom = ( dzf(k-1) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                dzf(k)  * ( ekm(i,j,k-1) + ekm(i-1,j,k-1) ) ) / &
                ( 4.0_field_r * dzh(k) )
        v_at_w_plus = 0.25_field_r * ( v0(i,j,k-1) + v0(i,j,k) + v0(i,j+1,k-1) + v0(i,j+1,k) )
        tau_wu_plus = log_wallaw(w0(i  ,j,k), v_at_w_plus, Cm_xwall)

        wp(i,j,k  ) = wp(i,j,k  ) + 0.5_field_r * emom * ( (w0(i,j,k  ) - w0(i-1,j,  k)) / dx ) / dx - 0.5_field_r * tau_wu_plus / dx
      end if

      ! for w(i,j,k+1):
      emop = ( dzf(k) * ( ekm(i,j,k+1)  + ekm(i-1,j,k+1)  )  + &
              dzf(k+1)  * ( ekm(i,j,k) + ekm(i-1,j,k) ) ) / &
              ( 4.0_field_r * dzh(k+1) )
      v_at_w_plus = 0.25_field_r * ( v0(i,j,k) + v0(i,j,k+1) + v0(i,j+1,k) + v0(i,j+1,k+1) )
      tau_wu_plus = log_wallaw(w0(i  ,j,k+1), v_at_w_plus, Cm_xwall)

      wp(i,j,k+1) = wp(i,j,k+1) + 0.5_field_r * emop * ( (w0(i,j,k+1) - w0(i-1,j,k+1) ) / dx ) / dx - 0.5_field_r * tau_wu_plus / dx

      ! Enforcing zero flux by correction of tendencies of temperature, moisture and other scalars (by negating diffusion term)

      ! scalars at (i,j,k) are "to the right" of the x-walls, so correct on s(i,j,k)
      ! "+" because we reflect it back, and 0.5_field_r comes from interpolation of ekh
      thlp(i,j,k) = thlp(i,j,k) + 0.5_field_r * ( ekh(i,j,k) + ekh(i-1,j,k) ) * ( thl0(i,j,k) - thl0(i-1,j,k) ) * dx2i
      qtp(i,j,k)  =  qtp(i,j,k) + 0.5_field_r * ( ekh(i,j,k) + ekh(i-1,j,k) ) * (  qt0(i,j,k) -  qt0(i-1,j,k) ) * dx2i

      do nc=1,nsv
        svp(i,j,k,nc) = svp(i,j,k,nc) + 0.5_field_r * ( ekh(i,j,k) + ekh(i-1,j,k) ) * ( sv0(i,j,k,nc) - sv0(i-1,j,k,nc) ) * dx2i
      end do
      !call xwalle12(i,j,k) ! correction is ignored assuming u,v,w,subgrid TKE are near zero inside buildings

      ! Finally, set correct heat flux from wall to fluid
      uspeed = 0.5_field_r * ( ( v0(i,j,k) + v0(i,j+1,k) )**2 + ( w0(i,j,k) + w0(i,j,k+1) )**2 )**0.5_field_r
      thlp(i,j,k) = thlp(i,j,k) + Cd_xwall * uspeed * (thlwall - thl0(i,j,k)) / dx
      qtp (i,j,k) = qtp (i,j,k) + Cd_xwall * uspeed * ( qtwall -  qt0(i,j,k)) / dx

    end do

    ! Correct tendencies for walls in negative x-direction
    !$acc parallel loop gang vector default(present)
!!$omp target teams loop defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
    do nn = 1,Nxwalls_min
      i = ixw_m(nn,1)
      j = ixw_m(nn,2)
      k = ixw_m(nn,3)

      ! First set the normal velocity at the wall (correct for any drift)
      ! note: these should actually already be corrected for in loop over iobst (due to staggered grid arrangement)
      up(i+1,j,k)   = -u0(i+1,j,k)*rk3coefi

      ! Remove "wrong" diffusive tendencies and replace with wall drag
      ! for v(i,j,k):
      emmo = 0.25_field_r * ( ekm(i+1,j,k) + ekm(i+1,j-1,k) + ekm(i,j,k) + ekm(i,j-1,k) )
      w_at_v_min  = 0.25_field_r * ( w0(i,j,k) + w0(i,j,k+1) + w0(i,j-1,k) + w0(i,j-1,k+1) )  !at v(i,j,k)
      tau_vu_min = log_wallaw(v0(i,j,k), w_at_v_min , Cm_xwall) !if v0 > 0, tau > 0, minus sign in tendency enforces opposing friction

      vp(i,j  ,k) = vp(i,j  ,k) - 0.5_field_r * emmo * ( (v0(i+1,j,k) - v0(i,j,k) ) / dx) / dx  - 0.5_field_r * tau_vu_min / dx ! factor 0.5 originates to avoid double correction

      ! for v(i,j+1,k):
      empo = 0.25_field_r * ( ekm(i+1,j,k) + ekm(i+1,j+1,k) + ekm(i,j,k) + ekm(i,j+1,k) )
      w_at_v_min  = 0.25_field_r * ( w0(i,j+1,k) + w0(i,j+1,k+1) + w0(i,j,k) + w0(i,j,k+1) )
      tau_vu_min = log_wallaw(v0(i,j+1,k), w_at_v_min , Cm_xwall)

      vp(i,j+1,k) = vp(i,j+1,k) - 0.5_field_r * empo * ( (v0(i+1,j+1,k) - v0(i,j+1,k) ) / dx) / dx - 0.5_field_r * tau_vu_min / dx

      ! for w(i,j,k):
      if (k /= 1) then ! not correctable when at surface (k = 1)
        emom = ( dzf(k-1) * ( ekm(i+1,j,k)  + ekm(i,j,k)  )  + &
                dzf(k)  * ( ekm(i+1,j,k-1) + ekm(i,j,k-1) ) ) / &
                ( 4.0_field_r * dzh(k) )
        v_at_w_min = 0.25_field_r * (v0(i,j,k-1)+v0(i,j,k)+v0(i,j+1,k-1)+v0(i,j+1,k) )
        tau_wu_min = log_wallaw(w0(i,j,k), v_at_w_min , Cm_xwall)

        wp(i,j,k  ) = wp(i,j,k  ) - 0.5_field_r * emom * ( (w0(i+1,j,k  ) - w0(i,j,  k)) / dx ) / dx - 0.5_field_r * tau_wu_min / dx
      end if

      ! for w(i,j,k+1):
      emop = ( dzf(k) * ( ekm(i+1,j,k+1)  + ekm(i,j,k+1)  )  + &
              dzf(k+1)  * ( ekm(i+1,j,k) + ekm(i,j,k) ) ) / &
              ( 4.0_field_r * dzh(k+1) )
      v_at_w_min  = 0.25_field_r * (v0(i,j,k)+v0(i,j,k+1)+v0(i,j+1,k)+v0(i,j+1,k+1) )
      tau_wu_min = log_wallaw(w0(i,j,k+1), v_at_w_min , Cm_xwall)

      wp(i,j,k+1) = wp(i,j,k+1) - 0.5_field_r * emop * ( (w0(i+1,j,k+1) - w0(i,j,k+1) ) / dx ) / dx - 0.5_field_r * tau_wu_min / dx

      ! Enforcing zero flux by correction of tendencies of temperature, moisture and other scalars (by negating diffusion term)

      ! scalars at (i+1,j,k) are "to the right" of the x-walls, so correct on s(i,j,k) which is on the left for minus x-walls
      ! "-" because we reflect it back, and 0.5_field_r comes from interpolation of ekh
      thlp(i,j,k) = thlp(i,j,k) - 0.5_field_r * ( ekh(i+1,j,k) + ekh(i,j,k) ) * ( thl0(i+1,j,k) - thl0(i,j,k) ) * dx2i
      qtp (i,j,k) = qtp (i,j,k) - 0.5_field_r * ( ekh(i+1,j,k) + ekh(i,j,k) ) * (  qt0(i+1,j,k) -  qt0(i,j,k) ) * dx2i

      do nc=1,nsv
        svp(i,j,k,nc) = svp(i,j,k,nc) - 0.5_field_r * ( ekh(i+1,j,k) + ekh(i,j,k) ) * ( sv0(i+1,j,k,nc) - sv0(i,j,k,nc) ) * dx2i
      end do

      ! Finally, set correct heat flux from wall to fluid
      uspeed = 0.5_field_r * ( ( v0(i,j,k) + v0(i,j+1,k) )**2 + ( w0(i,j,k) + w0(i,j,k+1) )**2 )**0.5_field_r
      thlp(i,j,k) = thlp(i,j,k) + Cd_xwall * uspeed * (thlwall - thl0(i,j,k)) / dx
      qtp (i,j,k) = qtp (i,j,k) + Cd_xwall * uspeed * ( qtwall -  qt0(i,j,k)) / dx ! check this one..

    end do

    ! Correct tendencies for walls in positive y-direction
    !$acc parallel loop gang vector default(present)
!!$omp target teams loop defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
    do nn = 1,Nywalls_plus
      i = iyw_p(nn,1)
      j = iyw_p(nn,2)
      k = iyw_p(nn,3)

      ! First set the normal velocity at the wall (correct for any drift)
      vp(i,j,k)   = -v0(i,j,k)*rk3coefi

      ! Remove "wrong" diffusive tendencies and replace with wall drag
      ! for u(i,j,k):
      emmo = 0.25_field_r * ( ekm(i,j,k) + ekm(i,j-1,k) + ekm(i-1,j,k) + ekm(i-1,j-1,k) )
      w_at_u_plus = 0.25_field_r * ( w0(i,j,k) + w0(i,j,k+1) + w0(i-1,j,k) + w0(i-1,j,k+1) )  !at u(i,j,k)
      tau_uv_plus = log_wallaw(u0(i,j,k), w_at_u_plus, Cm_ywall)

      up(i,j  ,k) = up(i,j  ,k) + 0.5_field_r * emmo * ( (u0(i,j,k) - u0(i,j-1,k) ) / dy) / dy - 0.5_field_r * tau_uv_plus / dy

      ! for u(i+1,j,k):
      empo = 0.25_field_r * ( ekm(i,j,k) + ekm(i+1,j+1,k) + ekm(i,j,k) + ekm(i,j+1,k) )
      w_at_u_plus = 0.25_field_r * (w0(i+1,j,k) + w0(i+1,j,k+1) + w0(i,j,k) + w0(i,j,k+1) )
      tau_uv_plus = log_wallaw(u0(i+1,j,k), w_at_u_plus, Cm_ywall)

      up(i+1,j,k) = up(i+1,j,k) + 0.5_field_r * empo * ( (u0(i+1,j,k) - u0(i,j,k) ) / dy) / dy - 0.5_field_r * tau_uv_plus / dy

      ! for w(i,j,k):
      if (k /= 1) then ! not correctable when at surface (k = 1)
        emom = ( dzf(k-1) * ( ekm(i,j,k)  + ekm(i,j-1,k)  )  + &
                dzf(k)  * ( ekm(i,j,k-1) + ekm(i,j-1,k-1) ) ) / &
                ( 4.0_field_r * dzh(k) )
        u_at_w_plus = 0.25_field_r * ( u0(i,j,k-1) + u0(i,j,k) + u0(i+1,j,k-1) + u0(i+1,j,k) )
        tau_wv_plus = log_wallaw(w0(i  ,j,k), u_at_w_plus, Cm_ywall)

        wp(i,j,k  ) = wp(i,j,k  ) + 0.5_field_r * emom * ( (w0(i,j,k  ) - w0(i,j-1,  k)) / dy ) / dy - 0.5_field_r * tau_wv_plus / dy
      end if

      ! for w(i,j,k+1):
      emop = ( dzf(k) * ( ekm(i,j,k+1)  + ekm(i,j-1,k+1)  )  + &
              dzf(k+1)  * ( ekm(i,j,k) + ekm(i,j-1,k) ) ) / &
              ( 4.0_field_r * dzh(k+1) )
      u_at_w_plus = 0.25_field_r * ( u0(i,j,k) + u0(i,j,k+1) + u0(i+1,j,k) + u0(i+1,j,k+1) )
      tau_wv_plus = log_wallaw(w0(i  ,j,k+1), u_at_w_plus, Cm_ywall)

      wp(i,j,k+1) = wp(i,j,k+1) + 0.5_field_r * emop * ( (w0(i,j,k+1) - w0(i,j-1,k+1) ) / dy ) / dy - 0.5_field_r * tau_wv_plus / dy

      ! Enforcing zero flux by correction of tendencies of temperature, moisture and other scalars (by negating diffusion term)

      ! scalars at (i,j,k) are "to the right" of the y-walls, so correct on s(i,j,k)
      ! "+" because we reflect it back, and 0.5_field_r comes from interpolation of ekh
      thlp(i,j,k) = thlp(i,j,k) + 0.5_field_r * ( ekh(i,j,k) + ekh(i,j-1,k) ) * ( thl0(i,j,k) - thl0(i,j-1,k) ) * dy2i
      qtp(i,j,k)  =  qtp(i,j,k) + 0.5_field_r * ( ekh(i,j,k) + ekh(i,j-1,k) ) * (  qt0(i,j,k) -  qt0(i,j-1,k) ) * dy2i

      do nc=1,nsv
        svp(i,j,k,nc) = svp(i,j,k,nc) + 0.5_field_r * ( ekh(i,j,k) + ekh(i,j-1,k) ) * ( sv0(i,j,k,nc) - sv0(i,j-1,k,nc) ) * dy2i
      end do
      !call xwalle12(i,j,k) ! correction is ignored assuming u,v,w,subgrid TKE are near zero inside buildings

      ! Finally, set correct heat flux from wall to fluid
      uspeed = 0.5_field_r * ( ( u0(i,j,k) + u0(i+1,j,k) )**2 + ( w0(i,j,k) + w0(i,j,k+1) )**2 )**0.5_field_r
      thlp(i,j,k) = thlp(i,j,k) + Cd_ywall * uspeed * (thlwall - thl0(i,j,k)) / dy
      qtp (i,j,k) = qtp (i,j,k) + Cd_ywall * uspeed * ( qtwall -  qt0(i,j,k)) / dy ! check this one..

    end do

    ! Correct tendencies for walls in negative y-direction
    !$acc parallel loop gang vector default(present)
!!$omp target teams loop defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
    do nn = 1,Nywalls_min
      i = iyw_m(nn,1)
      j = iyw_m(nn,2)
      k = iyw_m(nn,3)

      ! First set the normal velocity at the wall (correct for any drift)
      ! note: these should actually already be corrected for in loop over iobst (due to staggered grid arrangement)
      vp(i,j+1,k)   = -v0(i,j+1,k)*rk3coefi

      ! Remove "wrong" diffusive tendencies and replace with wall drag
      ! for u(i,j,k):
      emmo = 0.25_field_r * ( ekm(i,j+1,k) + ekm(i,j,k) + ekm(i-1,j+1,k) + ekm(i-1,j,k) )
      w_at_u_min  = 0.25_field_r * ( w0(i,j,k) + w0(i,j,k+1) + w0(i-1,j,k) + w0(i-1,j,k+1) )  !at u(i,j,k)
      tau_uv_min = log_wallaw(u0(i,j,k), w_at_u_min , Cm_ywall) !if v0 > 0, tau > 0, minus sign in tendency enforces opposing friction

      up(i,j,k) = up(i,j  ,k) - 0.5_field_r * emmo * ( (u0(i,j+1,k) - u0(i,j,k) ) / dy) / dy  - 0.5_field_r * tau_uv_min / dy ! factor 0.5 originates to avoid double correction

      ! for u(i+1,j,k):
      empo = 0.25_field_r * ( ekm(i,j+1,k) + ekm(i,j,k) + ekm(i+1,j+1,k) + ekm(i+1,j,k) )
      w_at_u_min = 0.25_field_r * ( w0(i+1,j,k) + w0(i+1,j,k+1) + w0(i,j,k) + w0(i,j,k+1) )
      tau_uv_min = log_wallaw(u0(i+1,j,k), w_at_u_min , Cm_ywall)

      up(i+1,j,k) = up(i+1,j,k) - 0.5_field_r * empo * ( (u0(i+1,j+1,k) - u0(i+1,j,k) ) / dy) / dy - 0.5_field_r * tau_uv_min / dy

      ! for w(i,j,k):
      if (k /= 1) then ! not correctable when at surface (k = 1)
        emom = ( dzf(k-1) * ( ekm(i,j+1,k)  + ekm(i,j,k)  )  + &
                dzf(k)  * ( ekm(i,j+1,k-1) + ekm(i,j,k-1) ) ) / &
                ( 4.0_field_r * dzh(k) )
        u_at_w_min = 0.25_field_r * ( u0(i,j,k-1) + u0(i,j,k) + u0(i+1,j,k-1) + u0(i+1,j,k) )
        tau_wv_min = log_wallaw(w0(i,j,k), u_at_w_min , Cm_ywall)

        wp(i,j,k  ) = wp(i,j,k  ) - 0.5_field_r * emom * ( (w0(i,j+1,k  ) - w0(i,j,  k)) / dy ) / dy - 0.5_field_r * tau_wv_min / dy
      end if

      ! for w(i,j,k+1):
      emop = ( dzf(k) * ( ekm(i,j+1,k+1)  + ekm(i,j,k+1)  )  + &
              dzf(k+1)  * ( ekm(i,j+1,k) + ekm(i,j,k) ) ) / &
              ( 4.0_field_r * dzh(k+1) )
      u_at_w_min = 0.25_field_r * ( u0(i,j,k) + u0(i,j,k+1) + u0(i+1,j,k) + u0(i+1,j,k+1) )
      tau_wv_min = log_wallaw(w0(i,j,k+1), u_at_w_min , Cm_ywall)

      wp(i,j,k+1) = wp(i,j,k+1) - 0.5_field_r * emop * ( (w0(i,j+1,k+1) - w0(i,j,k+1) ) / dy ) / dy - 0.5_field_r * tau_wv_min / dy

      ! Enforcing zero flux by correction of tendencies of temperature, moisture and other scalars (by negating diffusion term)

      ! scalars at (i,j,k) are "to the right" of the y-walls, so correct on s(i,j-1,k) which is on the left for minus y-walls
      ! "-" because we reflect it back, and 0.5_field_r comes from interpolation of ekh
      thlp(i,j,k) = thlp(i,j,k) - 0.5_field_r * ( ekh(i,j+1,k) + ekh(i,j,k) ) * ( thl0(i,j+1,k) - thl0(i,j,k) ) * dy2i
      qtp(i,j,k)  =  qtp(i,j,k) - 0.5_field_r * ( ekh(i,j+1,k) + ekh(i,j,k) ) * (  qt0(i,j+1,k) -  qt0(i,j,k) ) * dy2i

      do nc=1,nsv
        svp(i,j,k,nc) = svp(i,j,k,nc) - 0.5_field_r * ( ekh(i,j+1,k) + ekh(i,j,k) ) * ( sv0(i,j+1,k,nc) - sv0(i,j,k,nc) ) * dy2i
      end do

      !> Finally, set correct heat flux from wall to fluid
      uspeed = 0.5_field_r * ( ( u0(i,j,k) + u0(i+1,j,k) )**2 + ( w0(i,j,k) + w0(i,j,k+1) )**2 )**0.5_field_r
      thlp(i,j,k) = thlp(i,j,k) + Cd_ywall * uspeed * (thlwall - thl0(i,j,k)) / dy
      qtp (i,j,k) = qtp (i,j,k) + Cd_ywall * uspeed * ( qtwall -  qt0(i,j,k)) / dy ! check this one..
    end do

    call timer_toc('modibm/applyibm')
    return
  end subroutine applyibm

  !> Force velocity at the boundaries to 0 for a better interaction with the poisson solver
  subroutine zerowallvelocity

    use modfields,      only : up, vp, wp, u0, v0, w0
    use modglobal,      only : rk3step, kmax, i1, j1, k1, ih, jh, rdt
    use modmpi,         only : excjs

    implicit none
    integer  :: i, j, k, nn
    real     :: rk3coef,rk3coefi

    rk3coef = rdt / (4. - dble(rk3step))
    rk3coefi = 1. / rk3coef

    ! Set tendencies inside obstacled (i.e., correct for any drift from previous integration step)
    do nn = 1,Nobst_wide
      i = iobst(nn,1)
      j = iobst(nn,2)
      k = iobst(nn,3)

      ! Correction of velocities also corrects walls w. normals in negative x,y-directions (due to staggered grid arrangement)
      up(i,j,k)   = -u0(i,j,k) * rk3coefi
      vp(i,j,k)   = -v0(i,j,k) * rk3coefi
      wp(i,j,k)   = -w0(i,j,k) * rk3coefi

      ! Do trick: moving one index up by 1 in each direction also corrects walls with normals pointing in positive direction
      up(i+1,j,k) = -u0(i+1,j,k) * rk3coefi
      vp(i,j+1,k) = -v0(i,j+1,k) * rk3coefi
      wp(i,j,k+1) = -w0(i,j,k+1) * rk3coefi

    end do

    return
  end subroutine zerowallvelocity

  !> Calculate drag using logarithmic law-of-wall.
  function log_wallaw(u1,u2,Cm_hor_wall) result(tau)
!!$omp declare target

    !$acc routine seq
    real(field_r), intent(in) :: u1,u2,Cm_hor_wall

    real(field_r)             :: tau

    tau  = Cm_hor_wall * sqrt(u1**2 + u2**2) * u1   ! minus sign included in subroutine above, where it checks direction of the wind
                                                    
  end function log_wallaw

  ! TODO: wall corection for e12 currenly neglected. To be tested in detail later.
  ! This allows for diffusive flux of e12 over the walls, while e12 and ekm/ekh should be insignificant inside buildings. Just outside ekm is probably too large to accurately determine such flux at all -> so we should consider proper cancellation.

  ! subroutine xwalle12(i,j,k)

  !   use modglobal,      only : dx2i, dx, dy, dzh
  !   use modsubgriddata, only : ekm
  !   use modfields,      only : e12p, e120, u0, v0, w0

  !   implicit none

  !   integer, intent(in)    :: i,j,k

  !   if(.not.(k==1)) then
  !     if(.not. (e12p(i,j,k)==0)) then
  !       e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
  !                                 + ekm(i,j,k)/(2*e120(i,j,k))* (&  !source terms
  !                                    -((w0(i,j,k+1)-w0(i-1,j,k+1))  / dx             + &
  !                                      (u0(i,j,k+1)-u0(i,j,k))      / dzh(k+1) )**2  + &
  !                                    +(2.*(w0(i,j,k+1))             / dx             + &
  !                                      (u0(i,j,k+1)-u0(i,j,k))      / dzh(k+1) )**2  + &

  !                                    -((w0(i,j,k)-w0(i-1,j,k))      / dx             + &
  !                                      (u0(i,j,k)-u0(i,j,k-1))      / dzh(k)   )**2  + &
  !                                    +(2.*(w0(i,j,k))               / dx             + &
  !                                      (u0(i,j,k)-u0(i,j,k-1))      / dzh(k)   )**2  + &

  !                                    -((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
  !                                      (v0(i,j+1,k)-v0(i-1,j+1,k))  / dx       )**2  + &
  !                                    +((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
  !                                      (2.*v0(i,j+1,k))             / dx       )**2  + &

  !                                    -((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
  !                                      (v0(i,j,k)-v0(i-1,j,k))      / dx       )**2  + &
  !                                    +((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
  !                                      (2.*v0(i,j,k))               / dx       )**2    &
  !                                   )
  !     elseif(.not. (e12p(i-1,j,k)==0)) then
  !       e12p(i-1,j,k) = e12p(i-1,j,k) + (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
  !                                 + ekm(i-1,j,k)/(2*e120(i-1,j,k))* (&  !source terms
  !                                      -((w0(i,j,k)-w0(i-1,j,k))    / dx             + &
  !                                        (u0(i,j,k)-u0(i,j,k-1))    / dzh(k)   )**2  + &
  !                                      +((-2.*w0(i-1,j,k))          / dx             + &
  !                                        (u0(i,j,k)-u0(i,j,k-1))    / dzh(k)   )**2  + &

  !                                      -((w0(i,j,k+1)-w0(i-1,j,k+1))/ dx             + &
  !                                        (u0(i,j,k+1)-u0(i,j,k))    / dzh(k+1) )**2  + &
  !                                      +((-2.*w0(i-1,j,k+1))        / dx             + &
  !                                        (u0(i,j,k+1)-u0(i,j,k))    / dzh(k+1) )**2  + &

  !                                      -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
  !                                      +((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
  !                                        (-2.*v0(i-1,j,k))          / dx       )**2  + &

  !                                      -((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
  !                                        (v0(i,j+1,k)-v0(i-1,j+1,k))/ dx       )**2  + &
  !                                      +((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
  !                                        (-2.*v0(i-1,j+1,k))        / dx       )**2    &
  !                                   )
  !   end if
  !   else !Special treatment for the lowest full level: k=1
  !     if(.not. (e12p(i,j,k)==0)) then
  !       e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
  !                                 + ekm(i,j,k)/(2*e120(i,j,k))* (&  !source terms
  !                                    -((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
  !                                      (v0(i,j+1,k)-v0(i-1,j+1,k))  / dx       )**2  + &
  !                                    +((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
  !                                      (2.*v0(i,j+1,k))             / dx       )**2  + &

  !                                    -((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
  !                                      (v0(i,j,k)-v0(i-1,j,k))      / dx       )**2  + &
  !                                    +((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
  !                                      (2.*v0(i,j,k))               / dx       )**2    &
  !                                   )
  !     elseif(.not. (e12p(i-1,j,k)==0)) then
  !       e12p(i-1,j,k) = e12p(i-1,j,k) + (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
  !                                 + ekm(i-1,j,k)/(2*e120(i-1,j,k))* (&  !source terms

  !                                      -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
  !                                      +((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
  !                                        (-2.*v0(i-1,j,k))          / dx       )**2  + &

  !                                      -((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
  !                                        (v0(i,j+1,k)-v0(i-1,j+1,k))/ dx       )**2  + &
  !                                      +((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
  !                                        (-2.*v0(i-1,j+1,k))        / dx       )**2    &
  !                                   )
  !     end if
  !   end if
  ! end subroutine xwalle12

  ! subroutine ywalle12(i,j,k)

  !   use modglobal,      only : dy2i, dx, dy, dzh
  !   use modsubgriddata, only : ekm
  !   use modfields,      only : e12p, e120, u0, v0, w0

  !   implicit none

  !   integer, intent(in)    :: i,j,k

  !   if(.not.(k==1)) then
  !     if(.not. (e12p(i,j,k)==0)) then
  !       e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
  !                                 + ekm(i,j,k)/(2.*e120(i,j,k))* (&  !source terms
  !                                      -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
  !                                      +((2.*u0(i,j,k))             / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

  !                                      -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
  !                                        (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
  !                                      +((2.*u0(i+1,j,k))           / dy             + &
  !                                        (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &

  !                                      -((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
  !                                        (w0(i,j,k+1)-w0(i,j-1,k+1))/ dy       )**2  + &
  !                                      +((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
  !                                        (2.*w0(i,j,k+1))           / dy       )**2  + &

  !                                      -((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
  !                                        (w0(i,j,k)-w0(i,j-1,k))    / dy       )**2  + &
  !                                      +((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
  !                                        (2.*w0(i,j,k))             / dy       )**2    &
  !                                   )
  !     elseif(.not. (e12p(i,j-1,k)==0)) then
  !       e12p(i,j-1,k) = e12p(i,j-1,k) + (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
  !                                 + ekm(i,j-1,k)/(2.*e120(i,j-1,k))* (&  !source terms
  !                                      -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
  !                                      +((-2.*u0(i,j-1,k))          / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

  !                                      -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
  !                                        (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
  !                                      +((-2.*u0(i+1,j-1,k))        / dy             + &
  !                                        (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &

  !                                      -((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
  !                                        (w0(i,j,k)-w0(i,j-1,k))    / dy       )**2  + &
  !                                      +((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
  !                                        (-2.*w0(i,j-1,k))          / dy       )**2  + &

  !                                      -((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
  !                                        (w0(i,j,k+1)-w0(i,j-1,k+1))/ dy       )**2  + &
  !                                      +((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
  !                                        (-2.*w0(i,j-1,k+1))        / dy       )**2    &
  !                                   )
  !     end if
  !   else !Special treatment for the lowest full level: k=1
  !     if(.not. (e12p(i,j,k)==0)) then
  !       e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
  !                                 + ekm(i,j,k)/(2.*e120(i,j,k))* (&  !source terms
  !                                      -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
  !                                      +((2.*u0(i,j,k))             / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

  !                                      -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
  !                                        (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
  !                                      +((2.*u0(i+1,j,k))           / dy             + &
  !                                        (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2    &
  !                                   )

  !     elseif(.not. (e12p(i,j-1,k)==0)) then
  !       e12p(i,j-1,k) = e12p(i,j-1,k) + (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
  !                                 + ekm(i,j-1,k)/(2.*e120(i,j-1,k))* (&  !source terms
  !                                      -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
  !                                      +((-2.*u0(i,j-1,k))          / dy             + &
  !                                        (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

  !                                      -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
  !                                        (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
  !                                      +((-2.*u0(i+1,j-1,k))        / dy             + &
  !                                        (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2    &
  !                                   )
  !     end if
  !   end if
  ! end subroutine ywalle12

end module modibm
