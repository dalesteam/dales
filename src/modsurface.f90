!> \file modsurface.f90
!!  Surface parameterization
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

!>
!! Surface routine including a land-surface scheme
!>
!! This module provides an interactive surface parameterization
!!
!! \par Revision list
!! \par Chiel van Heerwaarden
!! \todo documentation
!! \todo implement water reservoir at land surface for dew and interception water
!! \todo add moisture transport between soil layers
!!  \deprecated Modsurface replaces the old modsurf.f90
!
!   TODO: apply heterogeneity on pressure and temperature fields throughout the DALES code (e.g. thermodynamics) (HGO)
!
!
!Able to handle heterogeneous surfaces using the switch lhetero
!In case of heterogeneity an input file is needed
!EXAMPLE of surface.inp.xxx:

!#Surface input file - the standard land cover should be listed below. It is marked by typenr = 0
!#typenr    name       z0mav  z0hav    thls   ps    ustin  wtsurf  wqsurf  wsvsurf(01)  wsvsurf(02)  wsvsurf(03)  wsvsurf(04)
!    0   "standard  "  0.035  0.035   300.0  1.0e5   0.1    0.15   0.1e-3          1.0          0.0          0.0       0.0005
!    1   "forest    "  0.500  0.500   300.0  1.0e5   0.1    0.15   0.2e-3          1.0          0.0          0.0       0.0005
!    2   "grass     "  0.035  0.035   300.0  1.0e5   0.1    0.30   0.1e-3          1.0          0.0          0.0       0.0005

module modsurface
  use modsurfdata
  implicit none
  !public  :: initsurface, surface, exitsurface

save

contains
!> Reads the namelists and initialises the soil.
  subroutine initsurface

    use modglobal,  only : jmax, i1, i2, j1, j2, ih, jh, imax, jtot, cp, rlv, zf, nsv, ifnamopt, fname_options, ifinput, cexpnr
    use modraddata, only : iradiation
    use modfields,  only : thl0, qt0
    use modmpi,     only : myid, nprocs, comm3d, mpierr, my_real, mpi_logical, mpi_integer

    implicit none

    integer   :: i,j,k, landindex, ierr, defined_landtypes, landtype_0 = -1
    character(len=1500) :: readbuffer
    namelist/NAMSURFACE/ & !< Soil related variables
      isurf,tsoilav, tsoildeepav, phiwav, rootfav, &
      ! Land surface related variables
      lmostlocal, lsmoothflux, lneutral, z0mav, z0hav, rsisurf2, Cskinav, lambdaskinav, albedoav, Qnetav, cvegav, &
      ! Jarvis-Steward related variables
      rsminav, rssoilminav, LAIav, gDav, &
      ! Prescribed values for isurf 2, 3, 4
      z0, thls, ps, ustin, wtsurf, wqsurf, wsvsurf, &
      ! Heterogeneous variables
      lhetero, xpatches, ypatches, land_use

    ! 1    -   Initialize soil

    !if (isurf == 1) then

    ! 1.0  -   Read LSM-specific namelist

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMSURFACE,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMSURFACE'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMSURFACE'
      endif
      write(6 ,NAMSURFACE)
      close(ifnamopt)
    end if

    call MPI_BCAST(isurf        , 1       , MPI_INTEGER, 0, comm3d, mpierr)
    call MPI_BCAST(tsoilav      , ksoilmax, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(tsoildeepav  , 1       , MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(phiwav       , ksoilmax, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(rootfav      , ksoilmax, MY_REAL, 0, comm3d, mpierr)

    call MPI_BCAST(lmostlocal   , 1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(lsmoothflux  , 1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(lneutral     , 1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(z0mav        , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(z0hav        , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(rsisurf2     , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(Cskinav      , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(lambdaskinav , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(albedoav     , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(Qnetav       , 1, MY_REAL, 0, comm3d, mpierr)

    call MPI_BCAST(rsminav      , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(rssoilminav  , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(cvegav       , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(LAIav        , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(gDav         , 1, MY_REAL, 0, comm3d, mpierr)

    call MPI_BCAST(z0         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ustin      ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wtsurf     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wqsurf     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wsvsurf(1:nsv),nsv,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ps         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(thls       ,1,MY_REAL   ,0,comm3d,mpierr)

    call MPI_BCAST(lhetero                    ,            1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(xpatches                   ,            1, MPI_INTEGER, 0, comm3d, mpierr)
    call MPI_BCAST(ypatches                   ,            1, MPI_INTEGER, 0, comm3d, mpierr)
    call MPI_BCAST(land_use(1:mpatch,1:mpatch),mpatch*mpatch, MPI_INTEGER, 0, comm3d, mpierr)

    if(lhetero) then
      if(xpatches .gt. mpatch) then
        stop "NAMSURFACE: more xpatches defined than possible (change mpatch in modsurfdata to a higher value)"
      endif
      if(ypatches .gt. mpatch) then
        stop "NAMSURFACE: more ypatches defined than possible (change mpatch in modsurfdata to a higher value)"
      endif
      if (lsmoothflux .eqv. .true.) write(6,*) 'WARNING: You selected to use uniform heat fluxes (lsmoothflux) and heterogeneous surface conditions (lhetero) at the same time' 
      if (mod(imax,xpatches) .ne. 0) stop "NAMSURFACE: Not an integer amount of grid points per patch in the x-direction"
      if (mod(jtot,ypatches) .ne. 0) stop "NAMSURFACE: Not an integer amount of grid points per patch in the y-direction"

      allocate(z0mav_patch(xpatches,ypatches))
      allocate(z0hav_patch(xpatches,ypatches))
      allocate(thls_patch(xpatches,ypatches))
      allocate(qts_patch(xpatches,ypatches))
      allocate(thvs_patch(xpatches,ypatches))
      allocate(ps_patch(xpatches,ypatches))
      allocate(ustin_patch(xpatches,ypatches))
      allocate(wt_patch(xpatches,ypatches))
      allocate(wq_patch(xpatches,ypatches))
      allocate(wsv_patch(100,xpatches,ypatches))

      z0mav_patch = -1
      z0hav_patch = -1
      thls_patch  = -1
      qts_patch   = -1
      thvs_patch  = -1
      ps_patch    = -1
      ustin_patch = -1
      wt_patch    = -1
      wq_patch    = -1
      wsv_patch   = -1

      defined_landtypes = 0
      open (ifinput,file='surface.inp.'//cexpnr)
      ierr = 0
      do while (ierr == 0)  
        read(ifinput, '(A)', iostat=ierr) readbuffer
        if (ierr == 0) then                               !So no end of file is encountered
          if (readbuffer(1:1)=='#') then
            if (myid == 0)   print *,trim(readbuffer)
          else
            if (myid == 0)   print *,trim(readbuffer)
            defined_landtypes = defined_landtypes + 1
            i = defined_landtypes
            read(readbuffer, *, iostat=ierr) landtype(i), landname(i), z0mav_land(i), z0hav_land(i), thls_land(i), &
              ps_land(i), ustin_land(i), wt_land(i), wq_land(i), wsv_land(1:nsv,i) 

            if (ustin_land(i) .lt. 0) then
              if (myid == 0) stop "NAMSURFACE: A ustin value in the surface input file is negative"
            endif
            if(isurf .ne. 3) then
              if(z0mav_land(i) .lt. 0) then
                if (myid == 0) stop "NAMSURFACE: a z0mav value is not set or negative in the surface input file"
              end if
              if(z0hav_land(i) .lt. 0) then
                if (myid == 0) stop "NAMSURFACE: a z0hav value is not set or negative in the surface input file"
              end if
            end if

            if (landtype(i) .eq. 0) landtype_0 = i
            do j = 1, (i-1)
              if (landtype(i) .eq. landtype(j)) stop "NAMSURFACE: Two land types have the same type number"
            enddo
                
          endif
        endif
      enddo
      close(ifinput)

      if (myid == 0) then
        if (landtype_0 .eq. -1) then
          stop "NAMSURFACE: no standard land type (0) is defined"
        else
          print "(a,i2,a,i2)","There are ",defined_landtypes," land types defined in the surface input file. The standard land type is defined by line ",landtype_0
        endif
      endif

      thls   = 0
      ps     = 0
      ustin  = 0
      wtsurf = 0
      wqsurf = 0
      wsvsurf(1:nsv) = 0

      do i = 1, xpatches
        do j = 1, ypatches
          landindex = landtype_0
          do k = 1, defined_landtypes
            if (landtype(k) .eq. land_use(i,j)) then
              landindex = k
            endif
          enddo
        
          z0mav_patch(i,j)     = z0mav_land(landindex)
          z0hav_patch(i,j)     = z0hav_land(landindex)
          thls_patch(i,j)      = thls_land(landindex)
          ps_patch(i,j)        = ps_land(landindex)
          ustin_patch(i,j)     = ustin_land(landindex)
          wt_patch(i,j)        = wt_land(landindex)
          wq_patch(i,j)        = wq_land(landindex)
          wsv_patch(1:nsv,i,j) = wsv_land(1:nsv,landindex)

          thls   = thls   + ( thls_patch(i,j)  / ( xpatches * ypatches ) )
          ps     = ps     + ( ps_patch(i,j)    / ( xpatches * ypatches ) )
          ustin  = ustin  + ( ustin_patch(i,j) / ( xpatches * ypatches ) )
          wtsurf = wtsurf + ( wt_patch(i,j)    / ( xpatches * ypatches ) )
          wqsurf = wqsurf + ( wq_patch(i,j)    / ( xpatches * ypatches ) )
          wsvsurf(1:nsv) = wsvsurf(1:nsv) + ( wsv_patch(1:nsv,i,j) / ( xpatches * ypatches ) )
        enddo
      enddo
    else
      if((z0mav == -1 .and. z0hav == -1) .and. (z0 .ne. -1)) then
        z0mav = z0
        z0hav = z0
        write(6,*) "WARNING: z0m and z0h not defined, set equal to z0"
      end if

      if(isurf .ne. 3) then
        if(z0mav == -1) then
          stop "NAMSURFACE: z0mav is not set"
        end if
        if(z0hav == -1) then
          stop "NAMSURFACE: z0hav is not set"
        end if
      end if
  
    endif


    if(isurf == 1) then
      if(tsoilav(1) == -1 .or. tsoilav(2) == -1 .or. tsoilav(3) == -1 .or. tsoilav(4) == -1) then
        stop "NAMSURFACE: tsoil is not set"
      end if
      if(tsoildeepav == -1) then
        stop "NAMSURFACE: tsoildeep is not set"
      end if
      if(phiwav(1) == -1 .or. phiwav(2) == -1 .or. phiwav(3) == -1 .or. phiwav(4) == -1) then
        stop "NAMSURFACE: phiw is not set"
      end if
      if(rootfav(1) == -1 .or. rootfav(2) == -1 .or. rootfav(3) == -1 .or. rootfav(4) == -1) then
        stop "NAMSURFACE: rootf is not set"
      end if
      if(Cskinav == -1) then
        stop "NAMSURFACE: Cskinav is not set"
      end if
      if(lambdaskinav == -1) then
        stop "NAMSURFACE: lambdaskinav is not set"
      end if
      if(albedoav == -1) then
        stop "NAMSURFACE: albedoav is not set"
      end if
      if(Qnetav == -1) then
        stop "NAMSURFACE: Qnetav is not set"
      end if
      if(rsminav == -1) then
        stop "NAMSURFACE: rsminav is not set"
      end if
      if(LAIav == -1) then
        stop "NAMSURFACE: LAIav is not set"
      end if
      if(gDav == -1) then
        stop "NAMSURFACE: gDav is not set"
      end if

    end if

    ! 1.1  -   Allocate arrays
    if(isurf == 1) then
      allocate(zsoil(ksoilmax))
      allocate(dzsoil(ksoilmax))
      allocate(dzsoilh(ksoilmax))

      allocate(lambda(i2,j2,ksoilmax))
      allocate(lambdah(i2,j2,ksoilmax))
      allocate(Dh(i2,j2,ksoilmax))
      allocate(phiw(i2,j2,ksoilmax))
      allocate(pCs(i2,j2,ksoilmax))
      allocate(rootf(i2,j2,ksoilmax))
      allocate(tsoil(i2,j2,ksoilmax))
      allocate(tsoildeep(i2,j2))
      allocate(phitot(i2,j2))

      ! 1.2   -  Initialize arrays
      ! First test, pick ECMWF config
      dzsoil(1) = 0.07
      dzsoil(2) = 0.21
      dzsoil(3) = 0.72
      dzsoil(4) = 1.89

      !! 1.3   -  Calculate vertical layer properties
      zsoil(1)  = dzsoil(1)
      do k = 2, ksoilmax
        zsoil(k) = zsoil(k-1) + dzsoil(k)
      end do

      do k = 1, ksoilmax-1
        dzsoilh(k) = 0.5 * (dzsoil(k+1) + dzsoil(k))
      end do
      dzsoilh(ksoilmax) = 0.5 * dzsoil(ksoilmax)

      ! 1.4   -   Set evaporation related properties
      ! Set water content of soil - constant in this scheme
      phiw(:,:,1) = phiwav(1)
      phiw(:,:,2) = phiwav(2)
      phiw(:,:,3) = phiwav(3)
      phiw(:,:,4) = phiwav(4)

      phitot = 0.0

      do k = 1, ksoilmax
        phitot(:,:) = phitot(:,:) + phiw(:,:,k) * dzsoil(k)
      end do

      phitot(:,:) = phitot(:,:) / zsoil(ksoilmax)

      ! Set root fraction per layer for short grass
      rootf(:,:,1) = rootfav(1)
      rootf(:,:,2) = rootfav(2)
      rootf(:,:,3) = rootfav(3)
      rootf(:,:,4) = rootfav(4)

      ! Calculate conductivity saturated soil
      lambdasat = lambdasm ** (1. - phi) * lambdaw ** (phi)

      tsoil(:,:,1)   = tsoilav(1)
      tsoil(:,:,2)   = tsoilav(2)
      tsoil(:,:,3)   = tsoilav(3)
      tsoil(:,:,4)   = tsoilav(4)
      tsoildeep(:,:) = tsoildeepav

      ! Calculate the soil heat capacity and conductivity based on water content
      ! CvH - If we need prognostic soil moisture at one point, these 2 loops should move to the surface function
      do k = 1, ksoilmax
        do j = 2, j1
          do i = 2, i1
            pCs(i,j,k)    = (1. - phi) * pCm + phiw(i,j,k) * pCw
            Ke            = log10(phiw(i,j,k) / phi) + 1.
            lambda(i,j,k) = Ke * (lambdasat - lambdadry) + lambdadry
          end do
        end do
      end do

      do k = 1, ksoilmax-1
        do j = 2, j1
          do i = 2, i1
            lambdah(i,j,k) = (lambda(i,j,k) * dzsoil(k+1) + lambda(i,j,k+1) * dzsoil(k)) / dzsoilh(k)
          end do
        end do
      end do

      do j = 2, j1
        do i = 2, i1
          lambdah(i,j,ksoilmax) = lambda(i,j,ksoilmax)
        end do
      end do

    end if

    ! 2    -   Initialize land surface
    ! 2.1  -   Allocate arrays
    if(isurf == 1) then
      allocate(Qnet(i2,j2))
      allocate(LE(i2,j2))
      allocate(H(i2,j2))
      allocate(G0(i2,j2))

      Qnet = Qnetav
    end if

    if(isurf <= 2) then
      allocate(rs(i2,j2))
      allocate(rsveg(i2,j2))
      allocate(rsmin(i2,j2))
      allocate(rssoil(i2,j2))
      allocate(rssoilmin(i2,j2))
      allocate(cveg(i2,j2))
      allocate(cliq(i2,j2))
      allocate(ra(i2,j2))
      allocate(tendskin(i2,j2))
      allocate(tskinm(i2,j2))
      allocate(Cskin(i2,j2))
      allocate(lambdaskin(i2,j2))
      allocate(LAI(i2,j2))
      allocate(gD(i2,j2))

      Cskin      = Cskinav
      lambdaskin = lambdaskinav
      rsmin      = rsminav
      rssoilmin  = rsminav
      LAI        = LAIav
      gD         = gDav

      cveg       = cvegav
      cliq       = 0.
    end if

    allocate(albedo(i2,j2))
    allocate(z0m(i2,j2))
    allocate(z0h(i2,j2))
    allocate(obl(i2,j2))
    allocate(tskin(i2,j2))
    allocate(qskin(i2,j2))
    allocate(Cm(i2,j2))
    allocate(Cs(i2,j2))

    if(iradiation == 1) then
      allocate(swdavn(i2,j2,nradtime))
      allocate(swuavn(i2,j2,nradtime))
      allocate(lwdavn(i2,j2,nradtime))
      allocate(lwuavn(i2,j2,nradtime))
    end if

    albedo     = albedoav


    if(lhetero) then
      do j=1,j2
        do i=1,i2
          z0m(i,j)   = z0mav_patch(patchxnr(i),patchynr(j))
          z0h(i,j)   = z0hav_patch(patchxnr(i),patchynr(j))
        enddo
      enddo
    else
      z0m        = z0mav
      z0h        = z0hav
    endif

    ! 3. Initialize surface layer
    allocate(ustar   (i2,j2))

    allocate(dudz    (i2,j2))
    allocate(dvdz    (i2,j2))

    allocate(thlflux (i2,j2))
    allocate(qtflux  (i2,j2))
    allocate(dqtdz   (i2,j2))
    allocate(dthldz  (i2,j2))
    allocate(svflux  (i2,j2,nsv))
    allocate(svs(nsv))

    ! CvH set initial values for rs and ra to be able to compute qskin
    if(isurf <= 2) then
      ra = 50.
      if(isurf == 1) then
        rs = 100.
      end if
    end if

    return
  end subroutine initsurface

!> Calculates the interaction with the soil, the surface temperature and humidity, and finally the surface fluxes.
  subroutine surface
    use modglobal,  only : rdt, i1, i2, j1, j2, ih, jh, cp, rlv, fkar, zf, cu, cv, nsv, rk3step, timee, rslabs, pi, pref0, rd, rv, eps1, boltz
    use modraddata, only : iradiation, swu, swd, lwu, lwd, useMcICA
    use modfields,  only : thl0, qt0, u0, v0, rhof, ql0, exnf, presf, u0av, v0av
    use modmpi,     only : my_real, mpierr, comm3d, mpi_sum, myid, excj, excjs, mpi_integer
    use moduser,    only : surf_user
    implicit none

    real     :: f1, f2, f3, f4 ! Correction functions for Jarvis-Stewart
    integer  :: i, j, k, n, patchx, patchy
    real     :: upcu, vpcv, horv, horvav, horvpatch(xpatches,ypatches)
    real     :: upatch(xpatches,ypatches), vpatch(xpatches,ypatches)
    real     :: Supatch(xpatches,ypatches), Svpatch(xpatches,ypatches)
    integer  :: Npatch(xpatches,ypatches), SNpatch(xpatches,ypatches)
    real     :: lthls_patch(xpatches,ypatches)
    real     :: lqts_patch(xpatches,ypatches), qts_patch(xpatches,ypatches)
    real     :: phimzf, phihzf
    real     :: rk3coef, thlsl, qtsl

    real     :: ust
    real     :: ustl, wtsurfl, wqsurfl

    real     :: swdav, swuav, lwdav, lwuav
    real     :: exner, exnera, tsurfm, Tatm, e, esat, qsat, desatdT, dqsatdT, Acoef, Bcoef
    real     :: fH, fLE, fLEveg, fLEsoil, fLEpot

    if (isurf==10) then
      call surf_user
      return
    end if

    if(lhetero) then 
      upatch = 0
      vpatch = 0
      Npatch = 0

      do j = 2, j1
        do i = 2, i1
          patchx = patchxnr(i)
          patchy = patchynr(j) 
          upatch(patchx,patchy) = upatch(patchx,patchy) + 0.5 * (u0(i,j,1) + u0(i+1,j,1))
          vpatch(patchx,patchy) = vpatch(patchx,patchy) + 0.5 * (v0(i,j,1) + v0(i,j+1,1))
          Npatch(patchx,patchy) = Npatch(patchx,patchy) + 1
        enddo
      enddo

      call MPI_ALLREDUCE(upatch(1:xpatches,1:ypatches),Supatch(1:xpatches,1:ypatches),xpatches*ypatches,    MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(vpatch(1:xpatches,1:ypatches),Svpatch(1:xpatches,1:ypatches),xpatches*ypatches,    MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(Npatch(1:xpatches,1:ypatches),SNpatch(1:xpatches,1:ypatches),xpatches*ypatches,MPI_INTEGER,MPI_SUM, comm3d,mpierr)
          
      horvpatch = sqrt(((Supatch/SNpatch) + cu) **2. + ((Svpatch/SNpatch) + cv) ** 2.)
      horvpatch = max(horvpatch, 1.e-2)
    endif 

    ! 1     -   Calculate the surface temperature
    
    ! CvH start with computation of drag coefficients to allow for implicit solver
    if(isurf <= 2) then

      if(lneutral) then
        obl(:,:) = -1.e10
        oblav    = -1.e10
      else
        call getobl
      end if

      call MPI_BCAST(oblav ,1,MY_REAL ,0,comm3d,mpierr)

      do j = 2, j1
        do i = 2, i1
          if(isurf == 2) then
            rs(i,j) = rsisurf2
          else
            if(lhetero) then
              patchx = patchxnr(i)
              patchy = patchynr(j) 
            endif
            ! 2.1   -   Calculate the surface resistance 
            ! Stomatal opening as a function of incoming short wave radiation
            if (iradiation > 0) then
              f1  = 1. / min(1., (0.004 * max(0.,-swd(i,j,1)) + 0.05) / (0.81 * (0.004 * max(0.,-swd(i,j,1)) + 1.)))
            else
              f1  = 1.
            end if

            ! Soil moisture availability
            f2  = (phifc - phiwp) / (phitot(i,j) - phiwp)

            ! Response of stomata to vapor deficit of atmosphere
            esat = 0.611e3 * exp(17.2694 * (thl0(i,j,1) - 273.16) / (thl0(i,j,1) - 35.86))

            if(lhetero) then
              e    = qt0(i,j,1) * ps_patch(patchx,patchy) / 0.622
            else
              e    = qt0(i,j,1) * ps / 0.622
            endif

            f3   = 1. / exp(-gD(i,j) * (esat - e) / 100.)

            ! Response to temperature
            exnera  = (presf(1) / pref0) ** (rd/cp)
            Tatm    = exnera * thl0(i,j,1) + (rlv / cp) * ql0(i,j,1)
            f4      = 1./ (1. - 0.0016 * (298.0 - Tatm) ** 2.)

            rsveg(i,j)  = rsmin(i,j) / LAI(i,j) * f1 * f2 * f3 * f4

            ! 2.2   - Calculate soil resistance based on ECMWF method

            f2  = (phifc - phiwp) / (phiw(i,j,1) - phiwp)
            rssoil(i,j) = rssoilmin(i,j) * f2
            rssoil(i,j) = rssoilmin(i,j) * f2
          end if

          ! 3     -   Calculate the drag coefficient and aerodynamic resistance
          Cm(i,j) = fkar ** 2. / (log(zf(1) / z0m(i,j)) - psim(zf(1) / obl(i,j)) + psim(z0m(i,j) / obl(i,j))) ** 2.
          Cs(i,j) = fkar ** 2. / (log(zf(1) / z0m(i,j)) - psim(zf(1) / obl(i,j)) + psim(z0m(i,j) / obl(i,j))) / (log(zf(1) / z0h(i,j)) - psih(zf(1) / obl(i,j)) + psih(z0h(i,j) / obl(i,j)))

          if(lmostlocal) then
            upcu  = 0.5 * (u0(i,j,1) + u0(i+1,j,1)) + cu
            vpcv  = 0.5 * (v0(i,j,1) + v0(i,j+1,1)) + cv
            horv  = sqrt(upcu ** 2. + vpcv ** 2.)
            horv  = max(horv, 1.e-2)
            ra(i,j) = 1. / ( Cs(i,j) * horv )
          else
            if (lhetero) then
              ra(i,j) = 1. / ( Cs(i,j) * horvpatch(patchx,patchy) )
            else
              horvav  = sqrt(u0av(1) ** 2. + v0av(1) ** 2.)
              horvav  = max(horvav, 1.e-2)
              ra(i,j) = 1. / ( Cs(i,j) * horvav )
            endif
          end if

        end do
      end do
    end if

    ! Solve the surface energy balance and the heat and moisture transport in the soil
    if(isurf == 1) then
      thlsl = 0.0
      if(lhetero) then
        lthls_patch = 0.0
        Npatch      = 0
      endif
      do j = 2, j1
        do i = 2, i1
          if(lhetero) then
            patchx = patchxnr(i)
            patchy = patchynr(j) 
          endif
          ! 1.1   -   Calculate the heat transport properties of the soil
          ! CvH I put it in the init function, as we don't have prognostic soil moisture at this stage

          ! 1.2   -   Calculate the skin temperature as the top boundary conditions for heat transport
          if(iradiation > 0) then
            if(iradiation == 1 .and. useMcICA) then
              if(rk3step == 1) then
                swdavn(i,j,2:nradtime) = swdavn(i,j,1:nradtime-1)  
                swuavn(i,j,2:nradtime) = swuavn(i,j,1:nradtime-1)  
                lwdavn(i,j,2:nradtime) = lwdavn(i,j,1:nradtime-1)  
                lwuavn(i,j,2:nradtime) = lwuavn(i,j,1:nradtime-1)  

                swdavn(i,j,1) = swd(i,j,1) 
                swuavn(i,j,1) = swu(i,j,1) 
                lwdavn(i,j,1) = lwd(i,j,1) 
                lwuavn(i,j,1) = lwu(i,j,1) 

              end if
              
              swdav = sum(swdavn(i,j,:)) / nradtime
              swuav = sum(swuavn(i,j,:)) / nradtime
              lwdav = sum(lwdavn(i,j,:)) / nradtime
              lwuav = sum(lwuavn(i,j,:)) / nradtime

              Qnet(i,j) = -(swdav + swuav + lwdav + lwuav)
            else
              Qnet(i,j) = -(swd(i,j,1) + swu(i,j,1) + lwd(i,j,1) + lwu(i,j,1))
            end if
          end if

          ! CvH solve the surface temperature implicitly including variations in LWout
          if(rk3step == 1 .and. timee > 0.) then
            tskinm(i,j) = tskin(i,j)
          end if

          if(lhetero) then
            exner   = (ps_patch(patchx,patchy) / pref0) ** (rd/cp)
          else
            exner   = (ps / pref0) ** (rd/cp)
          endif
          tsurfm  = tskinm(i,j) * exner
          
          esat    = 0.611e3 * exp(17.2694 * (tsurfm - 273.16) / (tsurfm - 35.86))
          if(lhetero) then
            qsat    = 0.622 * esat / ps_patch(patchx,patchy)
          else
            qsat    = 0.622 * esat / ps 
          endif
          desatdT = esat * (17.2694 / (tsurfm - 35.86) - 17.2694 * (tsurfm - 273.16) / (tsurfm - 35.86)**2.)
          if(lhetero) then
            dqsatdT = 0.622 * desatdT / ps_patch(patchx,patchy)
          else
            dqsatdT = 0.622 * desatdT / ps 
          endif

          ! First, remove LWup from Qnet calculation
          Qnet(i,j) = Qnet(i,j) + boltz * tsurfm ** 4.

          ! Calculate coefficients for surface fluxes
          fH      = rhof(1) * cp / ra(i,j)

          ! Allow for dew fall
          if(qsat - qt0(i,j,1) < 0.) then
            rsveg(i,j)  = 0.
            rssoil(i,j) = 0.
          end if

          fLEveg  = (1. - cliq(i,j)) * cveg(i,j) * rhof(1) * rlv / (ra(i,j) + rsveg(i,j))
          fLEsoil = (1. - cveg(i,j))             * rhof(1) * rlv / (ra(i,j) + rssoil(i,j))
          fLEpot  = cliq(i,j) * cveg(i,j)        * rhof(1) * rlv /  ra(i,j)
          
          fLE     = fLEveg + fLEsoil + fLEpot

          exnera  = (presf(1) / pref0) ** (rd/cp)
          Tatm    = exnera * thl0(i,j,1) + (rlv / cp) * ql0(i,j,1)
          
          ! Calculate coefficients for surface fluxes
          fH      = rhof(1) * cp / ra(i,j)

          ! Allow for dew fall
          if(qsat - qt0(i,j,1) < 0.) then
            rsveg(i,j)  = 0.
            rssoil(i,j) = 0.
          end if

          fLEveg  = (1. - cliq(i,j)) * cveg(i,j) * rhof(1) * rlv / (ra(i,j) + rsveg(i,j))
          fLEsoil = (1. - cveg(i,j))             * rhof(1) * rlv / (ra(i,j) + rssoil(i,j))
          fLEpot  = cliq(i,j) * cveg(i,j)        * rhof(1) * rlv /  ra(i,j)
          
          fLE     = fLEveg + fLEsoil + fLEpot

          exnera  = (presf(1) / pref0) ** (rd/cp)
          Tatm    = exnera * thl0(i,j,1) + (rlv / cp) * ql0(i,j,1)
          
          rk3coef = rdt / (4. - dble(rk3step))
          
          Acoef   = Qnet(i,j) - boltz * tsurfm ** 4. + 4. * boltz * tsurfm ** 4. / rk3coef + fH * Tatm + fLE * (dqsatdT * tsurfm - qsat + qt0(i,j,1)) + lambdaskin(i,j) * tsoil(i,j,1)
          Bcoef   = 4. * boltz * tsurfm ** 3. / rk3coef + fH + fLE * dqsatdT + lambdaskin(i,j)
          
          Acoef   = Qnet(i,j) - boltz * tsurfm ** 4. + 4. * boltz * tsurfm ** 4. / rk3coef + fH * Tatm + fLE * (dqsatdT * tsurfm - qsat + qt0(i,j,1)) + lambdaskin(i,j) * tsoil(i,j,1)
          Bcoef   = 4. * boltz * tsurfm ** 3. / rk3coef + fH + fLE * dqsatdT + lambdaskin(i,j)

          if (Cskin(i,j) == 0.) then
            tskin(i,j) = Acoef * Bcoef ** (-1.) / exner
          else
            tskin(i,j) = (1. + rk3coef / Cskin(i,j) * Bcoef) ** (-1.) * (tsurfm + rk3coef / Cskin(i,j) * Acoef) / exner
          end if

          Qnet(i,j)     = Qnet(i,j) - (boltz * tsurfm ** 4. + 4. * boltz * tsurfm ** 3. * (tskin(i,j) * exner - tsurfm) / rk3coef)
          G0(i,j)       = lambdaskin(i,j) * ( tskin(i,j) * exner - tsoil(i,j,1) )
          LE(i,j)       = - fLE * ( qt0(i,j,1) - (dqsatdT * (tskin(i,j) * exner - tsurfm) + qsat))
          if(LE(i,j) == 0.) then
            rs(i,j)     = 1.e8
          else
            rs(i,j)     = - rhof(1) * rlv * (qt0(i,j,1) - (dqsatdT * (tskin(i,j) * exner - tsurfm) + qsat)) / LE(i,j) - ra(i,j)
          end if

          H(i,j)        = - fH  * ( Tatm - tskin(i,j) * exner ) 
          tendskin(i,j) = Cskin(i,j) * (tskin(i,j) - tskinm(i,j)) * exner / rk3coef

          ! 1.4   -   Solve the diffusion equation for the heat transport
          tsoil(i,j,1) = tsoil(i,j,1) + rdt / pCs(i,j,1) * ( lambdah(i,j,ksoilmax) * (tsoil(i,j,2) - tsoil(i,j,1)) / dzsoilh(1) + G0(i,j) ) / dzsoil(1)
          do k = 2, ksoilmax-1
            tsoil(i,j,k) = tsoil(i,j,k) + rdt / pCs(i,j,k) * ( lambdah(i,j,k) * (tsoil(i,j,k+1) - tsoil(i,j,k)) / dzsoilh(k) - lambdah(i,j,k-1) * (tsoil(i,j,k) - tsoil(i,j,k-1)) / dzsoilh(k-1) ) / dzsoil(k)
          end do
          tsoil(i,j,ksoilmax) = tsoil(i,j,ksoilmax) + rdt / pCs(i,j,ksoilmax) * ( lambda(i,j,ksoilmax) * (tsoildeep(i,j) - tsoil(i,j,ksoilmax)) / dzsoil(ksoilmax) - lambdah(i,j,ksoilmax-1) * (tsoil(i,j,ksoilmax) - tsoil(i,j,ksoilmax-1)) / dzsoil(ksoilmax-1) ) / dzsoil(ksoilmax)

          thlsl  = thlsl + tskin(i,j)
          if (lhetero) then
            lthls_patch(patchx,patchy) = lthls_patch(patchx,patchy) + tskin(i,j)
            Npatch(patchx,patchy)      = Npatch(patchx,patchy)      + 1
          endif
        end do
      end do

      call MPI_ALLREDUCE(thlsl      , thls      , 1                ,     MY_REAL, MPI_SUM, comm3d,mpierr)
      thls = thls / rslabs

      if (lhetero) then
        call MPI_ALLREDUCE(lthls_patch(1:xpatches,1:ypatches), thls_patch(1:xpatches,1:ypatches), xpatches*ypatches,     MY_REAL, MPI_SUM, comm3d,mpierr)
        call MPI_ALLREDUCE(Npatch(1:xpatches,1:ypatches)     , SNpatch(1:xpatches,1:ypatches)   , xpatches*ypatches, MPI_INTEGER ,MPI_SUM, comm3d,mpierr)
        thls_patch = thls_patch / SNpatch
      endif

      call qtsurf

    elseif(isurf == 2) then
      do j = 2, j1
        do i = 2, i1
          if(lhetero) then
            tskin(i,j) = thls_patch(patchxnr(i),patchynr(j))
          else
            tskin(i,j) = thls
          endif
        end do
      end do

      call qtsurf

    end if

    ! 2     -   Calculate the surface fluxes
    if(isurf <= 2) then
      do j = 2, j1
        do i = 2, i1
          upcu   = 0.5 * (u0(i,j,1) + u0(i+1,j,1)) + cu
          vpcv   = 0.5 * (v0(i,j,1) + v0(i,j+1,1)) + cv
          horv   = sqrt(upcu ** 2. + vpcv ** 2.)
          horv   = max(horv, 1.e-2)
          horvav = sqrt(u0av(1) ** 2. + v0av(1) ** 2.)
          horvav = max(horvav, 1.e-2)
            
          if(lhetero) then
            patchx = patchxnr(i)
            patchy = patchynr(j) 
          endif
          
          if(lmostlocal) then
            ustar  (i,j) = sqrt(Cm(i,j)) * horv 
          else
            if(lhetero) then
              ustar  (i,j) = sqrt(Cm(i,j)) * horvpatch(patchx,patchy) 
            else
              ustar  (i,j) = sqrt(Cm(i,j)) * horvav
            endif
          end if
          
          thlflux(i,j) = - ( thl0(i,j,1) - tskin(i,j) ) / ra(i,j) 
          qtflux(i,j)  = - (qt0(i,j,1)  - qskin(i,j)) / ra(i,j)
          
          if(lhetero) then
            do n=1,nsv
              svflux(i,j,n) = wsv_patch(n,patchx,patchy) 
            enddo
          else
            do n=1,nsv
              svflux(i,j,n) = wsvsurf(n) 
            enddo
          endif

          ! commented lines are classical Businger-Dyer functions
          if (obl(i,j) < 0.) then
            !phimzf = (1.-16.*zf(1)/obl)**(-0.25)
            phimzf = (1. + 3.6 * (-zf(1)/obl(i,j))**(2./3.))**(-0.5)
            !phihzf = (1.-16.*zf(1)/obl)**(-0.50)
            phihzf = (1. + 7.9 * (-zf(1)/obl(i,j))**(2./3.))**(-0.5)
          elseif (obl(i,j) > 0.) then
            phimzf = (1.+5.*zf(1)/obl(i,j))
            phihzf = (1.+5.*zf(1)/obl(i,j))
          else
            phimzf = 1.
            phihzf = 1.
          endif

          dudz  (i,j) = ustar(i,j) * phimzf / (fkar*zf(1))*(upcu/horv)
          dvdz  (i,j) = ustar(i,j) * phimzf / (fkar*zf(1))*(vpcv/horv)
          dthldz(i,j) = - thlflux(i,j) / ustar(i,j) * phihzf / (fkar*zf(1))
          dqtdz (i,j) = - qtflux(i,j)  / ustar(i,j) * phihzf / (fkar*zf(1))

        end do
      end do

      if(lsmoothflux) then

        ustl    = sum(ustar  (2:i1,2:j1))
        wtsurfl = sum(thlflux(2:i1,2:j1))
        wqsurfl = sum(qtflux (2:i1,2:j1))

        call MPI_ALLREDUCE(ustl, ust, 1,  MY_REAL,MPI_SUM, comm3d,mpierr)
        call MPI_ALLREDUCE(wtsurfl, wtsurf, 1,  MY_REAL,MPI_SUM, comm3d,mpierr)
        call MPI_ALLREDUCE(wqsurfl, wqsurf, 1,  MY_REAL,MPI_SUM, comm3d,mpierr)

        wtsurf = wtsurf / rslabs
        wqsurf = wqsurf / rslabs

        do j = 2, j1
          do i = 2, i1

            thlflux(i,j) = wtsurf 
            qtflux (i,j) = wqsurf 

            do n=1,nsv
              svflux(i,j,n) = wsvsurf(n)
            enddo

            if (obl(i,j) < 0.) then
              !phimzf = (1.-16.*zf(1)/obl)**(-0.25)
              phimzf = (1. + 3.6 * (-zf(1)/obl(i,j))**(2./3.))**(-0.5)
              !phihzf = (1.-16.*zf(1)/obl)**(-0.50)
              phihzf = (1. + 7.9 * (-zf(1)/obl(i,j))**(2./3.))**(-0.5)
            elseif (obl(i,j) > 0.) then
              phimzf = (1.+5.*zf(1)/obl(i,j))
              phihzf = (1.+5.*zf(1)/obl(i,j))
            else
              phimzf = 1.
              phihzf = 1.
            endif

            upcu  = 0.5 * (u0(i,j,1) + u0(i+1,j,1)) + cu
            vpcv  = 0.5 * (v0(i,j,1) + v0(i,j+1,1)) + cv
            horv  = sqrt(upcu ** 2. + vpcv ** 2.)
            horv  = max(horv, 1.e-2)

            dudz  (i,j) = ustar(i,j) * phimzf / (fkar*zf(1))*(upcu/horv)
            dvdz  (i,j) = ustar(i,j) * phimzf / (fkar*zf(1))*(vpcv/horv)
            dthldz(i,j) = - thlflux(i,j) / ustar(i,j) * phihzf / (fkar*zf(1))
            dqtdz (i,j) = - qtflux(i,j)  / ustar(i,j) * phihzf / (fkar*zf(1))
          end do
        end do

      end if

    else

      if(lneutral) then
        obl(:,:) = -1.e10
        oblav    = -1.e10
      else
        call getobl
      end if

      thlsl = 0.
      qtsl  = 0.
      
      if(lhetero) then
        lthls_patch = 0.0
        lqts_patch  = 0.0
        Npatch      = 0
      endif
      do j = 2, j1
        do i = 2, i1
          if(lhetero) then
            patchx = patchxnr(i)
            patchy = patchynr(j) 
          endif

          upcu   = 0.5 * (u0(i,j,1) + u0(i+1,j,1)) + cu
          vpcv   = 0.5 * (v0(i,j,1) + v0(i,j+1,1)) + cv
          horv   = sqrt(upcu ** 2. + vpcv ** 2.)
          horv   = max(horv, 1.e-2)
          horvav = sqrt(u0av(1) ** 2. + v0av(1) ** 2.)
          horvav = max(horvav, 1.e-2)

          if(lhetero) then
            patchx = patchxnr(i)
            patchy = patchynr(j) 
          endif
          
          if( isurf == 4) then
            if(lmostlocal) then
              ustar (i,j) = fkar * horv  / (log(zf(1) / z0m(i,j)) - psim(zf(1) / obl(i,j)) + psim(z0m(i,j) / obl(i,j)))
            else
              if(lhetero) then
                ustar (i,j) = fkar * horvpatch(patchx,patchy) / (log(zf(1) / z0m(i,j)) - psim(zf(1) / obl(i,j)) + psim(z0m(i,j) / obl(i,j)))
              else
                ustar (i,j) = fkar * horvav / (log(zf(1) / z0m(i,j)) - psim(zf(1) / obl(i,j)) + psim(z0m(i,j) / obl(i,j)))
              endif
            end if
          else
            if(lhetero) then
              ustar (i,j) = ustin_patch(patchx,patchy)
            else
              ustar (i,j) = ustin
            endif
          end if

          ustar  (i,j) = max(ustar(i,j), 1.e-2)
          if(lhetero) then
            thlflux(i,j) = wt_patch(patchx,patchynr(j)) 
            qtflux (i,j) = wq_patch(patchx,patchynr(j))
          else
            thlflux(i,j) = wtsurf
            qtflux (i,j) = wqsurf
          endif

          if(lhetero) then
            do n=1,nsv
              svflux(i,j,n) = wsv_patch(n,patchx,patchynr(j)) 
            enddo
          else
            do n=1,nsv
              svflux(i,j,n) = wsvsurf(n) 
            enddo
          endif

          if (obl(i,j) < 0.) then
            !phimzf = (1.-16.*zf(1)/obl)**(-0.25)
            phimzf = (1. + 3.6 * (-zf(1)/obl(i,j))**(2./3.))**(-0.5)
            !phihzf = (1.-16.*zf(1)/obl)**(-0.50)
            phihzf = (1. + 7.9 * (-zf(1)/obl(i,j))**(2./3.))**(-0.5)
          elseif (obl(i,j) > 0.) then
            phimzf = (1.+5.*zf(1)/obl(i,j))
            phihzf = (1.+5.*zf(1)/obl(i,j))
          else
            phimzf = 1.
            phihzf = 1.
          endif

          dudz  (i,j) = ustar(i,j) * phimzf / (fkar*zf(1))*(upcu/horv)
          dvdz  (i,j) = ustar(i,j) * phimzf / (fkar*zf(1))*(vpcv/horv)
          dthldz(i,j) = - thlflux(i,j) / ustar(i,j) * phihzf / (fkar*zf(1))
          dqtdz (i,j) = - qtflux(i,j)  / ustar(i,j) * phihzf / (fkar*zf(1))

          Cs(i,j) = fkar ** 2. / (log(zf(1) / z0m(i,j)) - psim(zf(1) / obl(i,j)) + psim(z0m(i,j) / obl(i,j))) / (log(zf(1) / z0h(i,j)) - psih(zf(1) / obl(i,j)) + psih(z0h(i,j) / obl(i,j)))

          if(lhetero) then
            tskin(i,j) = wt_patch(patchx,patchy) / (Cs(i,j) * horv) + thl0(i,j,1)
            qskin(i,j) = wq_patch(patchx,patchy) / (Cs(i,j) * horv) + qt0(i,j,1)
          else
            tskin(i,j) = wtsurf / (Cs(i,j) * horv) + thl0(i,j,1)
            qskin(i,j) = wqsurf / (Cs(i,j) * horv) + qt0(i,j,1)
          end if
          thlsl      = thlsl + tskin(i,j)
          qtsl       = qtsl  + qskin(i,j)
          if (lhetero) then
            lthls_patch(patchx,patchy) = lthls_patch(patchx,patchy) + tskin(i,j)
            lqts_patch(patchx,patchy)  = lqts_patch(patchx,patchy)  + qskin(i,j)
            Npatch(patchx,patchy)      = Npatch(patchx,patchy)      + 1
          endif
        end do
      end do

      call MPI_ALLREDUCE(thlsl, thls, 1,  MY_REAL, MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(qtsl, qts, 1,  MY_REAL, MPI_SUM, comm3d,mpierr)

      thls = thls / rslabs
      qts  = qtsl / rslabs
      thvs = thls * (1. + (rv/rd - 1.) * qts)

      if (lhetero) then
        call MPI_ALLREDUCE(lthls_patch(1:xpatches,1:ypatches), thls_patch(1:xpatches,1:ypatches), xpatches*ypatches,     MY_REAL, MPI_SUM, comm3d,mpierr)
        call MPI_ALLREDUCE(lqts_patch(1:xpatches,1:ypatches),  qts_patch(1:xpatches,1:ypatches),  xpatches*ypatches,     MY_REAL, MPI_SUM, comm3d,mpierr)
        call MPI_ALLREDUCE(Npatch(1:xpatches,1:ypatches)     , SNpatch(1:xpatches,1:ypatches)   , xpatches*ypatches, MPI_INTEGER ,MPI_SUM, comm3d,mpierr)
        thls_patch = thls_patch / SNpatch
        qts_patch  = qts_patch  / SNpatch
        thvs_patch = thls_patch * (1. + (rv/rd - 1.) * qts_patch)
      endif

    !if (lhetero) then
    !  thvs_patch = thls_patch * (1. + (rv/rd - 1.) * qts_patch)
    !endif
      !call qtsurf ! CvH check this!

    end if

    ! Transfer ustar to neighbouring cells
    do j=1,j2
      ustar(1,j)=ustar(i1,j)
      ustar(i2,j)=ustar(2,j)
    end do

    call excj( ustar  , 1, i2, 1, j2, 1,1)

    return

  end subroutine surface

!> Calculate the surface humidity assuming saturation.
  subroutine qtsurf
    use modglobal,   only : tmelt,bt,at,rd,rv,cp,es0,pref0,rslabs,i1,j1
    use modfields,   only : qt0
    use modsurfdata, only : rs, ra
    use modmpi,      only : my_real,mpierr,comm3d,mpi_sum,myid,mpi_integer

    implicit none
    real       :: exner, tsurf, qsatsurf, surfwet, es, qtsl
    integer    :: i,j, patchx, patchy
    integer    :: Npatch(xpatches,ypatches), SNpatch(xpatches,ypatches)
    real       :: lqts_patch(xpatches,ypatches)!, qts_patch(xpatches,ypatches)
    !real       :: tsurf_patch(xpatches,ypatches), exner_patch(xpatches,ypatches), es_patch(xpatches,ypatches)

    if(isurf <= 2) then
      qtsl = 0.
      do j = 2, j1
        do i = 2, i1
          exner      = (ps / pref0)**(rd/cp)
          tsurf      = tskin(i,j) * exner
          es         = es0 * exp(at*(tsurf-tmelt) / (tsurf-bt))
          qsatsurf   = rd / rv * es / ps
          surfwet    = ra(i,j) / (ra(i,j) + rs(i,j))
          qskin(i,j) = surfwet * qsatsurf + (1. - surfwet) * qt0(i,j,1)
          qtsl       = qtsl + qskin(i,j)
        end do
      end do

      call MPI_ALLREDUCE(qtsl, qts, 1,  MY_REAL, &
                         MPI_SUM, comm3d,mpierr)
      qts  = qts / rslabs

      if (lhetero) then 
        lqts_patch = 0.
        Npatch     = 0
        do j = 2, j1
          do i = 2, i1
            patchx     = patchxnr(i)
            patchy     = patchynr(j)
            exner      = (ps_patch(patchx,patchy) / pref0)**(rd/cp)
            tsurf      = tskin(i,j) * exner
            es         = es0 * exp(at*(tsurf-tmelt) / (tsurf-bt))
            qsatsurf   = rd / rv * es / ps_patch(patchx,patchy)
            surfwet    = ra(i,j) / (ra(i,j) + rs(i,j))
            qskin(i,j) = surfwet * qsatsurf + (1. - surfwet) * qt0(i,j,1)

            lqts_patch(patchx,patchy) = lqts_patch(patchx,patchy) + qskin(i,j)
            Npatch(patchx,patchy)     = Npatch(patchx,patchy)     + 1
          enddo
        enddo  
        call MPI_ALLREDUCE(lqts_patch(1:xpatches,1:ypatches), qts_patch(1:xpatches,1:ypatches), xpatches*ypatches,     MY_REAL,MPI_SUM, comm3d,mpierr)
        call MPI_ALLREDUCE(Npatch(1:xpatches,1:ypatches)    , SNpatch(1:xpatches,1:ypatches)  , xpatches*ypatches,MPI_INTEGER ,MPI_SUM, comm3d,mpierr)
        qts_patch = qts_patch / SNpatch
       
      endif
    !else
    !  if (lhetero) then
    !    exner_patch = (ps_patch/pref0)**(rd/cp)
    !    tsurf_patch = thls_patch * exner_patch
    !    es_patch    = es0*exp(at*(tsurf_patch-tmelt)/(tsurf_patch-bt))
    !    qts_patch   = rd/rv*es_patch/(ps_patch-(1-rd/rv)*es_patch)
    !    qts_patch   = max(qts_patch, 0.)

    !  endif
    !  exner = (ps/pref0)**(rd/cp)
    !  tsurf = thls*exner
    !  es    = es0*exp(at*(tsurf-tmelt)/(tsurf-bt))
    !  !qts   = rd/rv*es/(ps-(1-rd/rv)*es)
    !  qts   = rd/rv*es/ps ! specific hum instead of mr
    !  ! check to prevent collapse of spinup u* = 0 and  u = 0 free convection case
    !  !qts   = max(qts, 0.)
    !  qskin(:,:) = qts
    end if
    
    !thvs = thls * (1. + (rv/rd - 1.) * qts)
    !if (lhetero) then
    !  thvs_patch = thls_patch * (1. + (rv/rd - 1.) * qts_patch)
    !endif

    return

  end subroutine qtsurf

!> Calculates the Obuhkov length iteratively.
  subroutine getobl
    use modglobal, only : zf, rv, rd, grav, rslabs, i1, j1, i2, j2, timee, cu, cv
    use modfields, only : thl0av, qt0av, u0, v0, thl0, qt0, u0av, v0av
    use modmpi,    only : my_real,mpierr,comm3d,mpi_sum,myid,excj,mpi_integer
    implicit none

    integer             :: i,j,iter,patchx,patchy
    real                :: thv, thvsl, L, horv2, oblavl, thvpatch(xpatches,ypatches), horvpatch(xpatches,ypatches)
    real                :: Rib, Lstart, Lend, fx, fxdif, Lold
    real                :: upatch(xpatches,ypatches), vpatch(xpatches,ypatches)
    real                :: Supatch(xpatches,ypatches), Svpatch(xpatches,ypatches)
    integer             :: Npatch(xpatches,ypatches), SNpatch(xpatches,ypatches)
    real                :: lthlpatch(xpatches,ypatches), thlpatch(xpatches,ypatches), lqpatch(xpatches,ypatches), qpatch(xpatches,ypatches)
    real                :: loblpatch(xpatches,ypatches), oblpatch(xpatches,ypatches) 

    if(lmostlocal) then

      oblavl = 0.

      do i=2,i1
        do j=2,j1
          thv    = thl0(i,j,1)  * (1. + (rv/rd - 1.) * qt0(i,j,1))
          thvsl  = tskin(i,j)   * (1. + (rv/rd - 1.) * qskin(i,j))
          horv2 = u0(i,j,1)*u0(i,j,1) + v0(i,j,1)*v0(i,j,1)
          horv2 = max(horv2, 1.e-4)

          if(lhetero) then
            patchx = patchxnr(i)
            patchy = patchynr(j)
            Rib   = grav / thvs_patch(patchx,patchy) * zf(1) * (thv - thvs_patch(patchx,patchy)) / horv2
          else
            Rib   = grav / thvs * zf(1) * (thv - thvsl) / horv2
          endif

          iter = 0
          L = obl(i,j)

          if(Rib * L < 0. .or. abs(L) == 1e5) then
            if(Rib > 0) L = 0.01
            if(Rib < 0) L = -0.01
          end if

          do while (.true.)
            iter    = iter + 1
            Lold    = L
            fx      = Rib - zf(1) / L * (log(zf(1) / z0h(i,j)) - psih(zf(1) / L) + psih(z0h(i,j) / L)) / (log(zf(1) / z0m(i,j)) - psim(zf(1) / L) + psim(z0m(i,j) / L)) ** 2.
            Lstart  = L - 0.001*L
            Lend    = L + 0.001*L
            fxdif   = ( (- zf(1) / Lstart * (log(zf(1) / z0h(i,j)) - psih(zf(1) / Lstart) + psih(z0h(i,j) / Lstart)) / (log(zf(1) / z0m(i,j)) - psim(zf(1) / Lstart) + psim(z0m(i,j) / Lstart)) ** 2.) - (-zf(1) / Lend * (log(zf(1) / z0h(i,j)) - psih(zf(1) / Lend) + psih(z0h(i,j) / Lend)) / (log(zf(1) / z0m(i,j)) - psim(zf(1) / Lend) + psim(z0m(i,j) / Lend)) ** 2.) ) / (Lstart - Lend)
            L       = L - fx / fxdif
            if(Rib * L < 0. .or. abs(L) == 1e5) then
              if(Rib > 0) L = 0.01
              if(Rib < 0) L = -0.01
            end if
            if(abs(L - Lold) < 0.0001) exit
          end do

          obl(i,j) = L

        end do
      end do
    
    elseif(lhetero) then 
      upatch    = 0
      vpatch    = 0
      Npatch    = 0
      lthlpatch = 0
      lqpatch   = 0
      loblpatch = 0 

      do j = 2, j1
        do i = 2, i1
          patchx = patchxnr(i)
          patchy = patchynr(j)
          upatch(patchx,patchy)    = upatch(patchx,patchy)    + 0.5 * (u0(i,j,1) + u0(i+1,j,1))
          vpatch(patchx,patchy)    = vpatch(patchx,patchy)    + 0.5 * (v0(i,j,1) + v0(i,j+1,1))
          Npatch(patchx,patchy)    = Npatch(patchx,patchy)    + 1
          lthlpatch(patchx,patchy) = lthlpatch(patchx,patchy) + thl0(i,j,1)
          lqpatch(patchx,patchy)   = lqpatch(patchx,patchy)   + qt0(i,j,1)
          loblpatch(patchx,patchy) = loblpatch(patchx,patchy) + obl(i,j)
        enddo
      enddo

      call MPI_ALLREDUCE(upatch(1:xpatches,1:ypatches)   ,Supatch(1:xpatches,1:ypatches) ,xpatches*ypatches,    MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(vpatch(1:xpatches,1:ypatches)   ,Svpatch(1:xpatches,1:ypatches) ,xpatches*ypatches,    MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(Npatch(1:xpatches,1:ypatches)   ,SNpatch(1:xpatches,1:ypatches) ,xpatches*ypatches,MPI_INTEGER,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(lthlpatch(1:xpatches,1:ypatches),thlpatch(1:xpatches,1:ypatches),xpatches*ypatches,    MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(lqpatch(1:xpatches,1:ypatches)  ,qpatch(1:xpatches,1:ypatches)  ,xpatches*ypatches,    MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(loblpatch(1:xpatches,1:ypatches),oblpatch(1:xpatches,1:ypatches),xpatches*ypatches,    MY_REAL,MPI_SUM, comm3d,mpierr)
          
      horvpatch = sqrt(((Supatch/SNpatch) + cu) **2. + ((Svpatch/SNpatch) + cv) ** 2.)
      horvpatch = max(horvpatch, 1.e-2)

      thlpatch  = thlpatch / SNpatch
      qpatch    = qpatch   / SNpatch
      oblpatch  = oblpatch / SNpatch
          
      thvpatch  = thlpatch * (1. + (rv/rd - 1.) * qpatch)

      do patchy = 1, ypatches
        do patchx = 1, xpatches
          Rib   = grav / thvs_patch(patchx,patchy) * zf(1) * (thvpatch(patchx,patchy) - thvs_patch(patchx,patchy)) / horvpatch(patchx,patchy)
          iter = 0
          L = oblpatch(patchx,patchy)

          if(Rib * L < 0. .or. abs(L) == 1e5) then
            if(Rib > 0) L = 0.01
            if(Rib < 0) L = -0.01
          end if

          do while (.true.)
            iter    = iter + 1
            Lold    = L
            fx      = Rib - zf(1) / L * (log(zf(1) / z0hav_patch(patchx,patchy)) - psih(zf(1) / L) + psih(z0hav_patch(patchx,patchy) / L)) / (log(zf(1) / z0mav_patch(patchx,patchy)) - psim(zf(1) / L) + psim(z0mav_patch(patchx,patchy) / L)) ** 2.
            Lstart  = L - 0.001*L
            Lend    = L + 0.001*L
            fxdif   = ( (- zf(1) / Lstart * (log(zf(1) / z0hav_patch(patchx,patchy)) - psih(zf(1) / Lstart) + psih(z0hav_patch(patchx,patchy) / Lstart)) / (log(zf(1) / z0mav_patch(patchx,patchy)) - psim(zf(1) / Lstart) + psim(z0mav_patch(patchx,patchy) / Lstart)) ** 2.) - (-zf(1) / Lend * (log(zf(1) / z0hav_patch(patchx,patchy)) - psih(zf(1) / Lend) + psih(z0hav_patch(patchx,patchy) / Lend)) / (log(zf(1) / z0mav_patch(patchx,patchy)) - psim(zf(1) / Lend) + psim(z0mav_patch(patchx,patchy) / Lend)) ** 2.) ) / (Lstart - Lend)
            L       = L - fx / fxdif
            if(Rib * L < 0. .or. abs(L) == 1e5) then
              if(Rib > 0) L = 0.01
              if(Rib < 0) L = -0.01
            end if
            if(abs(L - Lold) < 0.0001) exit
          end do

          oblpatch(patchx,patchy) = L
        enddo
      enddo
      do i=1,i2
        do j=1,j2
          obl(i,j) = oblpatch(patchxnr(i),patchynr(j))
        enddo
      enddo
    endif 
    
    !CvH also do a global evaluation if lmostlocal = .true. to get an appropriate local mean

    thv    = thl0av(1) * (1. + (rv/rd - 1.) * qt0av(1))

    horv2 = u0av(1)**2. + v0av(1)**2.
    horv2 = max(horv2, 1.e-4)

    Rib   = grav / thvs * zf(1) * (thv - thvs) / horv2

    iter = 0
    L = oblav

    if(Rib * L < 0. .or. abs(L) == 1e5) then
      if(Rib > 0) L = 0.01
      if(Rib < 0) L = -0.01
    end if

    do while (.true.)
      iter    = iter + 1
      Lold    = L
      fx      = Rib - zf(1) / L * (log(zf(1) / z0hav) - psih(zf(1) / L) + psih(z0hav / L)) / (log(zf(1) / z0mav) - psim(zf(1) / L) + psim(z0mav / L)) ** 2.
      Lstart  = L - 0.001*L
      Lend    = L + 0.001*L
      fxdif   = ( (- zf(1) / Lstart * (log(zf(1) / z0hav) - psih(zf(1) / Lstart) + psih(z0hav / Lstart)) / (log(zf(1) / z0mav) - psim(zf(1) / Lstart) + psim(z0mav / Lstart)) ** 2.) - (-zf(1) / Lend * (log(zf(1) / z0hav) - psih(zf(1) / Lend) + psih(z0hav / Lend)) / (log(zf(1) / z0mav) - psim(zf(1) / Lend) + psim(z0mav / Lend)) ** 2.) ) / (Lstart - Lend)
      L       = L - fx / fxdif
      if(Rib * L < 0. .or. abs(L) == 1e5) then
        if(Rib > 0) L = 0.01
        if(Rib < 0) L = -0.01
      end if
      if(abs(L - Lold) < 0.0001) exit
    end do

    if(.not. lmostlocal) then
      if(.not. lhetero) then 
        obl(:,:) = L
      endif
    end if
    oblav = L

    return

  end subroutine getobl

  function psim(zeta)
    implicit none

    real             :: psim
    real, intent(in) :: zeta
    real             :: x

    if(zeta <= 0) then
      !x     = (1. - 16. * zeta) ** (0.25)
      !psim  = 3.14159265 / 2. - 2. * atan(x) + log( (1.+x) ** 2. * (1. + x ** 2.) / 8.)
      ! CvH use Wilson, 2001 rather than Businger-Dyer for correct free convection limit
      x     = (1. + 3.6 * abs(zeta) ** (2./3.)) ** (-0.5)
      psim = 3. * log( (1. + 1. / x) / 2.)
    else
      psim  = -2./3. * (zeta - 5./0.35)*exp(-0.35 * zeta) - zeta - (10./3.) / 0.35
    end if

    return
  end function psim

  function psih(zeta)

    implicit none

    real             :: psih
    real, intent(in) :: zeta
    real             :: x

    if(zeta <= 0) then
      !x     = (1. - 16. * zeta) ** (0.25)
      !psih  = 2. * log( (1. + x ** 2.) / 2. )
      ! CvH use Wilson, 2001
      x     = (1. + 7.9 * abs(zeta) ** (2./3.)) ** (-0.5)
      psih  = 3. * log( (1. + 1. / x) / 2.)
    else
      psih  = -2./3. * (zeta - 5./0.35)*exp(-0.35 * zeta) - (1. + (2./3.) * zeta) ** (1.5) - (10./3.) / 0.35 + 1.
    end if

    return

  end function psih

  function patchxnr(xpos)
    use modglobal,  only : imax
    implicit none
    integer             :: patchxnr
    integer             :: positionx
    integer, intent(in) :: xpos

    positionx = xpos - 2                                   !First grid point lies at i = 2. This lines makes sure that position = 0 for first grid point 
    if (positionx .lt. 0)    positionx = positionx + imax  !To account for border grid points
    if (positionx .ge. imax) positionx = positionx - imax  !To account for border grid points
    patchxnr  = 1 + (positionx*xpatches)/imax              !Converts position to patch number
    return
  end function

  function patchynr(ypos)
    use modmpi,     only : myid
    use modglobal,  only : jmax,jtot
    implicit none
    integer             :: patchynr
    integer             :: yposreal, positiony
    integer, intent(in) :: ypos
    
    yposreal  = ypos + (myid * jmax)                       !Converting the j position to the real j position by taking the processor number into account
    positiony = yposreal - 2                               !First grid point lies at j = 2. This lines makes sure that position = 0 for first grid point
    if (positiony .lt. 0)    positiony = positiony + jtot  !To account for border grid points
    if (positiony .ge. jtot) positiony = positiony - jtot  !To account for border grid points
    patchynr  = 1 + (positiony*ypatches)/jtot              !Converts position to patch number
    return
  end function
!>
  subroutine exitsurface
    implicit none
    return
  end subroutine exitsurface

end module modsurface
