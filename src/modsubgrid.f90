!!> \file modsubgrid.f90
!!!  Calculates and applies the Sub Filter Scale diffusion
!
!>
!!  Calculates and applies the Sub Filter Scale diffusion
!>
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!!  \todo Documentation
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
! This routine contains the Sullivan subgrid model coded by Victor Vertregt (September 2015)
! dalesinclsullivanorgineel.tar -> Gabls1_Sullivan_orgclosure as described in
! Sullivan et al., 1994, A subgrid-scale model for large-eddy simulation of planetary boundary-layer flows. BLM.
! Victor has tested various variants of the Sullivan model which are stored at 
! /nfs/livedata/victor/Les/Les_versions/Dales4/Sullivan
!
! modifications in: 
! modfields
! modgenstat
! modsubgriddata
! modsubgrid
! modfields
! 

module modsubgrid
use modsubgriddata
use modtimestat, only: calcblheight !Victor: om te bepalen wanneer het sullivan model uit kan

implicit none
save
  public :: subgrid, initsubgrid, exitsubgrid, subgridnamelist

contains
  subroutine initsubgrid
    use modglobal, only : ih,i1,jh,j1,k1,delta,zf,fkar,pi
    use modmpi, only : myid

    implicit none

    integer   :: k

    real :: ceps, ch
    real :: mlen

    call subgridnamelist

    allocate(ekm(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(ekh(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(zlt(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(sbdiss(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(sbshr(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(sbbuo(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(csz(k1))

 ! Victor: De rate of strain tensor definieren (symmetrisch dus 6 elementen)
    allocate(suu(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(suv(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(suw(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(svv(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(svw(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(sww(2-ih:i1+ih,2-jh:j1+jh,k1))
    
 
  
  allocate(suum(k1))
  allocate(suvm(k1))
  allocate(suwm(k1))
  allocate(svvm(k1))
  allocate(svwm(k1))
  allocate(swwm(k1))
  
  
  
  allocate(svaruul(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(svarvvl(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(svarwwl(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(svaruvl(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(svaruwl(2-ih:i1+ih,2-jh:j1+jh,k1))
  allocate(svarvwl(2-ih:i1+ih,2-jh:j1+jh,k1))
  
  
  allocate(svaruu(k1))
  allocate(svarvv(k1))
  allocate(svarww(k1))
  allocate(svaruv(k1))
  allocate(svaruw(k1))
  allocate(svarvw(k1))
  
  allocate(sacc(k1))
  allocate(sslab(k1))
  allocate(gamma(k1))
  allocate(sdemean(2-ih:i1+ih,2-jh:j1+jh,k1))
  
  allocate(ekm2(k1))
  allocate(uwr(2-ih:i1+ih,2-jh:j1+jh))
  allocate(vwr(2-ih:i1+ih,2-jh:j1+jh))


    ! Initialize variables to avoid problems when not using subgrid scheme JvdD
    ekm=0.; ekh=0.; zlt=0.; sbdiss=0.; sbshr=0.; sbbuo=0.; csz=0.

    ! Victor: doe hetzelfde met rate of strain tensor
    suu=0.; suv=0.; suw=0.; svv=0.; svw=0.; sww=0.;
   
    suum=0.; suvm=0.; suwm=0.; svvm=0.; svwm=0.; swwm=0.;
 
    svaruu=0.; svaruv=0.; svaruw=0.; svarvw=0.; svarvv=0.; svarww=0.;
     svaruul=0.; svaruvl=0.; svaruwl=0.; svarvvl=0.; svarvwl=0.; svarwwl=0.;
     sdemean=0.; 

    cm = cf / (2. * pi) * (1.5*alpha_kolm)**(-1.5)

!     ch   = 2. * alpha_kolm / beta_kolm
    ch   = 1.0/Prandtl
    ch2  = ch-ch1

    ceps = 2. * pi / cf * (1.5*alpha_kolm)**(-1.5)
    ce1  = (cn**2)* (cm/Rigc - ch1*cm )
    ce2  = ceps - ce1

    if(cs == -1.) then
      csz(:)  = (cm**3/ceps)**0.25   !< Smagorinsky constant
    else
      csz(:)  = cs
    end if

    if(lmason) then
      do k = 1,k1
        mlen   = (1. / (csz(k) * delta(k))**nmason + 1. / (fkar * zf(k))**nmason)**(-1./nmason)
        csz(k) = mlen / delta(k)
      end do
    end if

    if (myid==0) then
      write (6,*) 'cf    = ',cf
      write (6,*) 'cm    = ',cm
      write (6,*) 'ch    = ',ch
      write (6,*) 'ch1   = ',ch1
      write (6,*) 'ch2   = ',ch2
      write (6,*) 'ceps  = ',ceps
      write (6,*) 'ceps1 = ',ce1
      write (6,*) 'ceps2 = ',ce2
      write (6,*) 'cs    = ',cs
      write (6,*) 'Rigc  = ',Rigc
    endif

  end subroutine initsubgrid

  subroutine subgridnamelist
    use modglobal, only : ifnamopt,fname_options
    use modmpi,    only : myid, comm3d, mpierr, my_real, mpi_logical

    implicit none

    integer :: ierr

    namelist/NAMSUBGRID/ &
        ldelta,lmason, cf,cn,Rigc,Prandtl,lsmagorinsky,cs,nmason,sgs_surface_fix

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMSUBGRID,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMSUBGRID'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMSUBGRID'
      endif
      write(6 ,NAMSUBGRID)
      close(ifnamopt)
    end if

    call MPI_BCAST(ldelta     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lmason     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(nmason     ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(lsmagorinsky,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(cs         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(cf         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(cn         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(Rigc       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(Prandtl    ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(sgs_surface_fix ,1,MPI_LOGICAL   ,0,comm3d,mpierr)


  end subroutine subgridnamelist

  subroutine subgrid

 ! Diffusion subroutines
! Thijs Heus, Chiel van Heerwaarden, 15 June 2007

    use modglobal, only : nsv, lmoist
    use modfields, only : up,vp,wp,e12p,thl0,thlp,qt0,qtp,sv0,svp
    use modsurfdata,only : thlflux,qtflux,svflux
    implicit none
    integer n

    call calcblheight !Victor: om te bepalen wanneer het Sullivan model uit kan

    call closure
    call closure2 !Victor
    call diffu(up)
    call diffv(vp)
    call diffw(wp)
    if (.not. lsmagorinsky) call diffe(e12p)
    call diffc(thl0,thlp,thlflux)
    if (lmoist) call diffc( qt0, qtp, qtflux)
    do n=1,nsv
      call diffc(sv0(:,:,:,n),svp(:,:,:,n),svflux(:,:,n))
    end do
    if (.not. lsmagorinsky) call sources
  end subroutine

  subroutine exitsubgrid
    implicit none
    deallocate(ekm,ekh,zlt,sbdiss,sbbuo,sbshr,csz)
  end subroutine exitsubgrid

   subroutine closure

!-----------------------------------------------------------------|
!                                                                 |
!*** *closure*  calculates K-coefficients                         |
!                                                                 |
!      Hans Cuijpers   I.M.A.U.   06/01/1995                      |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!     All the K-closure factors are calculated.                   |
!                                                                 |
!     ekm(i,j,k) = k sub m : for velocity-closure                 |
!     ekh(i,j,k) = k sub h : for temperture-closure               |
!     ekh(i,j,k) = k sub h = k sub c : for concentration-closure  |
!                                                                 |
!     We will use the next model for these factors:               |
!                                                                 |
!     k sub m = 0.12 * l * sqrt(E)                                |
!                                                                 |
!     k sub h = k sub c = ( 1 + (2*l)/D ) * k sub m               |
!                                                                 |
!           where : l = mixing length  ( in model = z2 )          |
!                   E = subgrid energy                            |
!                   D = grid-size distance                        |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!             *closure* is called from *program*.                 |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal,  only : i1, j1,kmax,k1,ih,jh,i2,j2,delta,ekmin,grav, zf, fkar, &
                         dxi,dyi,dzf,dzh
  use modfields,  only : dthvdz,e120,u0,v0,w0,thvf
  use modsurfdata, only : dudz,dvdz,z0m
  use modmpi,    only : excjs
  implicit none

  real    :: strain2,mlen
  integer :: i,j,k,kp,km,jp,jm

  if(lsmagorinsky) then
    do k = 1,kmax
      mlen        = csz(k) * delta(k)

      do i = 2,i1
        do j = 2,j1

          kp=k+1
          km=k-1
          jp=j+1
          jm=j-1

          if(k == 1) then
            strain2 =  ( &
              ((u0(i+1,j,k)-u0(i,j,k))   *dxi        )**2    + &
              ((v0(i,jp,k)-v0(i,j,k))    *dyi        )**2    + &
              ((w0(i,j,kp)-w0(i,j,k))    /dzf(k)     )**2    )

            strain2 = strain2 + 0.5 * ( &
              ( 0.25*(w0(i+1,j,kp)-w0(i-1,j,kp))*dxi + &
              dudz(i,j)   )**2 )

            strain2 = strain2 + 0.125 * ( &
              ((u0(i,jp,k)-u0(i,j,k))     *dyi     + &
              (v0(i,jp,k)-v0(i-1,jp,k))  *dxi        )**2    + &
              ((u0(i,j,k)-u0(i,jm,k))     *dyi     + &
              (v0(i,j,k)-v0(i-1,j,k))    *dxi        )**2    + &
              ((u0(i+1,j,k)-u0(i+1,jm,k)) *dyi     + &
              (v0(i+1,j,k)-v0(i,j,k))    *dxi        )**2    + &
              ((u0(i+1,jp,k)-u0(i+1,j,k)) *dyi     + &
              (v0(i+1,jp,k)-v0(i,jp,k))  *dxi        )**2    )

            strain2 = strain2 + 0.5 * ( &
              ( 0.25*(w0(i,jp,kp)-w0(i,jm,kp))*dyi + &
              dvdz(i,j)   )**2 )
      
          else

            strain2 =  ( &
              ((u0(i+1,j,k)-u0(i,j,k))   *dxi        )**2    + &
              ((v0(i,jp,k)-v0(i,j,k))    *dyi        )**2    + &
              ((w0(i,j,kp)-w0(i,j,k))    /dzf(k)     )**2    )

            strain2 = strain2 + 0.125 * ( &
              ((w0(i,j,kp)-w0(i-1,j,kp))  *dxi     + &
              (u0(i,j,kp)-u0(i,j,k))     / dzh(kp)  )**2    + &
              ((w0(i,j,k)-w0(i-1,j,k))    *dxi     + &
              (u0(i,j,k)-u0(i,j,km))     / dzh(k)   )**2    + &
              ((w0(i+1,j,k)-w0(i,j,k))    *dxi     + &
              (u0(i+1,j,k)-u0(i+1,j,km)) / dzh(k)   )**2    + &
              ((w0(i+1,j,kp)-w0(i,j,kp))  *dxi     + &
              (u0(i+1,j,kp)-u0(i+1,j,k)) / dzh(kp)  )**2    )

            strain2 = strain2 + 0.125 * ( &
              ((u0(i,jp,k)-u0(i,j,k))     *dyi     + &
              (v0(i,jp,k)-v0(i-1,jp,k))  *dxi        )**2    + &
              ((u0(i,j,k)-u0(i,jm,k))     *dyi     + &
              (v0(i,j,k)-v0(i-1,j,k))    *dxi        )**2    + &
              ((u0(i+1,j,k)-u0(i+1,jm,k)) *dyi     + &
              (v0(i+1,j,k)-v0(i,j,k))    *dxi        )**2    + &
              ((u0(i+1,jp,k)-u0(i+1,j,k)) *dyi     + &
              (v0(i+1,jp,k)-v0(i,jp,k))  *dxi        )**2    )

            strain2 = strain2 + 0.125 * ( &
              ((v0(i,j,kp)-v0(i,j,k))     / dzh(kp) + &
              (w0(i,j,kp)-w0(i,jm,kp))   *dyi        )**2    + &
              ((v0(i,j,k)-v0(i,j,km))     / dzh(k)+ &
              (w0(i,j,k)-w0(i,jm,k))     *dyi        )**2    + &
              ((v0(i,jp,k)-v0(i,jp,km))   / dzh(k)+ &
              (w0(i,jp,k)-w0(i,j,k))     *dyi        )**2    + &
              ((v0(i,jp,kp)-v0(i,jp,k))   / dzh(kp) + &
              (w0(i,jp,kp)-w0(i,j,kp))   *dyi        )**2    )
          end if

          ekm(i,j,k)  = mlen ** 2. * sqrt(2. * strain2)
          ekh(i,j,k)  = ekm(i,j,k) / Prandtl
        end do
      end do
    end do

  ! do TKE scheme
  else
    do k=1,kmax
      do j=2,j1
        do i=2,i1
          if (ldelta .or. (dthvdz(i,j,k)<=0)) then
            zlt(i,j,k) = delta(k)
            if (lmason) zlt(i,j,k) = (1. / zlt(i,j,k) ** nmason + 1. / ( fkar * (zf(k) + z0m(i,j)))**nmason) ** (-1./nmason)
            ekm(i,j,k) = cm * zlt(i,j,k) * e120(i,j,k)
            ekh(i,j,k) = (ch1 + ch2) * ekm(i,j,k)
          else
            zlt(i,j,k) = min(delta(k),cn*e120(i,j,k)/sqrt(grav/thvf(k)*abs(dthvdz(i,j,k))))
            if (lmason) zlt(i,j,k) = (1. / zlt(i,j,k) ** nmason + 1. / ( fkar * (zf(k) + z0m(i,j)))**nmason) ** (-1./nmason)

            ekm(i,j,k) = cm * zlt(i,j,k) * e120(i,j,k)
            ekh(i,j,k) = (ch1 + ch2 * zlt(i,j,k)/delta(k)) * ekm(i,j,k)
          endif
        end do
      end do
    end do
  end if

  do k=1,k1
    do j=2,j1
      do i=2,i1
        ekm(i,j,k) = max(ekm(i,j,k),ekmin)
        ekh(i,j,k) = max(ekh(i,j,k),ekmin)
      end do
    end do
  end do
!*************************************************************
!     Set cyclic boundary condition for K-closure factors.
!*************************************************************

  ekm(1, :,:) = ekm(i1,:,:)
  ekm(i2,:,:) = ekm(2, :,:)
  ekh(1, :,:) = ekh(i1,:,:)
  ekh(i2,:,:) = ekh(2, :,:)

  call excjs( ekm           , 2,i1,2,j1,1,k1,ih,jh)
  call excjs( ekh           , 2,i1,2,j1,1,k1,ih,jh)

  do j=1,j2
    do i=1,i2
      ekm(i,j,k1)  = ekm(i,j,kmax)
      ekh(i,j,k1)  = ekh(i,j,kmax)
    end do
  end do

  return
  end subroutine closure


  subroutine sources


!-----------------------------------------------------------------|
!                                                                 |
!*** *sources*                                                    |
!      calculates various terms from the subgrid TKE equation     |
!                                                                 |
!     Hans Cuijpers   I.M.A.U.     06/01/1995                     |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!      Subroutine sources calculates all other terms in the       |
!      subgrid energy equation, except for the diffusion terms.   |
!      These terms are calculated in subroutine diff.             |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!     *sources* is called from *program*.                         |
!                                                                 |
!-----------------------------------------------------------------|


  use modglobal,   only : i1,j1,kmax,delta,dx,dy,dxi,dyi,dzf,zf,dzh,grav, imax, jtot, rslabs, k1
  use modfields,   only : u0,v0,w0,e120,e12p,dthvdz,thvh,rhobf,thvf
  use modsurfdata,  only : dudz,dvdz
  use modmpi,    only : excjs, myid, nprocs, comm3d, mpierr, my_real, mpi_sum, slabsum
  use modsubgriddata, only: zisullivan

  implicit none

  real    tdef2, zinv
  integer i,j,k,jm,jp,km,kp, ip, im 
  
  !real, allocatable :: thetawflux(:)
  !real, allocatable :: thetawfluxl(:)
  
  !allocate (thetawflux(k1),thetawfluxl(k1))
 ! thetawflux=0.;thetawfluxl=0.;
  
  ! Onderstaande loop voor het inplementeren van het Sullivanmodel: strainij moet uitgerekend worden, en slabgemiddeld nemen met MPI_Reduce
  
  !---------------------------------------------------------------------------------------------------------------------------------------------------
  
  !Victor: bij zowel de smagorinsky, als het orginele model, en ook de subroutine sources is de rate of strain tensor Sij nodig, deze wordt daarom hieronder uitgerekend en kan daarna gebruikt worden door vervolgende modules.
  
  

  do k=1,kmax 
  
    do i=2,i1
      do j=2,j1
   
          kp=k+1
          km=k-1
          jp=j+1
          jm=j-1
          ip=i+1
          im=i-1
       
          suu(i,j,k) = (u0(ip,j,k)-u0(i,j,k)) * dxi
          svv(i,j,k) = (v0(i,jp,k)-v0(i,j,k)) * dyi
          sww(i,j,k) = (w0(i,j,kp)-w0(i,j,k)) / dzf(k)
          
          suv(i,j,k) = 0.125 * &
                       (((u0(i,jp,k)-u0(i,j,k)) * dyi + (v0(i,jp,k)-v0(im,jp,k)) * dxi ) + & 
                        ((u0(i,j,k)-u0(i,jm,k)) * dyi + (v0(i,j,k)-v0(im,j,k)) *dxi ) + &
                        ((u0(ip,j,k)-u0(ip,jm,k)) * dyi + (v0(ip,j,k)-v0(i,j,k)) *dxi ) + &
                        ((u0(ip,jp,k)-u0(ip,j,k)) *dyi + (v0(ip,jp,k)-v0(i,jp,k)) *dxi ) )

          if (k == 1) then 
          
          suw(i,j,k) = 0.25 * ( ( w0(i+1,j,kp)-w0(i-1,j,kp) )*dxi + dudz(i,j)) 
          
          svw(i,j,k) = 0.25 * ( (w0(i,jp,kp)-w0(i,jm,kp))*dyi + dvdz(i,j)) 
    
         else
          
          suw(i,j,k) = 0.125 * &
                       (((w0(i,j,kp)-w0(i-1,j,kp)) * dxi +  (u0(i,j,kp)-u0(i,j,k)) / dzh(kp) )  + &
                        ((w0(i,j,k)-w0(i-1,j,k)) *dxi + (u0(i,j,k)-u0(i,j,km)) / dzh(k) ) + &
                        ((w0(i+1,j,k)-w0(i,j,k)) *dxi + (u0(i+1,j,k)-u0(i+1,j,km)) / dzh(k) ) + &
                        ((w0(i+1,j,kp)-w0(i,j,kp))  *dxi + (u0(i+1,j,kp)-u0(i+1,j,k)) / dzh(kp) ) )

          svw(i,j,k) = 0.125 * &
                       (((v0(i,j,kp)-v0(i,j,k)) / dzh(kp) + (w0(i,j,kp)-w0(i,jm,kp)) * dyi ) + &
                        ((v0(i,j,k)-v0(i,j,km)) / dzh(k) + (w0(i,j,k)-w0(i,jm,k)) * dyi ) + &
                        ((v0(i,jp,k)-v0(i,jp,km)) / dzh(k) + (w0(i,jp,k)-w0(i,j,k)) * dyi ) + &
                        ((v0(i,jp,kp)-v0(i,jp,k)) / dzh(kp) + (w0(i,jp,kp)-w0(i,j,kp)) * dyi ) )
          
          end if
          
     
          
          
      end do
    end do
    end do 
    
    call slabsum(suum ,1,kmax,suu  ,2,i1,2,j1,1,kmax,2,i1,2,j1,1,kmax)
    call slabsum(suvm ,1,kmax,suv  ,2,i1,2,j1,1,kmax,2,i1,2,j1,1,kmax)
    call slabsum(suwm ,1,kmax,suw  ,2,i1,2,j1,1,kmax,2,i1,2,j1,1,kmax)
    call slabsum(svvm ,1,kmax,svv  ,2,i1,2,j1,1,kmax,2,i1,2,j1,1,kmax)
    call slabsum(svwm ,1,kmax,svw  ,2,i1,2,j1,1,kmax,2,i1,2,j1,1,kmax)
    call slabsum(swwm ,1,kmax,sww  ,2,i1,2,j1,1,kmax,2,i1,2,j1,1,kmax)
    

    
    suum=suum/rslabs ! Victor: is dit inderdaad de totale domeingrootte?
    suvm=suvm/rslabs
    suwm=suwm/rslabs
    svvm=svvm/rslabs
    svwm=svwm/rslabs
    swwm=swwm/rslabs
    

    
   

    
    !http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance gebruiken! Hieronder "lompe" methode, daarna een sneller algoritme inbouwen
    
    do k=1,kmax
      do i=2,i1
      do j=2,j1
      
      svaruul(i,j,k)=(suu(i,j,k)-suum(k))**2
      svaruvl(i,j,k)=(suv(i,j,k)-suvm(k))**2
      svaruwl(i,j,k)=(suw(i,j,k)-suwm(k))**2
      svarvvl(i,j,k)=(svv(i,j,k)-svvm(k))**2
      svarvwl(i,j,k)=(svw(i,j,k)-svwm(k))**2
      svarwwl(i,j,k)=(sww(i,j,k)-swwm(k))**2
      
      end do
      end do
      end do
      
      call slabsum(svaruu ,1,kmax,svaruul ,2,i1,2,j1,1,kmax,2,i1,2,j1,1,kmax)
      call slabsum(svaruv ,1,kmax,svaruvl ,2,i1,2,j1,1,kmax,2,i1,2,j1,1,kmax)
      call slabsum(svaruw ,1,kmax,svaruwl ,2,i1,2,j1,1,kmax,2,i1,2,j1,1,kmax)
      call slabsum(svarvv ,1,kmax,svarvvl ,2,i1,2,j1,1,kmax,2,i1,2,j1,1,kmax)
      call slabsum(svarvw ,1,kmax,svarvwl ,2,i1,2,j1,1,kmax,2,i1,2,j1,1,kmax)
      call slabsum(svarww ,1,kmax,svarwwl ,2,i1,2,j1,1,kmax,2,i1,2,j1,1,kmax)
      
      !call MPI_ALLREDUCE(svar11l(k),svar11(k),1,MY_REAL,MPI_SUM,comm3d,mpierr)
      !call MPI_ALLREDUCE(svar12l(k),svar12(k),1,MY_REAL,MPI_SUM,comm3d,mpierr)
      !call MPI_ALLREDUCE(svar13l(k),svar13(k),1,MY_REAL,MPI_SUM,comm3d,mpierr)
      !call MPI_ALLREDUCE(svar22l(k),svar22(k),1,MY_REAL,MPI_SUM,comm3d,mpierr)
      !call MPI_ALLREDUCE(svar23l(k),svar23(k),1,MY_REAL,MPI_SUM,comm3d,mpierr)
      !call MPI_ALLREDUCE(svar33l(k),svar33(k),1,MY_REAL,MPI_SUM,comm3d,mpierr)
      
      svaruu=svaruu/(rslabs)
      svaruv=svaruv/(rslabs)
      svaruw=svaruw/(rslabs)
      svarvv=svarvv/(rslabs)
      svarvw=svarvw/(rslabs)
      svarww=svarww/(rslabs)
      
      do k=1,kmax
      sacc(k)=sqrt(2*svaruu(k)+2*svarvv(k)+2*svarww(k)+4*svaruv(k)+4*svaruw(k)+4*svarvw(k))
             
    
              
    sslab(k)=sqrt(2*(suum(k))**2+2*(svvm(k))**2+2*(swwm(k))**2+ & 
               4*(suvm(k))**2+4*(suwm(k))**2+4*(svwm(k))**2)
      
    if ((sacc(k)+sslab(k)).gt.0) then
        gamma(k)=sacc(k)/(sacc(k)+sslab(k)) 
    else
        gamma(k)=0
    end if
    
    if (gamma(k).lt.0.01) then
       gamma(k)=0.01
    end if
 
 !Victor: meer boundary layer implementation
 
! do i=2,i1
 ! do j=2,j1
 ! thetawfluxl(k)=thetawfluxl(k)+ekh(i,j,k)*dthvdz(i,j,k)  
 ! end do
! end do
 
!call MPI_ALLREDUCE(thetawflux(k),thetawfluxl(k),1,MY_REAL,MPI_SUM,comm3d,mpierr)

!thetawflux(k) = thetawflux(k)/rslabs
 
 
 !Victor: Boundary layer implementation
    if ((zf(k).gt.(zisullivan/2))) then !.OR.(thetawflux(k)< 0.01*thetawflux(1)) !Victor: ook het onderste level helemaal Sullivan vrij houden? .OR.(k==1)
       zinv=zf(k)
       gamma(k)=1
       suum(k)=0
       suvm(k)=0
       suwm(k)=0
       svvm(k)=0
       svwm(k)=0
       swwm(k)=0
   end if
   
   
  
 
      
    do i=2,i1
      do j=2,j1
      
      
      
      sdemean(i,j,k)=(suu(i,j,k)-suum(k))**2+(svv(i,j,k)-svvm(k))**2+(sww(i,j,k)-swwm(k))**2 + &
              2*(suv(i,j,k)-suvm(k))**2+2*(suw(i,j,k)-suwm(k))**2+2*(svw(i,j,k)-svwm(k))**2
      
      
      sbshr(i,j,k)  = (ekm(i,j,k)/rhobf(k))*gamma(k)*sdemean(i,j,k)/ ( 2*e120(i,j,k))
      sbbuo(i,j,k)  = -(ekh(i,j,k)/rhobf(k))*grav/thvf(k)*dthvdz(i,j,k)/ ( 2*e120(i,j,k))
      sbdiss(i,j,k) = - (ce1 + ce2*zlt(i,j,k)/delta(k)) * e120(i,j,k)**2 /(2.*zlt(i,j,k))
      
      
      
      
      end do
    end do
    
    
  end do
  
   
  
  e12p(2:i1,2:j1,1:kmax) = e12p(2:i1,2:j1,1:kmax)+ &
            sbshr(2:i1,2:j1,1:kmax)+sbbuo(2:i1,2:j1,1:kmax)+sbdiss(2:i1,2:j1,1:kmax)
  
    
  !deallocate (thetawflux,thetawfluxl)
  return
  end subroutine sources


  subroutine diffc (putin,putout,flux)

    use modglobal, only : i1,ih,i2,j1,jh,j2,k1,kmax,dx2i,dzf,dy2i,dzh
    use modfields, only : rhobf,rhobh
    implicit none

    real, intent(in)    :: putin(2-ih:i1+ih,2-jh:j1+jh,k1)
    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real, intent(in)    :: flux (i2,j2)

    integer i,j,k,jm,jp,km,kp

    do k=2,kmax
      kp=k+1
      km=k-1

      do j=2,j1
        jp=j+1
        jm=j-1

        do i=2,i1
          putout(i,j,k) = putout(i,j,k) &
                    +  0.5 * ( &
                  ( (ekh(i+1,j,k)+ekh(i,j,k))*(putin(i+1,j,k)-putin(i,j,k)) &
                    -(ekh(i,j,k)+ekh(i-1,j,k))*(putin(i,j,k)-putin(i-1,j,k)))*dx2i &
                    + &
                  ( (ekh(i,jp,k)+ekh(i,j,k)) *(putin(i,jp,k)-putin(i,j,k)) &
                    -(ekh(i,j,k)+ekh(i,jm,k)) *(putin(i,j,k)-putin(i,jm,k)) )*dy2i &
                  + &
                  ( rhobh(kp)/rhobf(k) * (dzf(kp)*ekh(i,j,k) + dzf(k)*ekh(i,j,kp)) &
                    *  (putin(i,j,kp)-putin(i,j,k)) / dzh(kp)**2 &
                    - &
                    rhobh(k)/rhobf(k) * (dzf(km)*ekh(i,j,k) + dzf(k)*ekh(i,j,km)) &
                    *  (putin(i,j,k)-putin(i,j,km)) / dzh(k)**2           )/dzf(k) &
                            )

        end do
      end do
    end do

    do j=2,j1
      do i=2,i1

        putout(i,j,1) = putout(i,j,1) &
                  + 0.5 * ( &
                ( (ekh(i+1,j,1)+ekh(i,j,1))*(putin(i+1,j,1)-putin(i,j,1)) &
                  -(ekh(i,j,1)+ekh(i-1,j,1))*(putin(i,j,1)-putin(i-1,j,1)) )*dx2i &
                  + &
                ( (ekh(i,j+1,1)+ekh(i,j,1))*(putin(i,j+1,1)-putin(i,j,1)) &
                  -(ekh(i,j,1)+ekh(i,j-1,1))*(putin(i,j,1)-putin(i,j-1,1)) )*dy2i &
                  + &
                ( rhobh(2)/rhobf(1) * (dzf(2)*ekh(i,j,1) + dzf(1)*ekh(i,j,2)) &
                  *  (putin(i,j,2)-putin(i,j,1)) / dzh(2)**2 &
                  + rhobh(1)/rhobf(1)*flux(i,j) *2.                        )/dzf(1) &
                          )

      end do
    end do

  end subroutine diffc



  subroutine diffe(putout)

    use modglobal, only : i1,ih,j1,jh,k1,kmax,dx2i,dzf,dy2i,dzh
    use modfields, only : e120,rhobf,rhobh
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    integer             :: i,j,k,jm,jp,km,kp

    do k=2,kmax
      kp=k+1
      km=k-1

      do j=2,j1
        jp=j+1
        jm=j-1

        do i=2,i1

          putout(i,j,k) = putout(i,j,k) &
                  +  ( &
              ((ekm(i+1,j,k)+ekm(i,j,k))*(e120(i+1,j,k)-e120(i,j,k)) &
              -(ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k)))*dx2i &
                  + &
              ((ekm(i,jp,k)+ekm(i,j,k)) *(e120(i,jp,k)-e120(i,j,k)) &
              -(ekm(i,j,k)+ekm(i,jm,k)) *(e120(i,j,k)-e120(i,jm,k)) )*dy2i &
                  + &
              (rhobh(kp)/rhobf(k) * (dzf(kp)*ekm(i,j,k) + dzf(k)*ekm(i,j,kp)) &
              *(e120(i,j,kp)-e120(i,j,k)) / dzh(kp)**2 &
              - rhobh(k)/rhobf(k) * (dzf(km)*ekm(i,j,k) + dzf(k)*ekm(i,j,km)) &
              *(e120(i,j,k)-e120(i,j,km)) / dzh(k)**2        )/dzf(k) &
                            )

        end do
      end do
    end do

  !     --------------------------------------------
  !     special treatment for lowest full level: k=1
  !     --------------------------------------------

    do j=2,j1
      do i=2,i1

        putout(i,j,1) = putout(i,j,1) + &
            ( (ekm(i+1,j,1)+ekm(i,j,1))*(e120(i+1,j,1)-e120(i,j,1)) &
              -(ekm(i,j,1)+ekm(i-1,j,1))*(e120(i,j,1)-e120(i-1,j,1)) )*dx2i &
            + &
            ( (ekm(i,j+1,1)+ekm(i,j,1))*(e120(i,j+1,1)-e120(i,j,1)) &
              -(ekm(i,j,1)+ekm(i,j-1,1))*(e120(i,j,1)-e120(i,j-1,1)) )*dy2i &
            + &
              ( rhobh(2)/rhobf(1) * (dzf(2)*ekm(i,j,1) + dzf(1)*ekm(i,j,2)) &
              *  (e120(i,j,2)-e120(i,j,1)) / dzh(2)**2              )/dzf(1)

      end do
    end do

  end subroutine diffe


  
  !-----------------------------------------------------------------------------------------
  ! Victor: hier verder met het Sullivan model bouwen. Grote uitdaging: ekm2 --> zie closure2
  
  
  

  subroutine diffu (putout)

    use modglobal, only : i1,ih,i2,j1,jh,j2,k1,kmax,dxi,dx2i,dzf,dy,dyi,dy2i,dzh, cu,cv
    use modfields, only : u0,v0,w0,rhobf,rhobh
    use modsurfdata,only : ustar
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: emmo,emom,emop,empo
    real                :: fu
    real                :: ucu, upcu
    real                :: gammam, suwmm, ekm2m, gammap, suwmp, ekm2p !Victor: geinterpoleerde versies van de orginiele variabele
    integer             :: i,j,k,jm,jp,km,kp

    do k=2,kmax
      kp=k+1
      km=k-1

      do j=2,j1
        jp=j+1
        jm=j-1

        do i=2,i1

          emom = ( dzf(km) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i-1,j,km) ) ) / &
                    ( 4.   * dzh(k) )

          emop = ( dzf(kp) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,kp) + ekm(i-1,j,kp) ) ) / &
                    ( 4.   * dzh(kp) )

          empo = 0.25 * ( &
                  ekm(i,j,k)+ekm(i,jp,k)+ekm(i-1,jp,k)+ekm(i-1,j,k)  )

          emmo = 0.25 * ( &
                  ekm(i,j,k)+ekm(i,jm,k)+ekm(i-1,jm,k)+ekm(i-1,j,k)  )
                  
         gammam = 0.5 * (dzf(k-1)*gamma(k) + dzf(k)*gamma(k-1) )/dzh(k)
         
         suwmm = 0.5 * (dzf(k-1)*suwm(k) + dzf(k)*suwm(k-1) )/dzh(k)
         
         ekm2m = 0.5 * ( dzf(k-1)*ekm2(k) + dzf(k)*ekm2(k-1) )/dzh(k)
         
         gammap = 0.5 * (dzf(k)*gamma(k+1) + dzf(k+1)*gamma(k) )/dzh(k+1)
         
         suwmp = 0.5 * (dzf(k+1)*suwm(k) + dzf(k)*suwm(k+1) )/dzh(k+1)
         
         ekm2p = 0.5 * ( dzf(k+1)*ekm2(k) + dzf(k)*ekm2(k+1) )/dzh(k+1)

         !Victor: pas op met ekm2: gedefinieerd in hetzelfde vlak als w!!!!!!

          putout(i,j,k) = putout(i,j,k) &
                  + (1./rhobf(k)) * &
                  ( ekm(i,j,k) *gamma(k) * (u0(i+1,j,k)-u0(i,j,k)  ) + ekm2p*suum(k) &
                    -ekm(i-1,j,k)*gamma(k)* (u0(i,j,k)-u0(i-1,j,k)) - ekm2p*suum(k ) ) * 2. * dx2i &   ! victor: omdat suum en ekm2 gelijk zijn in het vlak, verandert er hier niets!
                  + (1./rhobf(k)) * &
                  ( empo * gamma(k)* ( (u0(i,jp,k)-u0(i,j,k))   *dyi &
                            +(v0(i,jp,k)-v0(i-1,jp,k))*dxi) +ekm2p*suvm(k) &
                    -emmo * gamma(k) * ( (u0(i,j,k)-u0(i,jm,k))   *dyi &
                            +(v0(i,j,k)-v0(i-1,j,k))  *dxi)  -ekm2p*suvm(k)  ) / dy &
                  + (1./rhobf(k)) * &
                  ( emop * gammap * ( (u0(i,j,kp)-u0(i,j,k))   /dzh(kp) &
                            +(w0(i,j,kp)-w0(i-1,j,kp))*dxi) + ekm2(kp)*suwmp &
                    -emom * gammam* ( (u0(i,j,k)-u0(i,j,km))   /dzh(k) &
                            +(w0(i,j,k)-w0(i-1,j,k))  *dxi) -ekm2(k)*suwmm  ) /dzf(k)

        end do
      end do
    end do

  !     --------------------------------------------
  !     special treatment for lowest full level: k=1
  !     --------------------------------------------

    do j=2,j1
      jp = j+1
      jm = j-1

      do i=2,i1

        empo = 0.25 * ( &
              ekm(i,j,1)+ekm(i,jp,1)+ekm(i-1,jp,1)+ekm(i-1,j,1)  )

        emmo = 0.25 * ( &
              ekm(i,j,1)+ekm(i,jm,1)+ekm(i-1,jm,1)+ekm(i-1,j,1)  )

        emop = ( dzf(2) * ( ekm(i,j,1) + ekm(i-1,j,1) )  + &
                    dzf(1) * ( ekm(i,j,2) + ekm(i-1,j,2) ) ) / &
                  ( 4.   * dzh(2) )
                  
         suwmp = 0.5 * (dzf(2)*suwm(1) + dzf(1)*suwm(2) )/dzh(2)
         
         ekm2p = 0.5 * ( dzf(2)*ekm2(1) + dzf(1)*ekm2(2) )/dzh(2)
         
         gammap=0.5 *(dzf(2)*gamma(1) + dzf(1)*gamma(2) )/dzh(2)


        ucu   = 0.5*(u0(i,j,1)+u0(i+1,j,1))+cu

        if(ucu >= 0.) then
          upcu  = max(ucu,1.e-10)
        else
          upcu  = min(ucu,-1.e-10)
        end if


        fu = ( 0.5*( ustar(i,j)+ustar(i-1,j) ) )**2  * &
                upcu/sqrt(upcu**2  + &
                ((v0(i,j,1)+v0(i-1,j,1)+v0(i,jp,1)+v0(i-1,jp,1))/4.+cv)**2)

        putout(i,j,1) = putout(i,j,1) &
                + (1./rhobf(1))*&
              ( ekm(i,j,1) *gamma(1) * (u0(i+1,j,1)-u0(i,j,1))&
              -ekm(i-1,j,k)*gamma(1)* (u0(i,j,1)-u0(i-1,j,1))  ) * 2. * dx2i &
                + (1./rhobf(1)) * &
              ( empo * gamma(1)* ( (u0(i,jp,1)-u0(i,j,1))   *dyi &
                        +(v0(i,jp,1)-v0(i-1,jp,1))*dxi)  &
              -emmo  *gamma(1)*( (u0(i,j,1)-u0(i,jm,1))   *dyi &
                        +(v0(i,j,1)-v0(i-1,j,1))  *dxi)  ) / dy &
               + (1./rhobf(1))* &
              ( emop * gammap* ( (u0(i,j,2)-u0(i,j,1))    /dzh(2) &
                        +(w0(i,j,2)-w0(i-1,j,2))  *dxi) + ekm2(2)*suwmp & !Victor: ook hier weer een wijziging!!
                -rhobh(1)*fu  ) / dzf(1)  !Victor:  Zitten op deze plekken de grote fout!?!?!

      end do
    end do

  end subroutine diffu


  subroutine diffv (putout)

    use modglobal, only : i1,ih,i2,j1,jh,j2,k1,kmax,dx,dxi,dx2i,dzf,dyi,dy2i,dzh, cu,cv
    use modfields, only : u0,v0,w0,rhobf,rhobh
    use modsurfdata,only : ustar

    implicit none
    
    ! Hier even alleen de vw component inbouwen: rest wordt toch nul!

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: emmo, eomm,eomp,epmo
    real                :: fv, vcv,vpcv
    integer             :: i,j,k,jm,jp,km,kp
    real                :: gammam, svwmm, ekm2m, gammap, svwmp, ekm2p !Victor: geinterpoleerde versies van de orginiele variabele

    do k=2,kmax
      kp=k+1
      km=k-1

      do j=2,j1
        jp=j+1
        jm=j-1

        do i=2,i1

          eomm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i,jm,km) ) ) / &
                    ( 4.   * dzh(k) )

          eomp = ( dzf(kp) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                      dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) / &
                    ( 4.   * dzh(kp) )

          emmo = 0.25  * ( &
                ekm(i,j,k)+ekm(i,jm,k)+ekm(i-1,jm,k)+ekm(i-1,j,k)  )

          epmo = 0.25  * ( &
                ekm(i,j,k)+ekm(i,jm,k)+ekm(i+1,jm,k)+ekm(i+1,j,k)  )
                
         gammam = 0.5 * (dzf(k-1)*gamma(k) + dzf(k)*gamma(k-1) )/dzh(k)
         
         svwmm = 0.5 * (dzf(k-1)*svwm(k) + dzf(k)*svwm(k-1) )/dzh(k)
         
         ekm2m = 0.5 * ( dzf(k-1)*ekm2(k) + dzf(k)*ekm2(k-1) )/dzh(k)
         
         gammap = 0.5 * (dzf(k)*gamma(k+1) + dzf(k+1)*gamma(k) )/dzh(k+1)
         
         svwmp = 0.5 * (dzf(k+1)*svwm(k) + dzf(k)*svwm(k+1) )/dzh(k+1)
         
         ekm2p = 0.5 * ( dzf(k+1)*ekm2(k) + dzf(k)*ekm2(k+1) )/dzh(k+1)

         ! Hier alleen de vw contributie aangepast, want voor de rest valt de - K_M*<Sij> bijdrage toch weg

        putout(i,j,k) = putout(i,j,k) &
                + (1./rhobf(k))* &
              ( epmo *gamma(k)* ( (v0(i+1,j,k)-v0(i,j,k))   *dxi &
                        +(u0(i+1,j,k)-u0(i+1,jm,k))*dyi) &
                -emmo *gamma(k)* ( (v0(i,j,k)-v0(i-1,j,k))   *dxi &
                        +(u0(i,j,k)-u0(i,jm,k))    *dyi)   ) / dx &
                + (1./rhobf(k))* &
              (ekm(i,j,k) *gamma(k)* (v0(i,jp,k)-v0(i,j,k)) &
              -ekm(i,jm,k)*gamma(k)* (v0(i,j,k)-v0(i,jm,k))  ) * 2. * dy2i &
                + (1./rhobf(k))* &
              ( eomp *gammap* ( (v0(i,j,kp)-v0(i,j,k))    /dzh(kp) &
                        +(w0(i,j,kp)-w0(i,jm,kp))  *dyi) + ekm2(kp)* svwmp &
                -eomm *gammam* ( (v0(i,j,k)-v0(i,j,km))    /dzh(k) &
                        +(w0(i,j,k)-w0(i,jm,k))    *dyi) - ekm2(k)*svwmm  ) / dzf(k)

        end do
      end do
    end do

  !     --------------------------------------------
  !     special treatment for lowest full level: k=1
  !     --------------------------------------------

    do j=2,j1
      jp = j+1
      jm = j-1
      do i=2,i1

        emmo = 0.25 * ( &
              ekm(i,j,1)+ekm(i,jm,1)+ekm(i-1,jm,1)+ekm(i-1,j,1)  )

        epmo = 0.25  * ( &
              ekm(i,j,1)+ekm(i,jm,1)+ekm(i+1,jm,1)+ekm(i+1,j,1)  )

        eomp = ( dzf(2) * ( ekm(i,j,1) + ekm(i,jm,1)  )  + &
                    dzf(1) * ( ekm(i,j,2) + ekm(i,jm,2) ) ) / &
                  ( 4.   * dzh(2) )

        vcv   = 0.5*(v0(i,j,1)+v0(i,j+1,1))+cv
        if(vcv >= 0.) then
          vpcv  = max(vcv,1.e-10)
        else
          vpcv  = min(vcv,-1.e-10)
        end if
        
         svwmp = 0.5 * (dzf(2)*svwm(1) + dzf(1)*svwm(2) )/dzh(2)
         
         ekm2p = 0.5 * ( dzf(2)*ekm2(1) + dzf(1)*ekm2(2) )/dzh(2)
         
         gammap=0.5 *(dzf(2)*gamma(1) + dzf(1)*gamma(2) )/dzh(2)


        fv    = ( 0.5*( ustar(i,j)+ustar(i,j-1) ) )**2  * &
                    vpcv/sqrt(vpcv**2  + &
                ((u0(i,j,1)+u0(i+1,j,1)+u0(i,jm,1)+u0(i+1,jm,1))/4.+cu)**2)

        putout(i,j,1) = putout(i,j,1) &
                  + (1./rhobf(1))*&
                  ( epmo *gamma(1)* ( (v0(i+1,j,1)-v0(i,j,1))   *dxi &
                            +(u0(i+1,j,1)-u0(i+1,jm,1))*dyi) &
                    -emmo * gamma(1)* ( (v0(i,j,1)-v0(i-1,j,1))   *dxi &
                            +(u0(i,j,1)-u0(i,jm,1)) *dyi)   ) / dx &
                  + (1./rhobf(1))*&
                ( ekm(i,j,1) *gamma(1)* (v0(i,jp,1)-v0(i,j,1)) &
                  -ekm(i,jm,1)*gamma(1)* (v0(i,j,1)-v0(i,jm,1))  ) * 2. * dy2i &
                  + (1./rhobf(1))*&
                ( eomp * gammap*( (v0(i,j,2)-v0(i,j,1))     /dzh(2) &
                          +(w0(i,j,2)-w0(i,jm,2))    *dyi) +ekm2(2)*svwmp  & !Victor: Hier weer een aanpassing gedaan
                  -rhobh(1)*fv  ) / dzf(1)          !Victor:  Zitten op deze plekken de grote fout!?!?!

      end do
    end do

  end subroutine diffv



  subroutine diffw(putout)
 

    use modglobal, only : i1,ih,i2,j1,jh,j2,k1,kmax,dx,dxi,dx2i,dy,dyi,dy2i,dzf,dzh
    use modfields, only : u0,v0,w0,rhobh
    implicit none

  !*****************************************************************

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: emom, eomm, eopm, epom
    integer             :: i,j,k,jm,jp,km,kp
    real                :: gammam, svwmm, ekm2m, ekm2p, suwmm !Victor: geinterpoleerde versies van de orginiele variabele
    
    
    do k=2,kmax
      kp=k+1
      km=k-1
      do j=2,j1
        jp=j+1
        jm=j-1
        do i=2,i1

          emom = ( dzf(km) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i-1,j,km) ) ) / &
                    ( 4.   * dzh(k) )

          eomm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i,jm,km) ) ) / &
                    ( 4.   * dzh(k) )

          eopm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jp,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i,jp,km) ) ) / &
                    ( 4.   * dzh(k) )

          epom = ( dzf(km) * ( ekm(i,j,k)  + ekm(i+1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i+1,j,km) ) ) / &
                    ( 4.   * dzh(k) )
                    
         gammam = 0.5 * (dzf(k-1)*gamma(k) + dzf(k)*gamma(k-1) )/dzh(k)
         
         svwmm = 0.5 * (dzf(k-1)*svwm(k) + dzf(k)*svwm(k-1) )/dzh(k)
         
         suwmm =0.5 * (dzf(k-1)*suwm(k) + dzf(k)*suwm(k-1) )/dzh(k)
         
         ekm2m = 0.5 * ( dzf(k-1)*ekm2(k) + dzf(k)*ekm2(k-1) )/dzh(k)
         
         ekm2p = 0.5 * ( dzf(k+1)*ekm2(k) + dzf(k)*ekm2(k+1) )/dzh(k+1)        
         
         
         if (k==2) then !Victor: ook hier een special treatment for the lowest level ==> k=2
         
         putout(i,j,k) = putout(i,j,k) &
                + (1./rhobh(k))* &
                  ( epom *gammam* ( (w0(i+1,j,k)-w0(i,j,k))    *dxi &
                            +(u0(i+1,j,k)-u0(i+1,j,km)) /dzh(k)) &
                    -emom *gammam *( (w0(i,j,k)-w0(i-1,j,k))    *dxi &
                            +(u0(i,j,k)-u0(i,j,km))     /dzh(k) ))/dx &
                + (1./rhobh(k))* &
                  ( eopm * gammam* ( (w0(i,jp,k)-w0(i,j,k))     *dyi &
                            +(v0(i,jp,k)-v0(i,jp,km))   /dzh(k) ) &
                    -eomm * gammam *( (w0(i,j,k)-w0(i,jm,k))     *dyi &
                            +(v0(i,j,k)-v0(i,j,km))     /dzh(k) ))/dy &
                + (1./rhobh(k))*&
                  ( ekm(i,j,k) *gamma(k)* (w0(i,j,kp)-w0(i,j,k)) /dzf(k) + ekm2p*swwm(k) &
                  -ekm(i,j,km)* (w0(i,j,k)-w0(i,j,km)) /dzf(km) ) * 2. &   !Victor: hier de afhankelijkheid van gamma en ekm2 weggehaald bij S_ww op k==2 --> daarvoor is ekm2 op z=0 nodig (== ekm2(1)), en die is niet gedefinieerd!
                                                              / dzh(k)
         
         else 

          putout(i,j,k) = putout(i,j,k) &
                + (1./rhobh(k))* &
                  ( epom *gammam* ( (w0(i+1,j,k)-w0(i,j,k))    *dxi &
                            +(u0(i+1,j,k)-u0(i+1,j,km)) /dzh(k)) &
                    -emom *gammam *( (w0(i,j,k)-w0(i-1,j,k))    *dxi &
                            +(u0(i,j,k)-u0(i,j,km))     /dzh(k) ))/dx &
                + (1./rhobh(k))* &
                  ( eopm * gammam* ( (w0(i,jp,k)-w0(i,j,k))     *dyi &
                            +(v0(i,jp,k)-v0(i,jp,km))   /dzh(k) ) &
                    -eomm * gammam *( (w0(i,j,k)-w0(i,jm,k))     *dyi &
                            +(v0(i,j,k)-v0(i,j,km))     /dzh(k) ))/dy &
                + (1./rhobh(k))*&
                  ( ekm(i,j,k) *gamma(k)* (w0(i,j,kp)-w0(i,j,k)) /dzf(k) + ekm2p*swwm(k) &
                  -ekm(i,j,km)*gamma(km)* (w0(i,j,k)-w0(i,j,km)) /dzf(km) -ekm2m*swwm(km) ) * 2. &
                                                              / dzh(k)    
         end if

        end do
      end do
    end do

  end subroutine diffw
  
  
  subroutine closure2 !Victor
    use modglobal, only : i1,ih,i2,j1,jh,j2,k1,kmax,dx,dxi,dx2i,dy,dyi,dy2i,dzf,dzh, rslabs, fkar, cu, cv, zh, zf
    use modfields, only : u0,v0,w0,rhobh, ustarsullivan, oblsullivan
    use modmpi,    only : excjs, myid, nprocs, comm3d, mpierr, my_real, mpi_sum, slabsum
    use modsubgriddata, only : zisullivan
    implicit none
    
    real   ekmgammal
    real  ekmgamma
    real  uwresl
    real uwres
    real  vwresl
    real vwres
    real ustar2l
    real phimzf2l, phimzf
    integer :: i,j,k,kp,km,jp,jm
    
    
    ! Per definitie is deze ekm2 bepaald op de cell faces, 
    
    do j=2,j1
      do i=2,i1
        uwr(i,j)     = (w0(i,j,2)+w0(i-1,j,2)) &
                          *((u0(i,j,1)+cu)*dzf(2)+(u0(i,j,2)+cu)*dzf(1))/(4*dzh(2))
                          
        vwr(i,j)     = (w0(i,j,2)+w0(i,j-1,2)) &
                        *((v0(i,j,1)+cv)*dzf(2)+(v0(i,j,2)+cv)*dzf(1))/(4*dzh(2))
                        
        uwresl = uwresl + uwr(i,j)
        vwresl = vwresl + vwr(i,j)
        ustar2l= ustar2l + ustarsullivan(i,j)
        
        !Victor: onderstaand stuk komt uit modsurface, alleen zf(1) --> zh(2) om op het eerste w-punt te komen
        if (oblsullivan(i,j) < 0.) then
            phimzf = (1.-16.*zh(2)/oblsullivan(i,j))**(-0.25)
            !phimzf = (1. + 3.6 * (-zf(1)/obl(i,j))**(2./3.))**(-0.5)
            !phihzf = (1.-16.*zf(1)/obl(i,j))**(-0.50)
            !phihzf = (1. + 7.9 * (-zf(1)/obl(i,j))**(2./3.))**(-0.5)
          elseif (oblsullivan(i,j) > 0.) then
            phimzf = (1.+5.*zh(2)/oblsullivan(i,j))
            !phihzf = (1.+5.*zf(1)/obl(i,j))
          else
            phimzf = 1.
            !phihzf = 1.
          endif
        
        phimzf2l = phimzf2l +phimzf
        
        ekmgammal= ekmgammal + (ekm(i,j,1)+ekm(i,j,2))*(gamma(1)+gamma(2))/4 !Victor: gamma en ekm zijn beide op het full level gedefinieerd, maar ik heb ze op het half level nodig
        
         
       end do
    end do
    
    call MPI_ALLREDUCE(ekmgammal, ekmgamma, 1,    MY_REAL, MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(uwresl, uwres, 1,    MY_REAL, MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(vwresl, vwres, 1,    MY_REAL, MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(ustar2l, ustar2, 1,    MY_REAL, MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(phimzf2l, phimzf2, 1,    MY_REAL, MPI_SUM, comm3d,mpierr)
    
    uwres    = uwres    /rslabs
    vwres    = vwres    /rslabs
    ekmgamma = ekmgamma / rslabs
    phimzf2  = phimzf2 / rslabs
    ustar2   = ustar2 / rslabs
    
      
    do k=1,kmax
          
       if (k == 1) then 
        ekm2(1)=0 !Victor: Op k=1 (ground) is ekm2 niet gedefinieerd, gaat pas een rol spelen op 2 en hoger!!
        else if (k==2) then
        ekm2(2)=ustar2*fkar*zh(2)/phimzf2-ekmgamma-fkar*zh(2)/(ustar2*phimzf2)*sqrt((uwres)**2+(vwres)**2)
        else
        ekm2(k)=ekm2(2)*fkar*zh(2)/(ustar2*phimzf2)*(dzf(k)*sslab(k+1)+dzf(k+1)*sslab(k))/dzf(k)
        end if      
      
      !Victor Rekening houden met de maximale hoogte van de boundary layer
      if (zf(k).gt.zisullivan) then
      ekm2(k)=0
      gamma(k)=1
      end if
      
      
    end do
    
  end subroutine closure2

end module
