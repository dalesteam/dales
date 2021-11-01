!> \file modradiation.f90
!!  Calculates the radiative sources

!>
!!  Calculates the radiative sources
!>
!!  \author Stephan de Roode,TU Delft
!!  \author Thijs Heus, MPI-M
!!  \todo Documentation
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

module modradiation
use modraddata
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initradiation
    use modglobal,    only : i1,ih,j1,jh,k1,nsv,ih,jh,btime,tres,dt_lim,ifnamopt,fname_options,checknamelisterror
    use modmpi,       only : myid,comm3d,mpi_logical,mpi_integer,D_MPI_BCAST
    implicit none

    integer :: ierr

    namelist/NAMDE/ &
      SSA,iDE,laero

    namelist/NAMRADIATION/ &
      lCnstZenith, cnstZenith, lCnstAlbedo, ioverlap, &
      inflglw, iceflglw, liqflglw, inflgsw, iceflgsw, liqflgsw, &
      ocean, usero3, co2factor, doperpetual, doseasons, iyear

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMDE,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMDE')
      write(6 ,NAMDE)

      rewind(ifnamopt)

      read (ifnamopt,NAMRADIATION,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMRADIATION')
      write(6 ,NAMRADIATION)

      close(ifnamopt)
    end if

    call D_MPI_BCAST(SSA,1,0,comm3d,ierr)
    call D_MPI_BCAST(iDE,1,0,comm3d,ierr)
    call D_MPI_BCAST(laero,1,0,comm3d,ierr)

    call D_MPI_BCAST(lCnstZenith,1,0,comm3d,ierr)
    call D_MPI_BCAST(cnstZenith, 1, 0,comm3d,ierr)
    call D_MPI_BCAST(lCnstAlbedo,1,0,comm3d,ierr)
    call D_MPI_BCAST(ioverlap,   1,0,comm3d,ierr)
    call D_MPI_BCAST(inflglw,    1,0,comm3d,ierr)
    call D_MPI_BCAST(iceflglw,   1,0,comm3d,ierr)
    call D_MPI_BCAST(liqflglw,   1,0,comm3d,ierr)
    call D_MPI_BCAST(inflgsw,    1,0,comm3d,ierr)
    call D_MPI_BCAST(iceflgsw,   1,0,comm3d,ierr)
    call D_MPI_BCAST(liqflgsw,   1,0,comm3d,ierr)
    call D_MPI_BCAST(ocean,      1,0,comm3d,ierr)
    call D_MPI_BCAST(usero3,     1,0,comm3d,ierr)
    call D_MPI_BCAST(co2factor,  1, 0,comm3d,ierr)
    call D_MPI_BCAST(doperpetual,1,0,comm3d,ierr)
    call D_MPI_BCAST(doseasons,  1,0,comm3d,ierr)
    call D_MPI_BCAST(iyear,      1,0,comm3d,ierr)

    allocate(thlprad   (2-ih:i1+ih,2-jh:j1+jh,k1) )
    allocate(swd       (2-ih:i1+ih,2-jh:j1+jh,k1) )
    allocate(swu       (2-ih:i1+ih,2-jh:j1+jh,k1) )
    allocate(lwd       (2-ih:i1+ih,2-jh:j1+jh,k1) )
    allocate(lwu       (2-ih:i1+ih,2-jh:j1+jh,k1) )

    allocate(swdca     (2-ih:i1+ih,2-jh:j1+jh,k1) )
    allocate(swuca     (2-ih:i1+ih,2-jh:j1+jh,k1) )
    allocate(lwdca     (2-ih:i1+ih,2-jh:j1+jh,k1) )
    allocate(lwuca     (2-ih:i1+ih,2-jh:j1+jh,k1) )

    allocate(swdir     (2-ih:i1+ih,2-jh:j1+jh,k1) )
    allocate(swdif     (2-ih:i1+ih,2-jh:j1+jh,k1) )
    allocate(lwc       (2-ih:i1+ih,2-jh:j1+jh,k1) )

    allocate(SW_up_TOA (2-ih:i1+ih,2-jh:j1+jh)    )
    allocate(SW_dn_TOA (2-ih:i1+ih,2-jh:j1+jh)    )
    allocate(LW_up_TOA (2-ih:i1+ih,2-jh:j1+jh)    )
    allocate(LW_dn_TOA (2-ih:i1+ih,2-jh:j1+jh)    )

    allocate(SW_up_ca_TOA(2-ih:i1+ih,2-jh:j1+jh)  )
    allocate(SW_dn_ca_TOA(2-ih:i1+ih,2-jh:j1+jh)  )
    allocate(LW_up_ca_TOA(2-ih:i1+ih,2-jh:j1+jh)  )
    allocate(LW_dn_ca_TOA(2-ih:i1+ih,2-jh:j1+jh)  )

    thlprad = 0.

    swd = 0.
    swu = 0.
    lwd = 0.
    lwu = 0.

    swdca = 0.
    swuca = 0.
    lwdca = 0.
    lwuca = 0.

    swdir = 0.
    swdif = 0.
    lwc   = 0.

    SW_up_TOA=0;SW_dn_TOA=0;LW_up_TOA=0;LW_dn_TOA=0
    SW_up_ca_TOA = 0. ;SW_dn_ca_TOA=0    ;LW_up_ca_TOA=0    ;LW_dn_ca_TOA=0

    if (irad/=-1) then
      if (myid==0) write (*,*) 'WARNING: The use of irad is deprecated. Please use the iradiation switch'
      select case (irad)
      case (0)
        iradiation = 0
      case (1)
        iradiation = 2
        rad_ls     = .true.
        rad_longw  = .false.
        rad_shortw = .false.
        rad_smoke  = .false.
      case (2)
        iradiation = 2
        rad_ls     = .false.
        rad_longw  = .true.
        rad_shortw = .false.
        rad_smoke  = .false.
      case (3)
        iradiation = 1
      case (4)
        iradiation = 2
        rad_ls     = .false.
        rad_longw  = .true.
        rad_shortw = .true.
        rad_smoke  = .false.
      case (10)
        iradiation = 2
        rad_ls     = .false.
        rad_longw  = .false.
        rad_shortw = .false.
        rad_smoke  = .true.
      end select
    end if
    if(iradiation==0 .or. iradiation==10) then
      rad_shortw = .false.
      rad_longw  = .false.
      rad_smoke  = .false.
    end if

    if (iradiation == 0) return

    itimerad = floor(timerad/tres)
    !setting tnext is done in modstartup, after btime has been set
    !tnext = itimerad+btime
    !dt_lim = min(dt_lim,tnext)

    if (rad_smoke.and.isvsmoke>nsv) then
      if (rad_shortw) then
         stop 'you want to compute solar radiative transfer through a smoke cloud'
      endif
      stop 'Smoke radiation with wrong (non-existent?) scalar field'
    endif

  end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine radiation
    use modmpi, only: myid, MPI_Wtime
    use modglobal, only : timee, dt_lim,rk3step
    use modfields, only : thlp
    use moduser,   only : rad_user
    use modradfull,only : radfull
    use modradrrtmg, only : radrrtmg
    implicit none
    real wtime

    if(timee<tnext .and. rk3step==3) then
      dt_lim = min(dt_lim,tnext-timee)
    end if
    if((itimerad==0 .or. timee>=tnext) .and. rk3step==1) then
      tnext = tnext+itimerad

      wtime = MPI_Wtime()
      thlprad = 0.0
      select case (iradiation)
          case (irad_none)
          case (irad_full)
            call radfull
          case (irad_par)
             if(rad_longw.or.rad_shortw) then
              call radpar
            endif
          case (irad_lsm)
            call radlsm
          case (irad_rrtmg)
            call radrrtmg
          case (irad_user)
! EWB: the if statement should came first because moduser uses a radpar variable
            if(rad_longw.or.rad_shortw) then
              call radpar
            endif

            call rad_user

      end select
      if (rad_ls) then
        call radprof
      endif

      if (myid == 0 .and. iradiation /= irad_none) then
         wtime = MPI_Wtime() - wtime
         write (*,*) 'Radiation', timee, 'Time spent:', wtime, 's'
      end if
    end if

    thlp = thlp + thlprad
  end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine exitradiation
    implicit none
    deallocate(thlprad,swd,swdir,swdif,swu,lwd,lwu,swdca,swuca,lwdca,lwuca,lwc)
    deallocate(SW_up_TOA, SW_dn_TOA,LW_up_TOA,LW_dn_TOA, &
               SW_up_ca_TOA,SW_dn_ca_TOA,LW_up_ca_TOA,LW_dn_ca_TOA)

  end subroutine exitradiation


!> calculates tendency due to parameterized radiation
subroutine radpar

  use modglobal,    only : i1,j1,kmax, k1,ih,jh,dzf,cp,xtime,rtimee,xday,xlat,xlon
  use modfields,    only : ql0, sv0, rhof,exnf
  use modsurfdata,  only : tauField, lrsAgs
  implicit none
  real, allocatable :: lwpt(:),lwpb(:)
  real, allocatable :: tau(:)
  real, allocatable :: absorber(:,:,:)  !  liquid water content or smoke

  real thlpld,thlplu,thlpsw
  real tauc
  integer :: i=0, j=0, k

  real :: rho_l= 1000.  !liquid water density (kg/m3)

  allocate(lwpt(k1),lwpb(k1))
  allocate(tau(k1))
  allocate(absorber(2-ih:i1+ih,2-jh:j1+jh,k1))
  absorber = 0.0
  tau  = 0.0
  lwpt = 0.0
  lwpb = 0.0


! MPI
! initialise local variables

   if (rad_longw) then
     if (rad_smoke) then
        absorber(2-ih:i1+ih,2-jh:j1+jh,1:k1) = sv0(2-ih:i1+ih,2-jh:j1+jh,1:k1,isvsmoke)
     else
        absorber(2-ih:i1+ih,2-jh:j1+jh,1:k1) = ql0(2-ih:i1+ih,2-jh:j1+jh,1:k1)
     endif

     do j=2,j1
     do i=2,i1
       lwpt = 0.
       lwpb = 0.

! **   Downward LWP

       do k=kmax,1,-1
        lwpt(k) = lwpt(k+1) + rhof(k)*absorber(i,j,k)*dzf(k)
       end do

! **   Upward LWP

       do k=2,k1
         lwpb(k) = lwpb(k-1) + rhof(k)*absorber(i,j,k-1)*dzf(k-1)
       end do

       do k=1,k1
        lwd(i,j,k) = dlwtop*exp(-rka*lwpt(k))
        lwu(i,j,k) = dlwbot*exp(-rka*lwpb(k))
      end do

       do k=1,kmax
         thlpld         = -(lwd(i,j,k+1)-lwd(i,j,k))/(rhof(k)*cp*exnf(k)*dzf(k))
         thlplu         = -(lwu(i,j,k+1)-lwu(i,j,k))/(rhof(k)*cp*exnf(k)*dzf(k))
         thlprad(i,j,k) =   thlprad(i,j,k) + thlpld+thlplu
       end do

    end do
    end do  ! end i,j loop

  endif  !end longwave loop

!----------------------------------------------------------------------

  if (rad_shortw) then

    !compute solar zenith angle


    swd = 0.0
    mu=zenith(xtime*3600 + rtimee,xday,xlat,xlon)
    do j=2,j1
    do i=2,i1

      if (mu > 0.035) then  !factor 0.035 needed for security
        tauc = 0.           ! column-integrated tau cloud
        if (laero .or. lcloudshading) then ! not sure if I have to define the use of lcldoushading before
          do k = 1,kmax        
            tau(k) = 0.      ! tau laagje dz
            if(laero) then ! there are aerosols
              tau(k) = sv0(i,j,k,iDE)
            else if (lcloudshading) then ! there are clouds
              if (ql0(i,j,k) > 1e-5)  tau(k)=1.5*ql0(i,j,k)*rhof(k)*dzf(k)/reff/rho_l
            end if
            tauc=tauc+tau(k)
          end do
        endif
        if (lrsAgs) tauField(i,j) = tauc
        call sunray(tau,tauc,i,j)
      end if

      do k=1,kmax
        thlpsw          = ((swd(i,j,k+1)-swu(i,j,k+1))-(swd(i,j,k)-swu(i,j,k)))/(rhof(k)*cp*exnf(k)*dzf(k))
        thlprad(i,j,k)  = thlprad(i,j,k) + thlpsw

      end do

    end do
    end do

  endif  !end shortwave loop

  deallocate(lwpt,lwpb,tau,absorber)



  return
  end subroutine radpar

!!   the sunray model is described by fouquart and  bonnel
!!   (1980, contr. atmos. phys.).

  subroutine sunray(tau,tauc,i,j)

  use modglobal, only :  k1,boltz
  use modsurfdata,  only : albedo,tskin
  use modfields,   only : thl0

  implicit none

  real, intent(inout), dimension (k1) :: tau
  integer, intent(in) :: i,j
  real,allocatable, dimension (:) :: taude
  real gcde,tauc,taucde &
           ,taupath,t1,t2,t3,c1,c2 &
           ,omega,omegade,ff,x1,x2,x3,rk,mu2,rp,alpha,beta,rtt &
           ,exmu0,expk,exmk,xp23p,xm23p,ap23b,Irr0,Irr1
  integer k
  allocate(taude(k1))

  taucde = 0.         ! tau' cloud
  taupath = 0.

  do k=1,k1
     taude(k) =0.     ! tau' laagje dz
  end do

  if(laero) then
    omega= SSA
    gc= 0.64 ! typical for aerosols
  else
    omega=1.-1.e-3*(0.9+2.75*(mu+1.)*exp(-0.09*tauc)) !fouquart (for clouds)
  endif


! the equations for the delta-eddington approximation are equal to those
! for the eddington approximation with transformed parameters g, omega
! and tau (joseph, wiscomb and weinman, 1976, j.a.s.).
! parameternames: x -> xde (delta-eddington)

  ff=gc*gc
  gcde=gc/(1.+gc)
  taucde=(1.0-omega*ff)*tauc

  do k =1,k1
    taude(k)=(1.e0-omega*ff)*tau(k)
  end do

  omegade=(1.0-ff)*omega/(1.e0-omega*ff)

! the solution of the eddington equations are given by shettle and weinman
! (1970, j.a.s.).

  x1=1.0-omegade*gcde
  x2=1.0-omegade
  rk=sqrt(3.0*x2*x1)
  mu2=mu*mu
  x3=4.0*(1.0-rk*rk*mu2)
  rp=sqrt(3.0*x2/x1)
  alpha=3.e0*omegade*mu2*(1.0+gcde*x2)/x3
  beta=3.*omegade*mu*(1.0+3.0*gcde*mu2*x2)/x3

  rtt=2.0/3.0
  exmu0= exp(-taucde/mu)
  expk=  exp(rk*taucde)
  exmk=1.0/expk
  xp23p=1.0+rtt*rp
  xm23p=1.0-rtt*rp
  ap23b=alpha+rtt*beta

  t1=1-albedo(i,j)-rtt*(1.+albedo(i,j))*rp
  t2=1-albedo(i,j)+rtt*(1.+albedo(i,j))*rp
  t3=(1-albedo(i,j))*alpha-rtt*(1+albedo(i,j))*beta+albedo(i,j)*mu
  c2=(xp23p*t3*exmu0-t1*ap23b*exmk)/(xp23p*t2*expk-xm23p*t1*exmk)
  c1=(ap23b-c2*xm23p)/xp23p

  do k = k1,1,-1
      taupath = taupath + taude(k)

      Irr0    = sw0*(    c1*exp(-rk*taupath) + c2*exp(rk*taupath) - alpha*exp(-taupath/mu)) ! Shettle & Weinmann JAS 1976
      Irr1    = sw0*(rp*(c1*exp(-rk*taupath) - c2*exp(rk*taupath))- beta *exp(-taupath/mu)) ! Shettle & Weinmann JAS 1976

      swd(i,j,k) = (Irr0 + (2./3.)*Irr1) + mu*sw0*exp(-taupath/mu) ! difuse down + direct down
      swu(i,j,k) = (Irr0 - (2./3.)*Irr1)                           ! diffuse up (lambertian)
      swdir(i,j,k) = mu*sw0*exp(-taupath/mu)
      swdif(i,j,k) = (Irr0 + (2./3.)*Irr1)
      lwd(i,j,1) =  0.8 * boltz * thl0(i,j,1) ** 4.
      lwu(i,j,1) =  1.0 * boltz * tskin(i,j) ** 4.
  end do

  deallocate(taude)

  return
  end subroutine sunray


!***********************************************************************
!***  In this subroutine the a precribed radiative tendency     ********
!***  is taken into account.                                    ********
!***********************************************************************

  subroutine radprof
  use modglobal,    only : i1,j1,kmax
  use modfields,    only : thlpcar
  implicit none
  integer k

  do k=1,kmax
    thlprad(2:i1,2:j1,k) = thlprad(2:i1,2:j1,k) + thlpcar(k)
  end do

  return
  end subroutine radprof


  subroutine radlsm
    use modsurfdata, only : albedo, tskin
    use modglobal,   only : i1, j1, rtimee, xtime, xday, xlat, xlon, boltz, dzf
    use modfields,   only : thl0, ql0, rhof
    implicit none
    integer        :: i,j
    real           :: Tr, sinlea, qlint, tau
    real,parameter :: S0 = 1376.
    real           :: rhow  = 1000.  ! density of water
    real           :: tauc  = 5.     ! Optically thin clouds (critical value)

    sinlea = zenith(xtime*3600 + rtimee, xday, xlat, xlon)

    Tr  = (0.6 + 0.2 * sinlea)

    do j=2,j1
      do i=2,i1
        qlint = sum(ql0(i,j,:) * rhof(:) * dzf(:))
        tau   = (3./2.)*(qlint/(rhow*reff))
        if(lcloudshading .and. (tau >= tauc)) then ! 'dense' cloud and cloud shading is active
          swd(i,j,1) = - S0 * Tr * sinlea * (5. - exp(-tau))/(4.+ 3.*tau*0.14)
        else
          swd(i,j,1) = - S0 * Tr * sinlea
        endif
        swu(i,j,1) = - albedo(i,j) * swd(i,j,1)
        lwd(i,j,1) = - 0.8 * boltz * thl0(i,j,1) ** 4.
        lwu(i,j,1) = boltz * tskin(i,j) ** 4.
      end do
    end do

    !write(6,*) "CvHrad", swd(2,2,1), swu(2,2,1), lwd(2,2,1), lwu(2,2,1), tskin(2,2)

  end subroutine radlsm


end module
