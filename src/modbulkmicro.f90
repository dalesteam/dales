!> \file modbulkmicro.f90

!>
!!  Bulk microphysics.
!>
!! Calculates bulk microphysics using a two moment scheme.
!! \see  Seifert and Beheng (Atm. Res., 2001)
!! \see  Seifert and Beheng (Met Atm Phys, 2006)
!! \see  Stevens and Seifert (J. Meteorol. Soc. Japan, 2008)  (rain sedim, mur param)
!! \see  Seifert (J. Atm Sc., 2008) (rain evap)
!! \see  Khairoutdinov and Kogan (2000) (drizzle param : auto, accr, sedim, evap)
!!  \author Olivier Geoffroy, K.N.M.I.
!!  \author Margreet van Zanten, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Marco de Bruine, UU/IMAU 
!!  \par Revision list
!! \todo documentation
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

! MdB TODO: Clean up: calculation of aer_tend only after mass conservation check

module modbulkmicro

!   Amount of liquid water is splitted into cloud water and precipitable
!   water (as such it is a two moment scheme).
!
!   same rhof value used for some diagnostics calculation (in modbulkmicrostat, modtimestat)
!
!   Cond. sampled timeav averaged profiles are weighted with fraction of condition,
!   similarly as is done in sampling.f90
!
!   bulkmicro is called from *modmicrophysics*
!
!*********************************************************************
  use modmicrodata

  implicit none
  real :: gamma25
  real :: gamma3
  real :: gamma35

  integer, parameter :: ncld = 10,     & ! In-cloud scav. lookuptable  for cloud
                        nrai = 5,      & ! Below-cloud scav. lookuptable 
                        naer_blc = 100, &
                        naer_inc = 60

  real, dimension(naer_inc) :: aerrad,     logaerrad
  real, dimension(naer_blc) :: aerrad_blc, logaerrad_blc
  real, dimension(ncld)     :: cldrad,     logcldrad
  real, dimension(nrai)     :: rainrate,   lograinrate

  real, dimension(ncld,naer_inc) :: incmass, incnumb
  real, dimension(nrai,naer_blc) :: blcmass, blcnumb
 
  contains

!> Initializes and allocates the arrays
  subroutine initbulkmicro
    use modglobal, only : ih,i1,jh,j1,k1
    use modmpi,    only : myid, comm3d, mpierr, my_real

    implicit none

    character(3) :: cldrad_char
    integer :: icld, iaer

    allocate( Nc       (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,Nr       (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,Nrp      (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,qltot    (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,qr       (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,qrp      (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,nuc      (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,thlpmcr  (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,qtpmcr   (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,sedc     (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,sedcn    (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,sed_qr   (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,sed_Nr   (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,Dvc      (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,xc       (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,Dvr      (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,xr       (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,mur      (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,lbdr     (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,au       (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,phi      (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,tau      (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,ac       (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,sc       (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,br       (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,evap     (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,Nevap    (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,qr_spl   (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,Nr_spl   (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,wfall_qr (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,wfall_Nr (2-ih:i1+ih,2-jh:j1+jh,k1))

    allocate( precep   (2-ih:i1+ih,2-jh:j1+jh,k1))

    allocate( qrmask   (2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,qcmask   (2-ih:i1+ih,2-jh:j1+jh,k1))

    allocate( sedcm    (2-ih:i1+ih,2-jh:j1+jh,k1,naer))      

    !Local aerosol fields
    allocate( aer_conc (2-ih:i1+ih,2-jh:j1+jh,k1,naer) & ! Concentration
             ,aer_tend (2-ih:i1+ih,2-jh:j1+jh,k1,naer) & ! Aggregated tendency
             ,aer_acti (2-ih:i1+ih,2-jh:j1+jh,k1,naer) & ! A<->C activation
             ,aer_scvc (2-ih:i1+ih,2-jh:j1+jh,k1,naer) & ! A<->C in-cloud scav 
             ,aer_evpc (2-ih:i1+ih,2-jh:j1+jh,k1,naer) & ! A<->C cloud evap
             ,aer_scvr (2-ih:i1+ih,2-jh:j1+jh,k1,naer) & ! A<->R in-rain scav
             ,aer_evpr (2-ih:i1+ih,2-jh:j1+jh,k1,naer) & ! A<->R rain evap
             ,aer_auto (2-ih:i1+ih,2-jh:j1+jh,k1,naer) & ! C<->R autoconversion
             ,aer_accr (2-ih:i1+ih,2-jh:j1+jh,k1,naer) & ! C<->R accretion
             ,aer_slfc (2-ih:i1+ih,2-jh:j1+jh,k1,naer) & ! C<->C clouddrop selfcol
             ,aer_slfr (2-ih:i1+ih,2-jh:j1+jh,k1,naer) & ! R<->R raindrop selfcol
             ,aer_sedr (2-ih:i1+ih,2-jh:j1+jh,k1,naer))  ! R<->R raindrop sedimen

  gamma25 = lacz_gamma(2.5)
  gamma3  = 2.
  gamma35 = lacz_gamma(3.5)
  
  !*********************************************************************
  ! Reading lookuptable for scavenging here
  ! Could also be placed somewhere else? What is the proper place?
  !*********************************************************************
  if (myid == 0 )  then 
      ! =============================================================
      ! IN-CLOUD SCAVENGING    
      ! =============================================================
      cldrad = (/5., 10., 15., 20., 25., 30., 35., 40., 45., 50./)

      open(100,file='rbar_aerosol')
        read(100,*) aerrad
      close(100)

      logaerrad = log(aerrad)
      logcldrad = log(cldrad)

      ! =============================================================

      do icld=1, ncld

         write(cldrad_char,'(i3.3)') int(cldrad(icld))

         open(101,file='CDNC1_cloudscavenging_'//cldrad_char//'.dat')
         open(102,file='CDNC1_ncloudscavenging_'//cldrad_char//'.dat')

         do iaer=1, naer_inc
            read(101,*) incmass(icld,iaer)
            read(102,*) incnumb(icld,iaer)
         enddo

         close(101)
         close(102)
      enddo
      
      ! =============================================================
      ! BELOW CLOUD SCAVENGING    
      ! =============================================================
      rainrate = (/.01, .1, 1., 10., 100./)

      open(100,file='rbar_aerosol_belowcloud.dat')
        read(100,*) aerrad_blc
      close(100)      
 
      logaerrad_blc = log(aerrad_blc)
      lograinrate   = log(rainrate)

      ! =============================================================

      open (101,file='belowcloud_m.dat')
      read (101,*) blcmass
      close(101)

      open (102,file='belowcloud_n.dat')
      read (102,*) blcnumb
      close(102)

  end if

  call MPI_BCAST(logaerrad,     naer_inc,      MY_REAL,0,comm3d,mpierr)
  call MPI_BCAST(logaerrad_blc, naer_blc,      MY_REAL,0,comm3d,mpierr)
  call MPI_BCAST(logcldrad,     ncld,          MY_REAL,0,comm3d,mpierr)
  call MPI_BCAST(lograinrate,   nrai,          MY_REAL,0,comm3d,mpierr)
  call MPI_BCAST(incmass,       naer_inc*ncld, MY_REAL,0,comm3d,mpierr)
  call MPI_BCAST(incnumb,       naer_inc*ncld, MY_REAL,0,comm3d,mpierr)
  call MPI_BCAST(blcmass,       naer_blc*nrai, MY_REAL,0,comm3d,mpierr)
  call MPI_BCAST(blcnumb,       naer_blc*nrai, MY_REAL,0,comm3d,mpierr)

  end subroutine initbulkmicro

!> Cleaning up after the run
  subroutine exitbulkmicro
  !*********************************************************************
  ! subroutine exitbulkmicro
  !*********************************************************************
    implicit none

    deallocate(Nc,Nr,Nrp,qltot,qr,qrp,nuc,thlpmcr,qtpmcr)

    deallocate(sedc,sedcn,sedcm,sed_qr,sed_Nr,Dvc,xc,Dvr,xr,mur,lbdr, &
               au,phi,tau,ac,sc,br,evap,Nevap,qr_spl,Nr_spl,wfall_qr,wfall_Nr)

    deallocate(qcmask, qrmask)

    deallocate(precep)

    deallocate(aer_conc, aer_tend)

    deallocate(aer_acti,aer_scvc,aer_evpc, & ! A<->C
                        aer_scvr,aer_evpr, & ! A<->R
               aer_auto,aer_accr,          & ! C<->R
               aer_slfc,aer_slfr,aer_sedr  ) ! C<->C / R<->R

  end subroutine exitbulkmicro

!> Calculates the microphysical source term.
  subroutine bulkmicro
    use modglobal,        only : i1,j1,k1,rdt,rk3step,timee,rlv,cp 
    use modfields,        only : sv0,svm,svp,qtp,thlp,ql0,exnf,rhof
    use modbulkmicrostat, only : bulkmicrotend, bulkaertend
    use modmpi,           only : myid
    implicit none
    integer :: i,j,k, taer, iaer, ispec, imod
    real :: qrtest,nrtest,nctest,mctest,aertest,fix, &
            sinktotal, f_acti, f_scvc, f_scvr, &
            f_evpc,  f_auto, f_accr, f_slfc, &
            f_acs, f_cos     

    Nrp     = 0.0 
    qrp     = 0.0
    thlpmcr = 0.0
    qtpmcr  = 0.0

    aer_tend = 0.0     
  
    aer_acti = 0.0
    aer_scvc = 0.0
    aer_evpc = 0.0
    aer_scvr = 0.0
    aer_evpr = 0.0
    aer_auto = 0.0
    aer_accr = 0.0           
    aer_slfc = 0.0
    aer_slfr = 0.0
    aer_sedr = 0.0
    
    do j=2,j1
    do i=2,i1
    do k=1,k1

       qr(i,j,k) = sv0(i,j,k,iqr) ! kg tracer / kg air
      
       ! For clarity use duplicate arrays for Nc and Nr
       Nc(i,j,k) = sv0(i,j,k,inc+iaer_offset)!For now assume tracer is in #/m3 anyway *rhof(k) ! Cloud drops [ #/m3 ]
       Nr(i,j,k) = sv0(i,j,k,inr+iaer_offset)!Idem *rhof(k) ! Rain drops  [ #/m3 ]
    
       do iaer=1,naer
          aer_conc(i,j,k,iaer) = sv0(i,j,k,iaer+iaer_offset)! Idem *rhof(k) ! kg/m3 OR #/m3
       enddo  

    enddo
    enddo
    enddo

    delt = rdt/ (4. - dble(rk3step))

    if ( timee .eq. 0. .and. rk3step .eq. 1 .and. myid .eq. 0) then
      write(*,*) 'l_lognormal',l_lognormal
      write(*,*) 'rhof(1)', rhof(1),' rhof(10)', rhof(10)
      write(*,*) 'l_mur_cst',l_mur_cst,' mur_cst',mur_cst
      write(*,*) 'nuc = param'
    endif

    !*********************************************************************
    ! remove neg. values of Nr and qr
    !*********************************************************************

    if (l_rain) then    
       if (sum(qr, qr<0.) > 0.000001*sum(qr)) then
         write(*,*)'amount of neg. qr and Nr thrown away is too high  ',timee,' sec'
       end if
       if (sum(Nr, Nr<0.) > 0.000001*sum(Nr)) then
          write(*,*)'amount of neg. qr and Nr thrown away is too high  ',timee,' sec'
       end if
   
       do j=2,j1
       do i=2,i1
       do k=1,k1
   
          if (qr(i,j,k) < 0.)  then
            qr(i,j,k) = 0.
          endif
          if (Nr(i,j,k) < 0.)  then
            Nr(i,j,k) = 0.
          endif
           
          if (qr(i,j,k) > qrmin .and. Nr(i,j,k) > 0)  then
             qrmask (i,j,k) = .true.
          else
             qrmask (i,j,k) = .false.
          endif

       enddo
       enddo
       enddo
    end if !l_rain 
 
  !*********************************************************************
  ! calculate qltot and remove negative values of Nc
  !*********************************************************************

    do j=2,j1
    do i=2,i1
    do k=1,k1

       ql0   (i,j,k) = ql0 (i,j,k)
       qltot (i,j,k) = ql0  (i,j,k) + qr (i,j,k)
       if (ql0(i,j,k) > qcmin)  then
          qcmask (i,j,k) = .true.
       else
          qcmask (i,j,k) = .false.
       end if

       if (Nc(i,j,k) < 0.)  then
         Nc(i,j,k) = 0.
       endif

  !*********************************************************************
  ! Remove negative values in the aerosol fields
  !*********************************************************************

       do iaer=1,naer
         if (aer_conc(i,j,k,iaer) < 0. ) then
            aer_conc(i,j,k,iaer) = 0.
         end if
       end do

    enddo
    enddo
    enddo

  !*********************************************************************
  ! calculate Rain DSD integral properties & parameters xr, Dvr, lbdr, mur
  !*********************************************************************
    if (l_rain) then

      xr   (2:i1,2:j1,1:k1) = 0.
      Dvr  (2:i1,2:j1,1:k1) = 0.
      mur  (2:i1,2:j1,1:k1) = 30.
      lbdr (2:i1,2:j1,1:k1) = 0.

      if (l_sb) then

        do j=2,j1
        do i=2,i1
        do k=1,k1
           if (qrmask(i,j,k)) then
             xr (i,j,k) = rhof(k)*qr(i,j,k)/(Nr(i,j,k)+eps0) ! JvdD Added eps0 to avoid floating point exception
             xr (i,j,k) = min(max(xr(i,j,k),xrmin),xrmax) ! to ensure xr is within borders
             Dvr(i,j,k) = (xr(i,j,k)/pirhow)**(1./3.)
           endif
        enddo
        enddo
        enddo


        if (l_mur_cst) then
        ! mur = cst
          do j=2,j1
          do i=2,i1
          do k=1,k1
            if (qrmask(i,j,k)) then
               mur(i,j,k) = mur_cst
               lbdr(i,j,k) = ((mur(i,j,k)+3.)*(mur(i,j,k)+2.)*(mur(i,j,k)+1.))**(1./3.)/Dvr(i,j,k)
            endif
          enddo
          enddo
          enddo
        else
        ! mur = f(Dv)
          do j=2,j1
          do i=2,i1
          do k=1,k1
            if (qrmask(i,j,k)) then
!
!             mur(2:i1,2:j1,1:k1) = 10. * (1+tanh(1200.*(Dvr(2:i1,2:j1,1:k1)-0.0014)))
!             Stevens & Seifert (2008) param
!
              mur(i,j,k) = min(30.,- 1. + 0.008/ (qr(i,j,k)*rhof(k))**0.6)  ! G09b
              lbdr(i,j,k) = ((mur(i,j,k)+3.)*(mur(i,j,k)+2.)*(mur(i,j,k)+1.))**(1./3.)/Dvr(i,j,k)
            endif
          enddo
          enddo
          enddo

        endif

      else
         do j=2,j1
         do i=2,i1
         do k=1,k1
            if (qrmask(i,j,k)) then
              xr (i,j,k) = rhof(k)*qr(i,j,k)/(Nr(i,j,k)+eps0) ! JvdD Added eps0 to avoid floating point exception
              xr (i,j,k) = min(xr(i,j,k),xrmaxkk)      ! to ensure x_pw is within borders
              Dvr(i,j,k) = (xr(i,j,k)/pirhow)**(1./3.)
            endif
         enddo
         enddo
         enddo

      endif ! l_sb
    endif   ! l_rain

  !*********************************************************************
  ! call microphysical processes subroutines
  !*********************************************************************

    call cloudactivation
    call cloudselfcollection

    if (l_sedc) then 
      call sedimentation_cloud
    end if

    if (l_rain) then
      call bulkmicrotend
      call autoconversion
      call bulkmicrotend
      call accretion
      call bulkmicrotend
      call evaporation
      call bulkmicrotend
      call sedimentation_rain
      call bulkmicrotend
    endif    

    call scavenging  
    call cloudevaporation

    do k=1,k1
    do j=2,j1
    do i=2,i1

      ! NOTE ----------------------------------------------------------------- !
      ! Before transferring the tendencies calculated here, values have to be  !
      ! CHECKED for mass conservation. Rationale: If total tendency would      !
      ! put a tracer to negative values, all local sinks are reduced by the    !
      ! same fraction so that the tracer will become 0. instead.               ! 
      !                                                                        !
      ! Tracers are checked/corrected in groups based on the aerosol 'state',  !
      ! i.e. free aerosol, in-cloud or in-rain as sinks processes are different!
      ! for each group.                                                        !
      !                                                                        !
      ! Restrictions:                                                          !
      ! - Tendencies cannot change sign                                        !
      ! - Only local processes are taken into account, i.e. no                 ! 
      ! advection/sedimentation                                                ! 
      !                                                                        !
      ! Possible issue:                                                        !
      ! Aerosol number and corresponding mass are checked independently. This  !
      ! could cause situations where these are corrected for unevenly.         !  
      ! -----------------------------------------------------------------------!

      ! Transform to svp units (e.g. [kg m-3 s-1]/[kga m-3] = [kg kga-1 s-1])  
!      aer_tend(i,j,k,:) = aer_tend(i,j,k,:)! For now assume tracer is in #/ms anyway /rhof(k)
!      aer_evpc(i,j,k,:) = aer_evpc(i,j,k,:)! Idem /rhof(k)

      do iaer = 1,naer

         ispec = nspec_type(iaer)
         imod  = nmod_type(iaer)

         aertest = svm(i,j,k,iaer+iaer_offset)/delt + svp(i,j,k,iaer+iaer_offset) + aer_tend(i,j,k,iaer) !(# or kg) kg-1 s-1

         if ( aertest < 0. ) then

            if ( imod <= 7 ) then ! Checking and correction for free aerosol sinks
 
               select case(ispec)
                  case(1); taer = inc
                  case(2); taer = iso4cld
                  case(3); taer = ibccld
                  case(4); taer = ipomcld
                  case(5); taer = isscld
                  case(6); taer = iducld
               end select  

               sinktotal = aer_acti(i,j,k,iaer) + aer_scvc(i,j,k,iaer) + aer_scvr(i,j,k,iaer) 
  
               if (sinktotal < 0.) then
                                
                  f_acti  = aer_acti(i,j,k,iaer)/sinktotal
                  f_scvc  = aer_scvc(i,j,k,iaer)/sinktotal
                  f_scvr  = aer_scvr(i,j,k,iaer)/sinktotal
                
                  fix = max(aertest, sinktotal)                      

                  aer_acti(i,j,k,iaer) = aer_acti(i,j,k,iaer) - fix*f_acti
                  aer_scvc(i,j,k,iaer) = aer_scvc(i,j,k,iaer) - fix*f_scvc
                  aer_scvr(i,j,k,iaer) = aer_scvr(i,j,k,iaer) - fix*f_scvr
                  aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) - fix 

                  ! Activation has both mass and number target
                  aer_acti(i,j,k,taer) = aer_acti(i,j,k,taer) + fix*f_acti
                  aer_tend(i,j,k,taer) = aer_tend(i,j,k,taer) + fix*f_acti

                  ! Scavenging has no number target, i.e. does only yield
                  ! in-cloud/rain mass increase
                  if (ispec /= 1) then
                     aer_scvc(i,j,k,taer)   = aer_scvc(i,j,k,taer)   + fix*f_scvc
                     aer_tend(i,j,k,taer)   = aer_tend(i,j,k,taer)   + fix*f_scvc
 
                     aer_scvr(i,j,k,taer+1) = aer_scvr(i,j,k,taer+1) + fix*f_scvr
                     aer_tend(i,j,k,taer+1) = aer_tend(i,j,k,taer+1) + fix*f_scvr
                  endif

               end if 
                 
            elseif ( imod == 8 ) then ! Checking and correction for cloud evporation
               
               select case ( ispec ) ! Select free aerosol target variable
                  case(1); taer = iacs_n  ! number
                  case(2); taer = iso4acs ! SO4
                  case(3); taer = ibcacs  ! BC
                  case(4); taer = ipomacs ! POM
                  case(5); taer = issacs  ! SS
                  case(6); taer = iduacs  ! DUST
               end select

               if (ispec /= 1) then
                  sinktotal = aer_evpc(i,j,k,iaer) + aer_auto(i,j,k,iaer) + aer_accr(i,j,k,iaer)

                  f_evpc  =  aer_evpc(i,j,k,iaer)/sinktotal
                  f_auto  =  aer_auto(i,j,k,iaer)/sinktotal
                  f_accr  =  aer_accr(i,j,k,iaer)/sinktotal
               else
                  sinktotal = aer_evpc(i,j,k,iaer) + &
                              aer_auto(i,j,k,iaer) + &
                              aer_accr(i,j,k,iaer) + &
                              aer_slfc(i,j,k,iaer)

                  f_evpc  =   aer_evpc(i,j,k,iaer)/sinktotal
                  f_auto  =   aer_auto(i,j,k,iaer)/sinktotal
                  f_accr  =   aer_accr(i,j,k,iaer)/sinktotal
                  f_slfc  =   aer_slfc(i,j,k,iaer)/sinktotal
               endif
               
               if (sinktotal < -eps0 ) then
                  
                  fix = max(aertest, sinktotal)

                  if (aer_evpc(i,j,k,iaer) < 0.) then
                     aer_tend(i,j,k,taer)   = aer_tend(i,j,k,taer)   + fix*f_evpc*( -aer_evpc(i,j,k,taer)/aer_evpc(i,j,k,iaer))
                     aer_evpc(i,j,k,taer)   = aer_evpc(i,j,k,taer)   + fix*f_evpc*( -aer_evpc(i,j,k,taer)/aer_evpc(i,j,k,iaer))

                     aer_tend(i,j,k,taer+1) = aer_tend(i,j,k,taer+1) + fix*f_evpc*(1+aer_evpc(i,j,k,taer)/aer_evpc(i,j,k,iaer))
                     aer_evpc(i,j,k,taer+1) = aer_evpc(i,j,k,taer+1) + fix*f_evpc*(1+aer_evpc(i,j,k,taer)/aer_evpc(i,j,k,iaer))

                     aer_evpc(i,j,k,iaer)   = aer_evpc(i,j,k,iaer) - fix*f_evpc
                     aer_tend(i,j,k,iaer)   = aer_tend(i,j,k,iaer) - fix*f_evpc
                  end if
     
                  if (ispec /= 1) then
                     aer_auto(i,j,k,iaer+1) = aer_auto(i,j,k,iaer+1) + fix*f_auto
                     aer_tend(i,j,k,iaer+1) = aer_tend(i,j,k,iaer+1) + fix*f_auto

                     aer_accr(i,j,k,iaer+1) = aer_accr(i,j,k,iaer+1) + fix*f_accr
                     aer_tend(i,j,k,iaer+1) = aer_tend(i,j,k,iaer+1) + fix*f_accr
                  else
                     ! Note:
                     ! 1. Accretion does not influence raindrop number concentration (Nr)
                     ! 2. Due to difference in size of cloud and raindrops, 
                     !    autoconversion Nr tendency has to be adjusted accordingly
                     ! 3. Selfcollection only applies to Nc tendency
                                         
                     aer_tend(i,j,k,iaer+1) = aer_tend(i,j,k,iaer+1) + fix*f_auto*(-aer_auto(i,j,k,iaer+1)/(aer_auto(i,j,k,iaer)-eps0)) 
                     aer_auto(i,j,k,iaer+1) = aer_auto(i,j,k,iaer+1) + fix*f_auto*(-aer_auto(i,j,k,iaer+1)/(aer_auto(i,j,k,iaer)-eps0)) 
 
                     aer_slfc(i,j,k,iaer)   = aer_slfc(i,j,k,iaer)   - fix*f_slfc
                     aer_tend(i,j,k,iaer)   = aer_tend(i,j,k,iaer)   - fix*f_slfc 
                  endif         

                  aer_auto(i,j,k,iaer)   = aer_auto(i,j,k,iaer) - fix*f_auto
                  aer_tend(i,j,k,iaer)   = aer_tend(i,j,k,iaer) - fix*f_auto
                     
                  aer_accr(i,j,k,iaer)   = aer_accr(i,j,k,iaer) - fix*f_accr
                  aer_tend(i,j,k,iaer)   = aer_tend(i,j,k,iaer) - fix*f_accr
               endif
 
            elseif ( imod == 9 ) then

               select case ( ispec ) ! Select free aerosol target variable
                  case(1); taer = iacs_n  ! number
                  case(2); taer = iso4acs ! SO4
                  case(3); taer = ibcacs  ! BC
                  case(4); taer = ipomacs ! POM
                  case(5); taer = issacs  ! SS
                  case(6); taer = iduacs  ! DUST
               end select

               if (aer_evpr(i,j,k,iaer) < -eps0) then
                   
                  fix = max(aertest,aer_evpr(i,j,k,iaer))     

                  !Distribution between acs and cos modes
                  f_acs = -aer_evpr(i,j,k,taer  )/aer_evpr(i,j,k,iaer)!(aer_evpr(i,j,k,taer)+aer_evpr(i,j,k,taer+1)) !ACS
                  f_cos = -aer_evpr(i,j,k,taer+1)/aer_evpr(i,j,k,iaer)!(aer_evpr(i,j,k,taer)+aer_evpr(i,j,k,taer+1)) !COS 

                  aer_tend(i,j,k,taer  ) = aer_tend(i,j,k,taer  ) + fix*f_acs!( -aer_evpr(i,j,k,taer)/aer_evpr(i,j,k,iaer))
                  aer_evpr(i,j,k,taer  ) = aer_evpr(i,j,k,taer  ) + fix*f_acs!( -aer_evpr(i,j,k,taer)/aer_evpr(i,j,k,iaer))

                  aer_tend(i,j,k,taer+1) = aer_tend(i,j,k,taer+1) + fix*f_cos!(1+aer_evpr(i,j,k,taer)/aer_evpr(i,j,k,iaer))
                  aer_evpr(i,j,k,taer+1) = aer_evpr(i,j,k,taer+1) + fix*f_cos!(1+aer_evpr(i,j,k,taer)/aer_evpr(i,j,k,iaer))

                  aer_tend(i,j,k,iaer  ) = aer_tend(i,j,k,iaer  ) - fix
                  aer_evpr(i,j,k,iaer  ) = aer_evpr(i,j,k,iaer  ) - fix

               end if  

            end if !imod

         endif !aertest < 0 etc.
      enddo !iaer

      !*********************************************************************
      ! remove negative values and non physical low values
      !*********************************************************************

      qrtest=svm(i,j,k,iqr)             + (svp(i,j,k,iqr)             + qrp(i,j,k))         *delt
      nrtest=svm(i,j,k,inr+iaer_offset) + (svp(i,j,k,inr+iaer_offset) + aer_tend(i,j,k,inr))*delt

      if ((qrtest < qrmin) .or. (nrtest < 0.) ) then ! correction, after Jerome's implementation in Gales
         qtp(i,j,k)     = qtp(i,j,k)  + qtpmcr(i,j,k)                      + svm(i,j,k,iqr)/delt + svp(i,j,k,iqr) + qrp(i,j,k)
         thlp(i,j,k)    = thlp(i,j,k) + thlpmcr(i,j,k) - (rlv/(cp*exnf(k)))*(svm(i,j,k,iqr)/delt + svp(i,j,k,iqr) + qrp(i,j,k))
         svp(i,j,k,iqr) = - svm(i,j,k,iqr)/delt

      else

         ! Transferring tendencies 
         svp(i,j,k,iqr)     = svp(i,j,k,iqr)     + qrp(i,j,k)
         thlp(i,j,k)        = thlp(i,j,k)        + thlpmcr(i,j,k)
         qtp(i,j,k)         = qtp(i,j,k)         + qtpmcr(i,j,k)

      end if

      do iaer= 1,naer
        svp(i,j,k,iaer+iaer_offset) = svp(i,j,k,iaer+iaer_offset) + aer_tend(i,j,k,iaer)        
      end do    

      Nrp(i,j,k) = aer_tend(i,j,k,inr) !For statistical purposes    

    enddo
    enddo
    enddo

  ! Transfer individual tendencies for statistical purposes

  if (l_aertend) then
     call bulkaertend
  end if
  
  end subroutine bulkmicro

  ! -----------------------------------------------------------------------
  ! SUBROUTINE CLOUDACTIVATION
  ! 
  ! Calculates the amount of activated aerosols that form cloud droplets.
  !  
  ! 2 Routines available, based on:
  ! (1) k-Kohler (Petters & Kreidenweis, 2007); PK07)
  ! (2) updraft  (Pousse-Nottelman, 2015); PN15)
  !
  ! Switch with l_kohler in namoptions. 
  ! l_kohler = .true.  -> PK07
  ! l_kohler = .false. -> PN15
  !
  ! -----------------------------------------------------------------------
  ! For PK07: Equation for r_crit derived from Eq. (10) PK07. 
  ! k values used for species taken from Pringle et al. (2010):
  ! SO4 0.88,   BC 0,   POM 0.1,   SS 1.28,   DU: 0.
  ! Procedure:
  ! 1. Calculation of mode total mass and average density and hygrospocity
  ! 2. Calculation of activated number/mass fraction per mode
  ! 3. Determine source and target variables and apply activation
  ! -----------------------------------------------------------------------
  ! For PN15:
  ! Parameterisation parameterized according to the approach by Lin
  ! and Leaitch (1997) following the work by Muhlbauer and Lohmann
  ! (2008) and Zubler et al. (2011a)
  !
  ! Procedure:
  ! 1. Determine number of aerosol particles >35nm (Nt)
  ! 2. Calculate dNc/dt
  ! 3. Calculate mass/number tendencies of affected tracers
  ! -----------------------------------------------------------------------

  subroutine cloudactivation
    use modglobal, only : i1,j1,k1, rlv, cp, rhow, pi
    use modfields, only : ql0, thl0, exnf, w0, svm, svp

    implicit none
    integer         :: i,j,k, taer, iaer, ispec, imod
    real            :: A, B, T, sigma, r_crit, d_crit, dn_mean, dm_mean, fac
    real            :: d_ais, f_ais, Nt, dNcdt  

    real, dimension(nmod-2) :: mod_numb, mod_mass, mod_rho, mod_k, fn, fm
 
    real, parameter ::    & 
        R     = 8.1314,   & ! Gas constant
        Mw    = 0.018,    & ! Molar mass of water (kg mol-1)
        sw    = 0.075       ! Surface tension of water 

    do j=2,j1
    do i=2,i1
    do k=1,k1    

       if ( qcmask(i,j,k) ) then

          mod_numb = 0.0
          mod_mass = 0.0
          mod_rho  = 0.0   
          mod_k    = 0.0

          fn = 0.0
          fm = 0.0
   
          if (l_kohler) then

          if ( Nc(i,j,k) < 1e6 ) then

             !Aggregate free aerosol modes
             do iaer = 1, naer !Aggregate modes

             imod  = nmod_type(iaer)
             ispec = nspec_type(iaer)

                if (imod > 7 ) cycle  ! Only free aerosol modes
                if (ispec == 1 ) then ! Mode number
                   mod_numb(imod) = aer_conc(i,j,k,iaer)
                else
                   mod_mass(imod) = mod_mass(imod) + aer_conc(i,j,k,iaer)
                   mod_rho (imod) = mod_rho (imod) + aer_conc(i,j,k,iaer)/spec_rho(ispec-1)
                   if (spec_k(ispec-1) > 0.) then
                      mod_k(imod) = mod_k   (imod) + aer_conc(i,j,k,iaer)/spec_k  (ispec-1)
                   endif
                end if
             end do   

             T = thl0(i,j,k) * exnf(k) + (rlv/cp) * ql0(i,j,k)
             A = 4.*sw*Mw/(R*T*rhow)   

             ! Calculate activated fraction --------------------------------
             do imod=1,nmod-2 ! Activation per mode, only free aerosol modes
                if (mod_mass(imod) < eps0 .or. mod_numb(imod) < eps0 .or. mod_k(imod) == 0.) cycle
                
                   mod_rho(imod) = mod_mass(imod)/mod_rho(imod)
                   mod_k (imod)  = mod_mass(imod)/mod_k (imod)             
  
                   d_crit = (4.*A**3./(27.*mod_k(imod)*log(Ssat/100.+1.)**2.))**(1./3.) 

                   sigma = sigma_lognormal(imod)
                   ! Number and mass mean diameter
                   dn_mean  = (6.*mod_mass(imod)/(pi*mod_numb(imod)*scalefac*mod_rho(imod)))**(1./3.)*exp(-(3./2.)*log(sigma)**2.)
                   dm_mean  = dn_mean * exp(3.*log(sigma)**2.)
                
                   ! Number and mass activation fraction
                   fn(imod) = 1. - .5*erfc(-log(d_crit/dn_mean)/sqrt(2.)*log(sigma))
                   fm(imod) = 1. - .5*erfc(-log(d_crit/dm_mean)/sqrt(2.)*log(sigma))

             enddo !imod

          end if !ncmin

          else !l_kohler
            
             !Aggregate free aerosol modes
             do iaer = 1, naer !Aggregate modes

                imod  = nmod_type(iaer)
                ispec = nspec_type(iaer)

                if ( imod == 1 .or. imod >= 5 ) cycle  ! Only AIS,ACS,COS
                if ( ispec == 1 ) then ! Mode number
                   mod_numb(imod) = aer_conc(i,j,k,iaer)
                else
                   mod_mass(imod) = mod_mass(imod) + aer_conc(i,j,k,iaer)
                   mod_rho (imod) = mod_rho (imod) + aer_conc(i,j,k,iaer)/spec_rho(ispec-1)
                end if
             end do

             if ( mod_numb(iais_n) + mod_numb(iacs_n) + mod_numb(icos_n) > 0. ) then   
                !Complete mode mean density calculation
                mod_mass(iais_n) = mod_mass(iais_n) + eps0
                mod_numb(iais_n) = mod_numb(iais_n) + eps0
                mod_rho(iais_n)  = mod_mass(iais_n)/mod_rho(iais_n)

                ! Determine fraction of AIS aerosol with radius > 35 nm and calculate Nt    
                sigma  = sigma_lognormal(iais_n)
                r_crit = 35.e-9 !m
            
                d_ais  = (6.*mod_mass(iais_n)/(pi*mod_numb(iais_n)*scalefac*mod_rho(iais_n)))**(1./3.)*exp(-(3./2.)*log(sigma)**2.)
                f_ais  = 1. - .5*erfc(-log(2.*r_crit/d_ais)/sqrt(2.)*log(sigma))
           
                Nt = 1e-6*(f_ais*mod_numb(iais_n) + mod_numb(iacs_n) + mod_numb(icos_n)) !Eq.(4) PN15
                dNcdt = max(1.e6/delt*( 0.1*(w0(i,j,k)*100.*Nt/(w0(i,j,k)*100.+.023*Nt))**1.27 - 1.e-6*Nc(i,j,k)), 0. ) !Eq.(2) PN15

                ! Compare activated particles with mode and calculate if mode is activated completely
                ! If yes, activate all mode mass and continue
                ! In not, calculated mass fraction
                ! Order: 1:COS, 2:ACS, 3:AIS
                
                do imod = 4, 2, -1 !Caution imod is not the proper iterator, but coincidentally iais_n, iacs_n and icos_n = 2,3,4
                   if ( dNcdt*delt > aer_conc(i,j,k,imod) ) then
                      fn(imod) = 1.0
                      fm(imod) = 1.0
                      dNcdt    = dNcdt - aer_conc(i,j,k,imod)/delt
                   elseif ( aer_conc(i,j,k,imod) > 0. ) then
                      fn(imod) = dNcdt*delt/aer_conc(i,j,k,imod)
                      fm(imod) = 1. - .5*erfc(xinverfc(2.*fn(imod))-3.*log(sigma_lognormal(imod))/sqrt(2.))
                      dNcdt    = 0.0 
                   endif
                enddo ! imod
             end if   ! Nt = 0
          end if      ! l_kohler

          ! Apply activation  --------------------------------------------

          do iaer= 1,naer

             ispec = nspec_type(iaer) ! Retrieve species
             imod  = nmod_type(iaer)  ! Retrieve source mode

             if (imod > 7) cycle      ! Exclude cloud and rain modes

             select case ( ispec )    ! Select cloud target variable
                case(1); taer = inc;     fac = fn(imod) ! Number
                case(2); taer = iso4cld; fac = fm(imod) ! SO4
                case(3); taer = ibccld;  fac = fm(imod) ! BC
                case(4); taer = ipomcld; fac = fm(imod) ! POM
                case(5); taer = isscld;  fac = fm(imod) ! SS
                case(6); taer = iducld;  fac = fm(imod) ! DUST
             end select

             ! Caution with sign, source (iaer) tend. negative, target (taer) tend. positive
             aer_acti(i,j,k,iaer) = - fac*aer_conc(i,j,k,iaer)/delt   
             aer_acti(i,j,k,taer) = aer_acti(i,j,k,taer) - aer_acti(i,j,k,iaer) !Double negative = positive, sum of multiple sources  

             aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) + aer_acti(i,j,k,iaer)
             aer_tend(i,j,k,taer) = aer_tend(i,j,k,taer) - aer_acti(i,j,k,iaer) ! No aer_acti(taer) because that carries multiple sources 

          enddo !iaer 

       end if !qcmask 
    end do
    end do
    end do

  end subroutine cloudactivation

  subroutine cloudevaporation  
    use modglobal, only : i1, j1, k1
    use modfields, only : ql0, svp, rhof, qlm, svm, sv0

    implicit none
    integer :: i,j,k, iaer, imod, ispec, taer
    real    :: evapn, evapr, evap, &
               Dn, Dm, Dc, Fn, Fm, Fx, sigma, &
               frac, corr
    real, dimension(nspec-1) :: evapm

    !Distribution width for evaporated aerosol
    sigma = 1.5 !Following Pousse-Nottelman et al. (2015)    

    do j=2,j1
    do i=2,i1
    do k=1,k1
        
    evapn = 0.
    evapm = 0.
    evapr = 0.

    if ( Nc(i,j,k) > ncmin ) then

       ! Calculate the resuspension fraction
       if ( qcmask(i,j,k) ) then
          if ( qlm(i,j,k) > ql0(i,j,k) ) then !partial evaporation
             frac = min((qlm(i,j,k)-ql0(i,j,k))/(qlm(i,j,k)), 1.) !min function to prevent overflow when qlm=0.
          else !no evaporation
             frac = 0.
          endif !qlm>ql0
       else ! complete evaporation
          frac = 1.
       endif !qcmask

       ! Correction factor between evaporated water and resuspended aerosol
       ! based on Gong et al. (2006)
       corr = (1. - exp(-2.*sqrt(frac))*(1. +2.*sqrt(frac) + 2.*frac + 4./3.*frac**(3./2.)) )*(1.-frac) + frac**2.

       ! Aggregate resuspended aerosol number and mass 
       if ( frac > 1e-3 ) then 
          do iaer = 1,naer  

             imod  =  nmod_type(iaer)
             ispec = nspec_type(iaer)
             
             if (imod == 8) then !Clouds only
                if (ispec == 1) then
                   evapn          = frac/delt*aer_conc(i,j,k,iaer)
                else
                   evapm(ispec-1) = corr*frac/delt*aer_conc(i,j,k,iaer)
                end if
             endif !imod == 8 Clouds only   

           enddo !iaer
        endif !frac

        if ( evapn > ncmin .and. SUM(evapm) > mcmin ) then

           ! Complete evaporated tracer value calculation
           evapr = SUM(evapm) / SUM(evapm/spec_rho)

           ! Calculate number/mass median diameter for evaporated aerosol 
           Dn = 1.e6 * ( (6.*SUM(evapm))/(3.141592654*evapn*evapr) )**(1./3.) * exp(-(3./2.)*log(sigma)**2.) !in micron
           Dm = Dn * exp(3.*log(sigma)**2.)
           Dc = 1. !Critical radius 0.5 micron, i.e. diameter 1 micron

           ! Calculate partitioning accumulation/coarse aerosol mode
           Fn = .5*erfc(- log(Dc/Dn)/( sqrt(2.)*log(sigma)) )
           Fm = .5*erfc(- log(Dc/Dm)/( sqrt(2.)*log(sigma)) )
    
           do iaer = 1,naer

              imod = nmod_type(iaer)
              ispec = nspec_type(iaer)

              if (imod == 8) then

                 select case ( ispec )
                    case(1); taer = iacs_n;  evap = evapn;          Fx = Fn ! Number
                    case(2); taer = iso4acs; evap = evapm(ispec-1); Fx = Fm ! SO4
                    case(3); taer = ibcacs;  evap = evapm(ispec-1); Fx = Fm ! BC
                    case(4); taer = ipomacs; evap = evapm(ispec-1); Fx = Fm ! POM
                    case(5); taer = issacs;  evap = evapm(ispec-1); Fx = Fm ! SS
                    case(6); taer = iduacs;  evap = evapm(ispec-1); Fx = Fm ! DU
                 end select 
                 
                 aer_evpc(i,j,k,iaer)   = -evap                      
                 aer_evpc(i,j,k,taer)   =  evap*Fx
                 aer_evpc(i,j,k,taer+1) =  evap*(1-Fx)
                
                 aer_tend(i,j,k,iaer)   = aer_tend(i,j,k,iaer)   - evap 
                 aer_tend(i,j,k,taer)   = aer_tend(i,j,k,taer)   + evap*Fx
                 aer_tend(i,j,k,taer+1) = aer_tend(i,j,k,taer+1) + evap*(1-Fx)
                                
              endif                

          enddo !iaer
       endif !evapn > ncmin, evapm > mcmin
    endif !Nc > ncmin

    end do
    end do
    end do        
  end subroutine cloudevaporation      

  subroutine cloudselfcollection
    use modglobal, only : i1,j1,k1
    use modfields, only : ql0, rhof
    
    implicit none
    integer i,j,k
    real    nuc    

    do j=2,j1
    do i=2,i1
    do k=1,k1

       nuc = 1580.*rhof(k)*ql0(i,j,k) - 0.28                                                        ! kg water/m3 ?
       aer_slfc(i,j,k,inc) = - k_c*(nuc + 2.)/(nuc +1.)*(ql0(i,j,k)*rhof(k))**2. * 1.225 / rhof(k)  ! #/m3/s       
       aer_tend(i,j,k,inc) = aer_tend(i,j,k,inc) + aer_slfc(i,j,k,inc) ! #/m3/s
 
    end do
    end do
    end do
    
  end subroutine cloudselfcollection      

  !> Determine scavenging rate and adjust free aerosol and in-cloud/in-rain aerosol accordingly
  !!   
  !! The scavenging rate is based on the work of Croft et al. (2010). Scavenging
  !! efficiences are stored in a lookup table which is constructed during
  !! initialization. This table contains mean scavenging coefficients (s^-1) as
  !! a function of geometric mean aerosol radius (aerrad) and mean cloud drop radius (cldrad). 
  !!  
  !! Bilinear interpolation is applied to the logarithm of aerrad and cldrad to
  !! cover the many orders of magnitude.
  !!
  !! Notes 
  !! 1. Scavenging is applied to all aerosol modes individually. 
  !! 2. Mass is aggregated per species in the cloud/rain water. 
  !! 3. Scavenged aerosol number is 'lost' as all aerosol mass is assumed to dissolve and 
  !!    number is now equal to cloud/rain droplet number.
  !! 
  subroutine scavenging
  
  use modglobal, only : i1,j1,k1,rhow, pi
  use modfields, only : ql0,sv0,rhof 

  implicit none

  integer :: i,j,k, iaer, imod, ispec, taerinc, taerblc
  real    :: meancld, meanaer, meanaer_inc, rainrate, scavinc, scavblc

  real, dimension(nmod-2) :: mod_mass, mod_numb, mod_rho, &
                             scavinc_m, scavinc_n, scavblc_m, scavblc_n

  ! ---------------------------------------------------------------------------

  do i=2,i1
  do j=2,j1
  do k=1,k1

     meancld =  min( max( 1e6*(3.*ql0(i,j,k)*rhof(k)/(4.*3.1415*Nc(i,j,k)*rhow))**(1/3.) ,5.001), 49.999) ! in micron
     rainrate = min( max(sed_qr(i,j,k)*3600., 0.01001) , 99.999)                                           !kg m-2 s-1 --> mm/hr
  
     mod_numb = 0.0
     mod_mass = 0.0        
     mod_rho  = 0.0   

     ! Calculate number and mass for each mode and volume mean density
     do iaer = 1,naer

        imod = nmod_type(iaer)
        ispec = nspec_type(iaer)

         if (imod > 7 ) cycle !Exclude cloud & rain modes

         if (ispec == 1 ) then  ! Mode number
           mod_numb(imod) = aer_conc(i,j,k,iaer)
         else
           mod_mass(imod) = mod_mass(imod) + aer_conc(i,j,k,iaer)
           mod_rho (imod) = mod_rho (imod) + aer_conc(i,j,k,iaer)/spec_rho(ispec-1)  
         end if

     end do !iaer

     scavinc_n = 0.0
     scavinc_m = 0.0
     scavblc_n = 0.0
     scavblc_m = 0.0

     ! ---------------------------------------------------------------------------
     do imod= 1,nmod-2 !Only free aerosol modes

        ! Complete volume mean density calculation     
        mod_rho(imod) = mod_mass(imod)/mod_rho(imod)     

        if ( (mod_mass(imod) > eps0) .and. (mod_numb(imod) > ncmin)) then
   
           ! Get geometric mean aerosol radius (in meter)
           meanaer = 0.5*((6.*mod_mass(imod))/(pi*mod_numb(imod)*scalefac*mod_rho(imod)))**(1./3.)*exp(-(3./2.)*log(sigma_lognormal(imod))**2.)

           ! IN-CLOUD SCAVENGING  
           ! ==================================================================================== !
           !                                                                                      !
           ! Apply the interpolation to find efficiency coefficient
           ! Calculate tendency and transfer to corresponding arrays
           !   
           ! Note that:
           ! Stated in pers. comm. from Betty Croft:
           ! "Data files containing in-cloud scavenging coefficients for a cloud
           ! droplet number concentration (CDNC) of 1 cm^-3." 
           ! Also:
           ! "Coefficients are given as a function of mode radius and can be used
           ! directly by multiplying by the CDNC"
           ! So we multiply by the CDNC in units of cm^-1, i.e. Nc(i,j,k)*1e-6
           ! ==================================================================================== !
           if ( qcmask(i,j,k) .and. ( Nc(i,j,k) > ncmin )) then

              ! Transform aerosol radius to cm and apply bounds of lookup table (80 um > r > 0.1 nm)  
              meanaer_inc = max(min(100.*meanaer,8.0E-03),1.0E-8) ! in cm

              scavinc_n(imod) = 1e-6*Nc(i,j,k)*lookup_interpolate(ncld,logcldrad,naer,logaerrad,incnumb,log(meancld),log(meanaer_inc)) 
              scavinc_m(imod) = 1e-6*Nc(i,j,k)*lookup_interpolate(ncld,logcldrad,naer,logaerrad,incmass,log(meancld),log(meanaer_inc))

              if ((scavinc_n(imod)*delt > 1.) .or. (scavinc_m(imod)*delt > 1.)) then   
                 scavinc_n(imod) = 1./delt
                 scavinc_m(imod) = 1./delt
              end if

           end if !qcmask, Nc>0

           ! BELOW-CLOUD SCAVENGING  
           ! ==================================================================================== !
           ! EXPLANATION                                                                          !
           !                                                                                      !
           ! Method taken from Croft et al. (2009)                                                !
           ! Eq (1) : dC/dt = C*f*(R*F)                                                           !
           ! Where C is the ambient mixing ratio of tracer, f the cloud fraction,                 !
           ! R the normalized scavenging coeff. and F the precipitation flux                      !
           !                                                                                      !
           ! However:                                                                             !
           ! 1) The values in the lookuptable (L) are applied as R*F as stated on page 4656       !
           !    of Croft et al., and                                                              !
           ! 2) in DALES a gridbox is either cover by cloud or not. So cloud fraction             !
           !    is not used, i.e. equals 1 if  method is applied                                  !
           !                                                                                      !
           ! This leads to dC/dt      = C*L           or in the variables used in the code below: ! 
           !     aer_tend(i,j,k,iaer) = aer_conc(i,j,k,iaer) * scav                               !
           ! ==================================================================================== !
           if (l_rain) then
              if ( qrmask(i,j,k) .and. (sed_qr(i,j,k)*3600. > .01) ) then   
                 
                 ! Transform aerosol radius to micron and apply bounds of lookup table ( 1 nm > r > 1000 um )
                 meanaer  = max(min(.9999e3,meanaer*1e6 ),1.001e-3) !in micron 
                 
                 scavblc_n(imod) = lookup_interpolate(nrai,lograinrate,naer,logaerrad_blc,blcnumb,log(rainrate),log(meanaer))
                 scavblc_m(imod) = lookup_interpolate(nrai,lograinrate,naer,logaerrad_blc,blcmass,log(rainrate),log(meanaer))
              end if !qrmask, sed_qr
           end if !l_rain
        end if !mod_mass, mod_numb

     enddo !imod

     ! ---------------------------------------------------------------------------
     do iaer= 1,naer

        ispec = nspec_type(iaer) ! Retrieve species
        imod  = nmod_type(iaer)  ! Retrieve source mode

        if ( imod > 7 ) cycle    ! Exclude cloud and rain modes

        select case ( ispec )    ! Select cloud/rain target variable
           case(1); scavinc = scavinc_n(imod); scavblc = scavblc_n(imod); taerinc = 0;       taerblc = 0       ! Number
           case(2); scavinc = scavinc_m(imod); scavblc = scavblc_m(imod); taerinc = iso4cld; taerblc = iso4rai ! SO4
           case(3); scavinc = scavinc_m(imod); scavblc = scavblc_m(imod); taerinc = ibccld;  taerblc = ibcrai  ! BC
           case(4); scavinc = scavinc_m(imod); scavblc = scavblc_m(imod); taerinc = ipomcld; taerblc = ipomrai ! POM
           case(5); scavinc = scavinc_m(imod); scavblc = scavblc_m(imod); taerinc = isscld;  taerblc = issrai  ! SS
           case(6); scavinc = scavinc_m(imod); scavblc = scavblc_m(imod); taerinc = iducld;  taerblc = idurai  ! DUST
        end select

        aer_scvc(i,j,k,iaer) = - scavinc*max(0.,aer_conc(i,j,k,iaer))
        aer_scvr(i,j,k,iaer) = - scavblc*max(0.,aer_conc(i,j,k,iaer))
        
        aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) + aer_scvc(i,j,k,iaer)
        aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) + aer_scvr(i,j,k,iaer)
        
        if (ispec /= 1 ) then !Scavenging does not increase Nc 

           aer_scvc(i,j,k,taerinc) = aer_scvc(i,j,k,taerinc) - aer_scvc(i,j,k,iaer)
           aer_scvr(i,j,k,taerblc) = aer_scvr(i,j,k,taerblc) - aer_scvr(i,j,k,iaer)

           aer_tend(i,j,k,taerinc) = aer_tend(i,j,k,taerinc) - aer_scvc(i,j,k,iaer)
           aer_tend(i,j,k,taerblc) = aer_tend(i,j,k,taerblc) - aer_scvr(i,j,k,iaer)

        end if

     enddo !iaer                   
           
  end do
  end do
  end do

  end subroutine scavenging
  
  !> Determine autoconversion rate and adjust qrp and Nrp accordingly
  !!
  !!   The autoconversion rate is formulated for f(x)=A*x**(nuc)*exp(-Bx),
  !!   decaying exponentially for droplet mass x.
  !!   It can easily be reformulated for f(x)=A*x**(nuc)*exp(-Bx**(mu)) and
  !!   by chosing mu=1/3 one would get a gamma distribution in drop diameter
  !!   -> faster rain formation. (Seifert)
  subroutine autoconversion
    use modglobal, only : i1,j1,k1,kmax,rlv,cp
    use modmpi,    only : myid
    use modfields, only : exnf,rhof,ql0

    implicit none
    integer i,j,k,iaer

    au = 0.

    if (l_sb) then
    !
    ! SB autoconversion
    !
      tau(2:i1,2:j1,1:k1) = 0.
      phi(2:i1,2:j1,1:k1) = 0.
      nuc(2:i1,2:j1,1:k1) = 0.
      k_au = k_c/(20.*x_s)

      do j=2,j1
      do i=2,i1
      do k=1,k1
         if ( qcmask(i,j,k) .and. ( Nc(i,j,k) > ncmin ) ) then
            nuc    (i,j,k) = 1.58*(rhof(k)*ql0(i,j,k)*1000.) + 0.72 - 1. ! G09a
            xc     (i,j,k) = rhof(k) * ql0(i,j,k) / Nc(i,j,k)            ! No eps0 necessary
            au     (i,j,k) = k_au * (nuc(i,j,k)+2.) * (nuc(i,j,k)+4.) / (nuc(i,j,k)+1.)**2. &
                             * (ql0(i,j,k) * xc(i,j,k))**2. * 1.225 ! *rho**2/rho/rho (= 1)
            tau    (i,j,k) = 1. - ql0(i,j,k) / qltot(i,j,k)
            phi    (i,j,k) = k_1 * tau(i,j,k)**k_2 * (1. - tau(i,j,k)**k_2)**3
            au     (i,j,k) = au(i,j,k) * (1. + phi(i,j,k)/(1. - tau(i,j,k))**2)

            !Cannot take more cloud water than available !MdB    
            au(i,j,k) =  min( ql0(i,j,k)/delt , au(i,j,k) )   

            !Water    
            qtpmcr (i,j,k) = qtpmcr (i,j,k) - au(i,j,k)
            qrp    (i,j,k) = qrp    (i,j,k) + au(i,j,k)

            !Particle number
            aer_auto(i,j,k,inc) = - au(i,j,k)/xc(i,j,k)*rhof(k) !Nc is in #/m3, so correct with rhof   
            aer_auto(i,j,k,inr) =   au(i,j,k)/x_s*rhof(k)        

            aer_tend(i,j,k,inc) = aer_tend(i,j,k,inc) + aer_auto(i,j,k,inc)   
            aer_tend(i,j,k,inr) = aer_tend(i,j,k,inr) + aer_auto(i,j,k,inr)

            !Aerosol mass
            do iaer= 1,naer
               if (nmod_type(iaer) == 8 .and. nspec_type(iaer) > 1) then

                  aer_auto(i,j,k,iaer)   = - au(i,j,k)/ql0(i,j,k)*aer_conc(i,j,k,iaer)
                  aer_auto(i,j,k,iaer+1) = + au(i,j,k)/ql0(i,j,k)*aer_conc(i,j,k,iaer)
      
                  aer_tend(i,j,k,iaer)   = aer_tend(i,j,k,iaer)   + aer_auto(i,j,k,iaer)   !- au(i,j,k)/ql0(i,j,k)*aer_conc(i,j,k,iaer) 
                  aer_tend(i,j,k,iaer+1) = aer_tend(i,j,k,iaer+1) + aer_auto(i,j,k,iaer+1) !  au(i,j,k)/ql0(i,j,k)*aer_conc(i,j,k,iaer) 

               endif
            enddo    
            !Associated heat
            thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*au(i,j,k)

         endif
      enddo
      enddo
      enddo

    else
    !
    ! KK00 autoconversion
    !
      do j=2,j1
      do i=2,i1
      do k=1,k1

         if (qcmask(i,j,k)) then
            au     (i,j,k) = 1350.0 * ql0(i,j,k)**(2.47) * (Nc(i,j,k)/1.0E6)**(-1.79)

            qrp    (i,j,k) = qrp    (i,j,k) + au(i,j,k)
            nrp    (i,j,k) = nrp    (i,j,k) + au(i,j,k) * rhof(k)/(pirhow*D0_kk**3.)
            qtpmcr (i,j,k) = qtpmcr (i,j,k) - au(i,j,k)
            thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*au(i,j,k)
         endif

      enddo
      enddo
      enddo

    end if !l_sb


    if (any(ql0(2:i1,2:j1,1:kmax)/delt - au(2:i1,2:j1,1:kmax) .lt. 0.)) then
      write(6,*)'au too large', count(ql0(2:i1,2:j1,1:kmax)/delt - au(2:i1,2:j1,1:kmax) .lt. 0.),myid
    end if

  end subroutine autoconversion

  subroutine accretion
  !*********************************************************************
  ! determine accr. + self coll. + br-up rate and adjust qrp and Nrp
  ! accordingly. Break-up : Seifert (2007)
  !*********************************************************************
    use modglobal, only : ih,i1,jh,j1,k1,kmax,rlv,cp
    use modfields, only : exnf,rhof,ql0
    use modmpi,    only : myid
    implicit none
    real , allocatable :: phi_br(:,:,:)
    integer :: i,j,k,iaer
    allocate (phi_br(2-ih:i1+ih,2-jh:j1+jh,k1))

    ac(2:i1,2:j1,1:k1)=0

    if (l_sb) then
    !
    ! SB accretion
    !

     do j=2,j1
     do i=2,i1
     do k=1,k1
        if (qrmask(i,j,k) .and. qcmask(i,j,k) .and. Nc(i,j,k) > ncmin) then
           xc     (i,j,k) = rhof(k) * ql0(i,j,k) / Nc(i,j,k)
           tau    (i,j,k) = 1.0 - ql0(i,j,k)/(qltot(i,j,k))
           phi    (i,j,k) = (tau(i,j,k)/(tau(i,j,k) + k_l))**4.
           ac     (i,j,k) = k_r *rhof(k)*ql0(i,j,k) * qr(i,j,k) * phi(i,j,k) * &
                            (1.225/rhof(k))**0.5

           !Water
           qtpmcr (i,j,k) = qtpmcr (i,j,k) - ac(i,j,k)
           qrp    (i,j,k) = qrp    (i,j,k) + ac(i,j,k)

           !Number
           aer_accr(i,j,k,inc) = - ac(i,j,k)/xc(i,j,k)*rhof(k)     
           aer_tend(i,j,k,inc) = aer_tend(i,j,k,inc) + aer_accr(i,j,k,inc)
           !Number of raindrops does not change

           !Aerosol mass
           do iaer = 1,naer
              if (nmod_type(iaer) == 8 .and. nspec_type(iaer) > 1) then

                  aer_accr(i,j,k,iaer)   = - ac(i,j,k)/ql0(i,j,k)*aer_conc(i,j,k,iaer)
                  aer_accr(i,j,k,iaer+1) = + ac(i,j,k)/ql0(i,j,k)*aer_conc(i,j,k,iaer)

                  aer_tend(i,j,k,iaer)   = aer_tend(i,j,k,iaer)   + aer_accr(i,j,k,iaer)
                  aer_tend(i,j,k,iaer+1) = aer_tend(i,j,k,iaer+1) + aer_accr(i,j,k,iaer+1)

              endif
           enddo

           !Associated heat
           thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*ac(i,j,k)
        endif
     enddo
     enddo
     enddo

    !
    ! SB self-collection & Break-up
    !
     sc(2:i1,2:j1,1:k1)=0
     br(2:i1,2:j1,1:k1)=0

     do j=2,j1
     do i=2,i1
     do k=1,k1
        if (qrmask(i,j,k)) then
           sc(i,j,k) = k_rr *rhof(k)* qr(i,j,k) * Nr(i,j,k)  &
                       * (1 + kappa_r/lbdr(i,j,k)*pirhow**(1./3.))**(-9.)* (1.225/rhof(k))**0.5
        endif
        if (Dvr(i,j,k) .gt. 0.30E-3 .and. qrmask(i,j,k)) then
           phi_br(i,j,k) = k_br * (Dvr(i,j,k)-D_eq)
           br(i,j,k) = (phi_br(i,j,k) + 1.) * sc(i,j,k)
        else
           br(i,j,k) = 0. ! (phi_br = -1)
        endif
     enddo
     enddo
     enddo

     aer_tend(2:i1,2:j1,1:k1,inr) = aer_tend(2:i1,2:j1,1:k1,inr) - sc(2:i1,2:j1,1:k1) + br(2:i1,2:j1,1:k1)

    else
    !
    ! KK00 accretion
    !
     do j=2,j1
     do i=2,i1
     do k=1,k1
        if (qrmask(i,j,k) .and. qcmask(i,j,k)) then
           ac     (i,j,k) = 67.0 * ( ql0(i,j,k) * qr(i,j,k) )**1.15
           qrp    (i,j,k) = qrp     (i,j,k) + ac(i,j,k)
           qtpmcr (i,j,k) = qtpmcr  (i,j,k) - ac(i,j,k)
           thlpmcr(i,j,k) = thlpmcr (i,j,k) + (rlv/(cp*exnf(k)))*ac(i,j,k)
        endif
     enddo
     enddo
     enddo

    end if !l_sb


   if (any(ql0(2:i1,2:j1,1:kmax)/delt - ac(2:i1,2:j1,1:kmax) .lt. 0.)) then
     write(6,*)'ac too large', count(ql0(2:i1,2:j1,1:kmax)/delt - ac(2:i1,2:j1,1:kmax) .lt. 0.),myid
   end if

   deallocate (phi_br)

  end subroutine accretion

!> Sedimentation of cloud water
!!
!! The sedimentation of cloud droplets assumes a lognormal DSD in which the
!! geometric std dev. is assumed to be fixed at 1.3.
!! sedimentation of cloud droplets
!! lognormal CDSD is assumed (1 free parameter : sig_g)
!! terminal velocity : Stokes velocity is assumed (v(D) ~ D^2)
!! flux is calc. anal.

  subroutine sedimentation_cloud

    use modglobal, only : i1,j1,k1,kmax,rlv,cp,dzf,pi
    use modfields, only : rhof,exnf,ql0
    implicit none
    integer :: i,j,k,iaer

    sedc (2:i1,2:j1,1:k1) = 0.
    sedcn(2:i1,2:j1,1:k1) = 0.
    sedcm(2:i1,2:j1,1:k1,1:nspec-1) = 0.

    csed = c_St*(3./(4.*pi*rhow))**(2./3.)*exp(5.*log(sig_g)**2.)

    do j=2,j1
    do i=2,i1
    do k=1,k1
       if (qcmask(i,j,k) .and. Nc(i,j,k) > ncmin) then

          sedc (i,j,k) = csed*(Nc(i,j,k))**(-2./3.)*(ql0(i,j,k)*rhof(k))**(5./3.)
          
          !Particle number      
          sedcn(i,j,k) = sedc(i,j,k)/(rhof(k)*ql0(i,j,k))*Nc(i,j,k)

          !In-cloud aerosol mass
          do iaer = 1,naer
             if (nmod_type(iaer) == 8 .and. nspec_type(iaer) > 1) then
                sedcm(i,j,k,nspec_type(iaer)-1) = sedc(i,j,k)/(rhof(k)*ql0(i,j,k))*aer_conc(i,j,k,iaer)
             end if
          enddo !iaer

       endif
    enddo
    enddo
    enddo

    do k=1,kmax
    do j=2,j1
    do i=2,i1

      qtpmcr(i,j,k) = qtpmcr(i,j,k) + (sedc(i,j,k+1)-sedc(i,j,k))/(dzf(k)*rhof(k))
      thlpmcr(i,j,k) = thlpmcr(i,j,k) - (rlv/(cp*exnf(k))) &
                       *(sedc(i,j,k+1)-sedc(i,j,k))/(dzf(k)*rhof(k))

      aer_tend(i,j,k,inc) = aer_tend(i,j,k,inc) + (sedcn(i,j,k+1)-sedcn(i,j,k))/dzf(k)

      do iaer = 1,naer  
         if (nmod_type(iaer) == 8 .and. nspec_type(iaer) > 1) then
            aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) + (sedcm(i,j,k+1,nspec_type(iaer)-1)-sedcm(i,j,k,nspec_type(iaer)-1))/dzf(k)
         end if
      enddo

    enddo
    enddo
    enddo
  end subroutine sedimentation_cloud


!> Sedimentaion of rain
!! sedimentation of drizzle water
!! - gen. gamma distr is assumed. Terminal velocities param according to
!!   Stevens & Seifert. Flux are calc. anal.
!! - l_lognormal =T : lognormal DSD is assumed with D_g and N known and
!!   sig_g assumed. Flux are calc. numerically with help of a
!!   polynomial function
  subroutine sedimentation_rain
    use modglobal, only : ih,i1,jh,j1,k1,kmax,eps1,dzf
    use modfields, only : rhof, svm
    use modmpi,    only : myid
    implicit none
    integer :: i,j,k,jn, iaer, ispec
    integer :: n_spl      !<  sedimentation time splitting loop
    real    :: pwcont
    real, allocatable :: wvar    (:,:,:) &  
                        ,xr_spl  (:,:,:) &
                        ,Dvr_spl (:,:,:) &
                        ,mur_spl (:,:,:) &
                        ,lbdr_spl(:,:,:) &
                        ,Dgr     (:,:,:) &
                        ,aer_spl (:,:,:,:)

    real,save :: dt_spl,wfallmax, aer_in, aer_out

    allocate( wvar     (2-ih:i1+ih,2-jh:j1+jh,k1) & !<  work variable
              ,xr_spl  (2-ih:i1+ih,2-jh:j1+jh,k1) & !<  for time splitting
              ,Dvr_spl (2-ih:i1+ih,2-jh:j1+jh,k1) & !<     -
              ,mur_spl (2-ih:i1+ih,2-jh:j1+jh,k1) & !<     -
              ,lbdr_spl(2-ih:i1+ih,2-jh:j1+jh,k1) & !<     -
              ,Dgr     (2-ih:i1+ih,2-jh:j1+jh,k1) & !<  lognormal geometric diameter
              ,aer_spl (2-ih:i1+ih,2-jh:j1+jh,k1,nspec-1) ) 
  
    qr_spl(2:i1,2:j1,1:k1) = qr(2:i1,2:j1,1:k1)
    Nr_spl(2:i1,2:j1,1:k1) = Nr(2:i1,2:j1,1:k1)
 
    ! Copy all species of in-rain aerosol mass 
    do iaer = 1,naer    
       if (nmod_type(iaer) == 9 .and. nspec_type(iaer) > 1) then
          aer_spl(2:i1,2:j1,1:k1,nspec_type(iaer)-1) = aer_conc(2:i1,2:j1,1:k1,iaer)
       end if
    enddo

    wfallmax = 9.9
    n_spl = ceiling(wfallmax*delt/(minval(dzf)))
    dt_spl = delt/real(n_spl)

    do jn = 1 , n_spl ! time splitting loop

      sed_qr(2:i1,2:j1,1:k1) = 0.
      sed_Nr(2:i1,2:j1,1:k1) = 0.

      if (l_sb ) then

       do j=2,j1
       do i=2,i1
       do k=1,k1
        if (qr_spl(i,j,k) > qrmin) then
          xr_spl (i,j,k) = rhof(k)*qr_spl(i,j,k)/(Nr_spl(i,j,k)+eps0) ! JvdD Added eps0 to avoid division by zero
          xr_spl (i,j,k) = min(max(xr_spl(i,j,k),xrmin),xrmax) ! to ensure xr is within borders
          Dvr_spl(i,j,k) = (xr_spl(i,j,k)/pirhow)**(1./3.)
        endif
       enddo
       enddo
       enddo


      if (l_lognormal) then
        do j = 2,j1
        do i = 2,i1
        do k = 1,kmax
          if (qr_spl(i,j,k) > qrmin) then
            Dgr(i,j,k) = (exp(4.5*(log(sig_gr))**2))**(-1./3.)*Dvr_spl(i,j,k) ! correction for width of DSD
            sed_qr(i,j,k) = 1.*sed_flux(Nr_spl(i,j,k),Dgr(i,j,k),log(sig_gr)**2,D_s,3)
            sed_Nr(i,j,k) = 1./pirhow*sed_flux(Nr_spl(i,j,k),Dgr(i,j,k) ,log(sig_gr)**2,D_s,0)
!        correction for the fact that pwcont .ne. qr_spl
!        actually in this way for every grid box a fall velocity is determined
            pwcont = liq_cont(Nr_spl(i,j,k),Dgr(i,j,k),log(sig_gr)**2,D_s,3)         ! note : kg m-3
            if (pwcont > eps1) then
              sed_qr(i,j,k) = (qr_spl(i,j,k)*rhof(k)/pwcont)*sed_qr(i,j,k)  ! or qr_spl*(sed_qr/pwcont) = qr_spl*fallvel.
            end if
          end if ! qr_spl threshold statement
        end do
        end do
        end do

      else
    !
    ! SB rain sedimentation
    !
        if (l_mur_cst) then
          mur_spl(2:i1,2:j1,1:k1) = mur_cst
        else
          do j=2,j1
          do i=2,i1
          do k=1,k1
            if (qr_spl(i,j,k) > qrmin) then
!             mur_spl(i,j,k) = 10. * (1+tanh(1200.*(Dvr_spl(i,j,k)-0.0014))) ! SS08
              mur_spl(i,j,k) = min(30.,- 1. + 0.008/ (qr_spl(i,j,k)*rhof(k))**0.6)  ! G09b
            endif
          enddo
          enddo
          enddo

        endif

        do j=2,j1
        do i=2,i1
        do k=1,k1
          if (qr_spl(i,j,k) > qrmin) then
              lbdr_spl(i,j,k) = ((mur_spl(i,j,k)+3.)*(mur_spl(i,j,k)+2.)* &
                                 (mur_spl(i,j,k)+1.))**(1./3.)/Dvr_spl(i,j,k)
              wfall_qr(i,j,k) = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr_spl(i,j,k))**(-1.*(mur_spl(i,j,k)+4.))))
              wfall_Nr(i,j,k) = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr_spl(i,j,k))**(-1.*(mur_spl(i,j,k)+1.))))
              sed_qr  (i,j,k) = wfall_qr(i,j,k)*qr_spl(i,j,k)*rhof(k)
              sed_Nr  (i,j,k) = wfall_Nr(i,j,k)*Nr_spl(i,j,k)
          endif
        enddo
        enddo
        enddo

      endif !l_lognormal

    else
    !
    ! KK00 rain sedimentation
    !
      do j=2,j1
      do i=2,i1
      do k=1,k1
        if (qr_spl(i,j,k) > qrmin) then
           xr_spl(i,j,k) = rhof(k)*qr_spl(i,j,k)/(Nr_spl(i,j,k)+eps0) !JvdD added eps0 to avoid division by zero
           xr_spl(i,j,k) = min(xr_spl(i,j,k),xrmaxkk) ! to ensure xr is within borders
           Dvr_spl(i,j,k) = (xr_spl(i,j,k)/pirhow)**(1./3.)
           sed_qr(i,j,k) = max(0., 0.006*1.0E6*Dvr_spl(i,j,k)- 0.2) * qr_spl(i,j,k)*rhof(k)
           sed_Nr(i,j,k) = max(0.,0.0035*1.0E6*Dvr_spl(i,j,k)- 0.1) * Nr_spl(i,j,k)
        endif
      enddo
      enddo
      enddo

    end if !l_sb
!
    do k = 1,kmax
      do j=2,j1
      do i=2,i1
        wvar(i,j,k) = qr_spl(i,j,k) + (sed_qr(i,j,k+1) - sed_qr(i,j,k))*dt_spl/(dzf(k)*rhof(k))
      enddo
      enddo

      if (any(wvar(2:i1,2:j1,k) .lt. 0.)) then
        write(6,*)'sed qr too large', count(wvar(2:i1,2:j1,k) .lt. 0.),myid, minval(wvar), minloc(wvar)
      end if

      do j=2,j1
      do i=2,i1
        
        Nr_spl(i,j,k) = Nr_spl(i,j,k) + (sed_Nr(i,j,k+1) - sed_Nr(i,j,k))*dt_spl/dzf(k)

        do ispec = 1,nspec-1
           aer_in  = 0.
           aer_out = 0.     
     
           if (qr_spl(i,j,k+1) > 0.) then
              aer_in  = min(rhof(k+1)*dzf(k+1)/dt_spl      *aer_spl(i,j,k+1,ispec), &
                            sed_qr(i,j,k+1)/qr_spl(i,j,k+1)*aer_spl(i,j,k+1,ispec))
           endif
           if (qr_spl(i,j,k  ) > 0.) then
              aer_out = min(rhof(k  )*dzf(k  )/dt_spl      *aer_spl(i,j,k  ,ispec), &
                            sed_qr(i,j,k  )/qr_spl(i,j,k  )*aer_spl(i,j,k  ,ispec))
           endif
           
           aer_spl(i,j,k,ispec) = aer_spl(i,j,k,ispec) + (aer_in - aer_out)*dt_spl/(dzf(k)*rhof(k))
           if (aer_spl(i,j,k,ispec) < 0.) then
              write(6,"(A13, E11.4)") 'SEDIMENTATION', aer_in/(dt_spl/(dzf(k+1)*rhof(k+1))), aer_out/(dt_spl/(dzf(k)*rhof(k)))
           end if          
        enddo

        qr_spl(i,j,k) = qr_spl(i,j,k) + (sed_qr(i,j,k+1) - sed_qr(i,j,k))*dt_spl/(dzf(k)*rhof(k))

      enddo
      enddo

      if ( jn == 1. ) then
      do j=2,j1
      do i=2,i1
        precep(i,j,k) =  sed_qr(i,j,k)/rhof(k)   ! kg kg-1 m s-1
      enddo
      enddo
      endif
    end do  ! second k loop

    enddo ! time splitting loop

    qrp     (2:i1,2:j1,1:k1)    =      qrp(2:i1,2:j1,1:k1)     + (qr_spl(2:i1,2:j1,1:k1) - qr(2:i1,2:j1,1:k1))/delt
        
    aer_tend(2:i1,2:j1,1:k1,inr)= aer_tend(2:i1,2:j1,1:k1,inr) + (Nr_spl(2:i1,2:j1,1:k1) - Nr(2:i1,2:j1,1:k1))/delt
    aer_sedr(2:i1,2:j1,1:k1,inr)=                                (Nr_spl(2:i1,2:j1,1:k1) - Nr(2:i1,2:j1,1:k1))/delt           

    do iaer = 1,naer    
       if (nmod_type(iaer) == 9 .and. nspec_type(iaer) > 1) then
          aer_tend(2:i1,2:j1,1:k1,iaer)= aer_tend(2:i1,2:j1,1:k1,iaer) + (aer_spl(2:i1,2:j1,1:k1,nspec_type(iaer)-1) - aer_conc(2:i1,2:j1,1:k1,iaer))/delt
          aer_sedr(2:i1,2:j1,1:k1,iaer)=                                 (aer_spl(2:i1,2:j1,1:k1,nspec_type(iaer)-1) - aer_conc(2:i1,2:j1,1:k1,iaer))/delt
       endif
    enddo

    deallocate (wvar, xr_spl,Dvr_spl,mur_spl,lbdr_spl,Dgr)
  end subroutine sedimentation_rain

  !*********************************************************************
  !*********************************************************************

  subroutine evaporation
  !*********************************************************************
  ! Evaporation of prec. : Seifert (2008)
  ! Cond. (S>0.) neglected (all water is condensed on cloud droplets)
  !*********************************************************************

    use modglobal, only : ih,i1,jh,j1,k1,rv,rlv,cp,pi,mygamma251,mygamma21,lacz_gamma
    use modfields, only : exnf,qt0,svm,qvsl,tmp0,ql0,esl,qvsl,rhof,exnf,svm, svp

    implicit none
    integer :: i, j, k, iaer, ispec, imod 
    real    :: evapt, evapn, evapr, &
               Dn, Dm, Dc, Fn, Fm, &
               frac, corr 
    real, parameter :: sigma = 1.5
    real, dimension(nspec-1) :: evapm
     
    real, allocatable :: F(:,:,:),S(:,:,:),G(:,:,:)
    integer :: numel

    allocate( F(2-ih:i1+ih,2-jh:j1+jh,k1)   & ! ventilation factor
             ,S(2-ih:i1+ih,2-jh:j1+jh,k1)   & ! super or undersaturation
             ,G(2-ih:i1+ih,2-jh:j1+jh,k1) )   ! cond/evap rate of a drop

    evap(2:i1,2:j1,1:k1)  = 0.
    Nevap(2:i1,2:j1,1:k1) = 0.

    do j=2,j1
    do i=2,i1
    do k=1,k1
      if (qrmask(i,j,k)) then
        S(i,j,k) = min(0.,(qt0(i,j,k)-ql0(i,j,k))/qvsl(i,j,k)- 1.)
        G(i,j,k) = (Rv * tmp0(i,j,k)) / (Dv*esl(i,j,k)) + rlv/(Kt*tmp0(i,j,k))*(rlv/(Rv*tmp0(i,j,k)) -1.)
        G(i,j,k) = 1./G(i,j,k)
      endif
    enddo
    enddo
    enddo

    if ( l_sb ) then
       do j=2,j1
       do i=2,i1
       do k=1,k1
         if (qrmask(i,j,k)) then
           numel=nint(mur(i,j,k)*100.)
           F(i,j,k) = avf * mygamma21(numel)*Dvr(i,j,k) +  &
              bvf*Sc_num**(1./3.)*(a_tvsb/nu_a)**0.5*mygamma251(numel)*Dvr(i,j,k)**(3./2.) * &
              (1.-(1./2.)  *(b_tvsb/a_tvsb)    *(lbdr(i,j,k)/(   c_tvsb+lbdr(i,j,k)))**(mur(i,j,k)+2.5) &
                 -(1./8.)  *(b_tvsb/a_tvsb)**2.*(lbdr(i,j,k)/(2.*c_tvsb+lbdr(i,j,k)))**(mur(i,j,k)+2.5) &
                 -(1./16.) *(b_tvsb/a_tvsb)**3.*(lbdr(i,j,k)/(3.*c_tvsb+lbdr(i,j,k)))**(mur(i,j,k)+2.5) &
                 -(5./128.)*(b_tvsb/a_tvsb)**4.*(lbdr(i,j,k)/(4.*c_tvsb+lbdr(i,j,k)))**(mur(i,j,k)+2.5) )
! *lbdr(i,j,k)**(mur(i,j,k)+1.)/f_gamma_1(i,j,k) factor moved to F
            evap(i,j,k) = 2*pi*Nr(i,j,k)*G(i,j,k)*F(i,j,k)*S(i,j,k)/rhof(k)
            Nevap(i,j,k) = c_Nevap*evap(i,j,k)*rhof(k)/xr(i,j,k)
         endif
       enddo
       enddo
       enddo

    else
       do j=2,j1
       do i=2,i1
       do k=1,k1
        if (qrmask(i,j,k)) then
           evap(i,j,k) = c_evapkk*2*pi*Dvr(i,j,k)*G(i,j,k)*S(i,j,k)*Nr(i,j,k)/rhof(k)
           Nevap(i,j,k) = evap(i,j,k)*rhof(k)/xr(i,j,k)
        endif
       enddo
       enddo
       enddo

    endif

    do j=2,j1
    do i=2,i1
    do k=1,k1
       if (evap(i,j,k) < -svm(i,j,k,iqr)/delt .and. qrmask(i,j,k)) then
          Nevap(i,j,k) = - svm(i,j,k,inr+iaer_offset)/delt
          evap (i,j,k) = - svm(i,j,k,iqr)/delt
       endif
    enddo
    enddo
    enddo

    do j=2,j1
    do i=2,i1
    do k=1,k1

       !Apply tendencies for water and associated heat
       qrp    (i,j,k) = qrp    (i,j,k) + evap(i,j,k)
       qtpmcr (i,j,k) = qtpmcr (i,j,k) - evap(i,j,k)
       thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*evap(i,j,k)

       evapn = 0. 
       evapm = 0.
       evapr = 0.

       do iaer = 1,naer

          evapt = 0.      

          imod  =  nmod_type(iaer)
          ispec = nspec_type(iaer)                

          if (imod == 9) then
             ! Calculate the evaporated amount of aerosol
             if ( qrmask(i,j,k) ) then
                if ( ispec == 1 ) then
                   evapt = -Nevap(i,j,k)  
                else

                ! Resuspended aerosol fraction =/= evaprated rain water
                ! fraction, we apply a correction factor based on the work of
                ! Gong et al., 2006
                !
                ! Resuspended fraction aerosol is based on the evaporated
                ! fraction rainwater in the COMPLETE timestep, i.e. evap/qr*delt
                ! Afterwards, when applying the correction factor, we divide 
                ! by delt again to get the correct resuspension fraction per second.

                   frac = min(-evap(i,j,k)/qr(i,j,k)*delt,1.)
                   corr = (1. - exp(-2*sqrt(frac))*(1. + 2.*sqrt(frac) + 2.*frac + 4./3.*frac**(3./2.)) )*(1.-frac) + frac**2.
                   evapt = corr*frac*aer_conc(i,j,k,iaer)/delt

                endif !ispec
             else
                ! Remove trailing in-rain aerosol when there is no rain in gridbox
                evapt = aer_conc(i,j,k,iaer)/delt
             endif !qrmask      

             evapt = min(evapt,max(0.,svm(i,j,k,iaer+iaer_offset)/delt + min(0.,svp(i,j,k,iaer+iaer_offset)) + min(0.,aer_tend(i,j,k,iaer))))   

             ! Reduce source tracer, i.e. in-rain aerosol
             aer_evpr(i,j,k,iaer) = - evapt    
             aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) - evapt

             ! Gather evaporated tracer values
             if ( ispec == 1 ) then
                evapn = max(0.,evapt)
             else
                evapm(ispec-1) = max(0.,evapt)
             end if

          endif !imod == 9, i.e. rain   
       enddo !iaer

       if ( SUM(evapm) > 0. ) then
          ! Complete evaporated tracer value calculation
          evapr = SUM(evapm) / SUM(evapm/spec_rho)

          ! Calculate number/mass median diameter for evaporated aerosol
          Dn = 1.e6 * ( (6.*SUM(evapm))/(pi*(evapn+eps0)*scalefac*evapr) )**(1./3.) * exp(-(3./2.)*log(sigma)**2.) !in micron
          Dm = Dn * exp(3.*log(sigma)**2.)
          Dc = 1. !Critical radius 0.5 micron, i.e. diameter 1 micron

          ! Calculate partitioning accumulation/coarse aerosol mode
          Fn = .5*erfc(- log(Dc/Dn)/( sqrt(2.)*log(sigma)) )
          Fm = .5*erfc(- log(Dc/Dm)/( sqrt(2.)*log(sigma)) )

          do iaer = 1,naer
  
             imod = nmod_type(iaer)
             ispec = nspec_type(iaer)

             if (imod == 3) then     !ACS mode
  
                if (ispec == 1) then !Number
                   aer_evpr(i,j,k,iaer) =                        Fn*evapn
                   aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) + Fn*evapn
                else                 !Mass
                   aer_evpr(i,j,k,iaer) =                        Fm*evapm(ispec-1)
                   aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) + Fm*evapm(ispec-1)
                endif !ispec

             elseif (imod == 4) then !COS mode

                if (ispec == 1) then !Number
                   aer_evpr(i,j,k,iaer) =                        (1.-Fn)*evapn      
                   aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) + (1.-Fn)*evapn
                else                 !Mass
                   aer_evpr(i,j,k,iaer) =                        (1.-Fm)*evapm(ispec-1)      
                   aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) + (1.-Fm)*evapm(ispec-1)
                end if !ispec

              endif !imod
          enddo !iaer

       end if !evapn > 0.

    enddo !k
    enddo !i
    enddo !j

    deallocate (F,S,G)


  end subroutine evaporation

  !*********************************************************************
  !*********************************************************************

  real function sed_flux(Nin,Din,sig2,Ddiv,nnn)

  !*********************************************************************
  ! Function to calculate numerically the analytical solution of the
  ! sedimentation flux between Dmin and Dmax based on
  ! Feingold et al 1986 eq 17 -20.
  ! fall velocity is determined by alfa* D^beta with alfa+ beta taken as
  ! specified in Rogers and Yau 1989 Note here we work in D and in SI
  ! (in Roger+Yau in cm units + radius)
  ! flux is multiplied outside sed_flux with 1/rho_air to get proper
  ! kg/kg m/s units
  !
  ! M.C. van Zanten    August 2005
  !*********************************************************************
    use modglobal, only : pi,rhow
    implicit none

    real, intent(in) :: Nin, Din, sig2, Ddiv
    integer, intent(in) :: nnn
    !para. def. lognormal DSD (sig2 = ln^2 sigma_g), D sep. droplets from drops
    !,power of of D in integral

    real, parameter ::   C = rhow*pi/6.     &
                        ,D_intmin = 1e-6    &
                        ,D_intmax = 4.3e-3

    real ::  alfa         & ! constant in fall velocity relation
            ,beta         & ! power in fall vel. rel.
            ,D_min        & ! min integration limit
            ,D_max        & ! max integration limit
            ,flux           ![kg m^-2 s^-1]

    integer :: k

    flux = 0.0

    if (Din < Ddiv) then
      alfa = 3.e5*100  ![1/ms]
      beta = 2
      D_min = D_intmin
      D_max = Ddiv
      flux = C*Nin*alfa*erfint(beta,Din,D_min,D_max,sig2,nnn)
    else
      do k = 1,3
        select case(k)
        case(1)        ! fall speed ~ D^2
          alfa = 3.e5*100 ![1/m 1/s]
          beta = 2
          D_min = Ddiv
          D_max = 133e-6
        case(2)        ! fall speed ~ D
          alfa = 4e3     ![1/s]
          beta = 1
          D_min = 133e-6
          D_max = 1.25e-3
        case default         ! fall speed ~ sqrt(D)
          alfa = 1.4e3 *0.1  ![m^.5 1/s]
          beta = .5
          D_min = 1.25e-3
          D_max = D_intmax
        end select
        flux = flux + C*Nin*alfa*erfint(beta,Din,D_min,D_max,sig2,nnn)
      end do
    end if
      sed_flux = flux
  end function sed_flux

  !*********************************************************************
  !*********************************************************************

  real function liq_cont(Nin,Din,sig2,Ddiv,nnn)

  !*********************************************************************
  ! Function to calculate numerically the analytical solution of the
  ! liq. water content between Dmin and Dmax based on
  ! Feingold et al 1986 eq 17 -20.
  !
  ! M.C. van Zanten    September 2005
  !*********************************************************************
    use modglobal, only : pi,rhow
    implicit none

    real, intent(in) :: Nin, Din, sig2, Ddiv
    integer, intent(in) :: nnn
    !para. def. lognormal DSD (sig2 = ln^2 sigma_g), D sep. droplets from drops
    !,power of of D in integral

    real, parameter :: beta = 0           &
                      ,C = pi/6.*rhow     &
                      ,D_intmin = 80e-6    &   ! value of start of rain D
                      ,D_intmax = 3e-3         !4.3e-3    !  value is now max value for sqrt fall speed rel.

    real ::  D_min        & ! min integration limit
            ,D_max          ! max integration limit

    if (Din < Ddiv) then
    D_min = D_intmin
    D_max = Ddiv
    else
    D_min = Ddiv
    D_max = D_intmax
    end if

    liq_cont = C*Nin*erfint(beta,Din,D_min,D_max,sig2,nnn)

  end function liq_cont

  !*********************************************************************
  !*********************************************************************

  real function erfint(beta, D, D_min, D_max, sig2,nnn )

  !*********************************************************************
  ! Function to calculate erf(x) approximated by a polynomial as
  ! specified in 7.1.27 in Abramowitz and Stegun
  ! NB phi(x) = 0.5(erf(0.707107*x)+1) but 1 disappears by substraction
  !
  !*********************************************************************
    implicit none
    real, intent(in) :: beta, D, D_min, D_max, sig2
    integer, intent(in) :: nnn

    real, parameter :: eps = 1e-10       &
                      ,a1 = 0.278393    & !a1 till a4 constants in polynomial fit to the error
                      ,a2 = 0.230389    & !function 7.1.27 in Abramowitz and Stegun
                      ,a3 = 0.000972    &
                      ,a4 = 0.078108
    real :: nn, ymin, ymax, erfymin, erfymax, D_inv

    D_inv = 1./(eps + D)
    nn = beta + nnn

    ymin = 0.707107*(log(D_min*D_inv) - nn*sig2)/(sqrt(sig2))
    ymax = 0.707107*(log(D_max*D_inv) - nn*sig2)/(sqrt(sig2))

    erfymin = 1.-1./((1.+a1*abs(ymin) + a2*abs(ymin)**2 + a3*abs(ymin)**3 +a4*abs(ymin)**4)**4)
    erfymax = 1.-1./((1.+a1*abs(ymax) + a2*abs(ymax)**2 + a3*abs(ymax)**3 +a4*abs(ymax)**4)**4)
    if (ymin < 0.) then
      erfymin = -1.*erfymin
    end if
    if (ymax < 0.) then
      erfymax = -1.*erfymax
    end if
    erfint = D**nn*exp(0.5*nn**2*sig2)*0.5*(erfymax-erfymin)
    if (erfint < 0.) erfint = 0.
  end function erfint

  function binarysearch(length, array, value, delta)
     ! Given an array and a value, returns the index of the element
     ! that
     ! is closest to, but less than, the given value.
     ! Uses a binary search algorithm.
     ! "delta" is the tolerance used to determine if two values are
     ! equal
     ! if ( abs(x1 - x2) <= delta) then
     !    assume x1 = x2
     ! endif

     implicit none
     integer, intent(in) :: length
     real, dimension(length), intent(in) :: array
     real, intent(in) :: value
     real, intent(in), optional :: delta

     integer :: binarysearch

     integer :: left, middle, right
     real :: d

     if (present(delta) .eqv. .true.) then
         d = delta
     else
         d = 1e-9
     endif

     left = 1
     right = length
     do
         if (left > right) then
             exit
         endif
         middle = nint((left+right) / 2.0)
         if ( abs(array(middle) - value) <= d) then
             binarySearch = middle
             return
         else if (array(middle) > value) then
             right = middle - 1
         else
             left = middle + 1
         end if
     end do
     binarysearch = right

    end function binarysearch

    real function lookup_interpolate(x_len, x_array, y_len, y_array, f, x, y, delta)
        ! This function uses bilinear interpolation to estimate the
        ! value of a function f at point (x,y)
        ! f is assumed to be sampled on a regular grid, with the grid x
        ! values specified by x_array and the grid y values specified by y_array
        ! Reference: http://en.wikipedia.org/wiki/Bilinear_interpolation
        implicit none
        integer, intent(in) :: x_len, y_len
        real, dimension(x_len), intent(in) :: x_array
        real, dimension(y_len), intent(in) :: y_array
        real, dimension(x_len, y_len), intent(in) :: f
        real, intent(in) :: x,y
        real, intent(in), optional :: delta

        real :: denom, x1, x2, y1, y2
        integer :: i,j

        i = binarysearch(x_len, x_array, x)
        j = binarysearch(y_len, y_array, y)

        x1 = x_array(i)
        x2 = x_array(i+1)

        y1 = y_array(j)
        y2 = y_array(j+1)

        denom = (x2 - x1)*(y2 - y1)

        lookup_interpolate = (f(i,j)*(x2-x)*(y2-y) + f(i+1,j)*(x-x1)*(y2-y) + &
            f(i,j+1)*(x2-x)*(y-y1) + f(i+1, j+1)*(x-x1)*(y-y1))/denom

    end function lookup_interpolate

    real function xinverfc(x)
        ! This function calculates the value of the 
        ! inverse complementary errorfunction using the 
        ! MKL library defined function vserfcinv designed for arrays
        
        implicit none
        include 'mkl_vml.f90'
        integer, parameter :: n=1
        real, intent(in)   :: x
        real a(1),v(1)
        
        a(1) = x        
        call vderfcinv(n, a, v)
        xinverfc = v(1)        

    end function xinverfc
 
end module modbulkmicro


