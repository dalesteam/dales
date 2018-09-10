!> \file modbulkmicro.f90

!>
!!  Bulk microphysics.
!>
!! Calculates bulk microphysics using a two moment scheme.
!! \see   Seifert and Beheng (Atm. Res., 2001)
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
module modbulkmicro

!   Amount of liquid water is splitted into cloud water and precipitable
!   water (as such it is a two moment scheme). Cloud droplet number conc. is
!   fixed in place and time.
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

    allocate( aer_conc (2-ih:i1+ih,2-jh:j1+jh,k1,naer) &
             ,aer_tend (2-ih:i1+ih,2-jh:j1+jh,k1,naer) &
             ,sedcm    (2-ih:i1+ih,2-jh:j1+jh,k1,naer))      

!Individual tendencies
    allocate( aer_evpc (2-ih:i1+ih,2-jh:j1+jh,k1,naer))
    allocate( aer_acti (2-ih:i1+ih,2-jh:j1+jh,k1,naer))
    allocate( aer_scvc (2-ih:i1+ih,2-jh:j1+jh,k1,naer))
    allocate( aer_slfc (2-ih:i1+ih,2-jh:j1+jh,k1,naer))
!
  gamma25=lacz_gamma(2.5)
  gamma3=2.
  gamma35=lacz_gamma(3.5)

!>>>MdB
  
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
!<<<MdB  
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

    deallocate(aer_evpc,aer_acti,aer_scvc,aer_slfc)

  end subroutine exitbulkmicro

!> Calculates the microphysical source term.
  subroutine bulkmicro
    use modglobal, only : i1,j1,k1,rdt,rk3step,timee,rlv,cp
    use modfields, only : sv0,svm,svp,qtp,thlp,ql0,exnf,rhof
    use modbulkmicrostat, only : bulkmicrotend
    use modmpi,    only : myid
    implicit none
    integer :: i,j,k, taer, iaer, ispec, imod
    real :: qrtest,nrtest,nctest,mctest,aertest,aertest_hold

    do j=2,j1
    do i=2,i1
    do k=1,k1

      qr(i,j,k) = sv0(i,j,k,iqr) ! kg tracer / kg air
      
      ! For clarity use duplicate arrays for Nc and Nr
      Nc(i,j,k) = sv0(i,j,k,inc+iaer_offset)!For now assume tracer is in #/ms anyway*rhof(k) ! Cloud drops [ #/m3 ]
      Nr(i,j,k) = sv0(i,j,k,inr+iaer_offset)!Idem *rhof(k) ! Rain drops  [ #/m3 ]
    
      do iaer=1,naer
         aer_conc(i,j,k,iaer) = sv0(i,j,k,iaer+iaer_offset)! Idem *rhof(k) ! kg/m3 OR #/m3
      enddo  

    enddo
    enddo
    enddo

    Nrp     = 0.0 
    qrp     = 0.0
    thlpmcr = 0.0
    qtpmcr  = 0.0

    aer_tend = 0.0     
  
    aer_evpc = 0.0
    aer_acti = 0.0
    aer_scvc = 0.0
    aer_slfc = 0.0
 
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
!TODO If qrmask is TRUE, then Nr should automatically be >0 since this is a
!condition for qrmask
            if (qrmask(i,j,k).and.Nr(i,j,k).gt.0.) then
              xr  (i,j,k) = rhof(k)*qr(i,j,k)/(Nr(i,j,k)+eps0) ! JvdD Added eps0 to avoid floating point exception
              xr  (i,j,k) = min(xr(i,j,k),xrmaxkk) ! to ensure x_pw is within borders
              Dvr (i,j,k) = (xr(i,j,k)/pirhow)**(1./3.)
            endif
         enddo
         enddo
         enddo

      endif ! l_sb
    endif   ! l_rain
    !write (6,*) 'parameters calculated'

  !*********************************************************************
  ! call microphysical processes subroutines
  !*********************************************************************

    call bulkmicrotend
    call cloudactivation
    call bulkmicrotend
    call cloudselfcollection

    if (l_sedc) then 
      call bulkmicrotend
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
    endif    

    call bulkmicrotend
    call scavenging  
    call bulkmicrotend
    call cloudevaporation
    call bulkmicrotend

    do k=1,k1
    do j=2,j1
    do i=2,i1

      ! NOTE ------------------------------------------------------------------
      ! Before transferring the tendencies calculated here, values have to be
      ! TRANSFORMED to the right UNITS and to be CHECKED for mass conservation.
      ! Problematic processes for mass conservation are:
      ! (1) In-cloud scavenging, because the efficiencies in the lookup table
      !     can be > 1.0.
      ! (2) Cloud evaporation, because this process removes all should
      !     completely empty the in-cloud reservoir, possibly colliding with
      !     non-zero fluxes from i.e. advection
      ! -----------------------------------------------------------------------

      ! Transform to svp units (e.g. [kg m-3 s-1]/[kga m-3] = [kg kga-1 s-1])  
!      aer_tend(i,j,k,:) = aer_tend(i,j,k,:)! For now assume tracer is in #/ms anyway /rhof(k)
!      aer_evpc(i,j,k,:) = aer_evpc(i,j,k,:)! Idem /rhof(k)

      do iaer = 1,naer

         ispec = nspec_type(iaer)
         imod  = nmod_type(iaer)

         aertest = svm(i,j,k,iaer+iaer_offset) + (svp(i,j,k,iaer+iaer_offset) + aer_tend(i,j,k,iaer))*delt

         if ( aertest < 0. .and. aer_tend(i,j,k,iaer) < 0.) then ! .and. aer_tend(i,j,k,iaer) < 0. .and. aer_tend(i,j,k,iaer) < aertest ) then

            if ( imod <= 7 ) then ! Checking and correction for activation
 
!               write(6,"(A7,I3,A4,E10.3,A4,E10.3,A4,E10.3,A4,E10.3,A4,E10.3,A4,E10.3,A4,E10.3,A4,E10.3)") &
!               'AERFREE', iaer, &
!               'sum', svm(i,j,k,iaer+iaer_offset) + svp(i,j,k,iaer+iaer_offset)*delt + aer_tend(i,j,k,iaer)*delt, &
!               'svm', svm     (i,j,k,iaer+iaer_offset), &
!               'cnc', aer_conc(i,j,k,iaer), &
!               'svp', svp     (i,j,k,iaer+iaer_offset)*delt, &
!               'tnd', aer_tend(i,j,k,iaer)*delt, &
!               'act', aer_acti(i,j,k,iaer)*delt, &
!               'scv', aer_scvc(i,j,k,iaer)*delt, &
!               'evp', aer_evpc(i,j,k,iaer)*delt

               select case(ispec)
               case(1) 
                  taer = inc
               case(2)
                  taer = iso4cld
               case(3)
                  taer = ibccld
               case(4)
                  taer = ipomcld
               case(5)
                  taer = isscld
               case(6)
                  taer = iducld
               end select  
               
               !In-cloud scavenging 
               if (aer_scvc(i,j,k,iaer) > 0.) then
                  aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) + aer_scvc(i,j,k,iaer)
                  if (ispec /= 1 ) then; aer_tend(i,j,k,taer) = aer_tend(i,j,k,taer) - aer_scvc(i,j,k,iaer); endif

                  aertest_hold = aertest + aer_scvc(i,j,k,iaer)*delt
                  aer_scvc(i,j,k,iaer) = max(0., aer_scvc(i,j,k,iaer) + aertest/delt)
                  aertest = aertest_hold - aer_scvc(i,j,k,iaer)*delt

                  aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) - aer_scvc(i,j,k,iaer)
                  if (ispec /= 1 ) then; aer_tend(i,j,k,taer) = aer_tend(i,j,k,taer) + aer_scvc(i,j,k,iaer); endif
               endif 

               !Activation  
               if (aer_acti(i,j,k,iaer) > 0.) then
                  aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) + aer_acti(i,j,k,iaer)
                  aer_tend(i,j,k,taer) = aer_tend(i,j,k,taer) - aer_acti(i,j,k,iaer)

                  aertest_hold = aertest + aer_acti(i,j,k,iaer)*delt
                  aer_acti(i,j,k,iaer) = max(0., aer_acti(i,j,k,iaer) + aertest/delt)
                  aertest = aertest_hold - aer_acti(i,j,k,iaer)*delt

                  aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) - aer_acti(i,j,k,iaer)
                  aer_tend(i,j,k,taer) = aer_tend(i,j,k,taer) + aer_acti(i,j,k,iaer)
               endif 

!               if (aertest < -eps0) then 
!               write(6,"(A7,I3,A4,E10.3,A4,E10.3,A4,E10.3,A4,E10.3,A4,E10.3,A4,E10.3,A4,E10.3,A4,E10.3)") &
!               'AERFRE1', iaer, &
!               'sum', svm(i,j,k,iaer+iaer_offset) + svp(i,j,k,iaer+iaer_offset)*delt + aer_tend(i,j,k,iaer)*delt, &
!               'svm', svm     (i,j,k,iaer+iaer_offset), &
!               'cnc', aer_conc(i,j,k,iaer), &
!               'svp', svp     (i,j,k,iaer+iaer_offset)*delt, &
!               'tnd', aer_tend(i,j,k,iaer)*delt, &
!               'act', aer_acti(i,j,k,iaer)*delt, &
!               'scv', aer_scvc(i,j,k,iaer)*delt, &
!               'evp', aer_evpc(i,j,k,iaer)*delt
!               endif
 
!TODO Now aerosol number and mass are checked and corrected independently. This
!could cause situations with one being zero, while some of the other remains.
!This might cause problems.

            elseif ( imod == 8 ) then ! Checking and correction for cloud evporation
               
               select case ( ispec ) ! Select free aerosol target variable
               case(1) ! Number
                  taer = iacs_n
               case(2) ! SO4
                  taer = iso4acs
               case(3) ! BC
                  taer = ibcacs
               case(4) ! POM
                  taer = ipomacs
               case(5) ! SS
                  taer = issacs
               case(6) ! DUST
                  taer = iduacs
               end select

               !Cloud evaporation
               if (aer_evpc(i,j,k,iaer) > 0.) then
                  aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) + aer_evpc(i,j,k,iaer)
                  aer_tend(i,j,k,taer) = aer_tend(i,j,k,taer) - aer_evpc(i,j,k,iaer)

                  aertest_hold = aertest + aer_evpc(i,j,k,iaer)*delt
                  aer_evpc(i,j,k,iaer) = max(0., aer_evpc(i,j,k,iaer) + aertest/delt)
                  aertest = aertest_hold - aer_evpc(i,j,k,iaer)*delt

                  aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) - aer_evpc(i,j,k,iaer)
                  aer_tend(i,j,k,taer) = aer_tend(i,j,k,taer) + aer_evpc(i,j,k,iaer)
               endif

               !Cloud selfcollection
               if (aer_slfc(i,j,k,iaer) > 0.) then
                  aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) + aer_slfc(i,j,k,iaer)
                  aer_tend(i,j,k,taer) = aer_tend(i,j,k,taer) - aer_slfc(i,j,k,iaer)

                  aertest_hold = aertest + aer_slfc(i,j,k,iaer)*delt
                  aer_slfc(i,j,k,iaer) = max(0., aer_slfc(i,j,k,iaer) + aertest/delt)
                  aertest = aertest_hold - aer_slfc(i,j,k,iaer)*delt

                  aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) - aer_slfc(i,j,k,iaer)
                  aer_tend(i,j,k,taer) = aer_tend(i,j,k,taer) + aer_slfc(i,j,k,iaer)
               endif

!               if (aertest < -eps0) then
!               write(6,"(A7,I3,A4,E10.3, A4,E10.3,A4,E10.3,A4,E10.3,A4,E10.3,A4,E10.3,A4,E10.3,A4,E10.3,A4,E10.3,A4,E10.3)") &
!               'AERCLD1', iaer, &
!               'tst', aertest, & 
!               'sum', svm(i,j,k,iaer+iaer_offset) + svp(i,j,k,iaer+iaer_offset)*delt + aer_tend(i,j,k,iaer)*delt, &
!               'svm', svm     (i,j,k,iaer+iaer_offset), &
!               'cnc', aer_conc(i,j,k,iaer), &
!               'svp', svp     (i,j,k,iaer+iaer_offset)*delt, &
!               'tnd', aer_tend(i,j,k,iaer)*delt, &
!               'act', aer_acti(i,j,k,iaer)*delt, &
!               'scv', aer_scvc(i,j,k,iaer)*delt, &
!               'slf', aer_slfc(i,j,k,iaer)*delt, &
!               'evp', aer_evpc(i,j,k,iaer)*delt
!               endif

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

        do iaer = 1,naer
           if (nmod_type(iaer) == 9) then

           select case(nspec_type(iaer))
           case(1) ! Number
              taer = iacs_n
              aertest = nrtest
           case(2) ! SO4
              taer = iso4acs
              aertest = svm(i,j,k,iaer+iaer_offset)+(svp(i,j,k,iaer+iaer_offset)+aer_tend(i,j,k,iaer))*delt
           case(3) ! BC
              taer = ibcacs
              aertest = svm(i,j,k,iaer+iaer_offset)+(svp(i,j,k,iaer+iaer_offset)+aer_tend(i,j,k,iaer))*delt
           case(4) ! POM
              taer = ipomacs
              aertest = svm(i,j,k,iaer+iaer_offset)+(svp(i,j,k,iaer+iaer_offset)+aer_tend(i,j,k,iaer))*delt
           case(5) ! SS
              taer = issacs
              aertest = svm(i,j,k,iaer+iaer_offset)+(svp(i,j,k,iaer+iaer_offset)+aer_tend(i,j,k,iaer))*delt
           case(6) ! DUST
              taer = iduacs
              aertest = svm(i,j,k,iaer+iaer_offset)+(svp(i,j,k,iaer+iaer_offset)+aer_tend(i,j,k,iaer))*delt
           end select
                
           aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) - max(aer_tend(i,j,k,iaer),aertest/delt)
           aer_tend(i,j,k,taer) = aer_tend(i,j,k,taer) + max(aer_tend(i,j,k,iaer),aertest/delt)

           end if
        enddo

      else

         ! Transferring tendencies 
         svp(i,j,k,iqr)     = svp(i,j,k,iqr)     + qrp(i,j,k)
         thlp(i,j,k)        = thlp(i,j,k)        + thlpmcr(i,j,k)
         qtp(i,j,k)         = qtp(i,j,k)         + qtpmcr(i,j,k)
      end if

!TODO Check if in NO case the inbalance is greater than flux.
!If this would be the cause, applying the correction would cause a flux in the
!wrong direction. 
      do iaer= 1,naer
        svp(i,j,k,iaer+iaer_offset) = svp(i,j,k,iaer+iaer_offset) + aer_tend(i,j,k,iaer)        
      end do    

    Nrp(i,j,k) = aer_tend(i,j,k,inr) !For statistical purposes    

    enddo
    enddo
    enddo
  end subroutine bulkmicro

  ! -----------------------------------------------------------------------
  !
  ! SUBROUTINE CLOUDACTIVATION
  ! 
  ! Calculates the amount of activated aerosols that form cloud droplets.
  ! Routine is based on 'standard' Kohler thoery and applies this to each
  ! individual mode. Procedure is as follows:
  ! 1. Calculation of mode total mass and average density, molar mass and
  ! hygrospocity
  ! 2. Calculation of activated number/mass fraction per mode
  ! 3. Determine source and target variables and apply activation
  !
  ! -----------------------------------------------------------------------

  subroutine cloudactivation
    use modglobal, only : i1,j1,k1, rlv, cp, rhow
    use modfields, only : ql0, thl0, exnf, w0, svm, svp

    implicit none
    integer         :: i,j,k, taer, iaer, ispec, imod
    real            :: A, B, T, sigma, r_crit, dn_mean, dm_mean, fac
    real            :: d_ais, f_ais, Nt, dNcdt, acti 

    real, dimension(nmod-2) :: mod_numb, mod_mass, mod_rho, mod_mol, mod_k, fn, fm
 
    real, parameter ::    & 
        pi    = 3.14159,  & !
        R     = 8.1314,   & ! Gas constant
        !TODO Change fixed supersaturation (Sc) to updraft dependant?
        Sc    = 0.2,      & ! Assumed supersaturation (%)
        Mw    = 0.018,    & ! Molar mass of water (kg mol-1)
        sw    = 0.075       ! Surface tension of water 

        !Activation based on k-Kohler theory. Equation for r_crit derived from
        !Eq. (10) of Petters & Kreidenweis (2007). k values used for species
        !taken from Pringle et al. (2010):
        !SO4 0.88,   BC 0,   POM 0.1,   SS 1.28,   DU: 0.

    do j=2,j1
    do i=2,i1
    do k=1,k1    

       if ( qcmask(i,j,k) ) then

          mod_numb = 0.0
          mod_mass = 0.0
          mod_rho  = 0.0   
          mod_mol  = 0.0
          mod_k    = 0.0

          fn = 0.0
          fm = 0.0
   
          if (l_kohler) then

          if ( Nc(i,j,k) < ncmin ) then

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
                   mod_mol (imod) = mod_mol (imod) + aer_conc(i,j,k,iaer)/spec_mol(ispec-1)
                   if (spec_k(ispec-1) > 0.) then
                      mod_k(imod) = mod_k   (imod) + aer_conc(i,j,k,iaer)/spec_k  (ispec-1)
                   endif
                end if
             end do   

             T = thl0(i,j,k) * exnf(k) + (rlv/cp) * ql0(i,j,k)
             A = 2*sw*Mw/(R*T*rhow)   

             ! Calculate activated fraction --------------------------------
             do imod=1,nmod-2 ! Activation per mode, only free aerosol modes
                if (mod_mass(imod) < 1.e-20 .or. mod_numb(imod) < 1e2 .or. mod_k(imod) == 0.) cycle
                
                   mod_rho(imod) = mod_mass(imod)/mod_rho(imod)
                   mod_mol(imod) = mod_mass(imod)/mod_mol(imod)
                   mod_k (imod)  = mod_mass(imod)/mod_k (imod)             
  
                   r_crit = (4.*A**3./(27.*mod_k(imod)*log(Sc/100.+1.)**2.))**(1./3.) 

                   sigma = sigma_lognormal(imod)
                   ! Number mean diameter and activated number fraction
                   dn_mean  = (6.*mod_mass(imod)/(pi*mod_numb(imod)*mod_rho(imod)))**(1/3.)*exp(-3.*log(sigma)**2./2.)
                
                   !Double check!!! TODO     
                   fn(imod) = 1. - .5*erfc(-log(2*r_crit/dn_mean)/sqrt(2.)*log(sigma))

                   ! Mass mean diameter and activation mass fraction
                   dm_mean  = dn_mean * exp(3.*log(sigma)**2.)

                   !Double check!!! TODO     
                   fm(imod) = 1. - .5*erfc(-log(2.*r_crit/dm_mean)/sqrt(2.)*log(sigma))

             enddo !imod

          end if !ncmin

          else !l_kohler
            
             ! Aerosol activation based on Pousse-Nottelman et al., 2015 (PN15)
             ! Parameterisation parameterized according to the approach by Lin
             ! and Leaitch (1997) following the work by Muhlbauer and Lohmann
             ! (2008) and Zubler et al. (2011a)    
        
             ! 1. Determine number of aerosol particles >35nm (Nt)    
             ! 2. Calculate dNc/dt
             ! 3. Calcualte mass/number tendencies of affected tracers

             !Aggregate free aerosol modes
             do iaer = 1, naer !Aggregate modes

                imod  = nmod_type(iaer)
                ispec = nspec_type(iaer)

                if ( imod == 1 .or. imod > 4 ) cycle  ! Only AIS,ACS,COS
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
            
                d_ais  = (6.*mod_mass(iais_n)/(pi*mod_numb(iais_n)*mod_rho(iais_n)))**(1/3.)*exp(-3.*log(sigma)**2./2.)
                f_ais  = 1. - .5*erfc(-log(2.*r_crit/d_ais)/sqrt(2.)*log(sigma))
           
                Nt = 1e-6*(f_ais*mod_numb(iais_n) + mod_numb(iacs_n) + mod_numb(icos_n)) !Eq.(4) PN15
                dNcdt = max(1.e6/delt*( 0.1*(w0(i,j,k)*1e2*Nt/(w0(i,j,k)*1e2+.023*Nt))**1.27 - 1e-6*Nc(i,j,k)), 0. ) !Eq.(2) PN15

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
                      fm(imod) = .5*erfc(xinverfc(2*fn(imod))-3.*log(sigma_lognormal(imod))/sqrt(2.))
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
                case(1) ! Number
                   taer = inc
                   fac  = fn(imod)
                case(2) ! SO4
                   taer = iso4cld
                   fac  = fm(imod)
                case(3) ! BC
                   taer = ibccld
                   fac  = fm(imod)
                case(4) ! POM
                   taer = ipomcld
                   fac  = fm(imod)
                case(5) ! SS
                   taer = isscld
                   fac  = fm(imod)
                case(6) ! DUST
                   taer = iducld
                   fac  = fm(imod)
             end select

             acti = min(fac*aer_conc(i,j,k,iaer)/delt , svm(i,j,k,iaer+iaer_offset)/delt+svp(i,j,k,iaer+iaer_offset))   

             aer_tend(i,j,k,taer) = aer_tend(i,j,k,taer) + acti
             aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) - acti
             aer_acti(i,j,k,iaer) = acti !fac*aer_conc(i,j,k,iaer)/delt   
 
          enddo !iaer 

       end if !qcmask 
    end do
    end do
    end do

  end subroutine cloudactivation

  subroutine cloudevaporation  
    use modglobal, only : i1,j1,k1
    use modfields, only : ql0, svp, rhof, qlm, svm, sv0

    implicit none
    integer :: i,j,k, iaer, imod, ispec
    real    :: evap, evapn, evapr, &
               Dn, Dm, Dc, Fn, Fm, sigma, &
               frac, corr, aertest
    real, dimension(nspec-1) :: evapm

    !Distribution width for evaporated aerosol
    sigma = 1.5 !Following Pousse-Nottelman et al. (2015)    

    do j=2,j1
    do i=2,i1
    do k=1,k1
        
    evapn = 0.
    evapm = 0.
    evapr = 0.

    do iaer = 1,naer  

       evap = 0.
       imod  =  nmod_type(iaer)
       ispec = nspec_type(iaer)

       if (imod == 8) then 

          ! Calculate the resuspension fraction 
          if ( qcmask(i,j,k) ) then 
             if ( qlm(i,j,k) > ql0(i,j,k) ) then !partial evaporation

                frac = min((qlm(i,j,k)-ql0(i,j,k))/(qlm(i,j,k)), 1.) !min function to prevent overflow when qlm=0.

                if ( ispec /= 1 ) then
                   ! Resuspended aerosol fraction =/= evaporated rain water fraction, 
                   ! we apply a correction factor based on the work of Gong et al., 2006
                   corr = (1. - exp(-2*sqrt(frac))*(1. + 2.*sqrt(frac) + 2.*frac + 4./3.*frac**(3./2.)) )*(1.-frac) + frac**2.
                   if (corr*frac > 1.) then; write(6,"(I2, A30, F7.4, A8, F7.4)") iaer, 'Resuspended fraction, number:', frac, ',mass: ', frac*corr; end if
                   frac = corr*frac
                endif

             else !no evaporation
                frac = 0.
             endif !qlm>ql0

          else ! complete evaporation
             frac = 1.                 
          endif !qcmask

          !Calculate how much can be removed from clouds
          if (frac /= 0. .and. aer_conc(i,j,k,iaer) > eps0) then
             evap = frac/delt*min(aer_conc(i,j,k,iaer),svm(i,j,k,iaer+iaer_offset) + min(0.,svp(i,j,k,iaer+iaer_offset)*delt) + min(0.,aer_tend(i,j,k,iaer)*delt))
             evap = max(0., evap)

!             if (svm(i,j,k,iaer+iaer_offset) + min(0.,svp(i,j,k,iaer+iaer_offset)*delt) + min(0.,aer_tend(i,j,k,iaer)*delt) - evap*delt < 0. .and. evap > 0.) then 
!             write(6,"(A7,I3,A4,E11.4,A4,E11.4,A4,E11.4,A4,E11.4,A4,E11.4,A4,E11.4)") &
!             'INC AER', iaer, &
!             'sum', svm(i,j,k,iaer+iaer_offset) + svp(i,j,k,iaer+iaer_offset)*delt + aer_tend(i,j,k,iaer)*delt - evap*delt, &
!             'svm', svm     (i,j,k,iaer+iaer_offset), &
!             'cnc', aer_conc(i,j,k,iaer), &
!             'svp', svp     (i,j,k,iaer+iaer_offset)*delt, &
!             'tnd', aer_tend(i,j,k,iaer)*delt, &
!             'evp', evap*delt
!             endif
          endif
        
          ! Gather evaporated tracer values
          if (ispec == 1) then
             evapn = evap
          else
             evapm(ispec-1) = evap
          end if

        end if !nmod_type
     enddo !iaer

     if ( frac /= 0. .and. evapn > eps0 .and. sum(evapm) > eps0 ) then

        ! Complete evaporated tracer value calculation
        evapr = SUM(evapm) / SUM(evapm/spec_rho)

        ! Calculate number/mass median diameter for evaporated aerosol 
        Dn = 1e6 * ( (6*SUM(evapm))/(3.141592654*evapn*evapr) )**(1./3.) * exp(-(3/2)*log(sigma)**2) !in micron
        Dm = Dn * exp(3*log(sigma)**2)
        Dc = 1. !Critical radius 0.5 micron, i.e. diameter 1 micron

        ! Calculate partitioning accumulation/coarse aerosol mode
        Fn = .5*erfc(- log(Dc/Dn)/( sqrt(2.)*log(sigma)) )
        Fm = .5*erfc(- log(Dc/Dm)/( sqrt(2.)*log(sigma)) )
    
!        write(6,"(6E11.4, 2F7.4)") evapn, evapm(1), evapm(2), evapm(3), evapm(4), evapm(5), Fn, Fm
 
        do iaer = 1,naer

           imod = nmod_type(iaer)
           ispec = nspec_type(iaer)
                
           if (imod == 3) then     !ACS mode

              if (ispec == 1) then !Number
                 aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) + Fn*evapn
              else                 !Mass
                 aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) + Fm*evapm(ispec-1)
              endif !ispec

           elseif (imod == 4) then !COS mode

              if (ispec == 1) then !Number
                 aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) + (1.-Fn)*evapn
              else                 !Mass
                 aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) + (1.-Fm)*evapm(ispec-1)
              end if !ispec

           elseif (imod == 8) then !In-cloud
        
              if (ispec == 1) then !Number
                 aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) - evapn
                 aer_evpc(i,j,k,iaer) = evapn
              else                 !Mass
                 aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) - evapm(ispec-1)
                 aer_evpc(i,j,k,iaer) = evapm(ispec-1)
              end if !ispec
     
           endif !imod
       enddo !iaer
    endif !evapn > 0

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

       nuc = 1580.*rhof(k)*ql0(i,j,k) - 0.28                                                                           ! kg water/m3
       aer_tend(i,j,k,inc) = aer_tend(i,j,k,inc) - k_c*(nuc + 2.)/(nuc +1.)*(ql0(i,j,k)*rhof(k))**2. * 1.225 / rhof(k) ! #/m3/s
       aer_slfc(i,j,k,inc) = k_c*(nuc + 2.)/(nuc +1.)*(ql0(i,j,k)*rhof(k))**2. * 1.225 / rhof(k)
 
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
  
  use modglobal, only : i1,j1,k1,rhow
  use modfields, only : ql0,sv0,rhof 

  implicit none

  integer :: i,j,k, iaer, imod, ispec, taerinc, taerblc
  real    :: meancld, meanaer, meanaer_inc, rainrate, scavinc, scavblc,scvc

  real, dimension(nmod-2) :: mod_mass, mod_numb, mod_rho, &
                             scavinc_m, scavinc_n, scavblc_m, scavblc_n

  ! ---------------------------------------------------------------------------

  do i=2,i1
  do j=2,j1
  do k=1,k1

     meancld =  min( max( 1e6*(3.*ql0(i,j,k)*rhof(k)/(4.*3.1415*Nc(i,j,k)*rhow))**(1/3.) ,5.001), 49.999) ! in micron
     rainrate = min( max(sed_qr(i,j,k)*3600, 0.01001) , 99.999)                                           !kg m-2 s-1 --> mm/hr
  
!     if ( qrmask(i,j,k) ) then 
!        write(6,"(A9, E11.4,A7,E11.4,A2)") 'rainrate:', rainrate, 'raw:(', sed_qr(i,j,k)*3600, ')'        
!     end if   
 
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

        if ( (mod_mass(imod) > eps0) .and. (mod_numb(imod) > 0.)) then
   
           ! Get geometric mean aerosol radius
           meanaer = ((6.*mod_mass(imod))/(3.1415*mod_numb(imod)*mod_rho(imod)))**(1/3.)*exp(-3.*log(sigma_lognormal(imod))**2/2.)           ! in meter

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
           if ( qcmask(i,j,k) .and. ( Nc(i,j,k) > 0. )) then

               meanaer_inc = meanaer 
               if (meanaer < 1.0e-8) then
!                  write(6,"(A41, I3, A9, E11.4)") 'Mode mean aerosol radius too small. imod=', imod, 'meanaer=', meanaer
                  meanaer_inc = 1.0e-8
               endif  
                
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
           ! This leads to dC/dt      = C         * L or in the variables used in the code below: ! 
           !     aer_tend(i,j,k,iaer) = aer_conc(i,j,k,iaer) * scav                               !
           ! ==================================================================================== !
           if (l_rain) then
              if ( qrmask(i,j,k) .and. (sed_qr(i,j,k)*3600> 1.1e-2) ) then   

!                 write(6,"(A9, E11.4,A7,E11.4,A2)") 'rainrate:', rainrate, 'raw:(', sed_qr(i,j,k)*3600, ')'
  
                 meanaer  = max(1.001e-3,min(.9999e3,meanaer*1e6 ),meanaer*1e6 ) !in micron 

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

        if (imod > 7) cycle      ! Exclude cloud and rain modes

        select case ( ispec )    ! Select cloud target variable
               case(1) ! Number
                  scavinc = scavinc_n(imod)
                  scavblc = scavblc_n(imod)
                  taerinc = 0      
                  taerblc = 0
               case(2) ! SO4
                  taerinc = iso4cld
                  taerblc = iso4rai      
                  scavinc = scavinc_m(imod)
                  scavblc = scavblc_m(imod)
               case(3) ! BC
                  taerinc = ibccld
                  taerblc = ibcrai      
                  scavinc = scavinc_m(imod)
                  scavblc = scavblc_m(imod)
               case(4) ! POM
                  taerinc = ipomcld
                  taerblc = ipomrai      
                  scavinc = scavinc_m(imod)
                  scavblc = scavblc_m(imod)
               case(5) ! SS
                  taerinc = isscld
                  taerblc = issrai      
                  scavinc = scavinc_m(imod)
                  scavblc = scavblc_m(imod)
               case(6) ! DUST
                  taerinc = iducld
                  taerblc = idurai      
                  scavinc = scavinc_m(imod)
                  scavblc = scavblc_m(imod)
               end select

        scvc = scavinc*aer_conc(i,j,k,iaer)
        aer_scvc(i,j,k,iaer)       = scvc
        
        aer_tend(i,j,k,iaer)       = aer_tend(i,j,k,iaer)    - scvc !scavinc*aer_conc(i,j,k,iaer)
        aer_tend(i,j,k,iaer)       = aer_tend(i,j,k,iaer)    - scavblc*aer_conc(i,j,k,iaer)
        
        if (ispec /= 1 ) then 
           aer_tend(i,j,k,taerinc) = aer_tend(i,j,k,taerinc) + scvc !scavinc*aer_conc(i,j,k,iaer)
           aer_tend(i,j,k,taerblc) = aer_tend(i,j,k,taerblc) + scavblc*aer_conc(i,j,k,iaer)
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

    if (l_sb ) then
    !
    ! SB autoconversion
    !
      tau(2:i1,2:j1,1:k1) = 0.
      phi(2:i1,2:j1,1:k1) = 0.
      nuc(2:i1,2:j1,1:k1) = 0.
      k_au = k_c/(20*x_s)

      do j=2,j1
      do i=2,i1
      do k=1,k1
         if (qcmask(i,j,k) .and. ( Nc(i,j,k) > ncmin )) then
            nuc    (i,j,k) = 1.58*(rhof(k)*ql0(i,j,k)*1000.) +0.72-1. !G09a
!           nuc    (i,j,k) = 0. !
            xc     (i,j,k) = rhof(k) * ql0(i,j,k) / Nc(i,j,k) ! No eps0 necessary
            au     (i,j,k) = k_au * (nuc(i,j,k)+2.) * (nuc(i,j,k)+4.) / (nuc(i,j,k)+1.)**2.    &
                    * (ql0(i,j,k) * xc(i,j,k))**2. * 1.225 ! *rho**2/rho/rho (= 1)
            tau    (i,j,k) = 1.0 - ql0(i,j,k) / qltot(i,j,k)
            phi    (i,j,k) = k_1 * tau(i,j,k)**k_2 * (1.0 -tau(i,j,k)**k_2)**3
            au     (i,j,k) = au(i,j,k) * (1.0 + phi(i,j,k)/(1.0 -tau(i,j,k))**2)

            !Water    
            qtpmcr (i,j,k) = qtpmcr (i,j,k) - au(i,j,k)
            qrp    (i,j,k) = qrp    (i,j,k) + au(i,j,k)

            !Particle number
            aer_tend(i,j,k,inc) = aer_tend(i,j,k,inc) - au(i,j,k)/xc(i,j,k)*rhof(k) !Nc is in #/m3, so correct with rhof   
            aer_tend(i,j,k,inr) = aer_tend(i,j,k,inr) + au(i,j,k)/x_s*rhof(k)

            !Aerosol mass
            do iaer= 1,naer
               if (nmod_type(iaer) == 8 .and. nspec_type(iaer) > 1) then 
                  aer_tend(i,j,k,iaer)   = aer_tend(i,j,k,iaer)   - au(i,j,k)/ql0(i,j,k)*aer_conc(i,j,k,iaer) 
                  aer_tend(i,j,k,iaer+1) = aer_tend(i,j,k,iaer+1) + au(i,j,k)/ql0(i,j,k)*aer_conc(i,j,k,iaer) 
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
           aer_tend(i,j,k,inc) = aer_tend(i,j,k,inc) - ac(i,j,k)/xc(i,j,k)*rhof(k)
           !Number of raindrops does not change

           !Aerosol mass
           do iaer = 1,naer
              if (nmod_type(iaer) == 8 .and. nspec_type(iaer) > 1) then
                  aer_tend(i,j,k,iaer)   = aer_tend(i,j,k,iaer)   - ac(i,j,k)/ql0(i,j,k)*aer_conc(i,j,k,iaer)
                  aer_tend(i,j,k,iaer+1) = aer_tend(i,j,k,iaer+1) + ac(i,j,k)/ql0(i,j,k)*aer_conc(i,j,k,iaer)
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
    use modfields, only : rhof
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
              aer_in  = sed_qr(i,j,k+1)/qr_spl(i,j,k+1)*aer_spl(i,j,k+1,ispec)
           endif
           if (qr_spl(i,j,k  ) > 0.) then
              aer_out = sed_qr(i,j,k  )/qr_spl(i,j,k  )*aer_spl(i,j,k  ,ispec)
           endif

           aer_spl(i,j,k,ispec) = aer_spl(i,j,k,ispec) + (aer_in - aer_out)*dt_spl/(dzf(k)*rhof(k))
        enddo

        qr_spl(i,j,k) = qr_spl(i,j,k) + (sed_qr(i,j,k+1) - sed_qr(i,j,k))*dt_spl/(dzf(k)*rhof(k))

!        do iaer = 1,naer 
!           if (nspec_type(iaer) == 9 .and. nmod_type(iaer) > 1) then
!              aer_in = 0.
!              aer_out = 0.
!
!              if (qrmask(i,j,k+1)) then 
!                 aer_in =  sed_qr(i,j,k+1)/qr_spl(i,j,k+1)/(dzf(k+1)*rhof(k+1)) * aer_spl(i,j,k+1,nspec_type(iaer)-1)
!              end if
!              if (qrmask(i,j,k)) then 
!                 aer_out = sed_qr(i,j,k  )/qr_spl(i,j,k  )/(dzf(k  )*rhof(k  )) * aer_spl(i,j,k,  nspec_type(iaer)-1)
!              end if
!      
!              aer_spl(i,j,k,nspec_type(iaer)-1) = aer_spl(i,j,k,nspec_type(iaer)-1) + (aer_in - aer_out)*dt_spl
!
!           end if
!        enddo

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

!       do k = 1, kmax
!       do j = 2, j1
!       do i = 2, i1
!          do ispec = 1, nspec-1
!             aer_in = 0.
!             aer_out = 0.
!
!             if (qrmask(i,j,k+1)) then
!                aer_in = min(sed_qr(i,j,k+1)/qr_spl(i,j,k+1)/(dzf(k+1)*rhof(k+1)),1.) * aer_spl(i,j,k+1,ispec)
!                if 
!             end if      
!
!             if (qrmask(i,j,k  )) then
!                aer_in = min(sed_qr(i,j,k  )/qr_spl(i,j,k  )/(dzf(k  )*rhof(k  )),1.) * aer_spl(i,j,k  ,ispec)
!             end if
!
!             aer_spl(i,j,k,ispec) = aer_spl(i,j,k,ispec) + (aer_in - aer_out)*dt_spl
!          enddo    
!       enddo
!       enddo
!       enddo    

    enddo ! time splitting loop

    qrp     (2:i1,2:j1,1:k1)    =      qrp(2:i1,2:j1,1:k1)     + (qr_spl(2:i1,2:j1,1:k1) - qr(2:i1,2:j1,1:k1))/delt
    aer_tend(2:i1,2:j1,1:k1,inr)= aer_tend(2:i1,2:j1,1:k1,inr) + (Nr_spl(2:i1,2:j1,1:k1) - Nr(2:i1,2:j1,1:k1))/delt

    do iaer = 1,naer    
       if (nmod_type(iaer) == 9 .and. nspec_type(iaer) > 1) then
          aer_tend(2:i1,2:j1,1:k1,iaer)= aer_tend(2:i1,2:j1,1:k1,iaer) + (aer_spl(2:i1,2:j1,1:k1,nspec_type(iaer)-1) - aer_conc(2:i1,2:j1,1:k1,iaer))/delt
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
             aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) - evapt

             ! Gather evaporated tracer values
             if ( ispec == 1 ) then
                evapn = evapt
             else
                evapm(ispec-1) = max(0.,evapt)
             end if

          endif !imod == 9, i.e. rain   
       enddo !iaer

       if ( (evapn > 0.) .and. (SUM(evapm) > 0.) ) then
          ! Complete evaporated tracer value calculation
          evapr = SUM(evapm) / SUM(evapm/spec_rho)

          ! Calculate number/mass median diameter for evaporated aerosol
          Dn = 1e6 * ( (6*SUM(evapm))/(3.141592654*evapn*evapr) )**(1./3.) * exp(-(3/2)*log(sigma)**2) !in micron
          Dm = Dn * exp(3*log(sigma)**2)
          Dc = 1. !Critical radius 0.5 micron, i.e. diameter 1 micron

          ! Calculate partitioning accumulation/coarse aerosol mode
          Fn = .5*erfc(- log(Dc/Dn)/( sqrt(2.)*log(sigma)) )
          Fm = .5*erfc(- log(Dc/Dm)/( sqrt(2.)*log(sigma)) )

!          write(6,"(A16, 2F5.2, 2E9.2, A6, F5.0, A6, E9.2, A6, 5E9.2)") 'Evap Fn,Fm,Dn,Dm',Fn,Fm,Dn,Dm,'evapr',evapr,'evapn',evapn,'evapm',evapm(1),evapm(2),evapm(3),evapm(4),evapm(5)

          do iaer = 1,naer
  
             imod = nmod_type(iaer)
             ispec = nspec_type(iaer)

             if (imod == 3) then     !ACS mode
  
                if (ispec == 1) then !Number
                   aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) + Fn*evapn
                else                 !Mass
                   aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) + Fm*evapm(ispec-1)
                endif !ispec

             elseif (imod == 4) then !COS mode

                if (ispec == 1) then !Number
                   aer_tend(i,j,k,iaer) = aer_tend(i,j,k,iaer) + (1.-Fn)*evapn
                else                 !Mass
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
  !  if (erfint < 0.) write(*,*)'erfint neg'
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


