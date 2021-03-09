!> \file modtestbed.f90
!!  Testbed continuous forcing & nudging
!>

!>
!!  Testbed continuous forcing & nudging
!>
!!  \author Roel Neggers, IGMK
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



module modtestbed

use netcdf

implicit none
PRIVATE
PUBLIC :: inittestbed, testbednudge, exittestbed, ltestbed,testbed_getinttime, ntnudge, nknudge, &
          tb_time,tb_ps,tb_qts,tb_thls,tb_wqs,tb_wts, tb_z0h, tb_z0m, tb_alb, tb_Qnet, &
          tb_u,tb_v,tb_w,tb_thl,tb_qt,tb_ug,tb_vg, &
          tb_dqtdxls,tb_dqtdyls, &
          tb_qtadv,tb_thladv,tb_uadv,tb_vadv, &
          tb_tsoilav,tb_phiwav, &
          tbrad_p, tbrad_ql, tbrad_qv, tbrad_t, tbrad_o3
SAVE
  real, dimension(:,:), allocatable :: tnudge,tb_u,tb_v,tb_w,tb_thl,tb_qt,tb_ug,tb_vg, &
                                       tb_dqtdxls,tb_dqtdyls, &
                                       tb_qtadv,tb_thladv,tb_uadv,tb_vadv, &
                                       tb_tsoilav,tb_phiwav, &
                                       tbrad_p, tbrad_t, tbrad_qv, tbrad_ql, tbrad_o3
  real, dimension(:)  , allocatable :: tb_time, tb_ps, tb_qts, tb_thls, tb_wqs, tb_wts, tb_z0h, tb_z0m, tb_alb, tb_Qnet
  real :: tb_taunudge = 10800.
  logical :: ltestbed = .false., &
             ltb_nudge = .false., &
             ltb_u,ltb_v,ltb_w,ltb_thl,ltb_qt
  integer :: nknudge,ntnudge

contains
  subroutine inittestbed
    use modmpi,   only :myid,mpierr,comm3d,mpi_logical,mpi_integer &
                       , D_MPI_BCAST
    use modglobal,only :ifnamopt,fname_options,k1,&
                        grav,rd,cp,pref0,rlv,zf,checknamelisterror
    use modsurfdata,only : ksoilmax
    use modforces, only : lforce_user

    implicit none

    real, dimension(:,:), allocatable :: dumomega,dumqv,dumql,dumqi,dumt,dumpf, dumo3,&
                                         dumheight,dumqt,dumthl,dumu,dumv,dumw, &
                                         dumug,dumvg,dumqtadv,dumthladv,dumuadv,dumvadv, &
                                         dumqadv,dumladv,dumiadv,dumtadv, &
                                         dumtsoilav,dumphiwav,dumswi,&
                                         dumlwnet,dumswnet

    real, dimension(:), allocatable :: dumheights

    real :: dumphifc,dumphiwp

    INTEGER  NCID, STATUS, VARID, timID
    INTEGER start2(2), count2(2)
    character(len = nf90_max_name) :: RecordDimName

    integer :: ierr,i,k,ik,nknudgep1,nknudges
    real tv,rho,iexner,fac

    namelist /NAMTESTBED/ &
       ltestbed, ltb_nudge, tb_taunudge

    if(myid==0)then

      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMTESTBED,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMTESTBED')
      write(6 ,NAMTESTBED)
      close(ifnamopt)

    end if
 
    call D_MPI_BCAST(ltestbed     , 1,0,comm3d,mpierr)
    call D_MPI_BCAST(ltb_nudge    , 1,0,comm3d,mpierr)
    
    if (.not. ltestbed) return
    
    lforce_user = .true.

    if(myid==0) then
        !--- open nc file ---
        STATUS = NF90_OPEN('scm_in.nc', nf90_nowrite, NCID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
          
        !--- get time & height dimensions ---
        status = nf90_inq_dimid(ncid, "time", timID)
        if (status /= nf90_noerr) call handle_err(status)
        status = nf90_inquire_dimension(NCID, timID, len=ntnudge, name=RecordDimName)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

!        write(6,'(a15,i10," ",a10)') 'scm_in time:',ntnudge,RecordDimName

        status = nf90_inq_dimid(ncid, "nlev", timID)
        if (status /= nf90_noerr) call handle_err(status)
        status = nf90_inquire_dimension(NCID, timID, len=nknudge, name=RecordDimName)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

        status = nf90_inq_dimid(ncid, "nlevp1", timID)
        if (status /= nf90_noerr) call handle_err(status)
        status = nf90_inquire_dimension(NCID, timID, len=nknudgep1, name=RecordDimName)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

        status = nf90_inq_dimid(ncid, "nlevs", timID)
        if (status /= nf90_noerr) call handle_err(status)
        status = nf90_inquire_dimension(NCID, timID, len=nknudges, name=RecordDimName)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
    end if

    call D_MPI_BCAST(ntnudge    , 1,0,comm3d,mpierr)
    call D_MPI_BCAST(nknudge    , 1,0,comm3d,mpierr)
    call D_MPI_BCAST(nknudgep1  , 1,0,comm3d,mpierr)
    call D_MPI_BCAST(nknudges   , 1,0,comm3d,mpierr)

    !--- allocate space for input variables & reset---
    allocate(    tnudge    (ntnudge,k1), &
                 tb_u      (ntnudge,k1), &
                 tb_v      (ntnudge,k1), &
                 tb_w      (ntnudge,k1), &
                 tb_thl    (ntnudge,k1), &
                 tb_qt     (ntnudge,k1), &
                 tb_ug     (ntnudge,k1), &
                 tb_vg     (ntnudge,k1), &
                 tb_dqtdxls(ntnudge,k1), &
                 tb_dqtdyls(ntnudge,k1), &
                 tb_qtadv  (ntnudge,k1), &
                 tb_thladv (ntnudge,k1), &
                 tb_uadv   (ntnudge,k1), &
                 tb_vadv   (ntnudge,k1), &
                 tb_time   (ntnudge), &
                 tb_ps     (ntnudge), &
                 tb_qts    (ntnudge), &
                 tb_thls   (ntnudge), &
                 tb_wts    (ntnudge), &
                 tb_wqs    (ntnudge), &
                 tb_z0m    (ntnudge), &
                 tb_z0h    (ntnudge), &
                 tb_alb    (ntnudge), &
                 tb_Qnet   (ntnudge), &
                 tb_tsoilav(ntnudge,ksoilmax), &
                 tb_phiwav (ntnudge,ksoilmax), &
                 tbrad_p    (ntnudge, nknudge), &
                 tbrad_t    (ntnudge, nknudge), &
                 tbrad_qv   (ntnudge, nknudge), &
                 tbrad_ql   (ntnudge, nknudge), &
                 tbrad_o3   (ntnudge, nknudge) &
                 )

     tnudge = tb_taunudge     !nudging timescale

        tb_time=0
        tb_ps=0
        tb_qts=0
        tb_thls=0
        tb_wts=0
        tb_wqs=0
        tb_z0m=0
        tb_z0h=0
        tb_alb=0
        tb_Qnet=0

        tb_u=0
        tb_v=0
        tb_w=0
        tb_thl=0
        tb_qt=0
        tb_ug=0
        tb_vg=0
        tb_dqtdxls=0
        tb_dqtdyls=0
        tb_qtadv=0
        tb_thladv=0
        tb_uadv=0
        tb_vadv=0

        tb_tsoilav=0
        tb_phiwav=0

    
    if(myid==0) then

        allocate(dumomega (nknudge,ntnudge), &
                 dumheight(nknudge,ntnudge), &
                 dumpf    (nknudge,ntnudge), &
                 dumqv    (nknudge,ntnudge), &
                 dumql    (nknudge,ntnudge), &
                 dumqi    (nknudge,ntnudge), &
                 dumo3    (nknudge,ntnudge), &
                 dumt     (nknudge,ntnudge), &
                 dumqt    (nknudge,ntnudge), &
                 dumthl   (nknudge,ntnudge), &
                 dumu     (nknudge,ntnudge), &
                 dumv     (nknudge,ntnudge), &
                 dumw     (nknudge,ntnudge), &
                 dumug    (nknudge,ntnudge), &
                 dumvg    (nknudge,ntnudge), &
                 dumqtadv (nknudge,ntnudge), &
                 dumthladv(nknudge,ntnudge), &
                 dumuadv  (nknudge,ntnudge), &
                 dumvadv  (nknudge,ntnudge), &
                 dumqadv  (nknudge,ntnudge), &
                 dumladv  (nknudge,ntnudge), &
                 dumiadv  (nknudge,ntnudge), & 
                 dumtadv  (nknudge,ntnudge), & 
                 dumswnet  (nknudgep1,ntnudge), & 
                 dumlwnet  (nknudgep1,ntnudge), & 
                 dumheights (nknudges), & 
                 dumtsoilav (nknudges,ntnudge), & 
                 dumphiwav  (nknudges,ntnudge), & 
                 dumswi     (nknudges,ntnudge) &
                 )


        !--- timeseries ---

        !  time
        STATUS = NF90_INQ_VARID(NCID, 'time', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_time, start=(/1/), count=(/ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 
!        write(6,'(a30,5f10.2)') 'inittestbed: tb_time:',&
!             tb_time(1),tb_time(2),tb_time(3),tb_time(ntnudge-1),tb_time(ntnudge)

        !  surface pressure
        STATUS = NF90_INQ_VARID(NCID, 'ps', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_ps, start=(/1/), count=(/ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 

        !  surface temperature
        STATUS = NF90_INQ_VARID(NCID, 't_skin', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_thls, start=(/1/), count=(/ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 

        !  surface humidity
!        STATUS = NF90_INQ_VARID(NCID, '', VARID)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
!        STATUS = NF90_GET_VAR (NCID, VARID, tb_qts, start=(/1/), count=(/ntnudge/) )
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 

        !  surface T flux
        STATUS = NF90_INQ_VARID(NCID, 'sfc_sens_flx', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_wts, start=(/1/), count=(/ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 

        !  surface q flux
        STATUS = NF90_INQ_VARID(NCID, 'sfc_lat_flx', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_wqs, start=(/1/), count=(/ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 

        do i=1,ntnudge

          rho = tb_ps(i) / (rd * tb_thls(i))
          tb_wts(i) = -tb_wts(i) / (cp  * rho)        !Change sign: upward = positive in LES, but by convention upward = negative in most GCMs.
          tb_wqs(i) = -tb_wqs(i) / (rlv * rho)

          iexner = (tb_ps(i)/pref0)**(-rd/cp)
          tb_thls(i) = iexner * tb_thls(i)

        end do

        !  roughness length for momentum
        STATUS = NF90_INQ_VARID(NCID, 'mom_rough', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_z0m, start=(/1/), count=(/ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 
          
        !  roughness length for heat and moisture
        STATUS = NF90_INQ_VARID(NCID, 'heat_rough', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_z0h, start=(/1/), count=(/ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 
          
        !  surface albedo, for radiation
        STATUS = NF90_INQ_VARID(NCID, 'albedo', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_alb, start=(/1/), count=(/ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 
          
       

        !--- profiles full levels ---
        start2 = (/ 1      , 1       /)
        count2 = (/ nknudge, ntnudge /)
          
        !  height
        STATUS = NF90_INQ_VARID(NCID, 'height_f', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumheight, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
!        write(6,'(a30,91f10.2)') 'inittestbed: heightnudge:',&
!     &       ( dumheight(k,1),k=1,nknudge )
          
        !  u
        STATUS = NF90_INQ_VARID(NCID, 'u', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumu, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
        !  v
        STATUS = NF90_INQ_VARID(NCID, 'v', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumv, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
        !  qt
        STATUS = NF90_INQ_VARID(NCID, 'q', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumqv, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        
        STATUS = NF90_INQ_VARID(NCID, 'ql', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumql, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        
        STATUS = NF90_INQ_VARID(NCID, 'qi', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumqi, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        
        do i=1,ntnudge
        do k=1,nknudge
          dumqt(k,i) = dumqv(k,i) + dumql(k,i) + dumqi(k,i)
        enddo
        enddo       
          
        !  thl
        STATUS = NF90_INQ_VARID(NCID, 't', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumt, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

        STATUS = NF90_INQ_VARID(NCID, 'pressure_f', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumpf, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        
        do i=1,ntnudge
        do k=1,nknudge
          iexner = (dumpf(k,i)/pref0)**(-rd/cp)
          dumthl(k,i) = dumt(k,i) * iexner - (rlv * dumql(k,i)) / cp
        enddo
        enddo

        !  Ozone
        STATUS = NF90_INQ_VARID(NCID, 'o3', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumo3, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

        !  w
        STATUS = NF90_INQ_VARID(NCID, 'omega', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumomega, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
                
        do i=1,ntnudge
        do k=1,nknudge
          tv  = dumt(k,i) * (1.+0.61*dumqv(k,i))
          rho = dumpf(k,i) / (rd*tv)
          dumw(k,i) = - dumomega(k,i) / ( rho * grav )     !convert from Pa/s to m/s
        enddo
        enddo

        !  ug
        STATUS = NF90_INQ_VARID(NCID, 'ug', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumug, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
        !  vg
        STATUS = NF90_INQ_VARID(NCID, 'vg', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumvg, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)          
         
        !  uadv
        STATUS = NF90_INQ_VARID(NCID, 'uadv', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumuadv, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
        !  vadv
        STATUS = NF90_INQ_VARID(NCID, 'vadv', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumvadv, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
        !  qtadv
        STATUS = NF90_INQ_VARID(NCID, 'qadv', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumqadv, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        
        STATUS = NF90_INQ_VARID(NCID, 'ladv', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumladv, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        
        STATUS = NF90_INQ_VARID(NCID, 'iadv', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumiadv, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        
        do i=1,ntnudge
        do k=1,nknudge
          dumqtadv(k,i) = dumqadv(k,i) + dumladv(k,i) + dumiadv(k,i)
        enddo
        enddo
          
        !  thladv
        STATUS = NF90_INQ_VARID(NCID, 'tadv', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumtadv, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        
        do i=1,ntnudge
        do k=1,nknudge
          iexner = (dumpf(k,i)/pref0)**(-rd/cp)
          dumthladv(k,i) = dumtadv(k,i) * iexner
        enddo
        enddo
          
       

        !--- profiles half levels ---
        start2 = (/ 1        , 1       /)
        count2 = (/ nknudgep1, ntnudge /)
          
        !  net SW downward flux
        STATUS = NF90_INQ_VARID(NCID, 'fradSWnet', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumswnet, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
        !  net LW downward flux
        STATUS = NF90_INQ_VARID(NCID, 'fradLWnet', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumlwnet, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
        do i=1,ntnudge
          tb_Qnet(i) = dumswnet(nknudgep1,i) + dumlwnet(nknudgep1,i)      !flux at surface is stored in lowest half level of profile
!          write(6,*) "modtestbed: qnet:",i,tb_Qnet(i),nknudge,nknudgep1
        enddo


        !--- soil profiles ---

        STATUS = NF90_INQ_VARID(NCID, 'h_soil', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumheights, start=(/1/), count=(/nknudges/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 

        start2 = (/ 1       , 1       /)
        count2 = (/ nknudges, ntnudge /)

        !  tsoilav
        STATUS = NF90_INQ_VARID(NCID, 't_soil', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumtsoilav, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

        !  phiwav
        STATUS = NF90_INQ_VARID(NCID, 'q_soil', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumphiwav, start=start2, count=count2)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        
        !  field capacity
        status = nf90_inquire_attribute(ncid, nf90_global, "field_capacity")
        if (status /= nf90_noerr) call handle_err(status)
        status = nf90_get_att(ncid, nf90_global, "field_capacity", dumphifc)
        if (status /= nf90_noerr) call handle_err(status)
        
        !  wilting point
        status = nf90_inquire_attribute(ncid, nf90_global, "wilting_point")
        if (status /= nf90_noerr) call handle_err(status)
        status = nf90_get_att(ncid, nf90_global, "wilting_point", dumphiwp)
        if (status /= nf90_noerr) call handle_err(status)

        dumswi = ( dumphiwav - dumphiwp ) / ( dumphifc - dumphiwp )   !soil wetness index, using input values for wilting point and field capacity
          

         
        !--- close nc file ---
        STATUS = NF90_CLOSE(NCID)  
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        


        do i=1,ntnudge


          !--- interpolate towards LES levels, reverse height-order, switch dimensions ---
          ik = nknudge
          do k=1,k1
            
            do while( zf(k).gt.dumheight(ik,i) .and. ik.gt.1)
              ik=ik-1
            enddo
            if ( ik.lt.nknudge ) then 
              ik=ik+1  
            endif
          
            fac = ( zf(k)-dumheight(ik,i) ) / ( dumheight(ik-1,i)-dumheight(ik,i) )
        
            tb_thl    (i,k)   = dumthl    (ik,i) + fac * ( dumthl    (ik-1,i)    - dumthl    (ik,i)    )
            tb_qt     (i,k)   = dumqt     (ik,i) + fac * ( dumqt     (ik-1,i)    - dumqt     (ik,i)    )
            tb_u      (i,k)   = dumu      (ik,i) + fac * ( dumu      (ik-1,i)    - dumu      (ik,i)    )
            tb_v      (i,k)   = dumv      (ik,i) + fac * ( dumv      (ik-1,i)    - dumv      (ik,i)    )
            tb_w      (i,k)   = dumw      (ik,i) + fac * ( dumw      (ik-1,i)    - dumw      (ik,i)    )
            tb_ug     (i,k)   = dumug     (ik,i) + fac * ( dumug     (ik-1,i)    - dumug     (ik,i)    )
            tb_vg     (i,k)   = dumvg     (ik,i) + fac * ( dumvg     (ik-1,i)    - dumvg     (ik,i)    )
            tb_uadv   (i,k)   = dumuadv   (ik,i) + fac * ( dumuadv   (ik-1,i)    - dumuadv   (ik,i)    )
            tb_vadv   (i,k)   = dumvadv   (ik,i) + fac * ( dumvadv   (ik-1,i)    - dumvadv   (ik,i)    )
            tb_qtadv  (i,k)   = dumqtadv  (ik,i) + fac * ( dumqtadv  (ik-1,i)    - dumqtadv  (ik,i)    )
            tb_thladv (i,k)   = dumthladv (ik,i) + fac * ( dumthladv (ik-1,i)    - dumthladv (ik,i)    )

            !if (i.eq.1) write(6,*)  k, zf(k), " : ", ik, dumheight(ik,i), ik-1, dumheight(ik-1,i)

          enddo
          

          !--- soil & surface properties ---
          tb_qts(i) = tb_qt(i,1)    !qts seems not really used anymore (see subr. timedepsurf in modtimedep.f90)
          
!          if (i.eq.1) then
!            do k=1,nknudges
!              write(6,*)  "modtestbed: soil: ", k, dumheights(k), dumtsoilav(k,i), dumphiwav(k,i)
!            end do
!          end if

          !dzsoil

          !tb_tsoilav(i,:) = dumtsoilav(:,i)

          !tb_phiwav (i,:) = phiwp + dumswi(:,i) * (phifc - phiwp )     !scale soil moisture using field capacity and wilting point
          do k = 1, nknudge
            tbrad_p(i,k)  = dumpf(k,i)
            tbrad_t(i,k)  = dumt(k,i)
            tbrad_qv(i,k) = dumqv(k,i)
            tbrad_ql(i,k) = dumql(k,i) + dumqi(k,i) 
            tbrad_o3(i,k) = dumo3(k,i)
          end do

        enddo

        
        !--- clean-up ---
        deallocate(dumomega)
        deallocate(dumqv)
        deallocate(dumql)
        deallocate(dumqi)
        deallocate(dumt)
        deallocate(dumpf)
        deallocate(dumo3)

        deallocate(dumheight)
        deallocate(dumqt)
        deallocate(dumthl)
        deallocate(dumu)
        deallocate(dumv)
        deallocate(dumw)
        deallocate(dumug)
        deallocate(dumvg)

        deallocate(dumuadv)
        deallocate(dumvadv)
        deallocate(dumqtadv)
        deallocate(dumthladv)
        deallocate(dumqadv)
        deallocate(dumladv)
        deallocate(dumiadv)
        deallocate(dumtadv)
        deallocate(dumswnet)
        deallocate(dumlwnet)

        deallocate(dumheights)
        deallocate(dumtsoilav)
        deallocate(dumphiwav)
        deallocate(dumswi)


        !--- do some output to screen ---
!        do i=1,2
!        !do i=1,ntnudge
!
!        write(6,'(a20,f10.2,a15,3f10.2)') 'modtestbed: scm_in time:',tb_time(i),' sfc pressure:',tb_ps(i),tb_thls(i),tb_qts(i)
!
!        write(6,*) ' zf       tnudge    tb_u    tb_v    tb_w    tb_thl    tb_qt    tb_ug    tb_vg'
!        do k=kmax,1,-1
!          write (6,'(f7.1,8e12.4)') &
!                zf          (k), &
!                tnudge      (i,k), &
!                tb_u        (i,k), &
!                tb_v        (i,k), &
!                tb_w        (i,k), &
!                tb_thl      (i,k), &
!                tb_qt       (i,k), &
!                tb_ug       (i,k), &
!                tb_vg       (i,k)
!        end do
!
!        end do

    end if

    call D_MPI_BCAST(ntnudge    , 1,0,comm3d,mpierr)

    call D_MPI_BCAST(tb_time    ,ntnudge   ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_ps      ,ntnudge   ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_qts     ,ntnudge   ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_thls    ,ntnudge   ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_wts     ,ntnudge   ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_wqs     ,ntnudge   ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_z0h     ,ntnudge   ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_z0m     ,ntnudge   ,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_Qnet    ,ntnudge   ,0,comm3d,mpierr)

    call D_MPI_BCAST(tnudge     ,ntnudge*k1,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_u       ,ntnudge*k1,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_v       ,ntnudge*k1,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_w       ,ntnudge*k1,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_thl     ,ntnudge*k1,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_qt      ,ntnudge*k1,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_ug      ,ntnudge*k1,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_vg      ,ntnudge*k1,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_uadv    ,ntnudge*k1,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_vadv    ,ntnudge*k1,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_qtadv   ,ntnudge*k1,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_thladv  ,ntnudge*k1,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_uadv    ,ntnudge*k1,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_vadv    ,ntnudge*k1,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_dqtdxls ,ntnudge*k1,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_dqtdyls ,ntnudge*k1,0,comm3d,mpierr)

    call D_MPI_BCAST(tb_tsoilav ,ntnudge*ksoilmax,0,comm3d,mpierr)
    call D_MPI_BCAST(tb_phiwav  ,ntnudge*ksoilmax,0,comm3d,mpierr)

    call D_MPI_BCAST(tbrad_p      ,ntnudge*nknudge,0,comm3d,mpierr)
    call D_MPI_BCAST(tbrad_qv     ,ntnudge*nknudge,0,comm3d,mpierr)
    call D_MPI_BCAST(tbrad_ql     ,ntnudge*nknudge,0,comm3d,mpierr)
    call D_MPI_BCAST(tbrad_t      ,ntnudge*nknudge,0,comm3d,mpierr)
    call D_MPI_BCAST(tbrad_o3     ,ntnudge*nknudge,0,comm3d,mpierr)

    ltb_u   = any(abs(tb_u)>1e-8)
    ltb_v   = any(abs(tb_v)>1e-8)
    ltb_w   = any(abs(tb_w)>1e-8)
    ltb_thl = any(abs(tb_thl)>1e-8)
    ltb_qt  = any(abs(tb_qt)>1e-8)


  end subroutine inittestbed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine testbednudge
    use modglobal, only : timee,rtimee,i1,j1,kmax,rdt
    use modfields, only : up,vp,wp,thlp, qtp,u0av,v0av,qt0av,thl0av
    implicit none

    integer k,t
    real :: dtm,dtp,currtnudge, qttnudge,qtthres

    if (.not.(ltestbed .and. ltb_nudge)) return

    if (timee==0) return

    t=1
    do while(rtimee>tb_time(t))
      t=t+1
    end do
    if (rtimee>tb_time(1)) then
      t=t-1
    end if

    dtm = ( rtimee-tb_time(t) ) / ( tb_time(t+1)-tb_time(t) )
    dtp = ( tb_time(t+1)-rtimee)/ ( tb_time(t+1)-tb_time(t) )

    qtthres = 1e-6
    do k=1,kmax

      currtnudge = max(rdt,tnudge(t,k)*dtp+tnudge(t+1,k)*dtm)

      if (ltb_u)   up(2:i1,2:j1,k) = up(2:i1,2:j1,k)     - &
          ( u0av(k)   - (tb_u(t,k)  *dtp + tb_u(t+1,k)  *dtm) ) / currtnudge

      if (ltb_v)   vp(2:i1,2:j1,k) = vp(2:i1,2:j1,k)     - &
          ( v0av(k)   - (tb_v(t,k)  *dtp + tb_v(t+1,k)  *dtm) ) / currtnudge

      if (ltb_w)   wp(2:i1,2:j1,k) = wp(2:i1,2:j1,k)     - &
          (           - (tb_w(t,k)  *dtp + tb_w(t+1,k)  *dtm) ) / currtnudge

      if (ltb_thl) thlp(2:i1,2:j1,k) = thlp(2:i1,2:j1,k) - &
          ( thl0av(k) - (tb_thl(t,k)*dtp + tb_thl(t+1,k)*dtm) ) / currtnudge

      if (ltb_qt)  then
        if (qt0av(k)< qtthres) then
          qttnudge = rdt
        else
          qttnudge = currtnudge
        end if
          qtp(2:i1,2:j1,k) = qtp(2:i1,2:j1,k)   - &
            ( qt0av(k)  - (tb_qt(t,k) *dtp + tb_qt(t+1,k) *dtm) ) / qttnudge
      end if
    end do

    !write(6,*) 'testbednudge:', rtimee, t, tb_time(t), tb_time(t+1), currtnudge, dtm, dtp, qt0av (1),tb_qt (t,1),tb_qt (t+1,1)
    !write(6,*) 'testbednudge:', rtimee, t, tb_time(t), tb_time(t+1), currtnudge, dtm, dtp, qt0av (kmax),tb_qt (t,kmax),tb_qt (t+1,kmax)
    
    !write(6,*) 'testbednudge:', rtimee, t, tb_time(t), tb_time(t+1), currtnudge, dtm, dtp, thl0av(1),tb_thl(t,1),tb_thl(t+1,1)
    !write(6,*) 'testbednudge:', rtimee, t, tb_time(t), tb_time(t+1), currtnudge, dtm, dtp, thl0av(kmax),tb_thl(t,kmax),tb_thl(t+1,kmax)

  end subroutine testbednudge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine testbed_getinttime(t, dtm, dtp)
     use modglobal, only : rtimee
!     use modfields, only : up,vp,wp,thlp, qtp,u0av,v0av,qt0av,thl0av
!     use modmpi,    only : myid
    implicit none
    integer, intent(out) :: t
    real, intent(out)    :: dtm, dtp


    t=1
    do while(rtimee>tb_time(t))
      t=t+1
    end do
    if (rtimee>tb_time(1)) then
      t=t-1
    end if

    dtm = ( rtimee-tb_time(t) ) / ( tb_time(t+1)-tb_time(t) )
    dtp = ( tb_time(t+1)-rtimee)/ ( tb_time(t+1)-tb_time(t) )


  end subroutine testbed_getinttime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exittestbed
  if (allocated(tb_time)) then
    deallocate(tb_time)
  end if
  end subroutine exittestbed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
  subroutine handle_err(errcode)
      
  implicit none

  integer errcode
     
  write(6,*) 'Error: ', nf90_strerror(errcode)
  stop 2
      
  end subroutine handle_err


end module modtestbed
