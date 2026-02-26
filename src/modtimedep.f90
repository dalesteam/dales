!> \file modtimedep.f90
!!  Prescribes surface values, fluxes and LS forcings at certain times

!>
!!  Prescribes surface values, fluxes and LS forcings at certain times
!>
!!  \author Roel Neggers, KNMI
!!  \author Thijs Heus,MPI-M
!!  \author Stephan de Roode, TU Delft
!!  \author Simon Axelsen, UU
!!  \par Revision list
!! \todo documentation
!  This file is part of DALES.
!
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

module modtimedep

  use modlogging, only: finish, warning
  use modprecision, only: field_r

implicit none
character(len=*), parameter :: modname = 'modtimedep'
private
public :: inittimedep, timedep,ltimedep,ltimedepuv,exittimedep

save
! switches for timedependent surface fluxes and large scale forcings
  logical       :: ltimedep     = .false. !< Overall switch, input in namoptions
  logical       :: ltimedepuv   = .false. !< Switch for time-dependent u,v forcings from ls_flux.inp
  logical       :: ltimedepz    = .true.  !< Switch for large scale forcings
  logical       :: ltimedepsurf = .true.  !< Switch for surface fluxes
  integer    :: kflux
  integer    :: kls

  real, allocatable     :: timeflux (:)   !< time points for surface fluxes [s]
  real, allocatable     :: wqsurft  (:)   !< time dependent kinematic moisture flux [kg/kg m/s]
  real, allocatable     :: wtsurft  (:)   !< time dependent kinematic temperature flux [K m/s]
  real, allocatable     :: thlst    (:)   !< time dependent surface thl [K]
  real, allocatable     :: qtst     (:)   !< time dependent surface qt [kg/kg]
  real, allocatable     :: pst      (:)   !< time dependent surface pressure [Pa]
  real, allocatable     :: Qnetavt  (:)   !< time dependent average net surface radiative energy flux [W/m^2]
  real, allocatable     :: timels  (:)    !< time points for large scale forcings (are the same for netcdf input as timeflux) [s]
  real, allocatable     :: ugt     (:,:)  !< time dependent geostrophic eastward wind [m/s]
  real, allocatable     :: vgt     (:,:)  !< time dependent geostrophic northward wind [m/s]
  real, allocatable     :: dpdxlt  (:,:)  !< time dependent large-scale eastward pressure gradient [Pa/m]
  real, allocatable     :: dpdylt  (:,:)  !< time dependent large-scale northward pressure gradient [Pa/m]
  real, allocatable     :: wflst   (:,:)  !< time dependent large-scale subsidence [m/s]
  real, allocatable     :: dqtdxlst(:,:)  !< time dependent eastward gradient of the qt due to large-scale forcing [kg/kg/m]
  real, allocatable     :: dqtdylst(:,:)  !< time dependent northward gradient of the qt due to large-scale forcing [kg/kg/m]
  real, allocatable     :: dqtdtlst(:,:)  !< time dependent tendency of the total water mixing ratio due to large-scale forcing [kg/kg/s]
  real, allocatable     :: dthldtlst(:,:) !< time dependent tendency of the liquid water potential temperature due to large-scale forcing [K/s]
  real, allocatable     :: thlpcart(:,:)  !< time dependent tendency of the liquid water potential temperature due to radiative forcing [K/s]
  real, allocatable     :: dudtlst (:,:)  !< time dependent tendency of the eastward velocity due to large-scale forcing [m/s^2]
  real, allocatable     :: dvdtlst (:,:)  !< time dependent tendency of the northward velocity due to large-scale forcing [m/s^2]
  real, allocatable     :: thlproft(:,:)  !< time dependent profile of the liquid water potential temperature [K]
  real, allocatable     :: qtproft (:,:)  !< time dependent profile of the total water mixing ratio [kg/kg]

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine inittimedep
    use modmpi,    only :myid,mpierr,comm3d,D_MPI_BCAST
    use modglobal, only :cexpnr,k1,kmax,ifinput,runtime,zf,ntimedep,lcoriol,iinput
    use modsurfdata,only :ps,qts,wqsurf,wtsurf,thls, Qnetav
    use modtimedepsv, only : inittimedepsv

    use modtestbed,        only : ltestbed,ntnudge,&
                                  tb_time,tb_ps,tb_qts,tb_thls,tb_wqs,tb_wts,&
                                  tb_w,tb_ug,tb_vg,&
                                  tb_uadv,tb_vadv,tb_qtadv,tb_thladv,tb_Qnet

    implicit none

    character(len=*), parameter :: routine = modname//'/inittimedep'


    character (80):: chmess
    character (1) :: chmess1
    integer :: k,t, ierr
    real :: dummyr
    real, allocatable, dimension (:) :: height
    if (.not. ltimedep) return

    if (ltestbed) then
      kflux = ntnudge
      kls   = ntnudge
    else
      kflux = ntimedep
      kls   = ntimedep
    end if

    allocate(height   (k1))

    allocate(timeflux (0:kflux))
    allocate(wqsurft  (kflux))
    allocate(wtsurft  (kflux))
    allocate(thlst    (kflux))
    allocate(qtst     (kflux))
    allocate(pst      (kflux))
    allocate(Qnetavt  (kflux))

    allocate(timels   (0:kls))
    allocate(ugt      (k1,kls))
    allocate(vgt      (k1,kls))
    allocate(dpdxlt   (k1,kls))
    allocate(dpdylt   (k1,kls))
    allocate(wflst    (k1,kls))

    allocate(dqtdxlst (k1,kls))
    allocate(dqtdylst (k1,kls))

    allocate(dqtdtlst (k1,kls))
    allocate(dthldtlst(k1,kls))
    allocate(dudtlst  (k1,kls))
    allocate(dvdtlst  (k1,kls))

    allocate(thlpcart (k1,kls))

    allocate(thlproft (k1,kls))
    allocate(qtproft  (k1,kls))

    timeflux = 0
    timels   = 0

    wqsurft  = wqsurf
    wtsurft  = wtsurf
    thlst    = thls
    qtst     = qts
    pst      = ps
    Qnetavt  = Qnetav

    ugt      = 0
    vgt      = 0
    dpdxlt   = 0
    dpdylt   = 0
    wflst    = 0

    dqtdxlst = 0
    dqtdylst = 0

    dqtdtlst = 0
    dthldtlst= 0
    dudtlst  = 0
    dvdtlst  = 0

    thlpcart = 0

    thlproft = 0
    qtproft  = 0

    if (myid==0) then

      !--- load lsforcings---

      timeflux = 0
      timels   = 0

      if (ltestbed) then

        write(*,*) 'inittimedep: testbed mode: data for time-dependent forcing obtained from scm_in.nc'

        timeflux(1:kflux) = tb_time
        timels  (1:kls  ) = tb_time

        pst      = tb_ps
        qtst     = tb_qts
        thlst    = tb_thls
        wqsurft  = tb_wqs
        wtsurft  = tb_wts
        Qnetavt  = tb_Qnet

        height  (:) = zf
        do t=1,kls
          ugt      (:,t) = tb_ug    (t,:)
          vgt      (:,t) = tb_vg    (t,:)
          wflst    (:,t) = tb_w     (t,:)
          dqtdxlst (:,t) = 0.
          dqtdylst (:,t) = 0.
          dqtdtlst (:,t) = tb_qtadv (t,:)
          dthldtlst(:,t) = tb_thladv(t,:)
          dudtlst  (:,t) = tb_uadv  (t,:)
          dvdtlst  (:,t) = tb_vadv  (t,:)
        end do

      else
        if (iinput == 2) then
        call init_timedep_from_netcdf('forcings.'//cexpnr//'.nc', height, timeflux,kflux,kmax)
        if (size(timeflux,dim=1) < ntimedep) then
          call finish(routine, "Number of time points in forcings."//cexpnr//".nc is smaller than ntimedep = ", ntimedep)
        end if
        if(timeflux(1)>runtime) then
          call warning(routine,'Time dependent forcings do not change before end of simulation. Disabling time dependent large scale forcings and surface fluxes')
          ltimedepsurf=.false.
          ltimedepz=.false.
          endif
        else
          open(ifinput,file='ls_flux.inp.'//cexpnr)
          read(ifinput,'(a80)') chmess
          write(6,*) chmess
          read(ifinput,'(a80)') chmess
          write(6,*) chmess
          read(ifinput,'(a80)') chmess
          write(6,*) chmess

          timeflux = 0
          timels   = 0


        !--- load fluxes---
        t    = 0
        ierr = 0
        do while (timeflux(t) < runtime)
          t=t+1
          if (t > kflux) then
             call finish(routine, "Too many time points in file ", 'ls_flux.inp.'//cexpnr, ", the limit is kflux = ", kflux)
          end if
          read(ifinput,*, iostat = ierr) timeflux(t), wtsurft(t), wqsurft(t),thlst(t),qtst(t),pst(t)
          write(*,'(i8,6e12.4)') t,timeflux(t), wtsurft(t), wqsurft(t),thlst(t),qtst(t),pst(t)
          if (ierr < 0) then
            call finish(routine, 'STOP: No time dependend data for end of run (surface fluxes)')
          end if
        end do
        if(timeflux(1)>runtime) then
         write(6,*) 'Time dependent surface variables do not change before end of'
         write(6,*) 'simulation. --> only large scale forcings'
         ltimedepsurf=.false.
        endif
        ! flush to the end of fluxlist
        do while (ierr ==0)
          read (ifinput,*,iostat=ierr) dummyr
        end do
        backspace (ifinput)


        !---load large scale forcings----
        t = 0
        do while (timels(t) < runtime)
          t = t + 1
          if (t > kls) then
             call finish(routine, "Too many time points in file ", 'nudge.inp.'//cexpnr, ", the limit is kls = ", kls)
          end if
          chmess1 = "#"
          ierr = 1 ! not zero
          do while (.not.(chmess1 == "#" .and. ierr ==0)) !search for the next line consisting of "# time", from there onwards the profiles will be read
            read(ifinput,*,iostat=ierr) chmess1,timels(t)
            if (ierr < 0) then
              call finish(routine, 'STOP: No time dependend data for end of run')
            end if
          end do


          if (ltimedepuv) then
             ! new, optional format with u,v in ls_flux.inp.*
             do k=1,kmax
                read (ifinput,*) &
                     height  (k)  , &
                     ugt     (k,t), &
                     vgt     (k,t), &
                     wflst   (k,t), &
                     dqtdxlst(k,t), &
                     dqtdylst(k,t), &
                     dqtdtlst(k,t), &
                     thlpcart(k,t), &
                     dudtlst (k,t), &
                     dvdtlst (k,t)
             end do
          else
            ! if lcoriol, read in 2nd and 3rd column as ug and vg
            if (lcoriol) then
              ! old format without u,v in ls_flux.inp.*  (default)
              do k=1,kmax
                read (ifinput,*) &
                      height  (k)  , &
                      ugt     (k,t), &
                      vgt     (k,t), &
                      wflst   (k,t), &
                      dqtdxlst(k,t), &
                      dqtdylst(k,t), &
                      dqtdtlst(k,t), &
                      thlpcart(k,t)
              end do
            else ! else read in same columns as dpdx and dpdy
                do k=1,kmax
                  read (ifinput,*) &
                        height  (k)  , &
                        dpdxlt  (k,t), &
                        dpdylt  (k,t), &
                        wflst   (k,t), &
                        dqtdxlst(k,t), &
                        dqtdylst(k,t), &
                        dqtdtlst(k,t), &
                        thlpcart(k,t)
                end do
            end if


            end if
          end do

          close(ifinput)
        end if

      end if   !ltestbed

!      do k=kmax,1,-1
!        write (6,'(3f7.1,5e12.4)') &
!            height  (k)  , &
!            ugt     (k,t), &
!            vgt     (k,t), &
!            wflst   (k,t), &
!            dqtdxlst(k,t), &
!            dqtdylst(k,t), &
!            dqtdtlst(k,t), &
!            thlpcart(k,t)
!      end do

      if(timeflux(1)>runtime) then
        write(6,*) 'Time dependent surface variables do not change before end of'
        write(6,*) 'simulation. --> only large scale forcings'
        ltimedepsurf=.false.
      endif

      if ((timels(1) > runtime) .or. (timeflux(1) > runtime)) then
        write(6,*) 'Time dependent large scale forcings sets in after end of simulation -->'
        write(6,*) '--> only time dependent surface variables'
        ltimedepz=.false.
      end if

      close(ifinput)

    end if

    call D_MPI_BCAST(timeflux(1:kflux),kflux   ,0,comm3d,mpierr)
    call D_MPI_BCAST(wtsurft          ,kflux   ,0,comm3d,mpierr)
    call D_MPI_BCAST(wqsurft          ,kflux   ,0,comm3d,mpierr)
    call D_MPI_BCAST(thlst            ,kflux   ,0,comm3d,mpierr)
    call D_MPI_BCAST(qtst             ,kflux   ,0,comm3d,mpierr)
    call D_MPI_BCAST(pst              ,kflux   ,0,comm3d,mpierr)
    call D_MPI_BCAST(Qnetavt          ,kflux   ,0,comm3d,mpierr)
    call D_MPI_BCAST(timels(1:kls)    ,kls     ,0,comm3d,mpierr)
    call D_MPI_BCAST(ugt              ,kmax*kls,0,comm3d,mpierr)
    call D_MPI_BCAST(vgt              ,kmax*kls,0,comm3d,mpierr)
    call D_MPI_BCAST(dpdxlt           ,kmax*kls,0,comm3d,mpierr)
    call D_MPI_BCAST(dpdylt           ,kmax*kls,0,comm3d,mpierr)
    call D_MPI_BCAST(wflst            ,kmax*kls,0,comm3d,mpierr)
    call D_MPI_BCAST(dqtdxlst,kmax*kls ,0,comm3d,mpierr)
    call D_MPI_BCAST(dqtdylst,kmax*kls ,0,comm3d,mpierr)
    call D_MPI_BCAST(dqtdtlst,kmax*kls ,0,comm3d,mpierr)
    call D_MPI_BCAST(dthldtlst,kmax*kls,0,comm3d,mpierr)
    call D_MPI_BCAST(dudtlst,kmax*kls  ,0,comm3d,mpierr)
    call D_MPI_BCAST(dvdtlst,kmax*kls  ,0,comm3d,mpierr)
    call D_MPI_BCAST(thlpcart,kmax*kls ,0,comm3d,mpierr)
    call D_MPI_BCAST(thlproft,kmax*kls ,0,comm3d,mpierr)
    call D_MPI_BCAST(qtproft ,kmax*kls ,0,comm3d,mpierr)

    call D_MPI_BCAST(ltimedepsurf ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(ltimedepz    ,1,0,comm3d,mpierr)

    call inittimedepsv
    call timedep

    deallocate(height)


  end subroutine inittimedep

  subroutine timedep

!-----------------------------------------------------------------|
!                                                                 |
!*** *timedep*  calculates ls forcings and surface forcings       |
!               case as a funtion of timee                        |
!                                                                 |
!      Roel Neggers    K.N.M.I.     01/05/2001                    |
!                                                                 |
!                                                                 |
!    calls                                                        |
!    * timedepz                                                   |
!      calculation of large scale advection, radiation and        |
!      surface fluxes by interpolation between prescribed         |
!      values at certain times                                    |
!                                                                 |
!    * timedepsurf                                                |
!      calculation  surface fluxes by interpolation               |
!      between prescribed values at certain times                 |
!                                                                 |
!                                                                 |
!-----------------------------------------------------------------|
    use modtimedepsv, only : timedepsv
    implicit none

    if (.not. ltimedep) return
    call timedepz
    call timedepsurf
    call timedepsv
  end subroutine timedep

  subroutine timedepz
    use modfields,   only : ug, vg, wfls,whls, &
                            dqtdtls,dqtdxls,dqtdyls, &
                            dthldtls,dthldxls,dthldyls,thlpcar, &
                            dudtls,dudxls,dudyls, &
                            dvdtls,dvdxls,dvdyls, &
                            dpdxl,dpdyl

    use modglobal,   only : rtimee,om23_gs,dzf,dzh,k1,kmax,llsadv,lcoriol

    use modmpi,      only : myid

    implicit none

    character(len=*), parameter :: routine = modname//'/timedepz'

    integer t,k
    real fac

    if(.not.(ltimedepz)) return

    !---- interpolate ----
    t=1
    do while(rtimee>timels(t+1))
       t=t+1
    end do
    ! timels(t) < rtimee <= timels(t+1)
    ! or t = 1 if rtimee < timels(1)

    fac = ( rtimee-timels(t) ) / ( timels(t+1)-timels(t) )
    ug       = ugt      (:,t) + fac * ( ugt      (:,t+1) - ugt      (:,t) )
    vg       = vgt      (:,t) + fac * ( vgt      (:,t+1) - vgt      (:,t) )
    wfls     = wflst    (:,t) + fac * ( wflst    (:,t+1) - wflst    (:,t) )
    dqtdxls  = dqtdxlst (:,t) + fac * ( dqtdxlst (:,t+1) - dqtdxlst (:,t) )
    dqtdyls  = dqtdylst (:,t) + fac * ( dqtdylst (:,t+1) - dqtdylst (:,t) )
    dqtdtls  = dqtdtlst (:,t) + fac * ( dqtdtlst (:,t+1) - dqtdtlst (:,t) )
    dthldtls = dthldtlst(:,t) + fac * ( dthldtlst(:,t+1) - dthldtlst(:,t) )
    dudtls   = dudtlst  (:,t) + fac * ( dudtlst  (:,t+1) - dudtlst  (:,t) )
    dvdtls   = dvdtlst  (:,t) + fac * ( dvdtlst  (:,t+1) - dvdtlst  (:,t) )
    thlpcar  = thlpcart (:,t) + fac * ( thlpcart (:,t+1) - thlpcart (:,t) )

    if (lcoriol) then
      do k=1,kmax
        dpdxl(k) =  om23_gs*vg(k)
        dpdyl(k) = -om23_gs*ug(k)
      end do
    else ! if not lcoriol, update regular pressure gradients
      dpdxl      = dpdxlt (:,t) + fac * ( dpdxlt (:,t+1) - dpdxlt (:,t) )
      dpdyl      = dpdylt (:,t) + fac * ( dpdylt (:,t+1) - dpdylt (:,t) )
    end if

    whls(1)  = 0.0
    do k=2,kmax
      whls(k) = ( wfls(k)*dzf(k-1) +  wfls(k-1)*dzf(k) )/(2*dzh(k))
    end do
    whls(k1) = (wfls(kmax)+0.5*dzf(kmax)*(wfls(kmax)-wfls(kmax-1)) &
                                                  /dzh(kmax))

  !******include rho if rho = rho(z) /= 1.0 ***********

    if (llsadv) then
      if (myid==0) call finish(routine, 'llsadv should not be used anymore. Large scale gradients were calculated in a non physical way (and lmomsubs had to be set to true to retain conservation of mass)')
    end if

    dudxls   = 0.0
    dudyls   = 0.0
    dvdxls   = 0.0
    dvdyls   = 0.0
    dthldxls = 0.0
    dthldyls = 0.0


    return
  end subroutine timedepz

  subroutine timedepsurf
    use modglobal,   only : rtimee
    use modthermodynamics, only: lmoist
    use modsurfdata, only : wtsurf,wqsurf,thls,qts,ps, Qnetav
    use modsurface,  only : qtsurf
    implicit none
    integer t
    real fac

    if(.not.(ltimedepsurf)) return
  !     --- interpolate! ----
    t=1
    do while(rtimee>timeflux(t+1))
      t=t+1
    end do
    ! timeflux(t) <= rtimee <= timeflux(t+1)
    ! or t = 1 if rtimee < timeflux(1)

    fac = ( rtimee-timeflux(t) ) / ( timeflux(t+1)-timeflux(t))
    wqsurf = wqsurft(t) + fac * ( wqsurft(t+1) - wqsurft(t)  )
    wtsurf = wtsurft(t) + fac * ( wtsurft(t+1) - wtsurft(t)  )
    thls   = thlst(t)   + fac * ( thlst(t+1)   - thlst(t)    )
    ps     = pst(t)     + fac * ( pst(t+1)   - pst(t)    )
    Qnetav = Qnetavt(t) + fac * ( Qnetavt(t+1) - Qnetavt(t)  )
!cstep: not necessary to provide qts in ls_flux file qts    = qtst(t)    + fac * ( qtst(t+1)    - qtst(t)     )
    if (lmoist) then
       call qtsurf
    else
       qts = 0.
    endif

    return
  end subroutine timedepsurf


  subroutine exittimedep
    use modtimedepsv, only : exittimedepsv
    implicit none
    if (.not. ltimedep) return
    deallocate(timels,ugt,vgt,wflst,dqtdxlst,dqtdylst,dqtdtlst,dthldtlst,dudtlst,dvdtlst,thlpcart)
    deallocate(timeflux, wtsurft,wqsurft,thlst,qtst,pst,Qnetavt)
    call exittimedepsv

  end subroutine


  !> \brief Read initial profiles from forcings.XXX.nc
  subroutine init_timedep_from_netcdf(filename, height, time,ntimedep,kmax)
    use modstat_nc, only : read_nc_field, nchandle_error
    use netcdf, only : NF90_NOWRITE, nf90_open, nf90_close
    implicit none
    character(*),   intent(in)  :: filename !< Path to the netCDF file to read from.
    real,           intent(out) :: height(:) !< Vertical levels.
    real,           intent(out) :: time(:) !< Time steps.
    integer,        intent(in)  :: ntimedep !< Number of time steps for which time-dependent forcings are provided.
    integer,        intent(in)  :: kmax !< Index of highest vertical level.

    integer :: ncid

    call nchandle_error(nf90_open(filename, NF90_NOWRITE, ncid))

    ! "Regular" prognostic fields
    call read_nc_field(ncid, "zh", height, start=1, count=kmax)
    call read_nc_field(ncid, "time", time, start=1, count=ntimedep)

    ! Large-scale forcings
    call read_nc_field(ncid, "ug_timedep", ugt, start=(/1,1/), count=(/ntimedep,kmax/), &
                       requirefill=.false.)
    call read_nc_field(ncid, "vg_timedep", vgt, start=(/1,1/), count=(/ntimedep,kmax/), &
                       requirefill=.false.)
    call read_nc_field(ncid, "dpdx_ls_timedep", dpdxlt, start=(/1,1/), count=(/ntimedep,kmax/), &
                       requirefill=.false.)
    call read_nc_field(ncid, "dpdy_ls_timedep", dpdylt, start=(/1,1/), count=(/ntimedep,kmax/), &
                       requirefill=.false.)
    call read_nc_field(ncid, "wf_ls_timedep", wflst, start=(/1,1/), count=(/ntimedep,kmax/), &
                       requirefill=.false.)
    call read_nc_field(ncid, "dqtdx_ls_timedep", dqtdxlst, start=(/1,1/), count=(/ntimedep,kmax/), &
                       requirefill=.false.)
    call read_nc_field(ncid, "dqtdy_ls_timedep", dqtdylst, start=(/1,1/), count=(/ntimedep,kmax/), &
                       requirefill=.false.)
    call read_nc_field(ncid, "dqtdt_ls_timedep", dqtdtlst, start=(/1,1/), count=(/ntimedep,kmax/), &
                       requirefill=.false.)
    call read_nc_field(ncid, "dthldt_ls_timedep", dthldtlst, start=(/1,1/), count=(/ntimedep,kmax/), &
                       requirefill=.false.)
    call read_nc_field(ncid, "dudt_ls_timedep", dudtlst, start=(/1,1/), count=(/ntimedep,kmax/), &
                       requirefill=.false.)
    call read_nc_field(ncid, "dvdt_ls_timedep", dvdtlst, start=(/1,1/), count=(/ntimedep,kmax/), &
                       requirefill=.false.)
    call read_nc_field(ncid, "dthl_rad_timedep", thlpcart, start=(/1,1/), count=(/ntimedep,kmax/), &
                       requirefill=.false.)
    call read_nc_field(ncid, "wtsurf_timedep", wtsurft, start=1, count=ntimedep, &
                       requirefill=.false.)
    call read_nc_field(ncid, "wqsurf_timedep", wqsurft, start=1, count=ntimedep, &
                       requirefill=.false.)
    call read_nc_field(ncid, "thlsurf_timedep", thlst, start=1, count=ntimedep, &
                       requirefill=.false.)
    call read_nc_field(ncid, "qtsurf_timedep", qtst, start=1, count=ntimedep, &
                       requirefill=.false.)
    call read_nc_field(ncid, "psurf_timedep", pst, start=1, count=ntimedep, &
                       requirefill=.false.)
    call read_nc_field(ncid, "qnetavsurf_timedep", Qnetavt, start=1, count=ntimedep, &
                       requirefill=.false.)

    call nchandle_error(nf90_close(ncid))

  end subroutine init_timedep_from_netcdf
end module modtimedep
