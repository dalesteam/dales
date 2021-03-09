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


implicit none
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

  real, allocatable     :: timeflux (:)
  real, allocatable     :: wqsurft  (:)
  real, allocatable     :: wtsurft  (:)
  real, allocatable     :: thlst    (:)
  real, allocatable     :: qtst     (:)
  real, allocatable     :: pst      (:)
  real, allocatable     :: Qnetavt  (:)

  real, allocatable     :: timels  (:)
  real, allocatable     :: ugt     (:,:)
  real, allocatable     :: vgt     (:,:)
  real, allocatable     :: wflst   (:,:)
  real, allocatable     :: dqtdxlst(:,:)
  real, allocatable     :: dqtdylst(:,:)
  real, allocatable     :: dqtdtlst(:,:)
  real, allocatable     :: dthldtlst(:,:)
  real, allocatable     :: thlpcart(:,:)
  real, allocatable     :: dudtlst (:,:)
  real, allocatable     :: dvdtlst (:,:)
  real, allocatable     :: thlproft(:,:)
  real, allocatable     :: qtproft (:,:)



contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine inittimedep
    use modmpi,    only :myid,mpi_logical,mpierr,comm3d,D_MPI_BCAST
    use modglobal, only :cexpnr,k1,kmax,ifinput,runtime,zf,ntimedep
    use modsurfdata,only :ps,qts,wqsurf,wtsurf,thls, Qnetav
    use modtimedepsv, only : inittimedepsv

    use modtestbed,        only : ltestbed,ntnudge,&
                                  tb_time,tb_ps,tb_qts,tb_thls,tb_wqs,tb_wts,&
                                  tb_w,tb_ug,tb_vg,&
                                  tb_uadv,tb_vadv,tb_qtadv,tb_thladv,tb_Qnet

    implicit none

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
             write (*,*) "Too many time points in file ", 'ls_flux.inp.'//cexpnr, ", the limit is kflux = ", kflux 
             stop
          end if
          read(ifinput,*, iostat = ierr) timeflux(t), wtsurft(t), wqsurft(t),thlst(t),qtst(t),pst(t)
          write(*,'(i8,6e12.4)') t,timeflux(t), wtsurft(t), wqsurft(t),thlst(t),qtst(t),pst(t)
          if (ierr < 0) then
            stop 'STOP: No time dependend data for end of run (surface fluxes)'
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
             write (*,*) "Too many time points in file ", 'nudge.inp.'//cexpnr, ", the limit is kls = ", kls
             stop
          end if
          chmess1 = "#"
          ierr = 1 ! not zero
          do while (.not.(chmess1 == "#" .and. ierr ==0)) !search for the next line consisting of "# time", from there onwards the profiles will be read
            read(ifinput,*,iostat=ierr) chmess1,timels(t)
            if (ierr < 0) then
              stop 'STOP: No time dependend data for end of run'
            end if
          end do
          write (*,*) 'timels = ',timels(t)
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
          end if
       end do

        close(ifinput)

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

    use modglobal,   only : rtimee,om23_gs,dzf,dzh,k1,kmax,llsadv

    use modmpi,      only : myid

    implicit none

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


    do k=1,kmax
      dpdxl(k) =  om23_gs*vg(k)
      dpdyl(k) = -om23_gs*ug(k)
    end do

    whls(1)  = 0.0
    do k=2,kmax
      whls(k) = ( wfls(k)*dzf(k-1) +  wfls(k-1)*dzf(k) )/(2*dzh(k))
    end do
    whls(k1) = (wfls(kmax)+0.5*dzf(kmax)*(wfls(kmax)-wfls(kmax-1)) &
                                                  /dzh(kmax))

  !******include rho if rho = rho(z) /= 1.0 ***********

    if (llsadv) then
      if (myid==0) stop 'llsadv should not be used anymore. Large scale gradients were calculated in a non physical way (and lmomsubs had to be set to true to retain conservation of mass)'
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
    use modglobal,   only : rtimee, lmoist
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

end module modtimedep
