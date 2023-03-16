!> \file modibm.f90
!! By Michael Koene, email: mich4all@live.nl, TU Delft, section Atmospheric Physics, October 8, 2019
!! TODO: Write comments to output file
!!       clean up code (write statements)
!!       Test performance

module modibm

  use modibmdata, only : lapply_ibm,lpoislast, lreadfile_obstacles, lwallheat, &
                         thlwall, thlroof,qtroof, thlibm,qtibm, &
                         libm,Nair, ibm_adv_mask,kibm_max,&
                         z0m_wall
  implicit none
  save
  public :: initibm, exitibm,applyibm, zerowallvelocity

  ! Fields

  !< Normal immersed boundary layers for incoming x,y,z-velocities
  logical, allocatable :: lnorm_x(:,:,:), lnorm_y(:,:,:), lnorm_z(:,:,:)

  real, allocatable    :: tempsvp(:,:,:,:)
  real, allocatable    :: tempthlp(:,:,:), tempqtp(:,:,:)
  real, allocatable    :: tempup(:,:,:), tempvp(:,:,:), tempwp(:,:,:)

  real :: dx_half,dy_half,Cm_xwall,Cm_ywall, Cd_xwall, Cd_ywall


contains
  subroutine initibm
    use modglobal,  only : zh,zf,itot, jtot, ih, i1, jh, j1, k1, imax, jmax, kmax, cexpnr, ifnamopt, ifinput, &
                           fname_options, nsv, cu, cv, ijtot, &
                           iadv_mom,iadv_tke,iadv_thl,iadv_qt,iadv_sv, &
                           iadv_cd2,iadv_5th,iadv_52,iadv_cd6,iadv_62,iadv_kappa,iadv_upw,iadv_hybrid,& 
                           iadv_hybrid_f,ibas_prf,&
                           dx,dy,fkar
    use modmpi,     only : myid, comm3d, mpierr,  myidx, myidy, d_mpi_bcast, excjs, D_MPI_ALLREDUCE, &
                           mpi_max, mpi_sum
    use modfields,   only : ksfc
    !use mpi
    implicit none

     !< Field for the immersed boundary height
    real, allocatable :: bc_height(:,:)     !< Height of immersed boundary at grid pos x,y
    integer, allocatable :: Nairl(:)
    integer       :: i, j, k, ierr,ii,jj,kk,n  !cstep , kmin
    integer       :: advarr(4)
    integer       :: ibm_adv_mask_imin, ibm_adv_mask_imax, & !< use 2nd order advection near obstacles
                     ibm_adv_mask_kmin, ibm_adv_mask_kmax
    integer          :: kibm_maxl
    character(100) :: readstring

    namelist/NAMIBM/ lapply_ibm, lreadfile_obstacles, &
                               lwallheat, &
                               thlwall, thlibm, thlroof, qtibm,lpoislast, z0m_wall


    if(myid==0) then    !first myid
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMIBM,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMIBM'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMIBM'
      endif
      write(6 ,NAMIBM)
      close(ifnamopt)
    endif

    !if(.not.(myid==0)) return

    call D_MPI_BCAST(lapply_ibm   ,    1,  0, comm3d, mpierr)

    if (.not. (lapply_ibm)) return

    call D_MPI_BCAST(lreadfile_obstacles,    1, 0, comm3d, mpierr)
    call D_MPI_BCAST(lwallheat          ,    1, 0, comm3d, mpierr)  !zero (Cd->0) or nonzero heat flux from wall
    call D_MPI_BCAST(thlwall            ,    1, 0, comm3d, mpierr)
    call D_MPI_BCAST(thlibm             ,    1, 0, comm3d, mpierr)  !cstep thl inside obstacle
    call D_MPI_BCAST(qtibm              ,    1, 0, comm3d, mpierr)
    call D_MPI_BCAST(thlroof            ,    1, 0, comm3d, mpierr)  !cstep not used yet
    call D_MPI_BCAST(lpoislast          ,    1, 0, comm3d, mpierr)
    call D_MPI_BCAST(z0m_wall           ,    1, 0, comm3d, mpierr)

    dx_half = 0.5 * dx
    dy_half = 0.5 * dy
    Cm_xwall = (fkar/(log(dx_half/z0m_wall)))**2  !cstep similar to neutral boundary layer with a log wind profile
    Cm_ywall = (fkar/(log(dy_half/z0m_wall)))**2
    if (lwallheat) then
      Cd_xwall = Cm_xwall   !cstep offer possibility to set it to zero (zero heat flux BC)
      Cd_ywall = Cm_ywall
    else
      Cd_xwall = 0.
      Cd_ywall = 0.
    endif 


    if (lapply_ibm) then 
       ibas_prf = 2
       if (myid.eq.0) then
          write (6,*) 'ibas_prf is set to 2 (Boussinesq constant density)'
          write (6,*) 'height dependent density gives problems with correction vertical advective'
          write (6,*) 'tendencies at the tops of obstacles'
       endif
    endif


    if (abs(cu)>1e-15 .or. abs(cv)>1e-15) then
      if(myid==0) print *, 'Problem in namoptions'
      if(myid==0) print *, 'cu or cv cannot be nonzero while using IBM'
      if(myid==0) print *, 'The buildings would move in that case'
      if(myid==0) print *, 'Set cu and cv to 0. to solve this problem or simulate without buildings'
      stop 'ERROR: Problem in namoptions NAMIBM with cu and cv'
    endif

    write(6,*) 'allocating fields in modibm'

    allocate(bc_height (itot+1,jtot+1))
    allocate(libm(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(ibm_adv_mask(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(Nair(k1))
    allocate(Nairl(k1))

    allocate(tempsvp(2-ih:i1+ih,2-jh:j1+jh,k1,nsv))
    allocate(tempthlp(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(tempqtp(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(tempup(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(tempvp(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(tempwp(2-ih:i1+ih,2-jh:j1+jh,k1))

    allocate(lnorm_x (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(lnorm_y (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(lnorm_z (2-ih:i1+ih,2-jh:j1+jh,k1))


    write(6,*) 'succesfully allocated fields in modibm'

    ibm_adv_mask (:,:,:) = 0.
    bc_height(:,:) = 0
    libm (:,:,:) = .false.
    

    ! Definition of obstacles
    if (myid==0) then
      if (lreadfile_obstacles) then  !< Profile prescribed by use in the file ibm.inp.<expnr>
        write(6,*) 'Reading inputfile in modibm'
        open (ifinput,file='ibm.inp.'//cexpnr)
          do k=1,7 
           read (ifinput,'(a100)') readstring
           !  read (ifinput,*) readstring
           write (6,*) readstring
          enddo

          do j=jtot+1,2,-1  !cstep
!            read (ifinput,'(a1000)') readstring
            !cstep     write (6,*) 'j1',j, readstring
            !< If the program is unable to read the full line of points increasing the length of the string (a400) might help

!            do while (readstring(1:1)=='#')  ! Skip the lines that are commented (like headers)
!              read (ifinput,'(a1000)') readstring
!              write (6,*) 'readstring', readstring
!            end do
!            read(readstring,*) (bc_height(i+1,j+1),i=1,itot)
!            write (6,*) 'j2',j,bc_height(1:itot+1,j+1)

            do i=2,itot+1
               read(ifinput,'(F6.1)') bc_height(i,j)
            enddo
          end do
        close(ifinput)

        bc_height(1,:)=bc_height(itot+1,:)
        bc_height(:,1)=bc_height(:,jtot+1)

        write(6,*) 'Succesfully read inputfile in modibm'
      else           !< Simple block in the middle of the domain
        write(6,*) 'Generating standard boundary in modibm, it assumes indices k for zh and obstacle heights yet' 
        bc_height(NINT(itot*0.5):(NINT(itot*0.5)+1),NINT(jtot*0.5):(NINT(jtot*0.5)+1))=NINT(kmax*0.5) 

        bc_height(1,:)=bc_height(itot+1,:)
        bc_height(:,1)=bc_height(:,jtot+1)

        write(6,*) 'Succesfully generated immersed boundary in modibm'
      endif
    endif  !myid==0



    call D_MPI_BCAST(bc_height,(itot+1)*(jtot+1),0,comm3d,mpierr)

    kibm_maxl = 0. !cstep  the index of the highest obstacle in the entire domain  
    do i=2,i1
      do j=2,j1
       ! do k=1,kmax  !skip zero values which are ground surface points
       !   if (zh(k+1).LE.bc_height(i+myidx*imax,j+myidy*jmax)) then  !cstep read in heights rather than indices
              !if zh is 1 mm above building height, height is set to dz below (maybe better use zf as a criterion for nicer rounding off)

       !     libm (i,j,k) = .true.
       !     ksfc (i,j)   = k + 1

         do k=1,kmax
            if (zf(k).LE.bc_height(i+myidx*imax,j+myidy*jmax)) then  !obstacle height is above mid point of vertical grid
              libm (i,j,k) = .true.
              ksfc (i,j)   = k + 1   !half (flux) level
              if (ksfc(i,j).gt.kibm_maxl) then
                kibm_maxl = k
              endif
              write (6,*) 'libm',i+myidx*imax,j+myidy*jmax,i,j,k,libm(i,j,k),bc_height(i+myidx*imax,j+myidy*jmax),zh(ksfc(i,j))
           endif
        end do

      end do
    end do

    call excjs( libm  , 2,i1,2,j1,1,k1,ih,jh)
    call D_MPI_ALLREDUCE(kibm_maxl,kibm_max,1,MPI_MAX,comm3d,mpierr)


    Nair(:) = 0.
    Nairl(:) = 0
    do i=2,i1
      do j=2,j1
        do k=ksfc(i,j),k1 
          Nairl(k) = Nairl(k)+1
        enddo
      enddo
    enddo
    !cstep write (6,*) 'Nairl ',Nairl
    call D_MPI_ALLREDUCE(Nairl, Nair, k1, MPI_SUM, comm3d,mpierr)
    !cstep write (6,*) 'Nair',Nair

    call constructboundarytypes


    advarr = (/iadv_mom,iadv_tke,iadv_thl,iadv_qt/)
    if     (any(advarr==iadv_cd6).or.any(iadv_sv(1:nsv)==iadv_cd6)) then
      ibm_adv_mask_imin = 3
      ibm_adv_mask_imax = 3
      ibm_adv_mask_kmin = 3
      ibm_adv_mask_kmax = 3
    elseif (any(advarr==iadv_62).or.any(iadv_sv(1:nsv)==iadv_62)) then
      ibm_adv_mask_imin = 3
      ibm_adv_mask_imax = 3
      ibm_adv_mask_kmin = 1
      ibm_adv_mask_kmax = 1
    elseif (any(advarr==iadv_5th).or.any(iadv_sv(1:nsv)==iadv_5th)) then
      ibm_adv_mask_imin = 3
      ibm_adv_mask_imax = 3
      ibm_adv_mask_kmin = 2
      ibm_adv_mask_kmax = 2
    elseif (any(advarr==iadv_52).or.any(iadv_sv(1:nsv)==iadv_52)) then
      ibm_adv_mask_imin = 3
      ibm_adv_mask_imax = 3
      ibm_adv_mask_kmin = 1
      ibm_adv_mask_kmax = 1
    elseif (any(advarr==iadv_hybrid).or.any(iadv_sv(1:nsv)==iadv_hybrid)) then
      ibm_adv_mask_imin = 3
      ibm_adv_mask_imax = 2
      ibm_adv_mask_kmin = 3
      ibm_adv_mask_kmax = 2
    elseif (any(advarr==iadv_hybrid_f).or.any(iadv_sv(1:nsv)==iadv_hybrid_f)) then
      ibm_adv_mask_imin = 3
      ibm_adv_mask_imax = 2
      ibm_adv_mask_kmin = 3
      ibm_adv_mask_kmax = 2
    elseif (any(advarr==iadv_kappa).or.any(iadv_sv(1:nsv)==iadv_kappa)) then
      ibm_adv_mask_imin = 2
      ibm_adv_mask_imax = 1
      ibm_adv_mask_kmin = 2
      ibm_adv_mask_kmax = 1
    elseif (any(advarr==iadv_cd2).or.any(iadv_sv(1:nsv)==iadv_cd2)) then
      ibm_adv_mask_imin = 1
      ibm_adv_mask_imax = 1
      ibm_adv_mask_kmin = 1
      ibm_adv_mask_kmax = 1
    end if

    do k=1,kmax-ibm_adv_mask_kmax
    do i=2,imax
    do j=2,jmax
       if (libm(i,j,k)) then
          do ii=i-ibm_adv_mask_imax,i+ibm_adv_mask_imin
             ibm_adv_mask(ii,j,k) = 1
          enddo
          do jj=j-ibm_adv_mask_imax,j+ibm_adv_mask_imin
             ibm_adv_mask(i,jj,k) = 1
          enddo
          do kk=k,k+ibm_adv_mask_kmin
             ibm_adv_mask(i,j,kk) = 1
          enddo
       endif
    enddo
    enddo
    enddo


    write(6,* ) 'start deallocate'
    deallocate (bc_height,Nairl)
    write(6,*) 'exit initibm'

    return
  end subroutine initibm

  subroutine constructboundarytypes   !< Calculate the positions of the different boundary layers in multiple directions

    use modglobal,  only : imax, i1, ih, jmax, j1, jh, kmax, k1
    use modmpi,     only : myid, excjs
    implicit none
    integer i,j,k,ipos,jpos


    if(myid==0) write(6,*) 'Starting constructboundarytypes in modibm, myid =',myid

    !< Fill the layer types with .false.
    lnorm_x(:,:,:)=.false.
    lnorm_y(:,:,:)=.false.
    lnorm_z(:,:,:)=.false.

!cstep   else  !do for any lwall parameterization
      !< Find normal layers in x-direction
      do k=1,kmax
        do j=2,j1
          do i=2,i1
         !libm   ipos=i+myidx*imax
         !libm   jpos=j+myidy*jmax
         !libm   if (.not. (limmersed_boundary(ipos,jpos,k)==limmersed_boundary(ipos-1,jpos,k))) then

         !cstep            0     X      X     X      0       ,building position
         !cstep           i-2   i-1     i    i+1    i+2
         !cstep  libm      F     T      T     T      F
         !cstep  lnorm_x   F     T      F     F      T  , the true points refer to u-positions on the 
         !                                                grid box (on its left)

            if (.not. (libm(i,j,k).eqv.libm(i-1,j,k))) then
              lnorm_x(i,j,k)=.true.  !cstep a wall at position i with its normal pointing in the x-direction
            endif
          end do
        end do
      end do
      if(myid==0) write(6,*) 'Succesfully found normal layers in x-direction'

      !< Find normal layers in y-direction
      do k=1,kmax
        do i=2,i1
          do j=2,j1
           !libm ipos=i+myidx*imax
           !libm jpos=j+myidy*jmax
            if (.not. (libm(i,j,k).eqv.libm(i,j-1,k))) then
              lnorm_y(i,j,k)=.true.
            endif
          end do
        end do
      end do

      !< Find normal layers in z-direction
       do i=2,i1
        do j=2,j1
          do k=2,kmax
            if (.not. (libm(i,j,k).eqv.libm(i,j,k-1))) then
              lnorm_z(i,j,k)=.true.
            endif
          end do
        end do
      end do


      if(myid==0) write(6,*) 'Succesfully found normal layers in all directions'

    call excjs( lnorm_x  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( lnorm_y  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( lnorm_z  , 2,i1,2,j1,1,k1,ih,jh)

    if(myid==0) write(6,*) 'finished constructboundarytypes'


  end subroutine constructboundarytypes


  subroutine exitibm
    use modmpi, only : myid
    implicit none


    if (.not. (lapply_ibm)) return
    deallocate(Nair)

    deallocate(tempsvp)
    deallocate(tempthlp)
    deallocate(tempqtp)
    deallocate(tempup)
    deallocate(tempvp)
    deallocate(tempwp)

    deallocate(lnorm_x)
    deallocate(lnorm_y)
    deallocate(lnorm_z)
    deallocate(libm)
    deallocate(ibm_adv_mask)
    return
  end subroutine exitibm

  subroutine applyibm(simid)   !< apply immersed boundary method
    use modfields,      only : um, vm, wm, thlm, qtm, e12m, svm, &
                               u0, v0, w0, thl0, qt0, e120, sv0, &
                               up, vp, wp, thlp, qtp, e12p, svp, &
                               thl0av
    use modglobal,      only : rk3step,   kmax, i1, j1, k1, ih, jh, rdt, timee, dx, dy, dzh, dzf, nsv, e12min
    use modsubgriddata, only : ekm
    use modmpi,         only : excjs 
    !clater use modnudgeboundary, only : Nsim

    implicit none
    integer  :: i, j, k,  nc
    real     :: rk3coef,rk3coefi
    real     :: emmo, emom, eomm, empo, epmo, emop, eomp
    real     :: w_at_v_min,w_at_v_plus,w_at_u_min,w_at_u_plus
    real     :: u_at_w_min,u_at_w_plus,v_at_w_min,v_at_w_plus
    real     :: us_at_scalar_min,us_at_scalar_plus
    real     :: tau_vu_plus,tau_vu_min,tau_wu_min,tau_wu_plus,tau_uv_min,tau_uv_plus,tau_wv_min,tau_wv_plus

    integer  :: maxlocx(3)
    integer, intent(in) :: simid

    if (.not. lapply_ibm) return
    !if (.not. (lapply_ibm .and. simid == Nsim)) return  !cstep ensures buildings are not applied for precursor run
               !cstep run if lapply_ibm = true AND simid=1 without precursor simulation (so Nsim=1)
               !      or  if lapplY_ibm = true AND simid=2 with    precursor simulation (Nsim=2)

    thlibm = thl0av(1) !assumes inside air has the same temperature as the air outside

    rk3coef = rdt / (4. - dble(rk3step))
    rk3coefi = 1. / rk3coef

    tempsvp(:,:,:,:)=0.
    tempthlp(:,:,:)=0.
    tempqtp(:,:,:)=0.
    tempup(:,:,:)=0.
    tempvp(:,:,:)=0.
    tempwp(:,:,:)=0.

       do i=2,i1  
        do j=2,j1 
          do k=2,kmax          !cstep special treatment for k=1

            if (lnorm_x(i,j,k)) then     !< Wall in x-direction
              
              emmo = 0.25  * (  ekm(i,j,k)+ekm(i,j-1,k)+ekm(i-1,j-1,k)+ekm(i-1,j,k)  )
              w_at_v_min  = 0.25*(w0(i-1,j,k)+w0(i-1,j,k+1)+w0(i-1,j-1,k)+w0(i-1,j-1,k+1))  !at v(i-1,j,k)
              w_at_v_plus = 0.25*(w0(i  ,j,k)+w0(i  ,j,k+1)+w0(i  ,j-1,k)+w0(i  ,j-1,k+1))  !at v(i,j,k)
              call log_wallaw(v0(i-1,j,k),w_at_v_min ,Cm_xwall,tau_vu_min)   !if v0 > 0, tau > 0
              call log_wallaw(v0(i,j,k)  ,w_at_v_plus,Cm_xwall,tau_vu_plus)  !minus sign in tendency enforces opposing friction
         
              tempvp(i-1,j,k) = tempvp(i-1,j,k) - 0.5 * emmo*((v0(i,j,k)-v0(i-1,j,k))/dx) / dx - 0.5 * tau_vu_min /dx
              tempvp(i  ,j,k) = tempvp(i  ,j,k) + 0.5 * emmo*((v0(i,j,k)-v0(i-1,j,k))/dx) / dx - 0.5 * tau_vu_plus/dx
                !  dv/dt = - d/dx uv = -uv|_xright + uv|_xleft
                !two times sign can be understood from uv = - K dv/dx. right of building dv/dx > 0, left of it dv/dx < 0
                ! for v>0 in between buildings this gives uv|_right > 0, uv|_left < 0. in that case both have a damping tendency on v  
          
              empo = 0.25  * ( ekm(i,j+1,k)+ekm(i,j,k)+ekm(i-1,j,k)+ekm(i-1,j+1,k)  )
              w_at_v_min  = 0.25*(w0(i-1,j+1,k)+w0(i-1,j+1,k+1)+w0(i-1,j,k)+w0(i-1,j,k+1))
              w_at_v_plus = 0.25*(w0(i  ,j+1,k)+w0(i  ,j+1,k+1)+w0(i  ,j,k)+w0(i  ,j,k+1))
              call log_wallaw(v0(i-1,j+1,k),w_at_v_min ,Cm_xwall,tau_vu_min)
              call log_wallaw(v0(i,j+1,k)  ,w_at_v_plus,Cm_xwall,tau_vu_plus)

              tempvp(i-1,j+1,k) = tempvp(i-1,j+1,k) - 0.5 * empo*((v0(i,j+1,k)-v0(i-1,j+1,k))/dx) / dx - 0.5 * tau_vu_min   /dx
              tempvp(i  ,j+1,k) = tempvp(i,j+1,k)   + 0.5 * empo*((v0(i,j+1,k)-v0(i-1,j+1,k))/dx) / dx - 0.5 * tau_vu_plus  /dx

              emom = ( dzf(k-1) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,k-1) + ekm(i-1,j,k-1) ) ) / &
                    ( 4.   * dzh(k) )

              v_at_w_min  = 0.25 * (v0(i-1,j,k-1)+v0(i-1,j,k)+v0(i-1,j+1,k-1)+v0(i-1,j+1,k) )
              v_at_w_plus = 0.25 * (v0(i  ,j,k-1)+v0(i  ,j,k)+v0(i  ,j+1,k-1)+v0(i  ,j+1,k) )
              call log_wallaw(w0(i-1,j,k),v_at_w_min ,Cm_xwall,tau_wu_min)
              call log_wallaw(w0(i  ,j,k),v_at_w_plus,Cm_xwall,tau_wu_plus)


              tempwp(i-1,j,k) = tempwp(i-1,j,k) - 0.5 * emom * ((w0(i,j,k)-w0(i-1,j,k))/dx)/dx - 0.5 * tau_wu_min/dx
              tempwp(i  ,j,k) = tempwp(i,j,k)   + 0.5 * emom * ((w0(i,j,k)-w0(i-1,j,k))/dx)/dx - 0.5 * tau_wu_plus /dx
             

              emop = ( dzf(k) * ( ekm(i,j,k+1)  + ekm(i-1,j,k+1)  )  + &
                      dzf(k+1)  * ( ekm(i,j,k) + ekm(i-1,j,k) ) ) / &
                    ( 4.   * dzh(k+1) )

              v_at_w_min  = 0.25 * (v0(i-1,j,k)+v0(i-1,j,k+1)+v0(i-1,j+1,k)+v0(i-1,j+1,k+1) )
              v_at_w_plus = 0.25 * (v0(i  ,j,k)+v0(i  ,j,k+1)+v0(i  ,j+1,k)+v0(i  ,j+1,k+1) )
              call log_wallaw(w0(i-1,j,k+1),v_at_w_min ,Cm_xwall,tau_wu_min)
              call log_wallaw(w0(i  ,j,k+1),v_at_w_plus,Cm_xwall,tau_wu_plus)
              tempwp(i-1,j,k+1) = tempwp(i-1,j,k+1) - 0.5 * emom * ((w0(i,j,k+1)-w0(i-1,j,k+1))/dx)/dx - 0.5 * tau_wu_min/dx
              tempwp(i  ,j,k+1) = tempwp(i,j,k+1)   + 0.5 * emom * ((w0(i,j,k+1)-w0(i-1,j,k+1))/dx)/dx - 0.5 * tau_wu_plus /dx
              

              call xwallscalar(i,j,k,thl0,tempthlp)  ! zero subgrid flux through the boundary is applied here by subtraction
              call xwallscalar(i,j,k,qt0 ,tempqtp)

              us_at_scalar_min  = 0.5 * ((v0(i-1,j,k) + v0(i-1,j+1,k))**2 + (w0(i-1,j,k)+w0(i-1,j,k+1))**2)**0.5
              us_at_scalar_plus = 0.5 * ((v0(i  ,j,k) + v0(i  ,j+1,k))**2 + (w0(i  ,j,k)+w0(i  ,j,k+1))**2)**0.5 
              !cstep bulk_wall_temp applies bulk relation for heat flux at both sides of the wall
              ! for the ibm point this is no problem as it is overruled later in the code
              call bulk_wall_temp(us_at_scalar_min ,thlm(i-1,j,k),Cd_xwall,dx,tempthlp(i-1,j,k))  !cst
              call bulk_wall_temp(us_at_scalar_plus,thlm(i  ,j,k),Cd_xwall,dx,tempthlp(i  ,j,k))
 
              do nc=1,nsv
                call xwallscalar(i,j,k,sv0(:,:,:,nc),tempsvp(:,:,:,nc))
              end do
              !call xwalle12(i,j,k)   ! correction is ignored assuming u,v,w,subgrid TKE are near zero inside buildings
            endif

            if (lnorm_y(i,j,k)) then     !< Wall in y-direction
              emmo = 0.25  * ( &
                ekm(i,j,k)+ekm(i,j-1,k)+ekm(i-1,j-1,k)+ekm(i-1,j,k)  )

              w_at_u_min  = 0.25*(w0(i,j-1,k)+w0(i,j-1,k+1)+w0(i-1,j-1,k)+w0(i-1,j-1,k+1))
              w_at_u_plus = 0.25*(w0(i,j  ,k)+w0(i,j  ,k+1)+w0(i-1,j  ,k)+w0(i-1,j  ,k+1))
              call log_wallaw(u0(i,j-1,k),w_at_u_min ,Cm_ywall,tau_uv_min)
              call log_wallaw(u0(i,j  ,k),w_at_u_plus,Cm_ywall,tau_uv_plus)
              tempup(i,j-1,k) = tempup(i,j-1,k) - 0.5 * emmo * ((u0(i,j,k)-u0(i,j-1,k))/dy)/dy - 0.5 * tau_uv_min/dy
              tempup(i,j  ,k) = tempup(i,j,k)   + 0.5 * emmo * ((u0(i,j,k)-u0(i,j-1,k))/dy)/dy - 0.5 * tau_uv_plus/dy


              epmo = 0.25  * ( &
                ekm(i+1,j,k)+ekm(i+1,j-1,k)+ekm(i,j-1,k)+ekm(i,j,k)  )

              w_at_u_min  = 0.25*(w0(i+1,j-1,k)+w0(i+1,j-1,k+1)+w0(i,j-1,k)+w0(i,j-1,k+1))
              w_at_u_plus = 0.25*(w0(i+1,j  ,k)+w0(i+1,j  ,k+1)+w0(i,j  ,k)+w0(i,j  ,k+1))
              call log_wallaw(u0(i+1,j-1,k),w_at_u_min ,Cm_ywall,tau_uv_min)
              call log_wallaw(u0(i+1,j  ,k),w_at_u_plus,Cm_ywall,tau_uv_plus)

              tempup(i+1,j-1,k) = tempup(i+1,j-1,k) - 0.5 * epmo * ((u0(i+1,j,k)-u0(i+1,j-1,k))/dy)/dy - 0.5 * tau_uv_min/dy
              tempup(i+1,j,k)   = tempup(i+1,j,k)   + 0.5 * epmo * ((u0(i+1,j,k)-u0(i+1,j-1,k))/dy)/dy - 0.5 * tau_uv_plus/dy

              eomm = ( dzf(k-1) * ( ekm(i,j,k)  + ekm(i,j-1,k)  )  + &
                dzf(k) * ( ekm(i,j,k-1) + ekm(i,j-1,k-1) ) ) / ( 4.  * dzh(k) )

              u_at_w_min  = 0.25 * (u0(i,j-1,k-1)+v0(i,j-1,k)+u0(i+1,j-1,k-1)+u0(i+1,j-1,k) )
              u_at_w_plus = 0.25 * (u0(i,j  ,k-1)+u0(i,j  ,k)+u0(i+1,j  ,k-1)+u0(i+1,j  ,k) )
              call log_wallaw(w0(i,j-1,k),u_at_w_min ,Cm_ywall,tau_wv_min)
              call log_wallaw(w0(i,j  ,k),u_at_w_plus,Cm_ywall,tau_wv_plus)
              tempwp(i,j-1,k) = tempwp(i,j-1,k) - 0.5 * eomm * ((w0(i,j,k)-w0(i,j-1,k))/dy)/dy - 0.5 * tau_wv_min/dy
              tempwp(i,j,k)   = tempwp(i,j,k)   + 0.5 * eomm * ((w0(i,j,k)-w0(i,j-1,k))/dy)/dy - 0.5 * tau_wv_plus/dy

              eomp = ( dzf(k) * ( ekm(i,j,k+1)  + ekm(i,j-1,k+1)  )  + &
                dzf(k+1) * ( ekm(i,j,k) + ekm(i,j-1,k) ) ) / ( 4.  * dzh(k+1) )

              u_at_w_min  = 0.25 * (u0(i,j-1,k)+v0(i,j-1,k+1)+u0(i+1,j-1,k)+u0(i+1,j-1,k+1) )
              u_at_w_plus = 0.25 * (u0(i,j  ,k)+u0(i,j  ,k+1)+u0(i+1,j  ,k)+u0(i+1,j  ,k+1) )
              call log_wallaw(w0(i,j-1,k+1),u_at_w_min ,Cm_ywall,tau_wv_min)
              call log_wallaw(w0(i,j  ,k+1),u_at_w_plus,Cm_ywall,tau_wv_plus)
              tempwp(i,j-1,k+1) = tempwp(i,j-1,k+1) - 0.5 * eomm * ((w0(i,j,k+1)-w0(i,j-1,k+1))/dy)/dy - 0.5 * tau_wv_min/dy
              tempwp(i,j  ,k+1) = tempwp(i,j  ,k+1) + 0.5 * eomm * ((w0(i,j,k+1)-w0(i,j-1,k+1))/dy)/dy - 0.5 * tau_wv_plus/dy

              call ywallscalar(i,j,k,thl0,tempthlp)
              call ywallscalar(i,j,k,qt0 ,tempqtp)

              us_at_scalar_min  = 0.5 * ((u0(i,j-1,k) + u0(i+1,j-1,k))**2 + (w0(i,j-1,k)+w0(i,j-1,k+1))**2)**0.5
              us_at_scalar_plus = 0.5 * ((u0(i,j  ,k) + u0(i+1,j  ,k))**2 + (w0(i,j  ,k)+w0(i,j  ,k+1))**2)**0.5
              call bulk_wall_temp(us_at_scalar_min ,thlm(i,j-1,k),Cd_ywall,dy,tempthlp(i,j-1,k))
              call bulk_wall_temp(us_at_scalar_plus,thlm(i,j  ,k),Cd_ywall,dy,tempthlp(i,j  ,k))

              do nc=1,nsv
                call ywallscalar(i,j,k,sv0(:,:,:,nc),tempsvp(:,:,:,nc))
              end do
              !call ywalle12(i,j,k)
            endif


          end do  !k
        end do    !j
      end do      !i

      do i=2,i1
        do j=2,j1
          k=1
          if (lnorm_x(i,j,1)) then     !< Wall in x-direction
            emmo = 0.25  * ( &
                ekm(i,j,1)+ekm(i,j-1,1)+ekm(i-1,j-1,1)+ekm(i-1,j,1)  )

            w_at_v_min  = 0.25*(w0(i-1,j,k+1)+w0(i-1,j-1,k+1))  !w (k=1) = 0
            w_at_v_plus = 0.25*(w0(i  ,j,k+1)+w0(i  ,j-1,k+1))
            call log_wallaw(v0(i-1,j,k),w_at_v_min ,Cm_xwall,tau_vu_min)
            call log_wallaw(v0(i,j,k)  ,w_at_v_plus,Cm_xwall,tau_vu_plus)
            tempvp(i-1,j,k) = tempvp(i-1,j,k) - 0.5 * emmo*((v0(i,j,k)-v0(i-1,j,k))/dx) / dx - 0.5 * tau_vu_min /dx
            tempvp(i  ,j,k) = tempvp(i  ,j,k) + 0.5 * emmo*((v0(i,j,k)-v0(i-1,j,k))/dx) / dx - 0.5 * tau_vu_plus/dx

            empo = 0.25  * ( &
                ekm(i,j+1,1)+ekm(i,j,1)+ekm(i-1,j,1)+ekm(i-1,j+1,1)  )

              w_at_v_min  = 0.25*(w0(i-1,j+1,k+1)+w0(i-1,j,k+1))
              w_at_v_plus = 0.25*(w0(i  ,j+1,k+1)+w0(i  ,j,k+1))
              call log_wallaw(v0(i-1,j+1,k),w_at_v_min ,Cm_xwall,tau_vu_min)
              call log_wallaw(v0(i,j+1,k)  ,w_at_v_plus,Cm_xwall,tau_vu_plus)
              tempvp(i-1,j+1,k) = tempvp(i-1,j+1,k) - 0.5 * empo*((v0(i,j+1,k)-v0(i-1,j+1,k))/dx) / dx - 0.5 * tau_vu_min/dx
              tempvp(i  ,j+1,k) = tempvp(i,j+1,k)   + 0.5 * empo*((v0(i,j+1,k)-v0(i-1,j+1,k))/dx) / dx - 0.5 * tau_vu_plus  /dx


            call xwallscalar(i,j,1,thl0,tempthlp) 
            call xwallscalar(i,j,1,qt0, tempqtp)

              us_at_scalar_min  = 0.5 * ((v0(i-1,j,k) + v0(i-1,j+1,k))**2 + w0(i-1,j,k+1)**2)**0.5
              us_at_scalar_plus = 0.5 * ((v0(i  ,j,k) + v0(i  ,j+1,k))**2 + w0(i  ,j,k+1)**2)**0.5
              call bulk_wall_temp(us_at_scalar_min ,thlm(i-1,j,k),Cd_xwall,dx,tempthlp(i-1,j,k))
              call bulk_wall_temp(us_at_scalar_plus,thlm(i  ,j,k),Cd_xwall,dx,tempthlp(i  ,j,k))

            do nc=1,nsv
              call xwallscalar(i,j,1,sv0(:,:,:,nc),tempsvp(:,:,:,nc))
            end do
            !call xwalle12(i,j,1)
          endif

          if (lnorm_y(i,j,1)) then     !< Wall in y-direction
            emmo = 0.25  * ( &
              ekm(i,j,1)+ekm(i,j-1,1)+ekm(i-1,j-1,1)+ekm(i-1,j,1)  )

            w_at_u_min  = 0.25*(w0(i,j-1,k+1)+w0(i-1,j-1,k+1))
            w_at_u_plus = 0.25*(w0(i,j  ,k+1)+w0(i-1,j  ,k+1))
            call log_wallaw(u0(i,j-1,k),w_at_u_min ,Cm_ywall,tau_uv_min)
            call log_wallaw(u0(i,j  ,k),w_at_u_plus,Cm_ywall,tau_uv_plus)
            tempup(i,j-1,k) = tempup(i,j-1,k) - 0.5 * emmo * ((u0(i,j,k)-u0(i,j-1,k))/dy)/dy - 0.5 * tau_uv_min/dy
            tempup(i,j  ,k) = tempup(i,j,k)   + 0.5 * emmo * ((u0(i,j,k)-u0(i,j-1,k))/dy)/dy - 0.5 * tau_uv_plus/dy

            epmo = 0.25  * ( &
              ekm(i+1,j,1)+ekm(i+1,j-1,1)+ekm(i,j-1,1)+ekm(i,j,1)  )

             w_at_u_min  = 0.25*(w0(i+1,j-1,k+1)+w0(i,j-1,k+1))
             w_at_u_plus = 0.25*(w0(i+1,j  ,k+1)+w0(i,j  ,k+1))
             call log_wallaw(u0(i+1,j-1,k),w_at_u_min ,Cm_ywall,tau_uv_min)
             call log_wallaw(u0(i+1,j  ,k),w_at_u_plus,Cm_ywall,tau_uv_plus)

             tempup(i+1,j-1,k) = tempup(i+1,j-1,k) - 0.5 * epmo * ((u0(i+1,j,k)-u0(i+1,j-1,k))/dy)/dy - 0.5 * tau_uv_min/dy
             tempup(i+1,j,k)   = tempup(i+1,j,k)   + 0.5 * epmo * ((u0(i+1,j,k)-u0(i+1,j-1,k))/dy)/dy - 0.5 * tau_uv_plus/dy


            call ywallscalar(i,j,1,thl0,tempthlp)
            call ywallscalar(i,j,1,qt0 ,tempqtp)

             us_at_scalar_min  = 0.5 * ((u0(i,j-1,k) + u0(i+1,j-1,k))**2 + w0(i,j-1,k+1)**2)**0.5
             us_at_scalar_plus = 0.5 * ((u0(i,j  ,k) + u0(i+1,j  ,k))**2 + w0(i,j  ,k+1)**2)**0.5
             call bulk_wall_temp(us_at_scalar_min ,thlm(i,j-1,k),Cd_ywall,dy,tempthlp(i,j-1,k))
             call bulk_wall_temp(us_at_scalar_plus,thlm(i,j  ,k),Cd_ywall,dy,tempthlp(i,j  ,k))

            do nc=1,nsv
              call ywallscalar(i,j,1,sv0(:,:,:,nc),tempsvp(:,:,:,nc))
            end do
            !call ywalle12(i,j,1)

          endif


        end do
      end do
      do i=1,i1
        do j=1,j1
          do k=1,kmax
            if (libm(i,j,k)) then
              up(i,j,k)=-um(i,j,k)*rk3coefi
              vp(i,j,k)=-vm(i,j,k)*rk3coefi
              wp(i,j,k)=-wm(i,j,k)*rk3coefi
              thlp(i,j,k)=(thlibm-thlm(i,j,k))*rk3coefi  
              qtp (i,j,k)=(qtibm -qtm(i,j,k))*rk3coefi
              e12p(i,j,k)=(e12min-e12m(i,j,k))*rk3coefi
              do nc=1,nsv
                  svp(i,j,k,nc) = - svm(i,j,k,nc)*rk3coefi  
              enddo
            else
              up(i,j,k)=up(i,j,k)+tempup(i,j,k)
              vp(i,j,k)=vp(i,j,k)+tempvp(i,j,k)
              wp(i,j,k)=wp(i,j,k)+tempwp(i,j,k)
              thlp(i,j,k)=thlp(i,j,k)+tempthlp(i,j,k)
              qtp (i,j,k)=qtp (i,j,k)+tempqtp (i,j,k)
              do nc=1,nsv
                svp (i,j,k,nc) = svp(i,j,k,nc) + tempsvp(i,j,k,nc)
              enddo 
            endif
            !< Now the normal velocities
            if (lnorm_x(i,j,k)) then
              up(i,j,k)=-um(i,j,k)*rk3coefi
            end if
            if (lnorm_y(i,j,k)) then
              vp(i,j,k)=-vm(i,j,k)*rk3coefi
            end if
            if (lnorm_z(i,j,k)) then
              wp(i,j,k)=-wm(i,j,k)*rk3coefi  
            end if
          end do
        end do
      end do

    call excjs( up  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( vp  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( wp  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( e12p  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( thlp  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( qtp   , 2,i1,2,j1,1,k1,ih,jh)
    do nc=1,nsv
     ! svp(:,:,:,nc)=svp(:,:,:,nc)+tempsvp(:,:,:,nc)  !done above
      call excjs( svp(:,:,:,nc)  , 2,i1,2,j1,1,k1,ih,jh)
    enddo


    if(maxval(sqrt(u0**2.+v0**2.+w0**2.))>16) then
    !  maxlocx=maxloc(u0**2.+v0**2.+w0**2.)
    !  write(6,*) 'maxlocx = ',maxlocx(1),maxlocx(2),maxlocx(3)
    !  write(6,*) 'ERROR: vel>16, maxloc = ',maxloc(sqrt(u0**2.+v0**2.+w0**2.)), 'u0 and um and up are here:',u0(maxlocx(1),maxlocx(2),maxlocx(3)),um(maxlocx(1),maxlocx(2),maxlocx(3)),up(maxlocx(1),maxlocx(2),maxlocx(3))
    !  write(6,*) 'v0 and vm and vp are here:',v0(maxlocx(1),maxlocx(2),maxlocx(3)),vm(maxlocx(1),maxlocx(2),maxlocx(3)),vp(maxlocx(1),maxlocx(2),maxlocx(3))
    !  write(6,*) 'w0 and wm and wp are here:',w0(maxlocx(1),maxlocx(2),maxlocx(3)),wm(maxlocx(1),maxlocx(2),maxlocx(3)),wp(maxlocx(1),maxlocx(2),maxlocx(3))
    !  write(6,*) 'timee = ',timee
    endif

    !write (6,*) 'exit applyibm'
    return
  end subroutine applyibm


  subroutine zerowallvelocity(simid) !<- MK: Forcd velocity at the immersed boundaries to 0 for a better interaction with the poissonsolver

    use modfields,      only : um, vm, wm, up, vp, wp
    use modglobal,      only : rk3step, imax, jmax, kmax, i1, j1, k1, ih, jh, rdt
    use modmpi,         only : excjs   !cstep ,myidx,myidy
    !clater use modnudgeboundary, only : Nsim

    implicit none
    integer  :: i, j, k,ipos,jpos
    real     :: rk3coef,rk3coefi
    integer, intent(in) :: simid

    !cstep if (.not. (lapply_ibm)) return
    !clater if (.not. (lapply_ibm .and. simid==Nsim )) return  !added by bart w

    rk3coef = rdt / (4. - dble(rk3step))
    rk3coefi = 1. / rk3coef

    do i=1,i1
      do j=1,j1
        do k=1,kmax
          !libm ipos=i+myidx*imax
          !!libm jpos=j+myidy*jmax
          !libm if (limmersed_boundary(ipos,jpos,k)) then
          if (libm(i,j,k)) then
            up(i,j,k)=-um(i,j,k)*rk3coefi
            vp(i,j,k)=-vm(i,j,k)*rk3coefi
            wp(i,j,k)=-wm(i,j,k)*rk3coefi
          endif
          if (lnorm_x(i,j,k)) then
            up(i,j,k)=-um(i,j,k)*rk3coefi
          end if
          if (lnorm_y(i,j,k)) then
            vp(i,j,k)=-vm(i,j,k)*rk3coefi
          end if
          if (lnorm_z(i,j,k)) then
            wp(i,j,k)=-wm(i,j,k)*rk3coefi !added by Bart W
          end if
        end do
      end do
    end do
    call excjs( up  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( vp  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( wp  , 2,i1,2,j1,1,k1,ih,jh)

    return
  end subroutine zerowallvelocity


  subroutine xwallscalar(i,j,k,putin,putout)

    use modglobal,      only : ih, i1, jh, j1, k1, dx2i
    use modsubgriddata, only : ekh
    use modfields, only : sv0
!cstep    use modmpi, only : myid

    implicit none

    integer, intent(in)    :: i,j,k
    real,    intent(in)    :: putin(2-ih:i1+ih,2-jh:j1+jh,k1)
    real,    intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)

    ! zero flux through the boundary is applied here
    if(.not. (putin(i,j,k)==0)) then
      putout(i,j,k)   = putout(i,j,k) &
                      + 0.5 * (ekh(i,j,k)+ekh(i-1,j,k))*(putin(i,j,k)-putin(i-1,j,k))*dx2i
    elseif(.not. (putin(i-1,j,k)==0)) then
      putout(i-1,j,k) = putout(i-1,j,k) &
                      - 0.5 * (ekh(i,j,k)+ekh(i-1,j,k))*(putin(i,j,k)-putin(i-1,j,k))*dx2i
    endif

    return
  end subroutine xwallscalar

  subroutine ywallscalar(i,j,k,putin,putout)

    use modglobal,      only : ih, i1, jh, j1, k1, dy2i
    use modsubgriddata, only : ekh

    implicit none

    integer, intent(in)    :: i,j,k
    real,    intent(in)    :: putin(2-ih:i1+ih,2-jh:j1+jh,k1)
    real,    intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)

    ! zero flux through the boundary is applied here
    if(.not. (putin(i,j,k)==0)) then
      putout(i,j,k)   = putout(i,j,k) &
                      + 0.5 * (ekh(i,j,k)+ekh(i,j-1,k)) *(putin(i,j,k)-putin(i,j-1,k))*dy2i
    elseif(.not. (putin(i,j-1,k)==0)) then
      putout(i,j-1,k) = putout(i,j-1,k) &
                      - 0.5 * (ekh(i,j,k)+ekh(i,j-1,k)) *(putin(i,j,k)-putin(i,j-1,k))*dy2i
    endif

    return
  end subroutine ywallscalar


  subroutine xwalle12(i,j,k)

    use modglobal,      only : ih, i1, jh, j1, k1, dx2i, dx, dy, dzh
    use modsubgriddata, only : ekm
    use modfields,      only : e12p, e120, u0, v0, w0

    implicit none

    integer, intent(in)    :: i,j,k

    if(.not.(k==1)) then
      if(.not. (e12p(i,j,k)==0)) then
        e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
                                  + ekm(i,j,k)/(2*e120(i,j,k))* (&  !source terms
                                     -((w0(i,j,k+1)-w0(i-1,j,k+1))  / dx             + &
                                       (u0(i,j,k+1)-u0(i,j,k))      / dzh(k+1) )**2  + &
                                     +(2.*(w0(i,j,k+1))             / dx             + &
                                       (u0(i,j,k+1)-u0(i,j,k))      / dzh(k+1) )**2  + &

                                     -((w0(i,j,k)-w0(i-1,j,k))      / dx             + &
                                       (u0(i,j,k)-u0(i,j,k-1))      / dzh(k)   )**2  + &
                                     +(2.*(w0(i,j,k))               / dx             + &
                                       (u0(i,j,k)-u0(i,j,k-1))      / dzh(k)   )**2  + &

                                     -((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
                                       (v0(i,j+1,k)-v0(i-1,j+1,k))  / dx       )**2  + &
                                     +((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
                                       (2.*v0(i,j+1,k))             / dx       )**2  + &

                                     -((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
                                       (v0(i,j,k)-v0(i-1,j,k))      / dx       )**2  + &
                                     +((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
                                       (2.*v0(i,j,k))               / dx       )**2    &
                                    )
      elseif(.not. (e12p(i-1,j,k)==0)) then
        e12p(i-1,j,k) = e12p(i-1,j,k) + (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
                                  + ekm(i-1,j,k)/(2*e120(i-1,j,k))* (&  !source terms
                                       -((w0(i,j,k)-w0(i-1,j,k))    / dx             + &
                                         (u0(i,j,k)-u0(i,j,k-1))    / dzh(k)   )**2  + &
                                       +((-2.*w0(i-1,j,k))          / dx             + &
                                         (u0(i,j,k)-u0(i,j,k-1))    / dzh(k)   )**2  + &

                                       -((w0(i,j,k+1)-w0(i-1,j,k+1))/ dx             + &
                                         (u0(i,j,k+1)-u0(i,j,k))    / dzh(k+1) )**2  + &
                                       +((-2.*w0(i-1,j,k+1))        / dx             + &
                                         (u0(i,j,k+1)-u0(i,j,k))    / dzh(k+1) )**2  + &

                                       -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
                                       +((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (-2.*v0(i-1,j,k))          / dx       )**2  + &

                                       -((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
                                         (v0(i,j+1,k)-v0(i-1,j+1,k))/ dx       )**2  + &
                                       +((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
                                         (-2.*v0(i-1,j+1,k))        / dx       )**2    &
                                    )
    endif
    else !Special treatment for the lowest full level: k=1
      if(.not. (e12p(i,j,k)==0)) then
        e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
                                  + ekm(i,j,k)/(2*e120(i,j,k))* (&  !source terms
                                     -((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
                                       (v0(i,j+1,k)-v0(i-1,j+1,k))  / dx       )**2  + &
                                     +((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
                                       (2.*v0(i,j+1,k))             / dx       )**2  + &

                                     -((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
                                       (v0(i,j,k)-v0(i-1,j,k))      / dx       )**2  + &
                                     +((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
                                       (2.*v0(i,j,k))               / dx       )**2    &
                                    )
      elseif(.not. (e12p(i-1,j,k)==0)) then
        e12p(i-1,j,k) = e12p(i-1,j,k) + (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
                                  + ekm(i-1,j,k)/(2*e120(i-1,j,k))* (&  !source terms

                                       -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
                                       +((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (-2.*v0(i-1,j,k))          / dx       )**2  + &

                                       -((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
                                         (v0(i,j+1,k)-v0(i-1,j+1,k))/ dx       )**2  + &
                                       +((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
                                         (-2.*v0(i-1,j+1,k))        / dx       )**2    &
                                    )
      endif
    endif
  end subroutine xwalle12


  subroutine ywalle12(i,j,k)

    use modglobal,      only : ih, i1, jh, j1, k1, dy2i, dx, dy, dzh
    use modsubgriddata, only : ekm
    use modfields,      only : e12p, e120, u0, v0, w0

    implicit none

    integer, intent(in)    :: i,j,k

    if(.not.(k==1)) then
      if(.not. (e12p(i,j,k)==0)) then
        e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
                                  + ekm(i,j,k)/(2.*e120(i,j,k))* (&  !source terms
                                       -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
                                       +((2.*u0(i,j,k))             / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

                                       -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
                                       +((2.*u0(i+1,j,k))           / dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &

                                       -((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
                                         (w0(i,j,k+1)-w0(i,j-1,k+1))/ dy       )**2  + &
                                       +((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
                                         (2.*w0(i,j,k+1))           / dy       )**2  + &

                                       -((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
                                         (w0(i,j,k)-w0(i,j-1,k))    / dy       )**2  + &
                                       +((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
                                         (2.*w0(i,j,k))             / dy       )**2    &
                                    )
      elseif(.not. (e12p(i,j-1,k)==0)) then
        e12p(i,j-1,k) = e12p(i,j-1,k) + (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
                                  + ekm(i,j-1,k)/(2.*e120(i,j-1,k))* (&  !source terms
                                       -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
                                       +((-2.*u0(i,j-1,k))          / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

                                       -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
                                       +((-2.*u0(i+1,j-1,k))        / dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &

                                       -((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
                                         (w0(i,j,k)-w0(i,j-1,k))    / dy       )**2  + &
                                       +((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
                                         (-2.*w0(i,j-1,k))          / dy       )**2  + &

                                       -((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
                                         (w0(i,j,k+1)-w0(i,j-1,k+1))/ dy       )**2  + &
                                       +((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
                                         (-2.*w0(i,j-1,k+1))        / dy       )**2    &
                                    )
      endif
    else !Special treatment for the lowest full level: k=1
      if(.not. (e12p(i,j,k)==0)) then
        e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
                                  + ekm(i,j,k)/(2.*e120(i,j,k))* (&  !source terms
                                       -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
                                       +((2.*u0(i,j,k))             / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

                                       -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
                                       +((2.*u0(i+1,j,k))           / dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2    &
                                    )

      elseif(.not. (e12p(i,j-1,k)==0)) then
        e12p(i,j-1,k) = e12p(i,j-1,k) + (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
                                  + ekm(i,j-1,k)/(2.*e120(i,j-1,k))* (&  !source terms
                                       -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
                                       +((-2.*u0(i,j-1,k))          / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

                                       -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
                                       +((-2.*u0(i+1,j-1,k))        / dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2    &
                                    )
      endif
    endif
  end subroutine ywalle12

  subroutine log_wallaw(u1,u2,Cm_hor_wall,tau)
    implicit none
    
    real,intent(in) :: u1,u2,Cm_hor_wall
    real,intent(out) :: tau
    real :: uspd
   
    uspd = sqrt(u1**2 + u2**2)
    tau  = Cm_hor_wall * uspd * u1   !not a minus sign here but in the subroutine above, where it ensures force the direction of the wind
                                     !similar as michael who states "give tau the same sign as utan"     
!    write (6,*) 'u1,u2,dx_hor_half,Cm_hor_wall,tau',u1,u2,dx_hor_half,Cm_hor_wall,tau
    return
  end subroutine log_wallaw

  subroutine bulk_wall_temp(us,thl,Cd,dx,thlp)
    implicit none
    
    real,intent(in) :: us,thl,Cd,dx
    real,intent(out) :: thlp

    thlp = Cd * us * (thlwall - thl) / dx
!    write(6,*) 'us,thl,Cd,dx,thlp',us,thl,Cd,dx,thlp
    
    return
  end subroutine bulk_wall_temp

end module modibm
