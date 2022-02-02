!> \file modsimpleice.f90

!>
!!  Ice microphysics.
!>
!! Calculates ice microphysics in a cheap scheme without prognostic nr
!!  simpleice is called from *modmicrophysics*
!! \see  Grabowski, 1998, JAS and Khairoutdinov and Randall, 2006, JAS
!!  \author Steef B\"oing, TU Delft
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


module modsimpleice
  implicit none
  private
  public initsimpleice, exitsimpleice, simpleice
  real :: gamb1r
  real :: gambd1r
  real :: gamb1s
  real :: gambd1s
  real :: gamb1g
  real :: gambd1g
  real :: gam2dr
  real :: gam2ds
  real :: gam2dg
  real :: gammaddr3
  real :: gammadds3
  real :: gammaddg3
  contains

!> Initializes and allocates the arrays
  subroutine initsimpleice
    use modmicrodata, only : qr, qrp, nr, nrp, thlpmcr, qtpmcr, sed_qr, qr_spl, &
                             ilratio, rsgratio, sgratio, &
                             lambdar, lambdas, lambdag, &
                             qrmask, qcmask, precep, &
                             ccrz, ccsz, ccgz, bbg, bbr, bbs, ddg, ddr, dds

    use modglobal, only : ih,i1,jh,j1,k1,lacz_gamma

    implicit none

    allocate (qr(2:i1,2:j1,k1)        & ! qr (total precipitation!) converted from a scalar variable
             ,qrp(2:i1,2:j1,k1)       & ! qr tendency due to microphysics only, for statistics
             ,nr(2:i1,2:j1,k1)        & ! qr (total precipitation!) converted from a scalar variable
             ,nrp(2:i1,2:j1,k1)       & ! qr tendency due to microphysics only, for statistics
             ,thlpmcr(2:i1,2:j1,k1)   & ! thl tendency due to microphysics only, for statistics
             ,qtpmcr(2:i1,2:j1,k1)    & ! qt tendency due to microphysics only, for statistics
             ,sed_qr(2:i1,2:j1,k1)    & ! sedimentation rain droplets mixing ratio
             ,qr_spl(2:i1,2:j1,k1)    & ! time-splitting substep qr
             ,ilratio(2:i1,2:j1,k1)   & ! partition ratio cloud water vs cloud ice
             ,rsgratio(2:i1,2:j1,k1)  & ! partition ratio rain vs. snow/graupel
             ,sgratio(2:i1,2:j1,k1)   & ! partition ratio snow vs graupel
             ,lambdar(2:i1,2:j1,k1)   & ! slope parameter for rain
             ,lambdas(2:i1,2:j1,k1)   & ! slope parameter for snow
             ,lambdag(2:i1,2:j1,k1))    ! slope parameter for graupel

    allocate (qrmask(2:i1,2:j1,k1)    & ! mask for rain water
             ,qcmask(2:i1,2:j1,k1))     ! mask for cloud water

    allocate(precep(2:i1,2:j1,k1))      ! precipitation for statistics

    allocate(ccrz(k1),ccsz(k1),ccgz(k1))

     gamb1r=lacz_gamma(bbr+1)
     gambd1r=lacz_gamma(bbr+ddr+1)
     gamb1s=lacz_gamma(bbs+1)
     gambd1s=lacz_gamma(bbs+dds+1)
     gamb1g=lacz_gamma(bbg+1)
     gambd1g=lacz_gamma(bbg+ddg+1)
     gam2dr=lacz_gamma(2.5+0.5*ddr)
     gam2ds=lacz_gamma(2.5+0.5*dds)
     gam2dg=lacz_gamma(2.5+0.5*ddg)
     gammaddr3=lacz_gamma(3.+ddr)
     gammadds3=lacz_gamma(3.+dds)
     gammaddg3=lacz_gamma(3.+ddg)

  end subroutine initsimpleice

!> Cleaning up after the run
  subroutine exitsimpleice
    use modmicrodata, only : nr,nrp,qr,qrp,thlpmcr,qtpmcr,sed_qr,qr_spl, &
                             ilratio,rsgratio,sgratio,lambdar,lambdas,lambdag, &
                             qrmask,qcmask, &
                             precep, &
                             ccrz,ccsz,ccgz
    implicit none
    deallocate(nr,nrp,qr,qrp,thlpmcr,qtpmcr,sed_qr,qr_spl,ilratio,rsgratio,sgratio,lambdar,lambdas,lambdag)
    deallocate(qrmask,qcmask)
    deallocate(precep)
    deallocate(ccrz,ccsz,ccgz)
  end subroutine exitsimpleice

!> Calculates the microphysical source term.
  subroutine simpleice
    use modglobal, only : i1,j1,k1,rdt,rk3step,timee,rlv,cp,tup,tdn
    use modfields, only : sv0,svm,svp,qtp,thlp,ql0,exnf,rhof,tmp0,rhobf
    use modsimpleicestat, only : simpleicetend
    use modmicrodata, only : nr, nrp, iqr, qrp, sed_qr, qr_spl, qtpmcr, thlpmcr, delt, &
                             qcmask, qcmin, qrmask, qrmin, qr, &
                             ilratio, rsgratio, sgratio, &
                             aag, aar, aas, bbg, bbr, bbs, ccg, ccr, ccs, &
                             n0rg, n0rr, n0rs, &
                             tuprsg, tupsg, tdnrsg, tdnsg, &
                             ccgz, ccrz, ccsz, &
                             lambdag, lambdar, lambdas, &
                             l_graupel, l_rain, l_warm
    implicit none
    integer:: i,j,k
    real:: qrsmall, qrsum,qrtest

    delt = rdt/ (4. - dble(rk3step))

    ! used to check on negative qr and nr
    qrsum=0.
    qrsmall=0.
    ! reset microphysics tendencies
    qrp=0.
    nrp=0.
    nr=0.
    thlpmcr=0.
    qtpmcr=0.

    ! Density corrected fall speed parameters, see Tomita 2008
    do k=1,k1
    ccrz(k)=ccr*(1.29/rhobf(k))**0.5
    ccsz(k)=ccs*(1.29/rhobf(k))**0.5
    ccgz(k)=ccg*(1.29/rhobf(k))**0.5
    end do

    do k=1,k1
    do j=2,j1
    do i=2,i1
      ! initialise qr
      qr(i,j,k)= sv0(i,j,k,iqr)
      ! initialise qc mask
      if (ql0(i,j,k) > qcmin) then
        qcmask(i,j,k) = .true.
      else
        qcmask(i,j,k) = .false.
      end if
      ! initialise qr mask and check if we are not throwing away too much rain
      if (l_rain) then
        qrsum = qrsum+qr(i,j,k)
        if (qr(i,j,k) <= qrmin) then
          qrmask(i,j,k) = .false.
          if(qr(i,j,k)<0.) then
          qrsmall = qrsmall-qr(i,j,k)
          qr(i,j,k)=0.
          end if
        else
          qrmask(i,j,k)=.true.
        endif
      endif
    enddo
    enddo
    enddo

    if (qrsmall > 0.000001*qrsum) then
      write(*,*)'amount of neg. qr thrown away is too high  ',timee,' sec'
    end if


    if(l_warm) then !partitioning and determination of intercept parameter
      do k=1,k1
      do j=2,j1
      do i=2,i1
        ilratio(i,j,k)=1.   ! cloud water vs cloud ice partitioning
      enddo
      enddo
      enddo
    else
      do k=1,k1
      do j=2,j1
      do i=2,i1
        ilratio(i,j,k)=max(0.,min(1.,(tmp0(i,j,k)-tdn)/(tup-tdn)))! cloud water vs cloud ice partitioning
      enddo
      enddo
      enddo
    end if

    if(l_warm) then !partitioning and determination of intercept parameter
      do k=1,k1
      do j=2,j1
      do i=2,i1
        if(qrmask(i,j,k).eqv..true.) then
          rsgratio(i,j,k)=1.   ! rain vs snow/graupel partitioning
          sgratio(i,j,k)=0.   ! snow versus graupel partitioning
          lambdar(i,j,k)=(aar*n0rr*gamb1r/(rhof(k)*(qr(i,j,k)+1.e-6)))**(1./(1.+bbr)) ! lambda rain
          lambdas(i,j,k)=lambdar(i,j,k) ! lambda snow
          lambdag(i,j,k)=lambdar(i,j,k) ! lambda graupel
        end if
      enddo
      enddo
      enddo
    elseif(l_graupel) then
      do k=1,k1
      do j=2,j1
      do i=2,i1
        if(qrmask(i,j,k).eqv..true.) then
          rsgratio(i,j,k)=max(0.,min(1.,(tmp0(i,j,k)-tdnrsg)/(tuprsg-tdnrsg))) ! rain vs snow/graupel partitioning
          sgratio(i,j,k)=max(0.,min(1.,(tmp0(i,j,k)-tdnsg)/(tupsg-tdnsg))) ! snow versus graupel partitioning
          lambdar(i,j,k)=(aar*n0rr*gamb1r/(rhof(k)*(qr(i,j,k)*rsgratio(i,j,k)+1.e-6)))**(1./(1.+bbr)) ! lambda rain
          lambdas(i,j,k)=(aas*n0rs*gamb1s/(rhof(k)*(qr(i,j,k)*(1.-rsgratio(i,j,k))*(1.-sgratio(i,j,k))+1.e-6)))**(1./(1.+bbs)) ! snow
          lambdag(i,j,k)=(aag*n0rg*gamb1g/(rhof(k)*(qr(i,j,k)*(1.-rsgratio(i,j,k))*sgratio(i,j,k)+1.e-6)))**(1./(1.+bbg)) ! graupel
        endif
      enddo
      enddo
      enddo
    else
      do k=1,k1
      do j=2,j1
      do i=2,i1
        if(qrmask(i,j,k).eqv..true.) then
          rsgratio(i,j,k)=max(0.,min(1.,(tmp0(i,j,k)-tdnrsg)/(tuprsg-tdnrsg)))   ! rain vs snow/graupel partitioning
          sgratio(i,j,k)=0.
          lambdar(i,j,k)=(aar*n0rr*gamb1r/(rhof(k)*(qr(i,j,k)*rsgratio(i,j,k)+1.e-6)))**(1./(1.+bbr)) ! lambda rain
          lambdas(i,j,k)=(aas*n0rs*gamb1s/(rhof(k)*(qr(i,j,k)*(1.-rsgratio(i,j,k))+1.e-6)))**(1./(1.+bbs)) ! lambda snow
          lambdag(i,j,k)=lambdas(i,j,k)
        end if
      enddo
      enddo
      enddo
    endif

    if (l_rain) then
      call simpleicetend
      call autoconvert
      call simpleicetend
      call accrete
      call simpleicetend
      call evapdep
      call simpleicetend
      call precipitate
    endif

    do k=1,k1
    do j=2,j1
    do i=2,i1
      qrtest=svm(i,j,k,iqr)+(svp(i,j,k,iqr)+qrp(i,j,k))*delt
      if (qrtest .lt. qrmin) then ! correction, after Jerome's implementation in Gales
        qtp(i,j,k) = qtp(i,j,k) + qtpmcr(i,j,k) + svm(i,j,k,iqr)/delt + svp(i,j,k,iqr) + qrp(i,j,k)
        thlp(i,j,k) = thlp(i,j,k) +thlpmcr(i,j,k) - (rlv/(cp*exnf(k)))*(svm(i,j,k,iqr)/delt + svp(i,j,k,iqr) + qrp(i,j,k))
        svp(i,j,k,iqr) = - svm(i,j,k,iqr)/delt
      else
      svp(i,j,k,iqr)=svp(i,j,k,iqr)+qrp(i,j,k)
      thlp(i,j,k)=thlp(i,j,k)+thlpmcr(i,j,k)
      qtp(i,j,k)=qtp(i,j,k)+qtpmcr(i,j,k)
      ! adjust negative qr tendencies at the end of the time-step
     end if
    enddo
    enddo
    enddo

    if (l_rain) then
      call simpleicetend !after corrections
    endif
  end subroutine simpleice

  subroutine autoconvert
    use modglobal, only : i1,j1,k1,rlv,cp,tmelt
    use modfields, only : ql0,exnf,rhof,tmp0
    use modmicrodata, only : betakessi, delt, l_berry, Nc_0, qli0, qll0, timekessl, &
                             qcmask, qrp, qtpmcr, thlpmcr, ilratio
    implicit none
    real :: qll,qli,ddisp,lwc,autl,tc,times,auti,aut
    integer:: i,j,k

    if(l_berry.eqv..true.) then ! Berry/Hsie autoconversion
    do k=1,k1
    do j=2,j1
    do i=2,i1
        if (qcmask(i,j,k).eqv..true.) then
          ! ql partitioning
          qll=ql0(i,j,k)*ilratio(i,j,k)
          qli=ql0(i,j,k)-qll
          ddisp=0.146-5.964e-2*log(Nc_0/2.e9) ! Relative dispersion coefficient for Berry autoconversion
          lwc=1.e3*rhof(k)*qll ! Liquid water content in g/kg
          autl=1./rhof(k)*1.67e-5*lwc*lwc/(5. + .0366*Nc_0/(1.e6*ddisp*(lwc+1.e-6)))
          tc=tmp0(i,j,k)-tmelt ! Temperature wrt melting point
          times=min(1.e3,(3.56*tc+106.7)*tc+1.e3) ! Time scale for ice autoconversion
          auti=qli/times
          aut = min(autl + auti,ql0(i,j,k)/delt)
          qrp(i,j,k) = qrp(i,j,k)+aut
          qtpmcr(i,j,k) = qtpmcr(i,j,k)-aut
          thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlv/(cp*exnf(k)))*aut
        endif
      enddo
      enddo
      enddo
    else ! Lin/Kessler autoconversion as in Khairoutdinov and Randall, 2006
      do k=1,k1
      do j=2,j1
      do i=2,i1
        if (qcmask(i,j,k).eqv..true.) then
          ! ql partitioning
          qll=ql0(i,j,k)*ilratio(i,j,k)
          qli=ql0(i,j,k)-qll
          autl=max(0.,timekessl*(qll-qll0))
          tc=tmp0(i,j,k)-tmelt
          auti=max(0.,betakessi*exp(0.025*tc)*(qli-qli0))
          aut = min(autl + auti,ql0(i,j,k)/delt)
          qrp(i,j,k) = qrp(i,j,k)+aut
          qtpmcr(i,j,k) = qtpmcr(i,j,k)-aut
          thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlv/(cp*exnf(k)))*aut
        endif
      enddo
      enddo
      enddo
    endif

  end subroutine autoconvert

  subroutine accrete
    use modglobal, only : i1,j1,k1,rlv,cp,pi
    use modfields, only : ql0,exnf,rhof
    use modmicrodata, only : ddg, ddr, dds, aag, aar, aas, bbg, bbr, bbs, delt, &
                             lambdag, lambdar, lambdas, ccgz, ccrz, ccsz, &
                             ceffgi, ceffgl, ceffri, ceffrl, ceffsi, ceffsl, betakessi, &
                             qrmask, qcmask, qr, qrp, qtpmcr, thlpmcr, &
                             ilratio, rsgratio, sgratio
    implicit none
    real :: qll,qli,qrr,qrs,qrg,&
            gaccrl,gaccsl,gaccgl,gaccri,gaccsi,gaccgi,accr,accs,accg,acc
    integer:: i,j,k

    do k=1,k1
    do j=2,j1
    do i=2,i1
      if (qrmask(i,j,k).eqv..true.) then
      if (qcmask(i,j,k).eqv..true.) then ! apply mask
        ! ql partitioning
        qll=ql0(i,j,k)*ilratio(i,j,k)
        qli=ql0(i,j,k)-qll
        ! qr partitioning
        qrr=qr(i,j,k)*rsgratio(i,j,k)
        qrs=qr(i,j,k)*(1.-rsgratio(i,j,k))*(1.-sgratio(i,j,k))
        qrg=qr(i,j,k)*(1.-rsgratio(i,j,k))*sgratio(i,j,k)
        ! collection of cloud water by rain etc.
        gaccrl=pi/4.*ccrz(k)*ceffrl*rhof(k)*qll*qrr*lambdar(i,j,k)**(bbr-2.-ddr)*gammaddr3/(aar*gamb1r)
        gaccsl=pi/4.*ccsz(k)*ceffsl*rhof(k)*qll*qrs*lambdas(i,j,k)**(bbs-2.-dds)*gammadds3/(aas*gamb1s)
        gaccgl=pi/4.*ccgz(k)*ceffgl*rhof(k)*qll*qrg*lambdag(i,j,k)**(bbg-2.-ddg)*gammaddg3/(aag*gamb1g)
        gaccri=pi/4.*ccrz(k)*ceffri*rhof(k)*qli*qrr*lambdar(i,j,k)**(bbr-2.-ddr)*gammaddr3/(aar*gamb1r)
        gaccsi=pi/4.*ccsz(k)*ceffsi*rhof(k)*qli*qrs*lambdas(i,j,k)**(bbs-2.-dds)*gammadds3/(aas*gamb1s)
        gaccgi=pi/4.*ccgz(k)*ceffgi*rhof(k)*qli*qrg*lambdag(i,j,k)**(bbg-2.-ddg)*gammaddg3/(aag*gamb1g)
        accr=(gaccrl+gaccri)*qrr/(qrr+1.e-9)
        accs=(gaccsl+gaccsi)*qrs/(qrs+1.e-9)
        accg=(gaccgl+gaccgi)*qrg/(qrg+1.e-9)
        acc= min(accr+accs+accg,ql0(i,j,k)/delt)  ! total growth by accretion
        qrp(i,j,k) = qrp(i,j,k)+acc
        qtpmcr(i,j,k) = qtpmcr(i,j,k)-acc
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlv/(cp*exnf(k)))*acc
      end if
      end if
    enddo
    enddo
    enddo

  end subroutine accrete

  subroutine evapdep
    use modglobal, only : i1,j1,k1,rlv,cp,pi
    use modfields, only : qt0,ql0,exnf,rhof,tmp0,qvsl,qvsi,esl
    use modmicrodata, only : betag, betar, betas, ddg, ddr, dds, delt, &
                             n0rg, n0rr, n0rs, aag, aar, aas, &
                             ccgz, ccrz, ccsz, lambdag, lambdar, lambdas, &
                             evapfactor, qrmask, qr, qrp, qtpmcr, thlpmcr, &
                             qrmask, qcmask
    implicit none

    real :: ssl,ssi,ventr,vents,ventg,&
            thfun,evapdepr,evapdeps,evapdepg,devap
    integer:: i,j,k

    do k=1,k1
    do j=2,j1
    do i=2,i1
      if (qrmask(i,j,k).eqv..true.) then
        ! saturation ratios
        ssl=(qt0(i,j,k)-ql0(i,j,k))/qvsl(i,j,k)
        ssi=(qt0(i,j,k)-ql0(i,j,k))/qvsi(i,j,k)
        !integration over ventilation factors and diameters, see e.g. seifert 2008
        ventr=.78*n0rr/lambdar(i,j,k)**2 + gam2dr*.27*n0rr*sqrt(ccrz(k)/2.e-5)*lambdar(i,j,k)**(-2.5-0.5*ddr)
        vents=.78*n0rs/lambdas(i,j,k)**2 + gam2ds*.27*n0rs*sqrt(ccsz(k)/2.e-5)*lambdas(i,j,k)**(-2.5-0.5*dds)
        ventg=.78*n0rg/lambdag(i,j,k)**2 + gam2dg*.27*n0rg*sqrt(ccgz(k)/2.e-5)*lambdag(i,j,k)**(-2.5-0.5*ddg)
        thfun=1.e-7/(2.2*tmp0(i,j,k)/esl(i,j,k)+2.2e2/tmp0(i,j,k))  ! thermodynamic function
        evapdepr=(4.*pi/(betar*rhof(k)))*(ssl-1.)*ventr*thfun
        evapdeps=(4.*pi/(betas*rhof(k)))*(ssi-1.)*vents*thfun
        evapdepg=(4.*pi/(betag*rhof(k)))*(ssi-1.)*ventg*thfun
        ! total growth by deposition and evaporation
        ! limit with qr and ql after accretion and autoconversion
        devap= max(min(evapfactor*(evapdepr+evapdeps+evapdepg),ql0(i,j,k)/delt+qrp(i,j,k)),-qr(i,j,k)/delt-qrp(i,j,k))
        qrp(i,j,k) = qrp(i,j,k)+devap
        qtpmcr(i,j,k) = qtpmcr(i,j,k)-devap
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlv/(cp*exnf(k)))*devap
      end if
    enddo
    enddo
    enddo

  end subroutine evapdep

  subroutine precipitate
    use modglobal, only : i1,j1,k1,kmax,dzf,dzh
    use modfields, only : rhof,rhobf
    use modmicrodata, only : qr_spl, sed_qr, precep, qr, qrp, &
                             aag, aas, aar, bbg, bbs, bbr, ddg, dds, ddr, n0rg, n0rs, n0rr, &
                             betag, betar, betas, qrmin, &
                             lambdag, lambdar, lambdas, &
                             ccgz, ccrz, ccsz, &
                             sgratio, rsgratio, &
                             courantp, delt, qrmask
    implicit none
    integer :: i,j,k,jn
    integer :: n_spl      !<  sedimentation time splitting loop
    real :: dt_spl,wfallmax,vtr,vts,vtg,vtf

    wfallmax = 9.9
    n_spl = ceiling(wfallmax*delt/(minval(dzf)*courantp))
    dt_spl = delt/real(n_spl) !fixed time step

    sed_qr = 0. ! reset sedimentation fluxes

    do k=1,k1
    do j=2,j1
    do i=2,i1
      qr_spl(i,j,k) = qr(i,j,k)
      if (qrmask(i,j,k).eqv..true.) then
        vtr=ccrz(k)*(gambd1r/gamb1r)/(lambdar(i,j,k)**ddr)  ! terminal velocity rain
        vts=ccsz(k)*(gambd1s/gamb1s)/(lambdas(i,j,k)**dds)  ! terminal velocity snow
        vtg=ccgz(k)*(gambd1g/gamb1g)/(lambdag(i,j,k)**ddg)  ! terminal velocity graupel
        vtf=rsgratio(i,j,k)*vtr+(1.-rsgratio(i,j,k))*(1.-sgratio(i,j,k))*vts+(1.-rsgratio(i,j,k))*sgratio(i,j,k)*vtg ! weighted
        vtf = min(wfallmax,vtf)
        precep(i,j,k) = vtf*qr_spl(i,j,k)
        sed_qr(i,j,k) = precep(i,j,k)*rhobf(k) ! convert to flux
      else
        precep(i,j,k) = 0.
        sed_qr(i,j,k) = 0.
      end if
    enddo
    enddo
    enddo

    !  advect precipitation using upwind scheme
    do k=1,kmax
    do j=2,j1
    do i=2,i1
      qr_spl(i,j,k) = qr_spl(i,j,k) + (sed_qr(i,j,k+1) - sed_qr(i,j,k))*dt_spl/(dzh(k+1)*rhobf(k))
    enddo
    enddo
    enddo

    ! begin time splitting loop
    IF (n_spl > 1) THEN
      DO jn = 2 , n_spl

        ! reset fluxes at each step of loop
        sed_qr = 0.
        do k=1,k1
        do j=2,j1
        do i=2,i1
          if (qr_spl(i,j,k) > qrmin) then
            ! re-evaluate lambda
            lambdar(i,j,k)=(aar*n0rr*gamb1r/(rhof(k)*(qr_spl(i,j,k)*rsgratio(i,j,k)+1.e-6)))**(1./(1.+bbr)) ! lambda rain
            lambdas(i,j,k)=(aas*n0rs*gamb1s/(rhof(k)*(qr_spl(i,j,k)*(1.-rsgratio(i,j,k))*(1.-sgratio(i,j,k))+1.e-6)))**(1./(1.+bbs)) ! lambda snow
            lambdag(i,j,k)=(aag*n0rg*gamb1g/(rhof(k)*(qr_spl(i,j,k)*(1.-rsgratio(i,j,k))*sgratio(i,j,k)+1.e-6)))**(1./(1.+bbg)) ! lambda graupel
            vtr=ccrz(k)*(gambd1r/gamb1r)/(lambdar(i,j,k)**ddr)  ! terminal velocity rain
            vts=ccsz(k)*(gambd1s/gamb1s)/(lambdas(i,j,k)**dds)  ! terminal velocity snow
            vtg=ccgz(k)*(gambd1g/gamb1g)/(lambdag(i,j,k)**ddg)  ! terminal velocity graupel
            vtf=rsgratio(i,j,k)*vtr+(1.-rsgratio(i,j,k))*(1.-sgratio(i,j,k))*vts+(1.-rsgratio(i,j,k))*sgratio(i,j,k)*vtg  ! mass-weighted terminal velocity
            vtf=min(wfallmax,vtf)
            sed_qr(i,j,k) = vtf*qr_spl(i,j,k)*rhobf(k)
          else
            sed_qr(i,j,k) = 0.
          endif
        enddo
        enddo
        enddo

        do k=1,kmax
        do j=2,j1
        do i=2,i1
          qr_spl(i,j,k) = qr_spl(i,j,k) + (sed_qr(i,j,k+1) - sed_qr(i,j,k))*dt_spl/(dzh(k+1)*rhobf(k))
        enddo
        enddo
        enddo

      ! end time splitting loop and if n>1
      ENDDO
    ENDIF

    ! no thl and qt tendencies build in, implying no heat transfer between precipitation and air
    do k=1,kmax
    do j=2,j1
    do i=2,i1
      qrp(i,j,k)= qrp(i,j,k) + (qr_spl(i,j,k) - qr(i,j,k))/delt
    enddo
    enddo
    enddo

  end subroutine precipitate

end module modsimpleice
