!> \file moduser.f90
!! A dummy file for cases where one wants additional forcings
!----------------------------------------------------------------------------
! This file is part of DALES.
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
! Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!----------------------------------------------------------------------------
!
!
module moduser
implicit none
real :: bflx0
contains
subroutine force_user
  implicit none
  end subroutine force_user

subroutine rad_user
  implicit none
end subroutine rad_user

subroutine micro_user
  implicit none
end subroutine micro_user

subroutine initsurf_user
    use modglobal,   only : tmelt,bt,at,rd,rv,cp,es0,pref0
    use modsurfdata, only : ps,qts,thvs,thls, wtsurf
    real       :: exner, tsurf, es
    bflx0      = wtsurf
    exner      = (ps / pref0)**(rd/cp)
    tsurf      = thls * exner
    es         = es0 * exp(at*(tsurf-tmelt) / (tsurf-bt))
    qts        = rd / rv * es / ps
    thvs = thls * (1. + (rv/rd - 1.) * qts)
end subroutine initsurf_user

subroutine surf_user
 use modglobal,  only : zf,i1,j1,i2,j2,grav,nsv,fkar,cv,cu,tmelt,bt,at,rd,rv,cp,es0,pref0,ep2,rslabs
 use modsurfdata,only : ustar,dudz,dvdz,dqtdz,dthldz,&
                          svs,z0,qts,thls,thvs,thlflux,qtflux,wsvsurf,svflux, ps, wtsurf, wqsurf
  use modfields, only : u0,v0,thl0,qt0,sv0,u0av,v0av,qt0av,thl0av, sv0av
  use modmpi,    only :  excj
  implicit none
  integer i, j, n
  real phihzf, phimzf
  real upcu, vpcv
  real horv2, horv, stab, obl
  real, parameter :: v_bulk = 0.01
  real       :: es,bflx

    thls = 289.
    es   = es0 * exp(at*(thls-tmelt) / (thls-bt))
    qts  = rd / rv * es / ps   
    bflx = grav/thl0av(1)*v_bulk*((thls - thl0av(1))+ep2*thl0av(1)*(qts-qt0av(1)))

    do while ((bflx0-bflx)/bflx0 > 1e-4)
      thls = thls + 0.001
      es   = es0 * exp(at*(thls-tmelt) / (thls-bt))
      qts  = rd / rv * es / ps
      bflx = grav/thl0av(1)*v_bulk*((thls - thl0av(1))+ep2*thl0av(1)*(qts-qt0av(1)))
    end do
    wtsurf = v_bulk * (thls - thl0av(1))
    wqsurf = v_bulk * (qts - qt0av(1))

  do j=2,j1
  do i=2,i1
    thlflux(i,j) = v_bulk * (thls - thl0(i,j,1))
    qtflux(i,j)  = v_bulk * (qts - qt0(i,j,1))
    bflx         = thlflux(i,j)*grav/thl0av(1) + grav*ep2*qtflux(i,j)

    upcu  = 0.5*(u0(i,j,1)+u0(i+1,j,1))+cu
    vpcv  = 0.5*(v0(i,j,1)+v0(i,j+1,1))+cv
    horv  = max(0.1,sqrt(upcu**2 + vpcv**2))
    horv2 = (upcu**2 + vpcv**2)
    ustar (i,j) = diag_ustar(zf(1),z0,bflx,horv)
    obl   = -ustar(i,j)**3/(fkar*(grav/thvs)*(thlflux(i,j)+0.61*thvs*qtflux(i,j)))
    
    if (bflx > 0.) then
       phimzf = (1.-16.*zf(1)/obl)**(-0.25)
       phihzf = (1.-16.*zf(1)/obl)**(-0.50)
    endif

    if (bflx == 0.) then
       phimzf = 1.
       phihzf = 1.
    endif

    if (bflx < 0.) then
       phimzf = (1.+5.*zf(1)/obl)
       phihzf = (1.+8.*zf(1)/obl)
    endif
    


    dudz(i,j)   = ustar(i,j)*(phimzf/(fkar*zf(1)))*(upcu/horv)
    dvdz(i,j)   = ustar(i,j)*(phimzf/(fkar*zf(1)))*(vpcv/horv)
    dthldz(i,j) = - thlflux(i,j) / ustar(i,j) * phihzf / (fkar*zf(1))
    dqtdz (i,j) = - qtflux(i,j)  / ustar(i,j) * phihzf / (fkar*zf(1))

  end do
  end do

  do j=2,j1
    ustar(1,j)=ustar(i1,j)
  end do

  call excj( ustar  , 1, i2, 1, j2, 1,1)
  do n=1,nsv
    svflux(:,:,n) = wsvsurf(n)
    svs(n)        = wsvsurf(n)/v_bulk+ sv0av(1,n)
  enddo

  contains
  !
  ! ----------------------------------------------------------------------
  ! FUNCTION GET_USTAR:  returns value of ustar using the below 
  ! similarity functions and a specified buoyancy flux (bflx) given in
  ! kinematic units
  !
  ! phi_m (zeta > 0) =  (1 + am * zeta)
  ! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
  !
  ! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
  !
  ! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface 
  ! Layer, in Workshop on Micormeteorology, pages 67-100.
  !
  ! Code writen March, 1999 by Bjorn Stevens
  !
  real function diag_ustar(z,z0,bflx,wnd)

    implicit none

    real, parameter      :: am   =  4.8   !   "          "         "
    real, parameter      :: bm   = 19.3   !   "          "         "
    real, parameter      :: eps  = 1.e-10 ! non-zero, small number
    real, parameter      :: vonk = 0.40
    real, intent (in)    :: z             ! height where u locates
    real, intent (in)    :: z0            ! momentum roughness height
    real, intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
    real, intent (in)    :: wnd           ! wind speed at z

    integer :: iterate
    real    :: lnz, klnz, c1, x, psi1, zeta, lmo, ustar

    lnz   = log(z/z0) 
    klnz  = vonk/lnz              
    c1    = 3.14159/2. - 3.*log(2.)

    ustar =  wnd*klnz
    if (bflx /= 0.0) then 
       do iterate=1,4
          lmo   = -(ustar**3)/(bflx*vonk + eps)
          zeta  = z/lmo
          if (zeta > 0.) then
             ustar =  vonk*wnd  /(lnz + am*zeta)
          else
             x     = sqrt( sqrt( 1.0 - bm*zeta ) )
             psi1  = 2.*log(1.0+x) + log(1.0+x*x) - 2.*atan(x) + c1
             ustar = wnd*vonk/(lnz - psi1)
          end if
       end do
    end if

    diag_ustar = ustar

    return
  end function diag_ustar
end subroutine surf_user
end module moduser
