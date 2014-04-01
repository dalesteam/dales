module modcanopy
  implicit none
  save

  ! Namoptions
  logical :: lcanopy = .false.       !< Switch to enable canopy representation
  integer :: ncanopy = 10            !< Amount of layers to represent the canopy
  real    :: cd      = 0.15          !< Drag coefficient in the canopy

  ! Fields
  real, allocatable :: padf(:)       !< plant area density field full level
  real, allocatable :: padh(:)       !< plant area density field full level

contains
!-----------------------------------------------------------------------------------------
  SUBROUTINE initcanopy
    use modmpi,    only : myid, mpi_logical, mpi_integer, my_real, comm3d, mpierr
    use modglobal, only : kmax,k1, ifnamopt, fname_options, ifoutput, cexpnr

    implicit none

    integer ierr
  
    namelist/NAMCANOPY/ lcanopy
  
    if(myid==0) then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMCANOPY,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMCANOPY'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMCANOPY'
      endif
      write(6 ,NAMCANOPY)
      close(ifnamopt)
  
      ncanopy = min(ncanopy,kmax)
    endif
  
  
    call MPI_BCAST(lcanopy   , 1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(ncanopy   , 1, mpi_integer , 0, comm3d, mpierr)
    call MPI_BCAST(cd        , 1, my_real     , 0, comm3d, mpierr)
  
    allocate(padf(ncanopy))
    allocate(padh(ncanopy+1))

    padf=!!1.0
    padh=!!1.0

    return
  end subroutine initcanopy
  
  subroutine canopy
    use modfields, only : up,vp,wp,e12p,thlp,qtp,sv0,svp
     use modglobal, only : nsv

    integer n
  
!!  Momentum affected by trees
    call canopyu(up)
    call canopyv(vp)
    call canopyw(wp)
!!  TKE affected by trees
    call canopye(e12p)
!!  Emissions of heat, moisture and scalars by trees
!    call canopyc(thl0,thlp,thlflux)
!    call canopyc(qt0 ,qtp ,qtflux )
!    do n=1,nsv
!      call canopyc(sv0(:,:,:,n),svp(:,:,:,n),svflux(:,:,n))
!    end do
  
    return
  end subroutine canopy
  
  subroutine exitcanopy
    implicit none
    deallocate(padf)
    deallocate(padh)
    return
  end subroutine exitcanopy
  
  subroutine canopyu (putout)
    use modglobal, only  : i1, i2, ih, j1, j2, jh, k1, cu, cv, dzh, dzf, imax, jmax
    use modfields, only  : u0, v0, w0, rhobf, rhobh
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ucor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: vcor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ftau  (imax,jmax)
    integer             :: k, kp

    ucor = u0 + cu
    vcor = v0 + cv

    do k=1,ncanopy
      kp   = k+1
      ftau = cd * padf(k) * sqrt(ucor(2:i1,2:j1,k)**2 +  &
                ((vcor(1:imax,2:j1,k)+vcor(2:i1,2:j1,k)+vcor(1:imax,3:j2,k)+vcor(2:i1,3:j2,k))/4)**2 + &
                ((dzh(kp)*(w0(1:imax,2:j1,k)+w0(2:i1,2:j1,k))+dzh(k)*(w0(1:imax,2:j1,kp)+w0(2:i1,2:j1,kp)))/(2*(dzh(k)+dzh(kp))))**2 &
                )
      
      putout(2:i1,2:j1,k) = putout(2:i1,2:j1,k) - ftau * ucor(2:i1,2:j1,k)
    end do

    return
  end subroutine canopyu
    
  subroutine canopyv (putout)
    use modglobal, only  : i1, i2, ih, j1, j2, jh, k1, cu, cv, dzh, dzf, imax, jmax
    use modfields, only  : u0, v0, w0, rhobf, rhobh
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ucor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: vcor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ftau  (imax,jmax)
    integer             :: k, kp

    ucor = u0 + cu
    vcor = v0 + cv

    do k=1,ncanopy
      kp   = k+1
      ftau = cd * padf(k) * sqrt(vcor(2:i1,2:j1,k)**2 +  &
                ((ucor(3:i2,2:j1,k)+ucor(2:i1,2:j1,k)+ucor(3:i2,1:jmax,k)+ucor(2:i1,1:jmax,k))/4)**2 + &
                ((dzh(kp)*(w0(2:i1,1:jmax,k)+w0(2:i1,2:j1,k))+dzh(k)*(w0(2:i1,1:jmax,kp)+w0(2:i1,2:j1,kp)))/(2*(dzh(k)+dzh(kp))))**2 &
                )
      
      putout(2:i1,2:j1,k) = putout(2:i1,2:j1,k) - ftau * vcor(2:i1,2:j1,k)
    end do

    return
  end subroutine canopyv
    
  subroutine canopyw (putout)
    use modglobal, only  : i1, i2, ih, j1, j2, jh, k1, cu, cv, dzh, dzf, imax, jmax
    use modfields, only  : u0, v0, w0, rhobf, rhobh
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ucor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: vcor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ftau  (imax,jmax)
    integer             :: k, km

    ucor = u0 + cu
    vcor = v0 + cv

    do k=2,(ncanopy+1)
      km   = k-1
      ftau = cd * padh(k) * sqrt(w0(2:i1,2:j1,k)**2 +  &
                ((dzf(km)*(ucor(2:i1,2:j1,k)+ucor(3:i2,2:j1,k))+dzf(k)*(ucor(2:i1,2:j1,km)+ucor(3:i2,2:j1,km)))/(4*dzh(k)))**2 + &
                ((dzf(km)*(vcor(2:i1,2:j1,k)+vcor(2:i1,3:j2,k))+dzf(k)*(vcor(2:i1,2:j1,km)+vcor(2:i1,3:j2,km)))/(4*dzh(k)))**2 &
                )
      
      putout(2:i1,2:j1,k) = putout(2:i1,2:j1,k) - ftau * w0(2:i1,2:j1,k)
    end do

    return
  end subroutine canopyw
  
  subroutine canopye (putout)
    use modglobal, only  : i1, i2, ih, j1, j2, jh, k1, cu, cv, dzh, dzf, imax, jmax
    use modfields, only  : u0, v0, w0, rhobf, rhobh, e120
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ucor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: vcor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ftau  (imax,jmax)
    integer             :: k, km

    ucor = u0 + cu
    vcor = v0 + cv

    do k=2,(ncanopy+1)
      km   = k-1
      ftau = cd * padh(k) * sqrt(w0(2:i1,2:j1,k)**2 +  &
                ((dzf(km)*(ucor(2:i1,2:j1,k)+ucor(3:i2,2:j1,k))+dzf(k)*(ucor(2:i1,2:j1,km)+ucor(3:i2,2:j1,km)))/(4*dzh(k)))**2 + &
                ((dzf(km)*(vcor(2:i1,2:j1,k)+vcor(2:i1,3:j2,k))+dzf(k)*(vcor(2:i1,2:j1,km)+vcor(2:i1,3:j2,km)))/(4*dzh(k)))**2 &
                )
      
      putout(2:i1,2:j1,k) = putout(2:i1,2:j1,k) - 2 * e120(2:i1,2:j1,k) * ftau
    end do

    return
  end subroutine canopye
  
! ======================================================================
! subroutine from Ned Patton
! ======================================================================
      subroutine spline(x,y,n,yp1,ypn,y2)
      integer n, nmax
      real yp1, ypn, x(n), y(n), y2(n)
      parameter (nmax=500000)
      integer i, k
      real p, qn, sig, un, u(nmax)
      if(yp1 > .99e30) then
         y2(1) = 0.0
         u(1)  = 0.0
      else
         y2(1) = -0.5
         u(1) = (3./(x(2) - x(1)))*((y(2) - y(1))/(x(2) - x(1)) - yp1)
      endif
      do i=2,n-1
         sig = (x(i) - x(i-1))/(x(i+1) - x(i-1))
         p = sig*y2(i-1) + 2.0
         y2(i) = (sig - 1.0)/p
         u(i) = (6.0*((y(i+1) - y(i))/(x(i+1) - x(i)) - (y(i) - y(i-1)) &
              /(x(i) - x(i-1)))/(x(i+1) - x(i-1)) - sig*u(i-1))/p
      enddo
      if(ypn > .99e+30) then
         qn = 0.0
         un = 0.0
      else
         qn = 0.5
         un = (3.0/(x(n) - x(n-1)))* (ypn - (y(n) - y(n-1))/(x(n) - x(n-1)))
      endif
      y2(n) = (un - qn*u(n-1))/(qn*y2(n-1) + 1.0)
      do k=n-1,1,-1
         y2(k) = y2(k)*y2(k+1) + u(k)
      enddo

      return
      end subroutine spline

! ======================================================================
! subroutine from Ned Patton
! ======================================================================

      subroutine splint(xa,ya,y2a,n,x,y)
      integer n
      real x,y, xa(n), y2a(n), ya(n)
      integer k,khi,klo
      real a,b,h
      klo = 1
      khi = n
    1 continue
        if(khi - klo > 1) then
           k = (khi + klo)/2
           if(xa(k) > x) then
              khi = k
           else
              klo = k
           endif
           go to 1
      endif
      h = xa(khi) - xa(klo)
      a = (xa(khi) - x)/h
      b = (x - xa(klo))/h
      y = a*ya(klo) + b*ya(khi) + ((a**3 - a)*y2a(klo) + (b**3 - b)*y2a(khi))*(h**2)/6.0

      return
      end subroutine splint

end module modcanopy
