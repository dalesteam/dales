!> \file modfft2d.f90
!>
!!  Solve the poisson equation using a FFT + tridiagonal solver. Uses FFTW.
!>
!!  \author Jisk Attema, Netherlands eScience Center
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
!  Copyright 2014 Netherlands eScience Center
!


module modfftw
! All use arguments must be outside of #ifdef because CMAKE doesnt read beyond
! it for dependency generation
use, intrinsic  :: iso_c_binding
use modprecision, only : pois_r
use modglobal, only : itot, jtot, imax, jmax, i1, j1, ih, jh, kmax, ijtot &
                    , dxi,dyi,pi, lperiodic
use modmpi, only : commcol, commrow, mpierr, nprocx, D_MPI_ALLTOALL, nprocy &
                 ,  myid, myidx, myidy
implicit none

#ifdef USE_FFTW

include 'fftw3.f03'

save
  integer                         :: method
  integer                         :: konx, kony
  integer                         :: iony, jonx
  real(pois_r), dimension(:), allocatable :: bufin, bufout
  interface fftw_plan_many_r2r_if
    procedure :: d_fftw_plan_many_r2r
    procedure :: d_fftwf_plan_many_r2r
  end interface
  interface fftw_plan_guru_r2r_if
    procedure :: d_fftw_plan_guru_r2r
    procedure :: d_fftwf_plan_guru_r2r
  end interface
  interface fftw_execute_r2r_if
    procedure :: fftw_execute_r2r
    procedure :: fftwf_execute_r2r
  end interface

  ! C pointer to the actual (aligned) memory for FFTW
  type (C_ptr)                    :: ptr

  ! Method 1:
  !   domain composition in two dimensions (nprocx, nprocy)
  !    - transpose [x,k]
  !    - fft over x
  !    - transpose [y,x]
  !    - fft over y
  !    - transpose [k,y]
  !    - apply tri-diagonal solver over k
  !    - and the whole fft/transpose in reverse + inverted order

  ! FFTW plan and memory pointers for 1D transforms
  ! NOTE: these pointers point to sub arrays of ptr
  type (C_ptr)                    :: planx, planxi, plany, planyi
  real(pois_r), pointer                   :: p210(:,:,:), p201(:,:,:)
  real(pois_r), pointer, contiguous       :: p210_flat(:), p201_flat(:)

  ! Method 2:
  !   no domain decomposition nprocx = nprocy = 1
  !    - fft over xy
  !    - apply tri-diagonal solver over k
  !    - inverse fft over xy

  ! FFTW plan and memory pointers for 2D transforms
  ! NOTE: these pointers point to sub arrays of ptr
  type (C_ptr)                    :: planxy, planxyi
  real(pois_r), pointer                   :: p_nohalo(:)

contains

  subroutine fftwinit(p, Fp, d, xyrt, ps,pe,qs,qe)

    implicit none

    real(pois_r), pointer              :: p(:,:,:)
    real(pois_r), pointer              :: Fp(:,:,:)
    real(pois_r), allocatable          :: d(:,:,:)
    real(pois_r), allocatable          :: xyrt(:,:)
    integer,intent(out)        :: ps,pe,qs,qe

    integer(kind=8)     :: sz
    real(pois_r), pointer,contiguous :: fptr(:)
    integer             :: embed(1), kinds(2)
    type (fftw_iodim)   :: dimij(2), dimk(1)

! setup the matrix transposes.
! konx and kony are the number of vertical (k) points per processor in the x and y direction.
! iony and jonx are the number of i (j) points per processor in the y (x) direction.
! it is of course best if the points can be distributed equally, but if not
! just let one row or column of processors do less points

    if (nprocx == 1 .and. nprocy == 1) then
      method = 2
    else
      method = 1
    endif

    konx = kmax / nprocx
    if ( mod(kmax, nprocx) > 0 ) then
      konx = konx + 1
    endif

    kony = kmax / nprocy
    if ( mod(kmax, nprocy) > 0 ) then
      kony = kony + 1
    endif

    iony = itot / nprocy
    if ( mod(itot, nprocy) > 0 ) then
      iony = iony + 1
    endif

    jonx = jtot / nprocx
    if ( mod(jtot, nprocx) > 0 ) then
      jonx = jonx + 1
    endif

    if (myid == 0) then
       write(*,*) "-- modFFTW sizes --"
       write(*,*) "itot, jtot", itot, jtot
       write(*,*) "imax, jmax", imax, jmax
       
       write(*,*) "konx", konx
       write(*,*) "kony", kony
       write(*,*) "iony", iony
       write(*,*) "jonx", jonx
    end if
! Allocate communication buffers for the transpose functions
    sz = max( imax * jmax * konx * nprocx, & ! transpose a1
              iony * jmax * konx * nprocy, & ! transpose a2
              iony * jonx * konx * nprocx  ) ! transpose a3

    allocate(bufin (sz))
    bufin = 0 ! because some elements may not be used, but will be FFTed
    allocate(bufout(sz))
    if (myid == 0) write(*,*) "bufin, bufout sz", sz

! Allocate temporary arrays to hold transposed matrix

    ! calculate memory size as an long int (size_t)
    sz = max( (imax+2*ih) * (jmax+2*jh) * kmax, & ! p (p012)
              itot*jmax*konx,                   & ! p210
              jtot*konx*iony,                   & ! p201
              iony*jonx*kmax                    ) ! Fp (p102)

    ! get aligned memory for FFTW
    ptr = fftw_alloc_real(sz)
    
    if (myid == 0) write(*,*) "ptr sz", sz

    ! convert it to a fortran pointer, or 1D array
    call c_f_pointer(ptr, fptr, (/sz/))

    p(2-ih:i1+ih,2-jh:j1+jh,1:kmax) => fptr(1:(imax+2*ih)*(jmax+2*jh)*kmax)

    if (method == 1) then

    ! convert it to our 3D arrays with custom ubounds/lbounds
    p210(1:itot,1:jmax,1:konx) => fptr(1:itot*jmax*konx)
    p210_flat(1:itot*jmax*konx)=> fptr(1:itot*jmax*konx)
    p201(1:jtot,1:konx,1:iony) => fptr(1:jtot*konx*iony)
    p201_flat(1:jtot*konx*iony)=> fptr(1:jtot*konx*iony)
    Fp(1:iony,1:jonx,1:kmax) => fptr(1:iony*jonx*kmax)

    p = 0           ! because some elements may not be used, but will be FFTed
    p201_flat = 0   ! don't know which is the largest...
    p210_flat = 0
        
    ! Prepare 1d FFT transforms
    ! TODO: in plan_many, skip part where k > kmax
    embed(1) = itot
    if (lperiodic(1)) then
        kinds(1) = FFTW_R2HC
     else
        kinds(1) = FFTW_REDFT10 ! DCT for Neumann BC for pressure
     end if

    planx = fftw_plan_many_r2r_if( &
      1, &           ! rank
      embed, &       ! n (size)  [array]
      jmax*konx, &   ! howmany
      p210_flat,  &       ! array; location of transform k is: in + k * idist
      embed, &       ! inembed: subrank (halo) [array]
      1, &           ! istride
      itot, &        ! idist
      p210_flat, &        ! fftw_double *out
      embed, &       ! onembed [array]
      1, &           ! ostride
      itot, &        ! odist
      kinds, &       ! kind
      FFTW_MEASURE & ! flags (FFTW_MEASURE or FFTW_ESTIMATE)
    )

    embed(1) = itot
    if (lperiodic(1)) then
        kinds(1) = FFTW_HC2R
     else
        kinds(1) = FFTW_REDFT01 ! Inverse DCT
     end if

    planxi = fftw_plan_many_r2r_if( &
      1, &           ! rank
      embed, &       ! n (size)  [array]
      jmax*konx, &   ! howmany
      p210_flat,  &       ! array; location of transform k is: in + k * idist
      embed, &       ! inembed: subrank (halo) [array]
      1, &           ! istride
      itot, &        ! idist
      p210_flat, &        ! fftw_double *out
      embed, &       ! onembed [array]
      1, &           ! ostride
      itot, &        ! odist
      kinds, &       ! kind
      FFTW_MEASURE & ! flags (FFTW_MEASURE or FFTW_ESTIMATE)
    )

    embed(1) = jtot
    if (lperiodic(3)) then
       kinds(1) = FFTW_R2HC
    else
       kinds(1) = FFTW_REDFT10
    end if

    plany = fftw_plan_many_r2r_if( &
      1, &           ! rank
      embed, &       ! n (size)  [array]
      konx*iony, &   ! howmany
      p201_flat,  &       ! array; location of transform k is: in + k * idist
      embed, &       ! inembed: subrank (halo) [array]
      1, &           ! istride
      jtot, &        ! idist
      p201_flat, &        ! fftw_double *out
      embed, &       ! onembed [array]
      1, &           ! ostride
      jtot, &        ! odist
      kinds, &       ! kind
      FFTW_MEASURE & ! flags (FFTW_MEASURE or FFTW_ESTIMATE)
    )

    embed(1) = jtot
     if (lperiodic(3)) then
        kinds(1) = FFTW_HC2R
     else
        kinds(1) = FFTW_REDFT01 ! Inverse DCT
     end if
    planyi = fftw_plan_many_r2r_if( &
      1, &           ! rank
      embed, &       ! n (size)  [array]
      konx*iony, &   ! howmany
      p201_flat,  &       ! array; location of transform k is: in + k * idist
      embed, &       ! inembed: subrank (halo) [array]
      1, &           ! istride
      jtot, &        ! idist
      p201_flat, &        ! fftw_double *out
      embed, &       ! onembed [array]
      1, &           ! ostride
      jtot, &        ! odist
      kinds, &       ! kind
      FFTW_MEASURE & ! flags (FFTW_MEASURE or FFTW_ESTIMATE)
    )

    if ((.not. c_associated(planx)) .or. (.not. c_associated(plany)) &
        .or. (.not. c_associated(planxi)) .or. (.not. c_associated(planyi))) then

       STOP "FFTW plan creation failed (method 1)."
    end if

    allocate(xyrt(iony,jonx))
    allocate(d(iony,jonx,kmax))
    ps = 1
    pe = iony
    qs = 1
    qe = jonx

    ! special case for final task in x or y, which may have fewer elements
    if (myidx == nprocx-1) qe = jtot - (nprocx-1)*jonx
    if (myidy == nprocy-1) pe = itot - (nprocy-1)*iony
    
    
    else if (method == 2) then

    ! Make an array with the halo zone sliced off on the top and the left
    ! it practically starts at p(2,2,1) but fortran doesnt compile that
    p_nohalo(1:sz-ih-jh*(itot+2*ih)) => fptr(1+ih+jh*(itot+2*ih):sz)

    ! Prepare 2d FFT transforms

    dimij(1)%n = itot
    dimij(1)%is = 1
    dimij(1)%os = 1

    dimij(2)%n = jtot
    dimij(2)%is = itot + 2*ih
    dimij(2)%os = itot + 2*ih

    dimk(1)%n = kmax
    dimk(1)%is = (jtot + 2*jh) * (itot + 2*ih)
    dimk(1)%os = (jtot + 2*jh) * (itot + 2*ih)

    kinds(1) = FFTW_R2HC
    kinds(2) = FFTW_R2HC
    if (.not. lperiodic(1)) kinds(1) = FFTW_REDFT10
    if (.not. lperiodic(3)) kinds(2) = FFTW_REDFT10

    planxy = fftw_plan_guru_r2r_if(&
      2, &             ! rank
      dimij, &         ! dims
      1, &             ! howmany_rank
      dimk, &          ! howmany_dims
      p_nohalo, &      ! fftw_double *in
      p_nohalo, &      ! fftw_double *out
      kinds, &         ! kind
      FFTW_MEASURE &   ! flags (FFTW_MEASURE or FFTW_ESTIMATE)
    )

    kinds(1) = FFTW_HC2R
    kinds(2) = FFTW_HC2R
    if (.not. lperiodic(1)) kinds(1) = FFTW_REDFT01
    if (.not. lperiodic(3)) kinds(2) = FFTW_REDFT01

    planxyi = fftw_plan_guru_r2r_if(&
      2, &             ! rank
      dimij, &         ! dims
      1, &             ! howmany_rank
      dimk, &          ! howmany_dims
      p_nohalo, &      ! fftw_double *in
      p_nohalo, &      ! fftw_double *out
      kinds, &         ! kind
      FFTW_MEASURE &   ! flags (FFTW_MEASURE or FFTW_ESTIMATE)
    )

    if ((.not. c_associated(planxy)) .or. (.not. c_associated(planxyi))) then
       STOP "FFTW plan creation failed (method 2)."
    end if

    allocate(xyrt(2-ih:i1+ih,2-jh:j1+jh))
    allocate(   d(2-ih:i1+ih,2-jh:j1+jh,kmax))
    Fp(2-ih:i1+ih,2-jh:j1+jh,1:kmax) => fptr(1:(imax+2*ih)*(jmax+2*jh)*kmax)

    ps = 2
    pe = i1
    qs = 2
    qe = j1

    else
      stop 'Illegal method in fftwinit.'
    endif

    call fftwinit_factors(xyrt)

 end subroutine

 subroutine fftwexit(p,Fp,d,xyrt)
   implicit none
   real(pois_r), pointer :: p(:,:,:)
   real(pois_r), pointer :: Fp(:,:,:)
   real(pois_r), allocatable :: d(:,:,:)
   real(pois_r), allocatable :: xyrt(:,:)

   if (method == 1) then
     call fftw_destroy_plan(planx)
     call fftw_destroy_plan(planxi)
     call fftw_destroy_plan(plany)
     call fftw_destroy_plan(planyi)
   else if (method == 2) then
     call fftw_destroy_plan(planxy)
     call fftw_destroy_plan(planxyi)
   else
     stop 'Illegal method in fftwexit.'
   endif

   ! ptr, planx, planxi, plany, planyi are C pointers,
   ! so Nullify() doesnt work with them
   Nullify(p, p210, p201, Fp, p_nohalo)

   call fftw_free(ptr)

   deallocate(xyrt, d)
   deallocate(bufin,bufout)

 end subroutine

! Transpose functions:
!   Data is stored such that a whole column, 1:kmax, is on one node (pencil layout).
!   These funtions redistribute the data over the MPI nodes so that a whole dimension
!   is on one processor.
!   To facilitate processing after the MPI transpose, data is further transposed locally
!   such that the complete dimension, one of itot, jtot, kmax, is consecutive (ie. fastest).
!
!   The figure below shows the transpose between dimension 0 and 1, or 'k' and 'i':
!
!
!         /-------------/|                                /-------------/|
!        /../          / |                               /             / |
!  kmax  |------------|  |                               |------------|  |
!        |..|         |  |            ==>                |            |  |
!        |..|         |  |                         konx  |____________|. |
!        |..|         | /    jmax                        |............|./    jmax
!   1    |--|--|--|---|/   1                        1    |------------|/   1
!
!        1  imax                                         1            itot
!
!  transpose_a1: dimensions 0 and 2:
!            p012(2-ih:i1+ih,2-jh:j1+jh,kmax) <=> p210(itot,jmax,konx)
!
!  transpose_a2: dimensions 1 and 2:
!            p210(      itot,      jmax,konx) <=> p201(jtot,konx,iony)
!
!  transpose_a3: dimensions 0 and 2:
!            p201(      jtot,      konx,iony) <=> p102(iony,jonx,kmax)
!
  subroutine transpose_a1(p,p210)
    implicit none

    !real, intent(in)    :: p(2-ih:i1+ih,2-jh:j1+jh,kmax)
    !real, intent(out)   :: p210(itot,jmax,konx)
    real(pois_r), pointer :: p(:,:,:)
    real(pois_r), pointer :: p210(:,:,:)

    integer :: n, i,j,k, ii

    ii = 0
    do n=0,nprocx-1
    do k=n*konx + 1, (n+1)*konx
    do j=2,j1
    do i=2,i1
      ii = ii + 1
      if (k <= kmax) then
        bufin(ii) = p(i,j,k)
      endif
    enddo
    enddo
    enddo
    enddo

    call D_MPI_ALLTOALL(bufin,   (imax*jmax*konx), &
                        bufout,  (imax*jmax*konx), &
                        commrow,mpierr)

    ii = 0
    do n=0,nprocx-1
    do k=1,konx
    do j=1,jmax
    do i=n*imax + 1, (n+1)*imax
        ii = ii + 1
        p210(i,j,k) = bufout(ii)
    enddo
    enddo
    enddo
    enddo

  end subroutine

  subroutine transpose_a1inv(p,p210)
    implicit none

    !real, intent(out)   :: p(2-ih:i1+ih,2-jh:j1+jh,kmax)
    !real, intent(in)    :: p210(itot,jmax,konx)
    real(pois_r), pointer :: p(:,:,:)
    real(pois_r), pointer :: p210(:,:,:)

    integer :: n, i,j,k, ii

    ii = 0
    do n=0,nprocx-1
    do k=1,konx
    do j=1,jmax
    do i=n*imax + 1, (n+1)*imax
      ii = ii + 1
      bufin(ii) = p210(i,j,k)
    enddo
    enddo
    enddo
    enddo

    call D_MPI_ALLTOALL(bufin,   (imax*jmax*konx), &
                        bufout,  (imax*jmax*konx), &
                        commrow,mpierr)

    ii = 0
    do n=0,nprocx-1
    do k=n*konx + 1,(n+1)*konx
    do j=2,j1
    do i=2,i1
        ii = ii + 1
        if (k <= kmax) then
          p(i,j,k) = bufout(ii)
        endif
    enddo
    enddo
    enddo
    enddo

  end subroutine

  subroutine transpose_a2(p210, p201)
    implicit none

    !real, intent(in)    :: p210(itot,jmax,konx)
    !real, intent(out)   :: p201(jtot,konx,iony)
    real(pois_r), pointer :: p210(:,:,:)
    real(pois_r), pointer :: p201(:,:,:)

    integer :: n, i,j,k, ii

    ii = 0
    do n=0,nprocy-1
    do k=1,konx
    do j=1,jmax
    do i=n*iony + 1,(n+1)*iony
      ii = ii + 1
      if (i <= itot) then
        bufin(ii) = p210(i,j,k)
      endif
    enddo
    enddo
    enddo
    enddo

    call D_MPI_ALLTOALL(bufin,   (iony*jmax*konx), &
                        bufout,  (iony*jmax*konx), &
                        commcol,mpierr)

    ii = 0
    do n=0,nprocy-1
    do k=1,konx
    do j=n*jmax+1,(n+1)*jmax
    do i=1,iony
        ii = ii + 1
        p201(j,k,i) = bufout(ii)
    enddo
    enddo
    enddo
    enddo

  end subroutine

  subroutine transpose_a2inv(p210, p201)
    implicit none

    !real, intent(out)  :: p210(itot,jmax,konx)
    !real, intent(in)   :: p201(jtot,konx,iony)
    real(pois_r), pointer :: p210(:,:,:)
    real(pois_r), pointer :: p201(:,:,:)

    integer :: n, i,j,k, ii

    ii = 0
    do n=0,nprocy-1
    do k=1,konx
    do j=n*jmax + 1,(n+1)*jmax
    do i=1,iony
      ii = ii + 1
      bufin(ii) = p201(j,k,i)
    enddo
    enddo
    enddo
    enddo

    call D_MPI_ALLTOALL(bufin,   (iony*jmax*konx), &
                        bufout,  (iony*jmax*konx), &
                        commcol,mpierr)

    ii = 0
    do n=0,nprocy-1
    do k=1,konx
    do j=1,jmax
    do i=n*iony+1,(n+1)*iony
      ii = ii + 1
      if (i <= itot) then
        p210(i,j,k) = bufout(ii)
      else
        p210(i,j,k) = 0
      endif
    enddo
    enddo
    enddo
    enddo

  end subroutine

  subroutine transpose_a3(p201, Fp)
    implicit none

    !real, intent(in)    :: p201(jtot,konx,iony)
    !real, intent(out)   :: Fp(iony,jonx,kmax)
    real(pois_r), pointer :: p201(:,:,:)
    real(pois_r), pointer :: Fp(:,:,:)

    integer :: n, i,j,k, ii

    ii = 0
    do n=0,nprocx-1
    do k=1,konx
    do j=n*jonx+1,(n+1)*jonx
    do i=1,iony
      ii = ii + 1
      if (j <= jtot) then
        bufin(ii) = p201(j,k,i)
      else
        bufin(ii) = 0
      endif
    enddo
    enddo
    enddo
    enddo

    call D_MPI_ALLTOALL(bufin,   (iony*jonx*konx), &
                        bufout,  (iony*jonx*konx), &
                        commrow,mpierr)

    ii = 0
    do n=0,nprocx-1
    do k=n*konx+1,(n+1)*konx
    do j=1,jonx
    do i=1,iony
        ii = ii + 1
        if (k <= kmax) then
          Fp(i,j,k) = bufout(ii)
        endif
    enddo
    enddo
    enddo
    enddo

  end subroutine

  subroutine transpose_a3inv(p201, Fp)
    implicit none

    !real, intent(out)   :: p201(jtot,konx,iony)
    !real, intent(in)    :: Fp(iony,jonx,kmax)
    real(pois_r), pointer :: p201(:,:,:)
    real(pois_r), pointer :: Fp(:,:,:)

    integer :: n, i,j,k, ii

    ii = 0
    do n=0,nprocx-1
    do k=n*konx+1,(n+1)*konx
    do j=1,jonx
    do i=1,iony
      ii = ii + 1
      if (k <= kmax) then
        bufin(ii) = Fp(i,j,k)
      endif
    enddo
    enddo
    enddo
    enddo

    call D_MPI_ALLTOALL(bufin,   (iony*jonx*konx), &
                        bufout,  (iony*jonx*konx), &
                        commrow,mpierr)

    ii = 0
    do n=0,nprocx-1
    do k=1,konx
    do j=n*jonx+1,(n+1)*jonx
    do i=1,iony
        ii = ii + 1
        if (j <= jtot) then
          p201(j,k,i) = bufout(ii)
        endif
    enddo
    enddo
    enddo
    enddo

  end subroutine

  subroutine fftwf(p, Fp)
    implicit none

    real(pois_r), pointer :: p(:,:,:)
    real(pois_r), pointer :: Fp(:,:,:)

    if (method == 1) then
      call transpose_a1(p, p210)
      call fftw_execute_r2r_if(planx, p210_flat, p210_flat)

      call transpose_a2(p210, p201)
      call fftw_execute_r2r_if(plany, p201_flat, p201_flat)

      call transpose_a3(p201, Fp)
    else if (method == 2) then
      call fftw_execute_r2r_if(planxy, p_nohalo, p_nohalo)
    else
      stop 'Illegal method in fftwsolver.'
    endif

    !Fp(:,:,:) = Fp(:,:,:) / sqrt(ijtot)  ! moving normalization to fftwb
  end subroutine

  subroutine fftwb(p, Fp)
    implicit none

    real(pois_r), pointer :: p(:,:,:)
    real(pois_r), pointer :: Fp(:,:,:)
    real :: norm

    !Fp(:,:,:) = Fp(:,:,:) / sqrt(ijtot)

    if (method == 1) then
      call transpose_a3inv(p201, Fp)

      call fftw_execute_r2r_if(planyi, p201_flat, p201_flat)
      call transpose_a2inv(p210, p201)

      call fftw_execute_r2r_if(planxi, p210_flat, p210_flat)
      call transpose_a1inv(p, p210)

    else if (method == 2) then
      call fftw_execute_r2r_if(planxyi, p_nohalo, p_nohalo)
    else
      stop 'Illegal method in fftwsolver.'
   endif

   norm = 1.0/ijtot
   if (.not. lperiodic(1)) norm = norm / 2 ! different normalization for the DCT
   if (.not. lperiodic(3)) norm = norm / 2
   p(:,:,:) = p(:,:,:) * norm ! do all normalization at once
  end subroutine

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fftwinit_factors(xyrt)
    use modglobal, only : i1,j1,kmax,imax,jmax,itot,jtot,dxi,dyi,pi,ih,jh,lperiodic
    use modmpi, only    : myidx, myidy

    implicit none

    real(pois_r), allocatable :: xyrt(:,:)

    integer      :: i,j,iv,jv
    real(pois_r) :: fac
    real(pois_r) :: xrt(max(itot,nprocy*iony)), yrt(max(jtot,nprocx*jonx))

    xrt = 0
    yrt = 0
    
  ! Generate Eigenvalues xrt and yrt resulting from d**2/dx**2 F = a**2 F

  ! We are using FFTW's half-complex format. From the FFTW documentation:
  !
  !            r0, r1, r2, ..., rn/2, i(n+1)/2-1, ..., i2, i1
  !
  ! Here, rk is the real part of the kth output, and ik is the imaginary part. (Division by 2 is rounded down.)
  ! For a halfcomplex array hc[n], the kth component thus has its real part in hc[k] and its imaginary part in hc[n-k], with the
  ! exception of k == 0 or n/2 (the latter only if n is even)—in these two cases, the imaginary part is zero due to symmetries of the
  ! real-input DFT, and is not stored. Thus, the r2hc transform of n real values is a halfcomplex array of length n, and vice versa
  ! for hc2r.

  ! Fortran i:  1   2   3   4   5   6   7   8
  ! n=8:       r0, r1, r2, r3, r4, i3, i2, i1
  ! n=7:       r0, r1, r2, r3, i3, i2, i1

  ! I --> direction
    fac = 1./(2.*itot)
    if (.not. lperiodic(1)) fac = 1./(4.*itot)
    do i=2,itot
      xrt(i)=-4.*dxi*dxi*(sin(float(2*(i-1))*pi*fac))**2
    end do
    xrt(1) = 0.

  ! J --> direction
    fac = 1./(2.*jtot)
    if (.not. lperiodic(3)) fac = 1./(4.*jtot)
    do j=2,jtot
      yrt(j)=-4.*dyi*dyi*(sin(float(2*(j-1))*pi*fac))**2
    end do
    yrt(1) = 0.

  ! Combine I and J directions
  ! Note that:
  ! 1. MPI starts counting at 0 so it should be myidy * jmax
  ! 2. real data, ie. no halo, starts at index 2 in the array xyrt(2,2) <-> xrt(1), yrt(1)

  ! the last tile in x or y may have fewer elements than the rest
  ! filling xyrt will then access xrt or yrt beyond itot or jtot. xrt and yrt are padded above.
    
    if (method == 1) then
      ! nprocx /= 1, nprocy /= 1
      ! we will be working on matrix Fp(1:iony,1:jonx,1:kmax)
      do j=1,jonx
      do i=1,iony
        iv = myidy * iony + i
        jv = myidx * jonx + j
        xyrt(i,j)=(xrt(iv)+yrt(jv))
      enddo
      enddo
    else if (method == 2) then
      ! nprocx = nprocy = 1
      ! we will be working on matrix p(2-ih:i1+ih,2-jh:j1+ih,1:kmax)
      xyrt = 0.
      do j=2,j1
      do i=2,i1
        iv = i - 1
        jv = j - 1
        xyrt(i,j)=(xrt(iv)+yrt(jv))
      enddo
      enddo
    else
      stop 'Illegal method in fftwinit_factors.'
    endif
  end subroutine fftwinit_factors

  subroutine D_fftw_execute_r2r(p, in, out)
      implicit none
      type(C_PTR) :: p
      real(C_DOUBLE), pointer, contiguous, intent(inout) :: in(:)
      real(C_DOUBLE), pointer, contiguous, intent(out)   :: out(:)
      call fftw_execute_r2r(p, in, out)
  end subroutine

  subroutine D_fftwf_execute_r2r(p, in, out)
      implicit none
      type(C_PTR) :: p
      real(C_FLOAT), pointer, contiguous, intent(inout) :: in(:)
      real(C_FLOAT), pointer, contiguous, intent(out)   :: out(:)
      call fftwf_execute_r2r(p, in, out)
  end subroutine

  type(c_ptr) function D_fftw_plan_many_r2r(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,kind,flags)
      implicit none
      integer(C_INT) :: rank
      integer(C_INT), intent(in) :: n(:)
      integer(C_INT) :: howmany
      real(C_DOUBLE), pointer, intent(out) :: in(:)
      integer(C_INT), intent(in) :: inembed(:)
      integer(C_INT) :: istride
      integer(C_INT) :: idist
      real(C_DOUBLE), pointer, intent(out) :: out(:)
      integer(C_INT), intent(in) :: onembed(:)
      integer(C_INT) :: ostride
      integer(C_INT) :: odist
      integer(C_FFTW_R2R_KIND), intent(in) :: kind(:)
      integer(C_INT) :: flags

      D_fftw_plan_many_r2r = fftw_plan_many_r2r(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,kind,flags)
  end function

  type(c_ptr) function D_fftwf_plan_many_r2r(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,kind,flags)
      implicit none
      integer(C_INT) :: rank
      integer(C_INT), intent(in) :: n(:)
      integer(C_INT) :: howmany
      real(C_FLOAT), pointer, intent(out) :: in(:)
      integer(C_INT), intent(in) :: inembed(:)
      integer(C_INT) :: istride
      integer(C_INT) :: idist
      real(C_FLOAT), pointer, intent(out) :: out(:)
      integer(C_INT), intent(in) :: onembed(:)
      integer(C_INT) :: ostride
      integer(C_INT) :: odist
      integer(C_FFTW_R2R_KIND), intent(in) :: kind(:)
      integer(C_INT) :: flags

      D_fftwf_plan_many_r2r = fftwf_plan_many_r2r(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,kind,flags)
  end function

    type(C_PTR) function d_fftw_plan_guru_r2r(rank,dims,howmany_rank,howmany_dims,in,out,kind,flags)
      integer(C_INT) :: rank
      type(fftw_iodim), intent(in) :: dims(:)
      integer(C_INT) :: howmany_rank
      type(fftw_iodim), intent(in) :: howmany_dims(:)
      real(C_DOUBLE), intent(out) :: in(:)
      real(C_DOUBLE), intent(out) :: out(:)
      integer(C_FFTW_R2R_KIND), intent(in) :: kind(:)
      integer(C_INT) :: flags
      d_fftw_plan_guru_r2r = fftw_plan_guru_r2r(rank,dims,howmany_rank,howmany_dims,in,out,kind,flags)
    end function

    type(C_PTR) function d_fftwf_plan_guru_r2r(rank,dims,howmany_rank,howmany_dims,in,out,kind,flags)
      integer(C_INT) :: rank
      type(fftw_iodim), intent(in) :: dims(:)
      integer(C_INT) :: howmany_rank
      type(fftw_iodim), intent(in) :: howmany_dims(:)
      real(C_FLOAT), intent(out) :: in(:)
      real(C_FLOAT), intent(out) :: out(:)
      integer(C_FFTW_R2R_KIND), intent(in) :: kind(:)
      integer(C_INT) :: flags

      !fftw_iodim and fftwf_iodim are exactly the same, fortran needs some convincing
      type(fftwf_iodim) :: dimsf(size(dims))
      type(fftwf_iodim) :: howmany_dimsf(size(howmany_dims))
      dimsf = transfer(dims,dimsf)
      howmany_dimsf = transfer(howmany_dims,dimsf)

      d_fftwf_plan_guru_r2r = fftwf_plan_guru_r2r(rank,dimsf,howmany_rank,howmany_dimsf,in,out,kind,flags)
    end function


#else

contains

  subroutine fftwinit(p, Fp, d, xyrt, ps,pe,qs,qe)
    real(pois_r), pointer      :: p(:,:,:)
    real(pois_r), pointer      :: Fp(:,:,:)
    real(pois_r), allocatable  :: d(:,:,:)
    real(pois_r), allocatable  :: xyrt(:,:)
    integer,intent(out)        :: ps,pe,qs,qe
    call error_and_exit()
    ps=0 ! suppress warnings about intent(out) variables not being assigned
    pe=0
    qs=0
    qe=0
 end subroutine

 subroutine fftwexit(p,Fp,d,xyrt)
    real(pois_r), pointer     :: p(:,:,:)
    real(pois_r), pointer     :: Fp(:,:,:)
    real(pois_r), allocatable :: d(:,:,:)
    real(pois_r), allocatable :: xyrt(:,:)
    call error_and_exit()
 end subroutine

  subroutine fftwf(p, Fp)
    real(pois_r), pointer :: p(:,:,:)
    real(pois_r), pointer :: Fp(:,:,:)
    call error_and_exit()
  end subroutine

  subroutine fftwb(p, Fp)
    real(pois_r), pointer :: p(:,:,:)
    real(pois_r), pointer :: Fp(:,:,:)
    call error_and_exit()
  end subroutine

  subroutine error_and_exit
    write (*,*) 'DALES was compiled without FFTW.'
    write (*,*) 'Use the default poisson solver (solver_id=0),'
    write (*,*) 'or recompile DALES with the option USE_FFTW.'

    call exit(-1)
  end subroutine

#endif

end module modfftw
