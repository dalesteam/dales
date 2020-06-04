program tests
  ! include minimum number of fields and function to get a poisson solver running
  use modmpi, only: initmpicomm, commwrld, initmpi, myid, nprocx, nprocy
  use modglobal, only : ifnamopt, fname_options, version
  use modglobal, only : dxi,dyi

  ! variables below are set/overwritten by this test program
  use modglobal, only : itot,jtot,ijtot,i1,j1,k1,imax,jmax,kmax,ih,jh,solver_id
  use modversion, only : git_version

  implicit none

  integer :: ierr
  logical :: verbose

  namelist/DOMAIN/ &
    itot,jtot,kmax

  namelist/RUN/ &
    nprocx,nprocy

  namelist/TEST/ &
    verbose, ih, jh

  !==========================
  !          SETUP
  !==========================

  if (command_argument_count() >=1) then
    call get_command_argument(1,fname_options)
  else
    stop 'ERROR:No namelist'
  end if

  open(ifnamopt,file=fname_options,status='old',iostat=ierr)
  if (ierr /= 0) then
    stop 'ERROR:Namoptions does not exist'
  end if

  read (ifnamopt,DOMAIN,iostat=ierr)
  rewind(ifnamopt)
  read (ifnamopt,RUN,iostat=ierr)
  rewind(ifnamopt)
  read (ifnamopt,TEST,iostat=ierr)
  close(ifnamopt)

  call initmpicomm
  call initmpi

  ijtot = real(itot*jtot)
  imax = itot/nprocx
  jmax = jtot/nprocy
  i1=imax+1
  j1=jmax+1
  k1=kmax+1
  dxi = 1.
  dyi = 1.

  if (myid == 0 .and. verbose) then
    write (*, *) trim(version)//' git: '//trim(git_version)
    write (*, RUN)
    write (*, DOMAIN)
    write (*, TEST)
    write (*,*)
    write (*,*)
  endif

  call compare_ffts

end program tests


subroutine compare_ffts
  use modglobal, only : solver_id, ih,jh,i1,j1,kmax
  use modpois, only : initpois, exitpois, poisson, p, Fp, ps, pe, qs, qe
  use modfftw, only : method, fftwf, fftwb
  use modfft2d,  only : fft2df, fft2db
  use mpi, only : MPI_Wtime
  use modmpi, only: commwrld, myid
  implicit none
  character(20) logfilename
  real :: largest_diff, t1, t0
  real :: dummy, m
  real, allocatable :: res(:,:,:)
  integer :: loc(3)
  integer :: ierr
  integer :: logfile

  integer :: i,j,k

  write(logfilename, '(a,i0,a)') 'output.', myid, '.txt'
  open (logfile,file=logfilename,status='new')

  allocate(res(2-ih:i1+ih,2-jh:j1+jh,1:kmax))

  ! ==========================================================================
  solver_id = 0
  call initpois
  call fill_halo_array

  t0 = MPI_Wtime()
  call fft2df(p, Fp)
  dummy = maxval(p) - maxval(Fp)
  call fft2db(p, Fp)
  t1 = MPI_Wtime()

  res = p

  call check_halo_array(dummy)
  write (*,*) 'existing poisson solver, largets diff', dummy

  if (myid == 0) then
    write (*,*) 'Testing the existing poisson solver, solver_id = ', solver_id
    write (*,*) 'call took', t1 - t0
  endif

  call exitpois
  ! ==========================================================================

  call MPI_Barrier(commwrld, ierr)
  write (*,*)
  write (*,*)

  ! ==========================================================================
  solver_id = 100
  call initpois
  call fill_halo_array

  t0 = MPI_Wtime()
  call fftwf(p, Fp)
  dummy = maxval(p) - maxval(Fp)
  call fftwb(p, Fp)
  t1 = MPI_Wtime()

  dummy = -1.

  do k=1,kmax
  do i=2-ih,i1+ih
  do j=2-jh,j1+jh
    if (i<2 .or. i>i1 .and. j<2 .and. j>j1) then
      res(i,j,k) = 0.
    else if (abs(res(i,j,k)) .gt. dummy) then
      dummy = abs(res(i,j,k))
      loc(1) = i
      loc(2) = j
      loc(3) = k
    endif

  enddo
  enddo
  enddo

  call check_halo_array(dummy)
  write (*,*) 'fftw method, largets diff', method, dummy

  if (myid == 0) then
    write (*,*) 'Testing the FFTW poisson solver, solver_id, method = ', solver_id, method
    write (*,*) 'call took', t1 - t0
  endif

  write (*,*) 'Between methods id, p, res', myid, p(loc(1), loc(2), loc(3)), res(loc(1), loc(2), loc(3)), abs(p(loc(1), loc(2), loc(3)) - res(loc(1), loc(2), loc(3)))

  call exitpois
  ! ==========================================================================

  close(logfile)
  call MPI_Barrier(commwrld, ierr)
end subroutine compare_ffts

real function indexarray(gi,gj,gk)
  implicit none
  integer, intent(in) :: gi, gj, gk
  indexarray = sin( 1.0*(gi * 1000 + gj) * 1000 + gk ) + gk
  return
endfunction indexarray

subroutine fill_halo_array
  use modpois, only : p
  use modglobal, only : itot, jtot, i1, j1, imax, jmax, kmax, ih, jh
  use modmpi, only: myidx, myidy, myid

  implicit none
  real :: indexarray

  integer :: i,j,k
  integer :: gi, gj, gk

  do i=2-ih,i1+ih
  do j=2-jh,j1+jh
  do k=1,kmax
    gi = myidx * imax + i - 1
    gj = myidy * jmax + j - 1
    gk = k

    p(i,j,k) = indexarray(gi, gj, gk)
  enddo
  enddo
  enddo
end subroutine fill_halo_array

subroutine check_halo_array(largest_diff)
  use modpois, only : p
  use modglobal, only : i1, j1, imax, jmax, kmax, ih, jh, itot, jtot
  use modmpi, only: myidx, myidy

  implicit none
  real, intent(out)   :: largest_diff

  integer :: i,j,k
  integer :: gi, gj, gk
  integer :: li, lj, lk
  real :: diff

  real :: indexarray

  largest_diff = 0.

  do i=2,i1
  do j=2,j1
  do k=1,kmax
    gi = myidx * imax + i - 1
    gj = myidy * jmax + j - 1
    gk = k

    if (gi >= 1 .and. gj >= 1 .and. gk >= 1 .and. gi <= itot .and. gj <= jtot .and. gk <= kmax) then
      diff = abs(p(i,j,k) - indexarray(gi, gj, gk))
      if (diff > largest_diff) then
        largest_diff = diff
      endif
    endif
  enddo
  enddo
  enddo

end subroutine check_halo_array

subroutine print_halo_array(logfile, array_name, a)
  use modmpi, only: myid, myidx, myidy, nprocs
  use modglobal, only : i1,j1,k1,imax,jmax,kmax,ih,jh
  implicit none
  integer, intent(in) :: logfile
  real, intent(in) :: a(2-ih:i1+ih, 2-jh:j1+jh, 1:kmax)
  character(len = 8) :: array_name

  integer :: i,j,k

  write (logfile,*) array_name, 'sizes', lbound(a), ubound(a)
  do i=2-ih,i1+ih
  do j=2-jh,j1+jh
  do k=1,kmax
    write (logfile,*) array_name, i, j, k, a(i,j,k)
  enddo
  enddo
  enddo
end subroutine print_halo_array


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine poisson
    use modglobal, only : solver_id
    use mpi, only : MPI_Wtime
    use modmpi, only : myid
    use modhypre, only : solve_hypre
    use modfftw, only : fftwf, fftwb
    use modfft2d,  only : fft2df, fft2db

    use modglobal, only : i1,j1,ih,jh,kmax, itot, jtot
    implicit none

    real, allocatable ::  p_fft2d(:,:,:)
    real, allocatable ::  p_fftw(:,:,:)
    real, allocatable :: filled(:,:,:)
    integer :: i,j,k
    integer :: mi,mj,mk
    real :: maxdiff

    allocate(p_fft2d(2-ih:i1+ih,2-jh:j1+jh,1:kmax))
    allocate( p_fftw(2-ih:i1+ih,2-jh:j1+jh,1:kmax))
    allocate( filled(2-ih:i1+ih,2-jh:j1+jh,1:kmax))

    call exitpois


    ! ==== Use the old method

    solver_id = 0
    call initpois

    call fillps
    filled = p

    ! Forward FFT
    call fft2df(p, Fp)

    call solmpj

    ! Backward FFT
    call fft2db(p, Fp)

    p_fft2d = p

    call exitpois

    ! ==== Use the new method
    solver_id = 100
    call initpois

    p = filled

    ! Forward FFT
    call fftwf(p, Fp)

    call solmpj

    ! Backward FFT
    call fftwb(p, Fp)
    p_fftw = p

    ! ==== Look at diff

    write (*,*) '=========='

    maxdiff = 0.
    mi = -1
    mj = -1
    mk = -1
    do i=2,i1
    do j=2,j1
    do k=1,kmax
      if (abs(p_fft2d(i,j,k) - p_fftw(i,j,k)) > maxdiff) then
        maxdiff = p_fft2d(i,j,k) - p_fftw(i,j,k)
        mi = i
        mj = j
        mk = k
      endif
    enddo
    enddo
    enddo

    write (*,*) 'myid=',myid, 'i=',mi, 'j=',mj, 'k=',mk, 'd=',maxdiff

    call tderive

    deallocate(p_fft2d)
    deallocate(p_fftw)

  end subroutine poisson
