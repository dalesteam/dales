!> \file modchecksim.f90
!!  Monitors Courant and Peclet numbers, and divergence.

!>
!!  Monitors Courant and Peclet numbers, and divergence.
!>
!!  These numbers are put out to screen either every tcheck seconds, or every time step (if tcheck=0).
!!  \author Thijs Heus,MPI-M
!!  \author Hans Cuijpers, KNMI
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
!
module modchecksim

  use, intrinsic :: iso_fortran_env, only: real64, real32
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  use modprecision,   only: field_r
  use modglobal,      only: longint, i1, j1,ih, jh, kmax, dtmax, dx, dy, dzf, dzh, &
                            dt_reason, ifnamopt, checknamelisterror, tres, btime, &
                            ladaptive, timee, rtimee, rk3step, rdt, fname_options
  use modfields,      only: u0, v0, w0, qt0, thl0, e120, qtp, thlp, rhobf, rhobh
  use modsubgriddata, only: ekm
  use modstringutils, only: number2string
  use modmpi,         only: myid, comm3d, mpierr, mpi_sum, mpi_max, D_MPI_ALLREDUCE, &
                            D_MPI_BCAST
  use modtimer

  implicit none

  private

  character(len=*), parameter :: modname = 'modchecksim'

  public :: checksim_read_namelist
  public :: initchecksim
  public :: exitchecksim
  public :: checksim
  public :: chkdiv

  public :: checktend
  public :: check_array
  public :: lchecktend

  interface check_array !< Check array for invalid values and/or values outside of a given range.
    module procedure :: check_array_1d_int
    module procedure :: check_array_1d_r4
    module procedure :: check_array_1d_r8
    module procedure :: check_array_2d_int
    module procedure :: check_array_2d_r4
    module procedure :: check_array_2d_r8
    module procedure :: check_array_3d_int
    module procedure :: check_array_3d_r4
    module procedure :: check_array_3d_r8
  end interface

  real(field_r) :: &
    tcheck = 0,    &
    dtmn = 0,      &
    ndt = 0

  integer(longint) :: &
    tnext = 3600,     &
    itcheck

  ! explanations for dt_limit, determined in tstep_update()
  character (len=15) :: dt_reasons(0:5) = [character(len=15) :: &
    "initial step", "timee", "dt_lim" , "idtmax", "velocity", "diffusion"]

  logical :: lchecktend
  logical :: lstop

  real(field_r), allocatable :: &
    courx(:),                   &
    coury(:),                   &
    courz(:),                   &
    courtot(:),                 &
    peclettot(:)

contains
  
  !> Read checksim namelist.
  subroutine checksim_read_namelist(nml_filename)

    character(len=*), intent(in) :: nml_filename

    integer :: ierr

    namelist /NAMCHECKSIM/ tcheck, lchecktend, lstop

    if (myid == 0) then
      open(ifnamopt, file=nml_filename, status='old', iostat=ierr)
      read(ifnamopt, NAMCHECKSIM, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMCHECKSIM')
      write(6, NAMCHECKSIM) ! Maybe write to separate file, for cleaner terminal
      close(ifnamopt)
    end if

    call D_MPI_BCAST(tcheck, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(lchecktend, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(lstop, 1, 0, comm3d, mpierr)

  end subroutine checksim_read_namelist

  !> Initialize checksim variables.
  subroutine initchecksim

    character(len=*), parameter :: routine = modname//'/initchecksim'

    integer :: ierr

    call timer_tic(routine, 0)

    if (.not. ladaptive .and. tcheck < dtmax) then
      tcheck = dtmax
    end if

    itcheck = floor(tcheck / tres)
    tnext = itcheck + btime

    allocate(courx(kmax), coury(kmax), courz(kmax), courtot(kmax), peclettot(kmax))

    !$acc enter data create(courx, coury, courz, courtot, peclettot)

    call timer_toc(routine)

  end subroutine initchecksim

  !> Deallocate checksim arrays.
  subroutine exitchecksim

    !$acc exit data delete(courx, coury, courz, courtot, peclettot)

    deallocate(courx, coury, courz, courtot, peclettot)

  end subroutine exitchecksim

  !> Run checksim. Timekeeping, and output
  subroutine checksim
    
    character(len=*), parameter :: routine = modname//'/checksim'

    character(len=20) :: timeday

    if (timee == 0) return
    if (rk3step /= 3) return

    dtmn = dtmn + rdt
    ndt = ndt + 1

    if (timee < tnext) return

    call timer_tic('modchecksim/checksim', 0)

    tnext = tnext+itcheck
    dtmn  = dtmn / ndt

    if (myid == 0) then
      call date_and_time(time=timeday)
      write (*,*) '================================================================='
      write (*,'(7A,F11.2,A,F9.4)') 'Time of Day: ', timeday(1:2), ':', &
        timeday(3:4), ':', timeday(5:10),' Time of Simulation: ', &
        rtimee, '    dt: ',dtmn
    end if

    call calccourantandpeclet
    call chkdiv

    dtmn  = 0.
    ndt   = 0.

    call timer_toc('modchecksim/checksim')

  end subroutine checksim

  !> Calculates the courant number as in max(w)*deltat/deltaz
  !! and peclet number as max(ekm) *deltat/deltax**2
  subroutine calccourantandpeclet

    integer       :: i, j, k
    real(field_r) :: &
      velx_max,      &
      vely_max,      &
      velz_max,      &
      velmag_max,    &
      ekm_max

    !$acc parallel loop gang default(present) &
    !$acc private(velx_max, vely_max, velz_max, velmag_max, ekm_max)
    do k = 1, kmax
      velx_max = 0
      vely_max = 0
      velz_max = 0
      velmag_max = 0
      ekm_max = 0
      !$acc loop collapse(2) &
      !$acc reduction(max:velx_max, vely_max, velz_max, velmag_max, ekm_max)
      do j = 2, j1
        do i = 2, i1
          velx_max = max(velx_max, abs(u0(i,j,k)))
          vely_max = max(vely_max, abs(v0(i,j,k)))
          velz_max = max(velz_max, abs(w0(i,j,k)))
          velmag_max = max(velmag_max, u0(i,j,k)*u0(i,j,k)/(dx*dx) + &
                                       v0(i,j,k)*v0(i,j,k)/(dy*dy) + &
                                       w0(i,j,k)*w0(i,j,k)/(dzh(k)*dzh(k)))
          ekm_max = max(ekm_max, ekm(i,j,k))
        enddo
      enddo
      courx(k)=velx_max*dtmn/dx
      coury(k)=vely_max*dtmn/dy
      courz(k)=velz_max*dtmn/dzh(k)
      courtot(k)=velmag_max*dtmn*dtmn
      peclettot(k)=ekm_max*dtmn/min(dzh(k),dx,dy)**2
    end do

    !$acc update self(courx, coury, courz, courtot, peclettot)

    call D_MPI_ALLREDUCE(courx, kmax, MPI_MAX, comm3d, mpierr)
    call D_MPI_ALLREDUCE(coury, kmax, MPI_MAX, comm3d, mpierr)
    call D_MPI_ALLREDUCE(courz, kmax, MPI_MAX, comm3d, mpierr)
    call D_MPI_ALLREDUCE(courtot, kmax, MPI_MAX, comm3d, mpierr)
    call D_MPI_ALLREDUCE(peclettot, kmax, MPI_MAX, comm3d, mpierr)

    if (myid == 0) then
      write(*,'(A,3ES10.2,I5,ES10.2,I5)') 'Courant numbers (x,y,z,tot):', &
        maxval(courx(:)), maxval(coury(:)), maxval(courz(:)), maxloc(courz(:)), &
        sqrt(maxval(courtot(:))), maxloc(courtot(:))
      write(6,'(A,ES10.2,I5)') 'Cell Peclet number:', &
        maxval(peclettot(:)), maxloc(peclettot(:))
    end if

  end subroutine calccourantandpeclet

  !> Checks local and total divergence.
  subroutine chkdiv

    integer       :: i, j, k
    real(field_r) :: &
      div,           &
      divmax,        &
      divtot

    divmax = 0.
    divtot = 0.

    !$acc parallel loop collapse(3) default(present) private(div) &
    !$acc reduction(max:divmax) reduction(+:divtot)
    do k=1,kmax
      do j=2,j1
        do i=2,i1
          div = rhobf(k) * (u0(i+1,j,k) - u0(i,j,k) )/dx + &
                rhobf(k) * (v0(i,j+1,k) - v0(i,j,k) )/dy + &
                (rhobh(k+1)*w0(i,j,k+1) - rhobh(k)*w0(i,j,k) )/dzf(k)
          divmax = max(divmax,abs(div))
          divtot = divtot + div*dx*dy*dzf(k)
        end do
      end do
    end do

    call D_MPI_ALLREDUCE(divtot, 1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(divmax, 1, MPI_MAX, comm3d,mpierr)

    if (myid == 0) then
      write(6 ,'(A,2ES11.2,A,A)')'divmax, divtot = ', divmax, divtot,  &
        '       dt limited by ', dt_reasons(dt_reason)
   end if

  end subroutine chkdiv

  !> Check tendencies of various prognostic variables.
  subroutine checktend(step)

    character(len=*), intent(in) :: step

    call check_array(qtp, "qtp", step, [-0.01_field_r, 0.01_field_r])
    call check_array(thlp, "thlp", step, [-20.0_field_r, 20.0_field_r])
  
  end subroutine checktend

  subroutine check_array_1d_int(array, name, step, threshold, lacc)

    integer,          intent(in) :: array(:), threshold(2)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: step

    logical, intent(in), optional :: lacc

    integer           :: i
    integer           :: val
    character(len=32) :: cloc
    character(len=11) :: cval

    do i = 1, size(array, dim=1)
      val = array(i)
      cval = number2string(val)
      if ((val < threshold(1) .or. val > threshold(2))) then
        call print_warning_out_of_range(name, step, [i], cval, &
                [number2string(threshold(1)), number2string(threshold(2))])
      end if
    end do

  end subroutine check_array_1d_int

  subroutine check_array_1d_r4(array, name, step, threshold, lacc)

    real(real32),     intent(in) :: array(:)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: step

    real(real32), intent(in), optional :: threshold(2)
    logical,      intent(in), optional :: lacc

    integer           :: i
    real(real32)      :: val
    character(len=32) :: cloc
    character(len=11) :: cval

    do i = 1, size(array, dim=1)
      val = array(i)
      cval = number2string(val)
      if (.not. ieee_is_finite(val)) then
        call print_warning_non_finite(name, step, [i], cval)
      else if (present(threshold)) then
        if (val < threshold(1) .or. val > threshold(2)) then
          call print_warning_out_of_range(name, step, [i], cval, &
                  [number2string(threshold(1)), number2string(threshold(2))])
        end if
      end if
    end do

  end subroutine check_array_1d_r4

  subroutine check_array_1d_r8(array, name, step, threshold, lacc)

    real(real64),     intent(in) :: array(:)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: step

    real(real64), intent(in), optional :: threshold(2)
    logical,      intent(in), optional :: lacc

    integer           :: i
    real(real64)      :: val
    character(len=32) :: cloc
    character(len=11) :: cval

    do i = 1, size(array, dim=1)
      val = array(i)
      cval = number2string(val)
      if (.not. ieee_is_finite(val)) then
        call print_warning_non_finite(name, step, [i], cval)
      else if (present(threshold)) then
        if (val < threshold(1) .or. val > threshold(2)) then
          call print_warning_out_of_range(name, step, [i], cval, &
                  [number2string(threshold(1)), number2string(threshold(2))])
        end if
      end if
    end do

  end subroutine check_array_1d_r8

  subroutine check_array_2d_int(array, name, step, threshold, lacc)

    integer,          intent(in) :: array(:,:), threshold(2)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: step

    logical, intent(in), optional :: lacc

    integer           :: i, j
    integer           :: val
    character(len=32) :: cloc
    character(len=11) :: cval

    do j = 1, size(array, dim=2)
      do i = 1, size(array, dim=1)
        val = array(i,j)
        cval = number2string(val)
        if ((val < threshold(1) .or. val > threshold(2))) then
          call print_warning_out_of_range(name, step, [i, j], cval, &
                  [number2string(threshold(1)), number2string(threshold(2))])
        end if
      end do
    end do

  end subroutine check_array_2d_int

  subroutine check_array_2d_r4(array, name, step, threshold, lacc)

    real(real32),     intent(in) :: array(:,:)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: step

    real(real32), intent(in), optional :: threshold(2)
    logical,      intent(in), optional :: lacc

    integer           :: i, j
    real(real32)      :: val
    character(len=32) :: cloc
    character(len=11) :: cval

    do j = 1, size(array, dim=2)
      do i = 1, size(array, dim=1)
        val = array(i,j)
        cval = number2string(val)
        if (.not. ieee_is_finite(val)) then
          call print_warning_non_finite(name, step, [i, j], cval)
        else if (present(threshold)) then
          if (val < threshold(1) .or. val > threshold(2)) then
            call print_warning_out_of_range(name, step, [i, j], cval, &
                    [number2string(threshold(1)), number2string(threshold(2))])
          end if
        end if
      end do
    end do

  end subroutine check_array_2d_r4

  subroutine check_array_2d_r8(array, name, step, threshold, lacc)

    real(real64),     intent(in) :: array(:,:)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: step

    real(real64), intent(in), optional :: threshold(2)
    logical,      intent(in), optional :: lacc

    integer           :: i, j
    real(real64)      :: val
    character(len=32) :: cloc
    character(len=11) :: cval

    do j = 1, size(array, dim=2)
      do i = 1, size(array, dim=1)
        val = array(i,j)
        cval = number2string(val)
        if (.not. ieee_is_finite(val)) then
          call print_warning_non_finite(name, step, [i, j], cval)
        else if (present(threshold)) then
          if (val < threshold(1) .or. val > threshold(2)) then
            call print_warning_out_of_range(name, step, [i, j], cval, &
                    [number2string(threshold(1)), number2string(threshold(2))])
          end if
        end if
      end do
    end do

  end subroutine check_array_2d_r8

  subroutine check_array_3d_int(array, name, step, threshold, lacc)

    integer,          intent(in) :: array(:,:,:), threshold(2)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: step

    logical, intent(in), optional :: lacc

    integer           :: i, j, k
    integer           :: val
    character(len=32) :: cloc
    character(len=11) :: cval

    do k = 1, size(array, dim=3)
      do j = 1, size(array, dim=2)
        do i = 1, size(array, dim=1)
          val = array(i,j,k)
          cval = number2string(val)
          if ((val < threshold(1) .or. val > threshold(2))) then
            call print_warning_out_of_range(name, step, [i, j, k], cval, &
                    [number2string(threshold(1)), number2string(threshold(2))])
          else
            cycle
          end if
          if (lstop) then
            call dump_state([i,j,k])
            error stop
          end if
        end do
      end do
    end do

  end subroutine check_array_3d_int

  subroutine check_array_3d_r4(array, name, step, threshold, lacc)

    real(real32),     intent(in) :: array(:,:,:)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: step

    real(real32), intent(in), optional :: threshold(2)
    logical,      intent(in), optional :: lacc

    integer           :: i, j, k
    real(real32)      :: val
    character(len=32) :: cloc
    character(len=11) :: cval

    do k = 1, size(array, dim=3)
      do j = 1, size(array, dim=2)
        do i = 1, size(array, dim=1)
          val = array(i,j,k)
          cval = number2string(val)
          if (.not. ieee_is_finite(val)) then
            call print_warning_non_finite(name, step, [i, j, k], cval)
          else if (present(threshold)) then
            if (val < threshold(1) .or. val > threshold(2)) then
              call print_warning_out_of_range(name, step, [i, j, k], cval, &
                      [number2string(threshold(1)), number2string(threshold(2))])
            end if
          else
            cycle
          end if
          if (lstop) then
            call dump_state([i-ih+1,j-jh+1,k])
            error stop
          end if
        end do
      end do
    end do

  end subroutine check_array_3d_r4

  subroutine check_array_3d_r8(array, name, step, threshold, lacc)

    real(real64),     intent(in) :: array(:,:,:)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: step

    real(real64), intent(in), optional :: threshold(2)
    logical,      intent(in), optional :: lacc

    integer      :: i, j, k
    real(real64) :: val
    character(len=32) :: cloc
    character(len=11) :: cval

    do k = 1, size(array, dim=3)
      do j = 1, size(array, dim=2)
        do i = 1, size(array, dim=1)
          val = array(i,j,k)
          cval = number2string(val)
          if (.not. ieee_is_finite(val)) then
            call print_warning_non_finite(name, step, [i, j, k], cval)
          else if (present(threshold)) then
            if (val < threshold(1) .or. val > threshold(2)) then
              call print_warning_out_of_range(name, step, [i, j, k], cval, &
                      [number2string(threshold(1)), number2string(threshold(2))])
            end if
          else
            cycle
          end if
          if (lstop) then
            call dump_state([i-ih+1,j-jh+1,k])
            error stop
          end if
        end do
      end do
    end do

  end subroutine check_array_3d_r8

  !> Prints non-finite warning (Inf, NaN) to stderr.
  subroutine print_warning_non_finite(name, step, loc, cval)

    character(len=*), intent(in) :: name, step, cval
    integer,          intent(in) :: loc(:)

    character(len=128) :: cloc

    write(cloc, "(*(g0,:,','))") loc
    write(0,'(*(a))') &
      'modchecktend: Invalid value found in array: ', trim(name), '  [step=', &
      trim(step), ']  [location=', trim(cloc), ']  [value=', trim(cval), ']'

  end subroutine print_warning_non_finite

  !> Prints out of range warning to stderr.
  subroutine print_warning_out_of_range(name, step, loc, cval, cthreshold)

    character(len=*), intent(in) :: name, step, cval, cthreshold(2)
    integer,          intent(in) :: loc(:)

    character(len=128) :: cloc

    write(cloc, "(*(g0,:,','))") loc
    write(0,'(13a)') &
      'modchecktend: Value outside of valid range found in array: ', trim(name), &
      '  [step=', trim(step), ']  [location=', trim(cloc), ']  [value=', trim(cval), ']&
      &  [min,max=', trim(cthreshold(1)), ',', trim(cthreshold(2)), ']'

  end subroutine print_warning_out_of_range

  !> Dumps values of prognostic variables at specified location.
  !!
  !! @param[in] loc (i,j,k) location to print variables for. 
  subroutine dump_state(loc)

    integer, intent(in) :: loc(3)

    integer :: i, j, k

    i = loc(1)
    j = loc(2)
    k = loc(3)

    write(0, '(7(a,/))')"Prognostic variables:", &
      & "u   = "//number2string(u0(i,j,k)), &
      & "v   = "//number2string(v0(i,j,k)), &
      & "w   = "//number2string(w0(i,j,k)), &
      & "qt  = "//number2string(qt0(i,j,k)), &
      & "thl = "//number2string(thl0(i,j,k)), &
      & "e12 = "//number2string(e120(i,j,k))

  end subroutine dump_state

end module modchecksim
