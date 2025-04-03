!> \file modnudge.f90
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
! Copyright 1993-2024 The DALES team.
!
!> Module for nudging prognostic fields to some provided profiles.
module modnudge
  use modprecision, only: field_r
  use modtimer,     only: timer_tic, timer_toc

  implicit none

  private

  save

  character(*), parameter :: modname = "modnudge"

  ! Switches for enabling/disabling nudging
  logical :: lnudge = .false.
  logical :: lunudge = .true.
  logical :: lvnudge = .true.
  logical :: lwnudge = .true.
  logical :: lthlnudge = .true.
  logical :: lqtnudge = .true.

  ! Nudging profiles
  real(field_r), allocatable :: tnudge(:,:)
  real(field_r), allocatable :: unudge(:,:)
  real(field_r), allocatable :: vnudge(:,:)
  real(field_r), allocatable :: wnudge(:,:)
  real(field_r), allocatable :: thlnudge(:,:)
  real(field_r), allocatable :: qtnudge(:,:)
  ! Nudging constants (only supported with DEPHY input)
  real(field_r), allocatable :: tunudge(:,:)
  real(field_r), allocatable :: tvnudge(:,:)
  real(field_r), allocatable :: twnudge(:,:)
  real(field_r), allocatable :: tthlnudge(:,:)
  real(field_r), allocatable :: tqtnudge(:,:)

  real(field_r), allocatable :: timenudge(:)

  ! Nudging timescale
  real(field_r) :: tnudgefac = 1.
  ! Number of nudging time steps
  integer       :: ntnudge = 10000

  public :: initnudge
  public :: nudge
  public :: exitnudge

contains

  !> Read nudging profiles and nudging constants from input, and allocate
  !! memory.
  subroutine initnudge
    use modmpi,     only: myid, mpierr, comm3d, D_MPI_BCAST
    use modglobal,  only: ifnamopt, fname_options, runtime, cexpnr, ifinput, &
                          k1, kmax, checknamelisterror, iinput, input_netcdf
    use modstat_nc

    character(*), parameter :: routine = modname//"::initnudge"

    integer      :: ierr, k, t
    integer      :: ncid, varid, dimid
    character(1) :: chmess1
    real, allocatable, dimension(:) :: height

    namelist /NAMNUDGE/ lnudge, lunudge, lvnudge, lwnudge, lthlnudge, &
                        lqtnudge, tnudgefac

    if (myid == 0) then
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read(ifnamopt, NAMNUDGE, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMNUDGE')
      write(6, NAMNUDGE)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(lnudge, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(lunudge, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(lvnudge, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(lwnudge, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(lthlnudge, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(lqtnudge, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(tnudgefac, 1, 0, comm3d, mpierr)

    if (.not. lnudge) return

    call timer_tic(routine, 0)

    if (iinput == input_netcdf) then
      if (myid == 0) then
        call nchandle_error(nf90_open("init."//cexpnr//".nc", NF90_NOWRITE, &
                            ncid))

        ! Get the number of forcing time steps
        call nchandle_error(nf90_inq_dimid(ncid, "time", dimid))
        call nchandle_error(nf90_inquire_dimension(ncid, dimid, len=ntnudge))
      end if

      call D_MPI_BCAST(ntnudge, 1, 0, comm3d, mpierr)
      allocate(timenudge(0:ntnudge))
      timenudge(:) = 0
      if (myid == 0) then
        call nchandle_error(nf90_inq_varid(ncid, "time", varid))
        call nchandle_error(nf90_get_var(ncid, varid, timenudge(1:)))
      end if
      call D_MPI_BCAST(timenudge, ntnudge + 1, 0, comm3d, mpierr)
      
      if (lunudge) then
        allocate(unudge(k1,ntnudge), tunudge(k1,ntnudge))
        if (myid == 0) then
          call nchandle_error(nf90_inq_varid(ncid, "ua_nud", varid))
          call nchandle_error(nf90_get_var(ncid, varid, unudge(1:kmax,:)))
          call nchandle_error(nf90_inq_varid(ncid, "nudging_constant_ua", &
                              varid))
          call nchandle_error(nf90_get_var(ncid, varid, tunudge(1:kmax,:)))
        end if
      end if

      if (lvnudge) then
        allocate(vnudge(k1,ntnudge), tvnudge(k1,ntnudge))
        if (myid == 0) then
          call nchandle_error(nf90_inq_varid(ncid, "va_nud", varid))
          call nchandle_error(nf90_get_var(ncid, varid, vnudge(1:kmax,:)))
          call nchandle_error(nf90_inq_varid(ncid, "nudging_constant_va", &
                              varid))
          call nchandle_error(nf90_get_var(ncid, varid, tvnudge(1:kmax,:)))
        end if
      end if

      if (lwnudge) then
        allocate(wnudge(k1,ntnudge), twnudge(k1,ntnudge))
        if (myid == 0) then
          call nchandle_error(nf90_inq_varid(ncid, "wa_nud", varid))
          call nchandle_error(nf90_get_var(ncid, varid, wnudge))
          call nchandle_error(nf90_inq_varid(ncid, "nudging_constant_wa", &
                              varid))
          call nchandle_error(nf90_get_var(ncid, varid, twnudge))
        end if
      end if

      if (lthlnudge) then
        allocate(thlnudge(k1,ntnudge), tthlnudge(k1,ntnudge))
        if (myid == 0) then
          call nchandle_error(nf90_inq_varid(ncid, "thetal_nud", varid))
          call nchandle_error(nf90_get_var(ncid, varid, thlnudge(1:kmax,:)))
          call nchandle_error(nf90_inq_varid(ncid, "nudging_constant_thetal", &
                              varid))
          call nchandle_error(nf90_get_var(ncid, varid, tthlnudge(1:kmax,:)))
        end if
      end if

      if (lqtnudge) then
        allocate(qtnudge(k1,ntnudge), tqtnudge(k1,ntnudge))
        if (myid == 0) then
          call nchandle_error(nf90_inq_varid(ncid, "qt_nud", varid))
          call nchandle_error(nf90_get_var(ncid, varid, qtnudge(1:kmax,:)))
          call nchandle_error(nf90_inq_varid(ncid, "nudging_constant_qt", &
                              varid))
          call nchandle_error(nf90_get_var(ncid, varid, tqtnudge(1:kmax,:)))
        end if
      end if

      if (myid == 0) then
        call nchandle_error(nf90_close(ncid))
      end if
    else
      allocate(tnudge(k1,ntnudge), unudge(k1,ntnudge), vnudge(k1,ntnudge), &
               wnudge(k1,ntnudge), thlnudge(k1,ntnudge), qtnudge(k1,ntnudge))
      allocate(tunudge(k1,ntnudge), tvnudge(k1,ntnudge), &
               twnudge(k1,ntnudge), tthlnudge(k1,ntnudge), &
               tqtnudge(k1,ntnudge))
      allocate(timenudge(0:ntnudge), height(k1))

      tnudge = 0
      unudge = 0
      vnudge = 0
      wnudge = 0
      thlnudge = 0
      qtnudge = 0
      timenudge = 0

      t = 0
      if (myid == 0) then
        open(ifinput, file='nudge.inp.'//cexpnr)

        do while (timenudge(t) < runtime)
          t = t + 1
          if (t > ntnudge) then
            write(*, *) "Too many time points in file ", 'nudge.inp.'//cexpnr, &
            &", the limit is ntnudge = ", ntnudge
            stop
          end if

          chmess1 = "#"
          ierr = 1 ! not zero
          !search for the next line consisting of "# time", from there onwards the profiles will be read
          do while (.not. (chmess1 == "#" .and. ierr == 0))
            read(ifinput, *, iostat=ierr) chmess1, timenudge(t)
            if (ierr < 0) then
              stop 'STOP: No time dependend nudging data for end of run'
            end if

          end do
          write(6, *) 'time', timenudge(t)
          write(6, *) ' height    t_nudge    u_nudge    v_nudge    w_nudge    &
          &thl_nudge    qt_nudge'
          do k = 1, kmax
            read(ifinput, *) &
              height(k), &
              tnudge(k,t), &
              unudge(k,t), &
              vnudge(k,t), &
              wnudge(k,t), &
              thlnudge(k,t), &
              qtnudge(k,t)
          end do

          do k = kmax, 1, -1
            write(6, '(f7.1,6e12.4)') &
              height(k), &
              tnudge(k,t), &
              unudge(k,t), &
              vnudge(k,t), &
              wnudge(k,t), &
              thlnudge(k,t), &
              qtnudge(k,t)
          end do
        end do
        close (ifinput)
      end if

      lunudge = any(abs(unudge) > 1e-8)
      lvnudge = any(abs(vnudge) > 1e-8)
      lwnudge = any(abs(wnudge) > 1e-8)
      lthlnudge = any(abs(thlnudge) > 1e-8)
      lqtnudge = any(abs(qtnudge) > 1e-8)

      call D_MPI_BCAST(lunudge, 1, 0, comm3d, mpierr)
      call D_MPI_BCAST(lvnudge, 1, 0, comm3d, mpierr)
      call D_MPI_BCAST(lwnudge, 1, 0, comm3d, mpierr)
      call D_MPI_BCAST(lthlnudge, 1, 0, comm3d, mpierr)
      call D_MPI_BCAST(lqtnudge, 1, 0, comm3d, mpierr)

      tnudge = tnudgefac * tnudge

      tunudge(:,:) = tnudge(:,:)
      tvnudge(:,:) = tnudge(:,:)
      twnudge(:,:) = tnudge(:,:)
      tthlnudge(:,:) = tnudge(:,:)
      tqtnudge(:,:) = tnudge(:,:)
    end if


    call D_MPI_BCAST(timenudge, ntnudge + 1, 0, comm3d, mpierr)
    if (lunudge) then
      call D_MPI_BCAST(unudge, k1 * ntnudge, 0, comm3d, mpierr)
      call D_MPI_BCAST(tunudge, k1 * ntnudge, 0, comm3d, mpierr)
    end if
    if (lvnudge) then
      call D_MPI_BCAST(vnudge, k1 * ntnudge, 0, comm3d, mpierr)
      call D_MPI_BCAST(tvnudge, k1 * ntnudge, 0, comm3d, mpierr)
    end if
    if (lwnudge) then
      call D_MPI_BCAST(wnudge, k1 * ntnudge, 0, comm3d, mpierr)
      call D_MPI_BCAST(twnudge, k1 * ntnudge, 0, comm3d, mpierr)
    end if
    if (lthlnudge) then
      call D_MPI_BCAST(thlnudge, k1 * ntnudge, 0, comm3d, mpierr)
      call D_MPI_BCAST(tthlnudge, k1 * ntnudge, 0, comm3d, mpierr)
    end if
    if (lqtnudge) then
      call D_MPI_BCAST(qtnudge, k1 * ntnudge, 0, comm3d, mpierr)
      call D_MPI_BCAST(tqtnudge, k1 * ntnudge, 0, comm3d, mpierr)
    end if

    !$acc enter data copyin(timenudge, unudge, vnudge, wnudge, thlnudge, &
    !$acc&                  qtnudge, tunudge, tvnudge, twnudge, tthlnudge, &
    !$acc&                  tqtnudge)

    call timer_toc(routine)
  end subroutine initnudge

  !> Perform nudging of velocities, temperature and humidity fields.
  subroutine nudge
    use modglobal, only: timee, rtimee, i1, j1, kmax, rdt
    use modfields, only: up, vp, wp, thlp, qtp, u0av, v0av, qt0av, thl0av

    character(*), parameter :: routine = modname//"::nudge"

    integer       :: i, j, k, t
    real(field_r) :: dtm, dtp, currtnudge

    if (.not. (lnudge)) return

    if (timee == 0) return

    call timer_tic(routine, 0)

    t = 1
    do while (rtimee > timenudge(t))
      t = t + 1
    end do
    if (rtimee > timenudge(1)) then
      t = t - 1
    end if

    dtm = (rtimee - timenudge(t)) / (timenudge(t + 1) - timenudge(t))
    dtp = (timenudge(t + 1) - rtimee) / (timenudge(t + 1) - timenudge(t))

    if (lunudge) then
      !$acc parallel loop collapse(3) private(currtnudge) default(present) async
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1
            currtnudge = max(1.0_field_r * rdt, &
                             tunudge(k,t) * dtp + tunudge(k,t + 1) * dtm)
            up(i,j,k) = up(i,j,k) - (u0av(k) - (unudge(k,t) * dtp + &
                        unudge(k,t + 1) * dtm)) / currtnudge
          end do
        end do
      end do
    end if

    if (lvnudge) then
      !$acc parallel loop collapse(3) default(present) private(currtnudge) async
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1
            currtnudge = max(1.0_field_r * rdt, &
                             tvnudge(k,t) * dtp + tvnudge(k,t + 1) * dtm)
            vp(i,j,k) = vp(i,j,k) - (v0av(k) - (vnudge(k,t) * dtp + &
                        vnudge(k,t + 1) * dtm)) / currtnudge
          end do
        end do
      end do
    end if

    if (lwnudge) then
      !$acc parallel loop collapse(3) default(present) private(currtnudge) async
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1
            currtnudge = max(1.0_field_r * rdt, &
                             twnudge(k,t) * dtp + twnudge(k,t + 1) * dtm)
            wp(i,j,k) = wp(i,j,k) - ((wnudge(k,t) * dtp + wnudge(k,t + 1) &
                        * dtm)) / currtnudge
          end do
        end do
      end do
    end if

    if (lthlnudge) then
      !$acc parallel loop collapse(3) default(present) private(currtnudge) async
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1
            currtnudge = max(1.0_field_r * rdt, &
                             tthlnudge(k,t) * dtp + tthlnudge(k,t + 1) * dtm)
            thlp(i,j,k) = thlp(i,j,k) - (thl0av(k) - (thlnudge(k,t) * dtp + &
                          thlnudge(k,t + 1) * dtm)) / currtnudge
          end do
        end do
      end do
    end if

    if (lqtnudge) then
      !$acc parallel loop collapse(3) default(present) private(currtnudge) async
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1
            currtnudge = max(1.0_field_r * rdt, &
                             tqtnudge(k,t) * dtp + tqtnudge(k,t + 1) * dtm)
            qtp(i,j,k) = qtp(i,j,k) - (qt0av(k) - (qtnudge(k,t) * dtp + &
                        qtnudge(k,t + 1) * dtm)) / currtnudge
          end do
        end do
      end do
    end if

    !$acc wait

    call timer_toc(routine)
  end subroutine nudge

  !> Deallocates memory.
  subroutine exitnudge
    if (.not. lnudge) return

    deallocate(timenudge)

    if (allocated(tnudge)) deallocate(tnudge)
    
    if (lunudge) then
      deallocate(unudge, tunudge)
    end if

    if (lvnudge) then
      deallocate(vnudge, tvnudge)
    end if

    if (lwnudge) then
      deallocate(wnudge, twnudge)
    end if

    if (lthlnudge) then
      deallocate(thlnudge, tthlnudge)
    end if

    if (lqtnudge) then
      deallocate(qtnudge, tqtnudge)
    end if
  end subroutine exitnudge

end module
