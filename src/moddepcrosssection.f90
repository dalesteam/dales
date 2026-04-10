!> \file moddepcrosssection.f90
!!  Dumps fields of deposition data of scalars
!!
!!  Dumps fields of deposition data of scalars to depo_*.myidx.myidy.expnr
!! If netcdf is true, this module leads the depcross.myidx.myidy.expnr.nc output
!!
!!  \author Leon Geers, TNO
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
!  Copyright 1993-2023 Delft University of Technology, Wageningen University, Utrecht University, KNMI, TNO
!
module moddepcrosssection
  use modlsm, only : llsm
  use moddrydeposition, only : ldrydep
  use modglobal, only : longint, nsv, itot, jtot
  use modnetcdf_file_t, only : cross_section_file_t
  use modprecision, only : field_r
  use modstat_nc_files, only : add_output_file
  use modtracers, only: tracer_prop
  use modlogging, only: finish

  implicit none
  character(len=*), parameter :: modname = 'moddepcrossection'
  private
  public :: initdepcrosssection, depcrosssection, exitdepcrosssection
  save

  type(cross_section_file_t) :: dep_file
  integer :: dep_file_id = 0
  logical :: dep_file_enabled = .false.

  ! integer :: nvar = 0  !< Number of variables (for now, equal to nsv, not using svskip)
  real :: dtav
  integer(kind=longint) :: idtav, tnext
  logical :: ldepcrosssection = .false.  !< Switch for doing depcrosssection (on/off)


contains
  !> Initializing depcrosssection. Read out the namelist, initializing the variables
  subroutine initdepcrosssection
    use modmpi, only : myid, comm3d, &
        mpierr, D_MPI_BCAST
    use modglobal, only : dtav_glob, ifnamopt, fname_options, &
        checknamelisterror, tres, btime, dt_lim, ladaptive, &
      dtmax
    use modstat_nc, only : lnetcdf
    use moddrydeposition, only : ndeptracers
    use fortran_support,  only : nnml_output

    implicit none

    character(len=*), parameter :: routine = modname//'/initdepcrosssection'

    integer :: ierr, isv
    character(80) :: varname, varlongname

    namelist/NAMDEPCROSSSECTION/ ldepcrosssection, dtav

    dtav = dtav_glob
    if (myid==0) then
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read(ifnamopt, NAMDEPCROSSSECTION, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMDEPCROSSSECTION')
      write(nnml_output, NAMDEPCROSSSECTION)
      close(ifnamopt)
    end if

    if (ldepcrosssection .and. (.not. llsm .or. .not. ldrydep .or. ndeptracers == 0)) then
      ldepcrosssection = .false.
      write (6, *) "Ignoring ldepcrosssection, since no dry deposition &
        & and/or land surface model defined"
   end if

    call D_MPI_BCAST(dtav,             1, 0, comm3d, mpierr)
    call D_MPI_BCAST(ldepcrosssection, 1, 0, comm3d, mpierr)

    idtav = int(dtav / tres, kind=kind(idtav))
    tnext   = idtav+btime
    if(.not. ldepcrosssection) return
    dt_lim = min(dt_lim, tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      call finish(routine, 'depcrosssection: dtav should be a integer multiple of dtmax')
    end if

    if (lnetcdf) then
      dep_file = cross_section_file_t('depcross', nx=itot, ny=jtot, lgpu=.false.)
      do isv = 1, nsv
        if (.not. tracer_prop(isv)%ldep) cycle
        write (varname, '(a,a)') 'drydep_', trim(tracer_prop(isv)%tracname)
        write (varlongname, '(a,a)')  'Dry deposition flux of ', trim(tracer_prop(isv)%tracname)
        call dep_file%add_var(varname, varlongname, 'kg / (m2 * s)', 'tt0t')
      end do
      call add_output_file(dep_file, dtav, dep_file_id)
      dep_file_enabled = .true.
    end if
  end subroutine initdepcrosssection

  !> Do crosssection. Collect data to truncated (2 byte) integers, and write them to file
  subroutine depcrosssection
    use modglobal, only : rk3step, timee, dt_lim
    implicit none

    if (.not. ldepcrosssection) return
    if (rk3step/=3) return
    if (timee<tnext) then
      dt_lim = min(dt_lim, tnext-timee)
      return
    end if
    tnext = tnext + idtav
    dt_lim = minval((/dt_lim, tnext-timee/))

    call wrtdrydepfields

  end subroutine depcrosssection

  subroutine wrtdrydepfields
    use moddrydeposition, only : depfield
    use modglobal, only : i1, j1, nsv
    use modfields, only : rhof
    use modstat_nc, only : lnetcdf
    use modtracers, only : tracer_prop
    implicit none

    real(field_r), pointer :: dep_ptr(:, :)
    integer :: isv, idt
    real    :: MW_air = 28.9644
    character(80) :: varname

    if (.not. (lnetcdf .and. dep_file_enabled)) return

    ! Store the flux as a positive number
    ! dep_ptr(1:imax, 1:jmax, 1:ndeptracers) = -depfield(2:i1, 2:j1, 1:ndeptracers) &
    !     * rhof(1) * 1e-6  ! to go from ug*m/(s*g) to kg/(m2*s)
    idt = 1
    do isv = 1, nsv
      if (.not. tracer_prop(isv)%ldep) cycle
      write(varname, '(a,a)') 'drydep_', trim(tracer_prop(isv)%tracname)
      call dep_file%get_pointer(trim(varname), dep_ptr)
      dep_ptr(:,:) = -depfield(2:i1, 2:j1, idt) * rhof(1) * &
          (tracer_prop(isv)%molar_mass / MW_air) * 1e-9  ! from ppb m s-1 to kg/(m2*s)
      idt = idt + 1
    end do
  end subroutine wrtdrydepfields

  !> Clean up when leaving the run
  subroutine exitdepcrosssection
    implicit none
  end subroutine exitdepcrosssection

end module moddepcrosssection
