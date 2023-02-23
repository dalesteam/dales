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
  use modglobal, only : longint, nsv
  use modtracers, only: tracer_prop

  implicit none
  private
  public :: initdepcrosssection, depcrosssection, exitdepcrosssection
  save
  ! NetCDF variables
  integer :: ncid, nrec = 0
  character(80) :: fname
  character(80), allocatable, dimension(:, :) :: ncname
  character(80), dimension(1, 4) :: tncname

  ! integer :: nvar = 0  !< Number of variables (for now, equal to nsv, not using svskip) 
  real :: dtav
  integer(kind=longint) :: idtav, tnext
  logical :: ldepcrosssection = .false.  !< Switch for doing depcrosssection (on/off)
  logical :: lbinary = .false.  !< Switch for binary output (on/off)


contains
  !> Initializing depcrosssection. Read out the namelist, initializing the variables
  subroutine initdepcrosssection
    use mpi
    use modmpi, only : myid, comm3d, MPI_LOGICAL, MY_REAL, myidx, myidy, &
        mpierr
    use modglobal, only : dtav_glob, ifnamopt, fname_options, &
        checknamelisterror, tres, btime, dt_lim, ladaptive, &
        dtmax, iexpnr, imax, jmax
    use modstat_nc, only : lnetcdf, open_nc, define_nc, ncinfo, &
        nctiminfo, writestat_dims_nc
    use moddrydeposition, only : ndeptracers

    implicit none

    integer :: ierr, isv, idt
    character(80) :: varname, varlongname 

    namelist/NAMDEPCROSSSECTION/ ldepcrosssection, lbinary, dtav

    dtav = dtav_glob
    if (myid==0) then
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read(ifnamopt, NAMDEPCROSSSECTION, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMDEPCROSSSECTION')
      write(6, NAMDEPCROSSSECTION)
      close(ifnamopt)
    end if

    if (ldepcrosssection .and. (.not. llsm .or. .not. ldrydep .or. ndeptracers == 0)) then
      ldepcrosssection = .false.
      write (6, *) "Ignoring ldepcrosssection, since no dry deposition &
        & and/or land surface model defined"
    end if

    call MPI_BCAST(dtav,             1, MY_REAL,     0, comm3d, mpierr)
    call MPI_BCAST(ldepcrosssection, 1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(lbinary,          1, MPI_LOGICAL, 0, comm3d, mpierr)

    idtav = dtav/tres
    tnext   = idtav+btime
    if(.not. ldepcrosssection) return
    dt_lim = min(dt_lim, tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'depcrosssection: dtav should be a integer multiple of dtmax'
    end if

    if ( lnetcdf ) then
      write(fname, '(a,i3.3,a,i3.3,a,i3.3,a)') 'depcross.', myidx, &
        '.', myidy, '.', iexpnr, '.nc'

      allocate(ncname(ndeptracers, 4))

      call nctiminfo(tncname(1, :))
      idt = 1
      do isv = 1, nsv
        if (.not. tracer_prop(isv)%ldep) cycle
        write (varname, '(a,a)') 'drydep_', trim(tracer_prop(isv)%tracname)
        write (varlongname, '(a,a)')  'Dry deposition flux of ', trim(tracer_prop(isv)%tracname)
        call ncinfo(ncname(idt, :), varname, varlongname, 'kg / (m2 * s)', 'tt0t')
        idt = idt+1
      end do
      call open_nc(fname, ncid, nrec, n1=imax, n2=jmax)
      if (nrec==0) then
        call define_nc(ncid, 1, tncname)
        call writestat_dims_nc(ncid)
        call define_nc(ncid, ndeptracers, ncname)
      end if
    end if
  end subroutine initdepcrosssection

  !> Do crosssection. Collect data to truncated (2 byte) integers, and write them to file
  subroutine depcrosssection
    use modglobal, only : rk3step, timee, dt_lim
    use modstat_nc, only : writestat_nc
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
    use modglobal, only : i1, j1, nsv, imax, jmax, ifoutput, iexpnr, &
                          rtimee
    use modfields, only : rhof
    use modstat_nc, only : lnetcdf, writestat_nc
    use moddrydeposition, only : ndeptracers

    implicit none

    real, allocatable :: depfield_massflux(:, :, :)
    character(80) :: txtfname
    integer :: isv, idt, i, j
    
    allocate(depfield_massflux(1:imax, 1:jmax, ndeptracers))

    ! Store the flux as a positive number
    depfield_massflux(1:imax, 1:jmax, 1:ndeptracers) = -depfield(2:i1, 2:j1, 1:ndeptracers) &
        * rhof(1) * 1e-6  ! to go from ug*m/(s*g) to kg/(m2*s)
    
    if (lbinary) then
      idt = 1
      do isv = 1, nsv
        if (.not. tracer_prop(isv)%ldep) cycle
        write(txtfname, '("movh_depsv",i3.3,"."i3.3)') isv, iexpnr
        open(ifoutput, file=txtfname, position='append', action='write')
        write(ifoutput, '(es12.5)') ((depfield_massflux(i, j, idt), i=1, imax), j=1, jmax)
        close(ifoutput)
        idt = idt + 1
      end do
    end if
    
    if (lnetcdf) then
      call writestat_nc(ncid, 1, tncname, (/rtimee/), nrec, .true.)
      call writestat_nc(ncid, ndeptracers, ncname, depfield_massflux, nrec, imax, jmax)
    end if

    deallocate(depfield_massflux)
  end subroutine wrtdrydepfields

  !> Clean up when leaving the run
  subroutine exitdepcrosssection
    use modstat_nc, only : exitstat_nc, lnetcdf
    implicit none

    if(ldepcrosssection .and. lnetcdf) then
      call exitstat_nc(ncid)
      deallocate(ncname)
    end if
  end subroutine exitdepcrosssection

end module moddepcrosssection
