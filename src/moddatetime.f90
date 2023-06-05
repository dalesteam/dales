!> \ file moddatetime.f90
!! Utility module for date and time calculations

!>
!!  \author Marco de Bruine, VU
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
!  Copyright 1993-2020 Delft University of Technology, Wageningen
!  University, Utrecht University, KNMI, Vrije Universiteit Amsterdam
!

module moddatetime

  implicit none

  save

  ! Namelist variables
  logical :: l_datetime = .false.
  integer :: startyear  = 0, &
             startmonth = 0, &
             startday   = 0, &
             timezone   = 0  

  integer               :: jday
  integer, dimension(6) :: datex = 0, prevday = 0, nextday = 0

  contains
  
  subroutine initdatetime
    ! Using the namelist, construct datetime relevant variables

    use modmpi,    only : myid, comm3d, mpi_logical, D_MPI_BCAST
    use modglobal, only : ifnamopt, fname_options, xtime
   
    implicit none
   
    ! Auxiliary variables
    integer :: ierr

    ! Read & broadcast namelist DATETIME -----------------------------------
    namelist/NAMDATETIME/ l_datetime, startyear, startmonth, startday, timezone

    if (myid == 0) then

      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read(ifnamopt, NAMDATETIME, iostat=ierr)

      if (ierr > 0) then
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMDATETIME'
      endif

      write(6, NAMDATETIME)
      close(ifnamopt)

    endif

    call d_mpi_bcast(l_datetime, 1, 0, comm3d, ierr)
    call d_mpi_bcast(startyear,  1, 0, comm3d, ierr)
    call d_mpi_bcast(startmonth, 1, 0, comm3d, ierr)
    call d_mpi_bcast(startday,   1, 0, comm3d, ierr)
    call d_mpi_bcast(timezone,   1, 0, comm3d, ierr)

    ! Initialize datetime array -----------------------------------------------
    if (l_datetime) then

      datex(1) = startyear        ! year
      datex(2) = startmonth       ! month

      if (xtime + timezone > 24) then
        datex(3) = startday + 1   ! day (LT is a day ahead of GMT)
      elseif (xtime + timezone < 0) then
        datex(3) = startday - 1   ! day (LT is a day behind GMT)
      else
        datex(3) = startday       ! day
      endif

      datex(4) = mod(int(xtime) + timezone + 24, 24) ! hour
      datex(5) = 0                ! minutes
      datex(6) = 0                ! seconds
  
      ! Initialize julian day number --------------------------------------------
      jday = julday(datex(1), datex(2), datex(3)) 

      call caldat(jday-1, prevday(1), prevday(2), prevday(3))
      call caldat(jday+1, nextday(1), nextday(2), nextday(3))

      ! Catch non specified starting date while l_datetime = True ---------------
      ! We check on date only, because time is filled from variables shared with
      ! other functionalities, that thus have non-zero values.  
      if (sum(datex(1:3)) == 0) stop 'ERROR: Trying to use datetime functionality without specified date.' 
    endif
  end subroutine initdatetime

  subroutine datetime
    ! --------------------------------------------------------------------------
    !  
    !  Keep track of the (local) date and time by updating the datex array.
    !
    ! --------------------------------------------------------------------------
   
    use modglobal, only : rtimee, xtime 
    ! TODO: Check if rtimee really is seconds of if it is scaled

    implicit none
    integer :: xhr, itimee

    itimee = int(rtimee)
    xhr = mod(itimee/3600 + int(xtime) + timezone + 24, 24)

    if ( xhr < datex(4) ) then ! New day
      jday = jday + 1
      call caldat(jday, datex(1), datex(2), datex(3))
    endif

    datex(4) = xhr                ! hours
    datex(5) = mod(itimee/60, 60) ! minutes
    datex(6) = mod(itimee,    60) ! seconds

  end subroutine datetime

  subroutine caldat(julian, iyyy, mm, id)
    !--------------------------------------------------------------------------
    !  caldat
    !
    !  purpose
    !  -------
    !  calculate date from given julian day
    !
    !  parameters
    !  ----------
    !  on input : julday contains the julian day
    !  on output: mm, id, iyyy contain month, day and year
    !
    !  reference
    !  ---------
    !  J. Meeus, "Astronomical formulea for calculators, 4th Edition" 1988
    !--------------------------------------------------------------------------
    implicit none

    ! input/output
    integer,intent(in)  :: julian
    integer,intent(out) :: mm
    integer,intent(out) :: id
    integer,intent(out) :: iyyy

    ! local
    integer,parameter   :: igreg=2299161
    integer             :: jalpha, ja, jb, jc, jd, je
    !
    ! handle gregorian and julian date
    !
    if ( julian >= igreg )then
       jalpha=int(((julian-1867216)-0.25)/36524.25)
       ja=julian+1+jalpha-int(0.25*jalpha)
    else
       ja=julian
    end if
    jb=ja+1524
    jc=int(6680.+((jb-2439870)-122.1)/365.25)
    jd=365*jc+int(0.25*jc)
    je=int((jb-jd)/30.6001)
    id=jb-jd-int(30.6001*je)
    mm=je-1
    if ( mm > 12 ) mm=mm-12
    iyyy=jc-4715
    if ( mm > 2 ) iyyy=iyyy-1
    !
    ! handle dates before 0 AD
    !
    if ( iyyy <= 0 ) iyyy=iyyy-1

  end subroutine caldat

  integer function julday(iy,mm,id)
    !-----------------------------------------------------------------------
    !**** julday
    !
    !  purpose
    !  -------
    !  calculate julian day from given date
    !
    !  parameters
    !  ----------
    !  on input : mm, id, iyyy contain month, day and year
    !  on output: julday contains the julian day
    !
    !  reference
    !  ---------
    !  J. Meeus, "Astronomical formulea for calculators, 4th Edition" 1988
    !-----------------------------------------------------------------------
    implicit none

    ! input, output
    integer,intent(in) :: mm  ! month
    integer,intent(in) :: id  ! day
    integer,intent(in) :: iy  ! year

    ! local
    integer,parameter  :: igreg=15+31*(10+12*1582)
    integer            :: jy, jm, ja, iyyy

    ! handle dates before 0 AD
    !
    iyyy=iy
    if ( iy == 0 ) then
       stop 'julday:  ERROR invalid year 0 AD'
    end if
    if ( iy < 0 ) then
       iyyy=iy+1
    end if
    !
    !calculate julian day from date in gregorian calendar
    !
    if ( mm > 2 ) then
       jy=iyyy
       jm=mm+1
    else
       jy=iyyy-1
       jm=mm+13
    end if
    julday=int(365.25*jy)+int(30.6001*jm)+id+1720995
    !
    !handle julian calender
    !
    if ( id+31*(mm+12*iyyy) >= igreg ) then
       ja=int(0.01*jy)
       julday=julday+2-ja+int(0.25*ja)
    end if

  end function julday

end module moddatetime
