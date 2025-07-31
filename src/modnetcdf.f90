!> \file modnetcdf.f90
!!  Convenience functions for working with NetCDF.
!>
!!  \author Caspar Jungbacker, Delft University of Technology
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
! Copyright 2024 Delft University of Technology
!
module modnetcdf

  use netcdf

  implicit none

  
  public :: check, readfile_dim_len, readfile1d, readfile2d

contains
 
  !> Checks return code of netcdf calls
  !!
  !! If an error occurs, stops the program and print information about
  !! the location of the error
  subroutine check(status, file, line)

    use modmpi, only: myid

    integer,      intent(in) :: status, line
    character(*), intent(in) :: file

    if (status /= nf90_noerr) then
      if (myid == 0) then
        write(*,*) "NetCDF error in: ", file, " on line: ", line
        write(*,*) trim(nf90_strerror(status))
      end if
      stop
    end if

  end subroutine check


  subroutine readfile1d(file,varname,out1d)

   use modglobal, only : handle_err
   character(*), intent(in) :: file
   character(*), intent(in) :: varname
   real, intent(inout)      :: out1d(:)


   integer ::  NCID, STATUS, varID, var_len
   var_len = size(out1d)


   STATUS = NF90_OPEN(file, nf90_nowrite, NCID)
   if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
  
   STATUS = NF90_INQ_DIMID(NCID, varname, varID)
   if (STATUS .ne. nf90_noerr) call handle_err(status)
  
   STATUS = NF90_GET_VAR (NCID, varID, out1d, start=(/1/), count=(/var_len/) )
   if (STATUS .ne. nf90_noerr) call handle_err(STATUS)


   STATUS = NF90_CLOSE(NCID)
   if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
  end subroutine readfile1d


  subroutine readfile_dim_len(file,dimname,dim_len)
      

   use modglobal, only : handle_err
   character(*), intent(in) :: file
   character(*), intent(in) :: dimname
   integer, intent(out)     :: dim_len

   
   integer ::  NCID, STATUS, dimID 
   character(len = nf90_max_name) :: RecordDimName

   STATUS = NF90_OPEN(file, nf90_nowrite, NCID)
   if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
  
   STATUS = NF90_INQ_DIMID(NCID, dimname, dimID)
   if (STATUS .ne. nf90_noerr) call handle_err(status)
  
   STATUS = nf90_INQUIRE_DIMENSION(NCID, dimID, len=dim_len)
   if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
 
   STATUS = NF90_CLOSE(NCID)
   if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
  
  end subroutine readfile_dim_len

  subroutine readfile2d(file,varname,out2d)

   use modmpi,    only : myidx, myidy
   use modglobal, only : handle_err, i1, j1,imax,jmax

   character(*), intent(in) :: file
   character(*), intent(in) :: varname
   real, intent(inout)      :: out2d(:,:,:)
  
   integer ::  NCID, STATUS, timeID, nt, varID 
   character(len = nf90_max_name) :: RecordDimName
   real, allocatable :: time(:)


   STATUS = NF90_OPEN(file, nf90_nowrite, NCID)
   if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
  
   STATUS = NF90_INQ_DIMID(NCID, "time", timeID)
   if (STATUS .ne. nf90_noerr) call handle_err(status)
  
   STATUS = nf90_INQUIRE_DIMENSION(NCID, timeID, len=nt, name=RecordDimName)
   if (STATUS .ne. nf90_noerr) call handle_err(STATUS)


   allocate(time(nt))
   STATUS = NF90_INQ_VARID(NCID, 'time', timeID)
   if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
   STATUS = NF90_GET_VAR (NCID, timeID, time, start=(/1/), count=(/nt/) )
   if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
   PRINT *, time

   STATUS = NF90_INQ_VARID(NCID,'tskin', VARID)
   if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
   STATUS = NF90_GET_VAR (NCID, VARID, out2d, start=(/myidx*imax+1,myidy*jmax+1,1/), &
                         & count=(/imax,jmax,nt/))
   if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
   STATUS = NF90_CLOSE(NCID)
   if (STATUS .ne. nf90_noerr) call handle_err(STATUS)


  end subroutine readfile2d

end module modnetcdf
