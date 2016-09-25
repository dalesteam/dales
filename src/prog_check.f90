!> \file prog_check.f90
!! Binary using the daleslib api for testing purposes.
!  This file is part of DALESLIB.
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

program dales_check

    use daleslib,  only: initialize,step,finalize,grid_shape,get_field_3d,&
                        &get_field_layer_avg,FIELDID_THL
    use modglobal, only: timeleft,rk3step
    use modmpi,    only: myid

    character(512)      :: fname_options
    integer             :: s(3), ierr
    real,allocatable    ::a(:,:,:)
    real,allocatable    ::b(:)

    if (command_argument_count() >=1) then
        call get_command_argument(1,fname_options)
    end if

    call initialize(fname_options)

    s=grid_shape()
    allocate(a(s(1),s(2),s(3)))
    allocate(b(s(3)))
    ierr=0

    do while (timeleft>0 .or. rk3step < 3)
        call step()
        !ierr=get_field_3d(FIELDID_THL,a)
        ierr=get_field_layer_avg(FIELDID_THL,b)
        if(myid==0) then
            write(*,*) "top layer theta liquid is ",b(1)
        endif
    enddo

    deallocate(a)

    call finalize

end program dales_check
    
