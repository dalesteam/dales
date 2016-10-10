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
                        &get_field_layer_avg,get_field_2d,FIELDID_QT
    use modglobal, only: timeleft,rk3step
    use modmpi,    only: myid,commwrld

    character(512)      :: fname_options
    integer             :: s(3),sx,sy,sz,ierr,i
    real,allocatable    :: a(:,:,:)
    real,allocatable    :: b(:,:)
    real,allocatable    :: c(:)
    real                :: ca,cb,abserra,abserrb

    integer, parameter  ::fid=FIELDID_QT

    if (command_argument_count() >=1) then
        call get_command_argument(1,fname_options)
    end if

    call initialize(fname_options)

    s=grid_shape()
    sx=s(1)
    sy=s(2)
    sz=s(3)
    allocate(a(sx,sy,sz))
    allocate(b(sx,sy))
    allocate(c(sz))
    ierr=0

    !call MPI_BARRIER(commwrld,ierr)
    
    do while (timeleft>0 .or. rk3step < 3)
        call step()
        ierr=get_field_3d(fid,a)
        ierr=get_field_layer_avg(fid,c)
        if(myid==0) then
            write(*,*) 'delta ca           delta cb'
        endif
        abserra=0
        abserrb=0
        do i=1,sz
            ierr=get_field_2d(fid,i,b)
            ca=sum(a(:,:,i))/(sx*sy)
            cb=sum(b)/(sx*sy)
            abserra=max(abserra,abs(c(i)-ca))
            abserrb=max(abserrb,abs(c(i)-cb))
        enddo
        if(myid==0) then
            write(*,*) abserra,abserrb
        endif
    enddo

    deallocate(a)

    call finalize

end program dales_check 
