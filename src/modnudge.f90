!> \file modnudge.f90
!!  Nudges theta_l and q_t profiles to the initial profiles on a timescale tnudgeT
!>

!>
!!  Nudges theta_l and q_t profiles to the initial profiles on a timescale tnudgeT
!>
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!!  \todo Documentation
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

!cstep netcdf input file option. note that I have not tested the ncdf reading part yet.


module modnudge


implicit none
PRIVATE
PUBLIC :: initnudge, nudge,exitnudge
SAVE
  real, dimension(:,:), allocatable :: tnudge,unudge,vnudge,wnudge,thlnudge,qtnudge
  real, dimension(:)  , allocatable :: timenudge
  real :: tnudgefac = 1.
  logical :: lnudge = .false.,lunudge,lvnudge,lwnudge,lthlnudge,lqtnudge
  integer :: ntnudge = 100

contains
  subroutine initnudge
    use netcdf
    use modmpi,   only :myid,my_real,mpierr,comm3d,mpi_logical
    use modglobal,only :ifnamopt,fname_options,runtime,cexpnr,ifinput,k1,kmax,&
                        nc_input !cstep
    implicit none

    integer :: ierr,k,t
    real,allocatable,dimension(:) :: height
    character(1) :: chmess1

    !cstep : NCDF vars
    integer :: ncid,varID,sts ! checks for presence ncdf input file
    integer ::  dimIDtnudge
    integer ::  Nt_nudge
    character(len=nf90_max_name) :: tmpName

    namelist /NAMNUDGE/ &
       lnudge,tnudgefac
    allocate(tnudge(k1,ntnudge),unudge(k1,ntnudge),vnudge(k1,ntnudge),wnudge(k1,ntnudge),thlnudge(k1,ntnudge),qtnudge(k1,ntnudge))
    allocate(timenudge(0:ntnudge), height(k1))
    tnudge = 0
    unudge=0
    vnudge=0
    wnudge=0
    thlnudge=0
    qtnudge=0
    timenudge=0

    if(myid==0)then

      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMNUDGE,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMNUDGE'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMNUDGE'
      endif
      write(6 ,NAMNUDGE)
      close(ifnamopt)
    end if
    call MPI_BCAST(lnudge    , 1,MPI_LOGICAL,0,comm3d,mpierr)

    if (.not. lnudge) return

    if(myid==0) then
      if (nc_input) then 
          sts        = nf90_open('case_setup.nc',nf90_nowrite,ncid)
          if (sts.ne.nf90_noerr) stop 'ERROR: ncdf input file not found!'

          sts          = nf90_inq_dimid(ncid,"Ntimenudge",dimIDtnudge)
          if (sts.ne.nf90_noerr) stop  'ERROR:  nudging info not present'
          sts          = nf90_inquire_dimension(ncid, dimIDtnudge, tmpName,Nt_nudge)

          sts          = nf90_inq_varid(ncid,"timenudge",varID)
          
          sts          = nf90_get_var(ncid, varID, timenudge(1:Nt_nudge))

          sts          = nf90_inq_varid(ncid,"zf",varID)
          sts          = nf90_get_var(ncid, varID, height(1:kmax))

          sts          = nf90_inq_varid(ncid,"tnudge",varID)
          sts          = nf90_get_var(ncid, varID, tnudge(1:kmax,1:Nt_nudge))
    
          sts          = nf90_inq_varid(ncid,"unudge",varID)
          sts          = nf90_get_var(ncid, varID, unudge(1:kmax,1:Nt_nudge))

          sts          = nf90_inq_varid(ncid,"vnudge",varID)
          sts          = nf90_get_var(ncid, varID, vnudge(1:kmax,1:Nt_nudge))
  
          sts          = nf90_inq_varid(ncid,"wnudge",varID)
          sts          = nf90_get_var(ncid, varID, wnudge(1:kmax,1:Nt_nudge))      

          sts          = nf90_inq_varid(ncid,"thlnudge",varID)
          sts          = nf90_get_var(ncid, varID, thlnudge(1:kmax,1:Nt_nudge))

          sts          = nf90_inq_varid(ncid,"qtnudge",varID)
          sts          = nf90_get_var(ncid, varID, qtnudge(1:kmax,1:Nt_nudge))

      else
      t = 0
      open (ifinput,file='nudge.inp.'//cexpnr)

      do while (timenudge(t) < runtime)
        t = t + 1
        chmess1 = "#"
        ierr = 1 ! not zero
        !search for the next line consisting of "# time", from there onwards the profiles will be read
        do while (.not.(chmess1 == "#" .and. ierr ==0))
          read(ifinput,*,iostat=ierr) chmess1,timenudge(t)

           end do
           do  k=1,kmax
              read (ifinput,*) &
                height (k), &
                tnudge (k,t), &
                unudge (k,t), &
                vnudge (k,t), &
                wnudge (k,t), &
                thlnudge(k,t), &
                qtnudge(k,t)
            end do

         end do
         close(ifinput)
      endif !nc_input
   
      if (maxval(timenudge).lt.runtime) &
         stop 'STOP: No time dependend nudging data for end of run'

      tnudge  = tnudgefac*tnudge
      t = 0
      do while (timenudge(t) < runtime)
        t = t + 1
        write(6,*) ' height    t_nudge    u_nudge    v_nudge    w_nudge    thl_nudge    qt_nudge'
        do k=kmax,1,-1
          write (6,'(f7.1,6e12.4)') &
                height (k), &
                tnudge (k,t), &
                unudge (k,t), &
                vnudge (k,t), &
                wnudge (k,t), &
                thlnudge(k,t), &
                qtnudge(k,t)
        end do
      enddo
    end if

    call MPI_BCAST(timenudge ,ntnudge+1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tnudge    ,k1*ntnudge,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(unudge    ,k1*ntnudge,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(vnudge    ,k1*ntnudge,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(wnudge    ,k1*ntnudge,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(thlnudge  ,k1*ntnudge,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(qtnudge   ,k1*ntnudge,MY_REAL    ,0,comm3d,mpierr)
    lunudge = any(abs(unudge)>1e-8)
    lvnudge = any(abs(vnudge)>1e-8)
    lwnudge = any(abs(wnudge)>1e-8)
    lthlnudge = any(abs(thlnudge)>1e-8)
    lqtnudge  = any(abs(qtnudge)>1e-8)

  end subroutine initnudge

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nudge
    use modglobal, only : timee,rtimee,i1,j1,kmax,rdt
    use modfields, only : up,vp,wp,thlp, qtp,u0av,v0av,qt0av,thl0av
    implicit none

    integer k,t
    real :: dtm,dtp,currtnudge

    if (.not.(lnudge)) return
!     if (rk3step/=3) return
    if (timee==0) return

    t=1
    do while(rtimee>timenudge(t))
      t=t+1
    end do
    if (rtimee>timenudge(1)) then
      t=t-1
    end if

    dtm = ( rtimee-timenudge(t) ) / ( timenudge(t+1)-timenudge(t) )
    dtp = ( timenudge(t+1)-rtimee)/ ( timenudge(t+1)-timenudge(t) )

    do k=1,kmax
      currtnudge = max(rdt,tnudge(k,t)*dtp+tnudge(k,t+1)*dtm)
      if(lunudge  ) up  (2:i1,2:j1,k)=up  (2:i1,2:j1,k)-&
          (u0av  (k)-(unudge  (k,t)*dtp+unudge  (k,t+1)*dtm))/currtnudge
      if(lvnudge  ) vp  (2:i1,2:j1,k)=vp  (2:i1,2:j1,k)-&
          (v0av  (k)-(vnudge  (k,t)*dtp+vnudge  (k,t+1)*dtm))/currtnudge
      if(lwnudge  ) wp  (2:i1,2:j1,k)=wp  (2:i1,2:j1,k)-&
          (         -(wnudge  (k,t)*dtp+wnudge  (k,t+1)*dtm))/currtnudge
      if(lthlnudge) thlp(2:i1,2:j1,k)=thlp(2:i1,2:j1,k)-&
          (thl0av(k)-(thlnudge(k,t)*dtp+thlnudge(k,t+1)*dtm))/currtnudge
      if(lqtnudge ) qtp (2:i1,2:j1,k)=qtp (2:i1,2:j1,k)-&
          (qt0av (k)-(qtnudge (k,t)*dtp+qtnudge (k,t+1)*dtm))/currtnudge
    end do
  end subroutine nudge

  subroutine exitnudge
  deallocate(timenudge)
  end subroutine exitnudge

end module
