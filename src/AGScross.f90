!> \file modAGScross.f90
!!   Dumps an instantenous AGScross of the field

!>
!! Dumps an instantenous AGScross of the field.
!>
!! AGScrosss in the yz-plane and in the xy-plane            |
    !        of u,v,w,thl,thv,qt,ql. Written to movv_*.expnr and movh_*.expnr
!! If netcdf is true, this module leads the cross.myid.expnr.nc output

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
module modAGScross


  use modglobal, only : longint,kmax

implicit none
private
PUBLIC :: initAGScross, AGScross,exitAGScross
save
!NetCDF variables
  integer,parameter :: nvar = 20
  integer :: ncidAGS = 123
  integer :: nrecAGS = 0
  character(80) :: fnameAGS = 'crossAGS.xxx.xxx.nc'
  character(80),dimension(nvar,4) :: ncnameAGS
  character(80),dimension(1,4) :: tncnameAGS

  real    :: dtav
  integer(kind=longint) :: idtav,tnext
  logical :: lAGScross = .false. !< switch for doing the AGScross (on/off)

contains
!> Initializing AGScross. Read out the namelist, initializing the variables
  subroutine initAGScross
    use modmpi,   only :myid,my_real,mpierr,comm3d,mpi_logical,mpi_integer,cmyid
    use modglobal,only :imax,jmax,ifnamopt,fname_options,dtmax,rk3step, dtav_glob,ladaptive,j1,kmax,i1,dt_lim,cexpnr,tres,btime
    use modstat_nc,only : open_nc, define_nc,ncinfo,writestat_dims_nc
    use modsurfdata, only : lrsAgs, ksoilmax
   implicit none

    integer :: ierr,k

    namelist/NAMAGScross/ &
    lAGScross, dtav

    dtav = dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMAGScross,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMAGScross'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMAGScross'
      endif
      write(6 ,NAMAGScross)
      close(ifnamopt)
    end if

    if (.not. lrsAgs) lAGScross = .false.

    call MPI_BCAST(dtav       ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(lAGScross     ,1,MPI_LOGICAL,0,comm3d,mpierr)

    idtav = dtav/tres
    tnext   = idtav+btime
    if(.not.(lAGScross)) return
    dt_lim = min(dt_lim,tnext)
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'AGScross: dtav should be a integer multiple of dtmax'
    end if
    if (ksoilmax /= 4) stop 'ksoilmax is not equal to 4... this can give problems with AGScross.f90... update this file as well'
    fnameAGS(10:12) = cmyid
    fnameAGS(14:16) = cexpnr
    call ncinfo(tncnameAGS(1,:),'time  ','Time','s','time')
    call ncinfo(ncnameAGS( 1,:),'An    ', 'xy AGScross of An          ','mg/m2/s','tt0t')
    call ncinfo(ncnameAGS( 2,:),'Resp  ', 'xy AGScross of Resp        ','mg/m2/s','tt0t')
    call ncinfo(ncnameAGS( 3,:),'wco2  ', 'xy AGScross of wco2        ','ppm m/s','tt0t')
    call ncinfo(ncnameAGS( 4,:),'rs    ', 'xy AGScross of diagnosed rs','s/m    ','tt0t')
    call ncinfo(ncnameAGS( 5,:),'ra    ', 'xy AGScross of ra          ','s/m    ','tt0t')
    call ncinfo(ncnameAGS( 6,:),'rsCO2 ', 'xy AGScross of rsCO2       ','s/m    ','tt0t')
    call ncinfo(ncnameAGS( 7,:),'rsveg ', 'xy AGScross of rsveg=rsAgs ','s/m    ','tt0t')
    call ncinfo(ncnameAGS( 8,:),'rssoil', 'xy AGScross of rssoil      ','s/m    ','tt0t')
    call ncinfo(ncnameAGS( 9,:),'fstr  ', 'xy AGScross of stress fnct.','-      ','tt0t')
    call ncinfo(ncnameAGS(10,:),'phiw1 ', 'xy AGScross of phiw top    ','-      ','tt0t')
    call ncinfo(ncnameAGS(11,:),'phiw2 ', 'xy AGScross of phiw level 2','-      ','tt0t')
    call ncinfo(ncnameAGS(12,:),'phiw3 ', 'xy AGScross of phiw level 3','-      ','tt0t')
    call ncinfo(ncnameAGS(13,:),'phiw4 ', 'xy AGScross of phiw level 4','-      ','tt0t')
    call ncinfo(ncnameAGS(14,:),'CO2   ', 'xy AGScross of CO2 (grid 1)','ppm    ','tt0t')
    call ncinfo(ncnameAGS(15,:),'tskin ', 'xy AGScross of curr. tskin ','K      ','tt0t')
    call ncinfo(ncnameAGS(16,:),'tskinm', 'xy AGScross of prev. tskin ','K      ','tt0t')
    call ncinfo(ncnameAGS(17,:),'tsoil1', 'xy AGScross of tsoil top   ','K      ','tt0t')
    call ncinfo(ncnameAGS(18,:),'tsoil2', 'xy AGScross of tsoil lvl 2 ','K      ','tt0t')
    call ncinfo(ncnameAGS(19,:),'tsoil3', 'xy AGScross of tsoil lvl 3 ','K      ','tt0t')
    call ncinfo(ncnameAGS(20,:),'tsoil4', 'xy AGScross of tsoil lvl 4 ','K      ','tt0t')
    call open_nc(fnameAGS,  ncidAGS,nrecAGS,n1=imax,n2=jmax)
    if (nrecAGS == 0) then
      call define_nc( ncidAGS, 1, tncnameAGS)
      call writestat_dims_nc(ncidAGS)
    end if
    call define_nc( ncidAGS, NVar, ncnameAGS)

  end subroutine initAGScross
!>Run AGScross. Mainly timekeeping
  subroutine AGScross
    use modglobal, only : rk3step,timee,rtimee,dt_lim
    use modstat_nc, only : writestat_nc
    implicit none


    if (.not. lAGScross) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    call AGShorz

  end subroutine AGScross


!> Do the xy AGScrosss and dump them to file
  subroutine AGShorz
    use modglobal, only : imax,jmax,i1,j1,rtimee
    use modstat_nc, only : writestat_nc
    use modsurfdata, only : AnField, RespField, wco2Field,phiw,fstrField, rs, ra, rsco2Field, rsveg, rssoil, indCO2, tskin, tskinm, tsoil
    use modfields, only   : svm
    implicit none


    ! LOCAL
    integer i,j,n
    character(40) :: name
    real, allocatable :: vars(:,:,:)

      allocate(vars(1:imax,1:jmax,nvar))
      vars=0.
      vars(:,:, 1) = AnField   (2:i1,2:j1)
      vars(:,:, 2) = RespField (2:i1,2:j1)
      vars(:,:, 3) = wco2Field (2:i1,2:j1)
      vars(:,:, 4) = rs        (2:i1,2:j1)
      vars(:,:, 5) = ra        (2:i1,2:j1)
      vars(:,:, 6) = rsco2Field(2:i1,2:j1)
      vars(:,:, 7) = rsveg     (2:i1,2:j1)    
      vars(:,:, 8) = rssoil    (2:i1,2:j1)
      vars(:,:, 9) = fstrField (2:i1,2:j1)
      vars(:,:,10) = phiw      (2:i1,2:j1,1)
      vars(:,:,11) = phiw      (2:i1,2:j1,2)
      vars(:,:,12) = phiw      (2:i1,2:j1,3)
      vars(:,:,13) = phiw      (2:i1,2:j1,4)
      vars(:,:,14) = svm       (2:i1,2:j1,1,indCO2) / 1000.0
      vars(:,:,15) = tskin     (2:i1,2:j1)
      vars(:,:,16) = tskinm    (2:i1,2:j1)
      vars(:,:,17) = tsoil     (2:i1,2:j1,1)
      vars(:,:,18) = tsoil     (2:i1,2:j1,2)
      vars(:,:,19) = tsoil     (2:i1,2:j1,3)
      vars(:,:,20) = tsoil     (2:i1,2:j1,4)
      call writestat_nc(ncidAGS,1,tncnameAGS,(/rtimee/),nrecAGS,.true.)
      call writestat_nc(ncidAGS,nvar,ncnameAGS(1:nvar,:),vars,nrecAGS,imax,jmax)
      deallocate(vars)

  end subroutine AGShorz


!> Clean up when leaving the run
  subroutine exitAGScross
    use modstat_nc, only : exitstat_nc
    implicit none

    if(lAGScross) then
    call exitstat_nc(ncidAGS)
    end if

  end subroutine exitAGScross

end module modAGScross
