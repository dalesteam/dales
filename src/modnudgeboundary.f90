!> \file modnudgeboundary.f90
!! By Pim van Dorp (PVD), TU Delft, section Atmospheric Physics, 10 dec 2015
!! Nudge boundary to prescribed values

!! 4 oct 2021 - cstep allocation of variables only done if lnudgeboundary is true

module modnudgeboundary 
  use modglobal, only : longint

  implicit none 
  logical :: lnudgeboundary = .false. 
  logical :: lstatref = .false. 


  real, allocatable :: fnudgeglob(:,:,:) ! global array of fnudge values
  real, allocatable :: fnudgeloc(:,:,:) ! local, cpu dependent array of fnudge values

  integer :: nudgedepthgr = 10 ! number of nudge grid points
  

contains
  subroutine initnudgeboundary
    use modglobal, only: itot, jtot, kmax, i1, i2, j1, j2, k1, ih, jh,ifnamopt, fname_options, tres, ladaptive, dtmax,btime, dx, dy , pi, dt
    use modmpi, only : myidx,myidy,myid, comm3d,mpierr, d_mpi_bcast
    use mpi

    implicit none
    integer ierr


    namelist/NUDGEBOUNDARY/ lnudgeboundary, lstatref, nudgedepthgr

    if (myid==0) then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
        read (ifnamopt,NUDGEBOUNDARY,iostat=ierr)
        if (ierr > 0) then
          print *, 'Problem in namoptions NUDGEBOUNDARY'
          print *, 'iostat error: ', ierr
          stop 'ERROR: Problem in namoptions NUDGEBOUNDARY'
        endif
        write(6 ,NUDGEBOUNDARY)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(lnudgeboundary     ,1,0,comm3d,mpierr) 
    call D_MPI_BCAST(nudgedepthgr       ,1,0,comm3d,mpierr) 


    if (lnudgeboundary) then 

      allocate(fnudgeglob(1-ih:itot+ih,1-jh:jtot+jh,1:k1))
      allocate(fnudgeloc(2-ih:i1+ih,2-jh:j1+jh,k1))

      call calcfnudge

      if (myid == 0 ) then
        write(*,*) 'nudgedepthgr = ', nudgedepthgr 
        write(*,*) 'Succesfully initialized modnudgeboundary'
      end if
  else
       if (myid == 0 ) then
        write(*,*) 'no lateral nudging boundary  used'
        write(*,*) 'no modnudgeboundary variables defined'
      end if
  endif

  end subroutine initnudgeboundary

  subroutine calcfnudge
    use modglobal, only : pi, itot, jtot, ih, jh, k1, j1, i1, kmax
!cstep    use modwindturbinedata, only : turhzgr
    use modmpi, only : myidx, myidy

    implicit none
    integer i,j,k
    real fnudge

    fnudgeglob = 0.
    fnudgeloc = 0.

    do i=1,nudgedepthgr 

      fnudge = 0.5 + 0.5*COS((pi/(nudgedepthgr-1))*(i-1))

      fnudgeglob(i,i:jtot-i+1,:) = fnudge
      fnudgeglob(itot-i+1,i:jtot-i+1,:) = fnudge
      fnudgeglob(i+1:(itot-i),i,:) = fnudge
      fnudgeglob(i+1:(itot-i),jtot-i+1,:) = fnudge

    end do

    do k=1,kmax
      do j=2,j1
        do i=2,i1
          fnudgeloc(i,j,k) = fnudgeglob(iglob(i,myidx),jglob(j,myidy),k)
        end do
      end do
    end do

  end subroutine calcfnudge


  subroutine nudgeboundary 
    use modglobal, only : kmax, i1, j1, rdt,&
                          nsv
    use modmpi, only : myid !cstep can be removed if write statement is not needed anymore
    use modfields, only : &
                          svp,sv0

    implicit none
    integer i,j,k


    !if (myid.eq.0) then
    !   write(6,*) 'scalars are forced to zero at the lateral boundaries'  !assumes local emissions
    !endif
    if (nsv.gt.0) then
      do k=1,kmax
       do j=2,j1
        do i=2,i1
          svp(i,j,k,1:nsv) = (1-fnudgeloc(i,j,k))*svp(i,j,k,1:nsv) + fnudgeloc(i,j,k)*(0-sv0(i,j,k,1:nsv))/rdt
        end do
      end do
     end do
    endif

  end subroutine nudgeboundary

  subroutine exitnudgeboundary
    if (lnudgeboundary) then
      deallocate(fnudgeglob,fnudgeloc)
    endif
  end subroutine exitnudgeboundary

  function iglob(iloc,myidxloc)
    use modglobal, only : imax

    implicit none
    integer iloc,iglob,myidxloc

    iglob = iloc + imax*myidxloc - 1

  end function iglob

  function jglob(jloc,myidyloc)
    use modglobal, only : jmax

    implicit none
    integer jloc,jglob,myidyloc

    jglob = jloc + jmax*myidyloc - 1

  end function jglob

end module modnudgeboundary 

