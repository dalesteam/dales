!> \file modspraying.f99
!! Stephan de Roode and Annelot Broerze



module modspraying
   use modprecision, only: field_r
   use modsprayingdata, only: i_glob_spray,j_glob_spray,k_glob_spray,&
                              i_loc_spray,j_loc_spray,k_loc_spray,&
                              water_spray_rate,salt_spray_rate,&
                              lwater_spraying,lsalt_spraying,salinity,&
                              isv_salt,tracer,lsalt_sponge
   implicit none

  ! for lateral sponge
  real(field_r), allocatable :: fnudgeloc(:,:) ! local, cpu dependent array of fnudge values
  integer :: nudgedepthgr = 10 ! number of lateral nudge grid points

contains
  subroutine initspraying
  use modglobal,    only : i1,j1,imax,jmax,kmax,ifnamopt,fname_options,checknamelisterror
  use modmpi,       only : myid,myidx,myidy,comm3d, mpierr, d_mpi_bcast
  use modtracers,   only: add_tracer
  !use modnudgeboundary, only : lnudgeboundary


  integer ierr

  namelist/NAMSPRAYING/ lwater_spraying,lsalt_spraying,&
                        i_glob_spray,j_glob_spray,k_glob_spray,&
                        water_spray_rate,salt_spray_rate,salinity,&
                        tracer,lsalt_sponge

  if(myid==0) then    !first myid
    open(ifnamopt,file=fname_options,status='old',iostat=ierr)
    read (ifnamopt,NAMSPRAYING,iostat=ierr)
    call checknamelisterror(ierr, ifnamopt, 'NAMSPRAYING')
    write(6 ,NAMSPRAYING)
    close(ifnamopt)
  endif

  call D_MPI_BCAST(lwater_spraying ,    1,  0, comm3d, mpierr)
  call D_MPI_BCAST(lsalt_spraying  ,    1,  0, comm3d, mpierr)
  call D_MPI_BCAST(i_glob_spray    ,    1,  0, comm3d, mpierr)
  call D_MPI_BCAST(j_glob_spray    ,    1,  0, comm3d, mpierr)
  call D_MPI_BCAST(k_glob_spray    ,    1,  0, comm3d, mpierr)
  call D_MPI_BCAST(water_spray_rate,    1,  0, comm3d, mpierr)
  call D_MPI_BCAST(salt_spray_rate,     1,  0, comm3d, mpierr)
  call D_MPI_BCAST(tracer,             20,  0, comm3d, mpierr)
  call D_MPI_BCAST(lsalt_sponge,        1,  0, comm3d, mpierr)
  
  if (lwater_spraying) then
     lsalt_spraying  = .true.
     salt_spray_rate = water_spray_rate * salinity   !directy coupled to water spray rate
  else
     water_spray_rate = 0.
  endif

  if (lsalt_spraying) then
     call add_tracer(trim(tracer), long_name=trim(tracer)//" mixing ratio", &
          unit="kg/kg", isv=isv_salt)
  endif

  !determine local position of spraying from global position
  i_loc_spray = i_glob_spray - myidx*imax
  j_loc_spray = j_glob_spray - myidy*jmax
  k_loc_spray = k_glob_spray

  ! are the local coordinates actually in the domain?
  if (i_loc_spray >= 2 .and. i_loc_spray <= i1 .and. &
       j_loc_spray >= 2 .and. j_loc_spray <= j1 .and. &
       k_loc_spray >= 1 .and. k_loc_spray <= kmax) then
     write(6,*) 'spraying point at myid = ',myid
     write(6,*) 'global locations ',i_glob_spray,j_glob_spray,k_glob_spray
     write(6,*) 'local locations ',i_loc_spray,j_loc_spray,k_loc_spray
  else  ! if not, there is no sprayer here
     i_loc_spray = -999
     j_loc_spray = -999
     k_loc_spray = -999
  endif

  if (myid==0) then
     write(6,*) 'Spraying data used: '
     write(6,*) 'lwater_spraying     ',lwater_spraying
     write(6,*) 'lsalt_spraying      ',lsalt_spraying
     write(6,*) 'i_glob_spray        ',i_glob_spray
     write(6,*) 'j_glob_spray        ',j_glob_spray
     write(6,*) 'k_glob_spray        ',k_glob_spray
     write(6,*) 'water_spray_rate    ',water_spray_rate
     write(6,*) 'salt_spray_rate     ',salt_spray_rate
     write(6,*) 'salt scalar number  ',isv_salt
     write(6,*)
  endif

  if (lsalt_spraying .and. lsalt_sponge) then
     call initlateralsponge ! TODO: move this to be independent of spraying
  end if

  end subroutine initspraying


  ! smooth nudging of a scalar field to 0 at the boundary - a lateral sponge
  ! modified from code by Pim van Dorp 2015
  subroutine initlateralsponge
    use modglobal, only: itot, jtot, imax, jmax, i1, j1, ih, jh, pi
    use modmpi, only : myidx,myidy
    real(field_r) :: fnudge
    integer i, j, iglob, jglob
    real(field_r), allocatable :: fnudgeglob(:,:) ! global array of fnudge values

    allocate(fnudgeglob(1-ih:itot+ih,1-jh:jtot+jh))
    allocate(fnudgeloc(2-ih:i1+ih,2-jh:j1+jh))
    fnudgeglob = 0
    fnudgeloc = 0

    ! construct a 2D field of nudging constants
    do i=nudgedepthgr,1,-1
      fnudge = 0.5 + 0.5*cos((pi/(nudgedepthgr-1))*(i-1))

      fnudgeglob(i,i:jtot-i+1) = fnudge
      fnudgeglob(itot-i+1,i:jtot-i+1) = fnudge
      fnudgeglob(i+1:(itot-i),i) = fnudge
      fnudgeglob(i+1:(itot-i),jtot-i+1) = fnudge
   end do

   ! cut out the part for this processor
   do j=2,j1
      do i=2,i1
         iglob = i + imax*myidx - 1
         jglob = j + jmax*myidy - 1
         fnudgeloc(i,j) = fnudgeglob(iglob,jglob)
      end do
   end do

   deallocate(fnudgeglob)

  end subroutine initlateralsponge

  subroutine lateralsponge
    use modglobal, only : kmax, i1, j1, rdt, nsv
    use modfields, only : svp,sv0 !lsv_nudge_at_boundary

    integer i,j,k,isv
    if (.not. lsalt_sponge) return

    !if (nsv.gt.0) then
       !do isv=1,nsv
       !if (lsv_nudge_at_boundary(isv)) then
    isv = isv_salt  ! TODO: for now only nudge the sprayed salt scalar to 0
       do k=1,kmax
          do j=2,j1
             do i=2,i1
                ! svp(i,j,k,isv) = (1-fnudgeloc(i,j,k))*svp(i,j,k,isv) + fnudgeloc(i,j,k)*(0-sv0(i,j,k,isv))/rdt
                ! the original nudged also the tendency towards 0

                svp(i,j,k,isv) = svp(i,j,k,isv) + fnudgeloc(i,j)*(0-sv0(i,j,k,isv))/rdt
             end do
          end do
       end do
    !endif
    !end do
    !endif

  end subroutine lateralsponge


  subroutine exitlateralsponge
    deallocate(fnudgeloc)
  end subroutine exitlateralsponge


end module modspraying
