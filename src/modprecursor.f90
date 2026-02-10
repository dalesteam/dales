!> \file modnudgeboundary.f90
!! By Pim van Dorp (PVD), TU Delft, section Atmospheric Physics, 10 dec 2015
!! Nudge boundary to prescribed values
!!  Modifications by Dion Hoeksema (DH), section Geoscience and Remote sensing, 21 October 2024
!!      Modification for compatibility with newest DALES v4.4: MPI_BCAST -> D_MPI_BCAST, MPI_ALLREDUCE -> D_MPI_ALLREDUCE

module modprecursor
  use modglobal, only : longint, nsv
  use modfields, only : sv0, svm, svp, sv0av
  use modtimer, only: timer_tic, timer_toc
  use modprecision, only: field_r

  implicit none 

  character(len=*), parameter :: modname = 'modprecursor'

  !> Pointer swap
  interface swap
    module procedure :: swap_4d
    module procedure :: swap_3d
    module procedure :: swap_2d
    module procedure :: swap_1d
  end interface swap

  logical :: lprecursor = .false. 
  logical :: lstatref = .false. 

  integer :: Nsim = 1 
  integer :: statid = 1
  integer :: refid = 1
  integer :: turid = 1

  real(field_r), allocatable :: fnudgeglob(:,:,:) ! global array of fnudge values
  real(field_r), allocatable :: fnudgeloc(:,:,:) ! local, cpu dependent array of fnudge values

  integer :: nudgedepthgr = 10 ! number of nudge grid points
  
  ! Prognostic variables; first dimension: 1=turbine simulation, 2=reference simulation
  real(field_r), pointer :: umsave(:,:,:)        !<   x-component of velocity at time step t-1
  real(field_r), pointer :: vmsave(:,:,:)        !<   y-component of velocity at time step t-1
  real(field_r), pointer :: wmsave(:,:,:)        !<   z-component of velocity at time step t-1
  real(field_r), pointer :: thlmsave(:,:,:)      !<   liq. water pot. temperature at time step t-1
  real(field_r), pointer :: e12msave(:,:,:)      !<   square root of turb. kin. energy at time step t-1
  real(field_r), pointer :: qtmsave(:,:,:)       !<   total specific humidity at time step t

  real(field_r), pointer :: u0save(:,:,:)        !<   x-component of velocity at time step t
  real(field_r), pointer :: v0save(:,:,:)        !<   y-component of velocity at time step t
  real(field_r), pointer :: w0save(:,:,:)        !<   z-component of velocity at time step t
  real(field_r), pointer :: thl0save(:,:,:)      !<   liq. water pot. temperature at time step t
  real(field_r), pointer :: qt0save(:,:,:)       !<   total specific humidity at time step t
  real(field_r), pointer :: ql0save(:,:,:)   
  real(field_r), pointer :: ql0hsave(:,:,:)  
  real(field_r), pointer :: e120save(:,:,:)      !<   square root of turb. kin. energy at time step t

  real(field_r), pointer :: dthvdzsave(:,:,:)  
  real(field_r), pointer :: ekmsave(:,:,:)  
  real(field_r), pointer :: tmp0save(:,:,:)  
  real(field_r), pointer :: eslsave(:,:,:)  
  real(field_r), pointer :: qvslsave(:,:,:)  
  real(field_r), pointer :: qvsisave(:,:,:)  

  real(field_r), pointer :: thv0hsave(:,:,:)  

  real(field_r), pointer :: presfsave(:)  
  real(field_r), pointer :: preshsave(:)  

  real(field_r), pointer :: thvhsave(:)  

  real(field_r), pointer :: u0avsave(:)
  real(field_r), pointer :: v0avsave(:)
  real(field_r), pointer :: thl0avsave(:)
  real(field_r), pointer :: qt0avsave(:)

  real(field_r), pointer :: svmsave(:,:,:,:)
  real(field_r), pointer :: sv0save(:,:,:,:)
  real(field_r), pointer :: sv0avsave(:,:)

contains
  subroutine init_precursor
    use modglobal, only: itot, jtot, kmax, i1, i2, j1, j2, k1, ih, jh,ifnamopt, fname_options, tres, ladaptive, dtmax,btime, dx, dy , pi, dt
    use modmpi, only : myidx,myidy,myid,MPI_INTEGER,D_MPI_BCAST, MPI_SUM,MPI_COMM_WORLD,MPI_LOGICAL,comm3d,mpierr
    use modfields, only : u0, v0, w0, e120, thl0, qt0, ql0, ql0h, tmp0, &
                          um, vm, wm, e12m, thlm, qtm, &
                          presf, presh, dthvdz, esl, qvsl, qvsi, thv0h, &
                          u0av, v0av, thl0av, qt0av, thvh
    use modsubgrid, only : ekm

    implicit none
    integer i,j,k,n,ierr, simid


    namelist/precursor/ lprecursor, lstatref, nudgedepthgr

    if (myid==0) then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
        read (ifnamopt,precursor,iostat=ierr)
        if (ierr > 0) then
          print *, 'Problem in namoptions precursor'
          print *, 'iostat error: ', ierr
          stop 'ERROR: Problem in namoptions precursor'
        endif
        write(6 ,precursor)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(lprecursor ,1,0,MPI_COMM_WORLD,mpierr) !DH
    call D_MPI_BCAST(lstatref,1,0,MPI_COMM_WORLD,mpierr) !DH
    call D_MPI_BCAST(nudgedepthgr       ,1,0,MPI_COMM_WORLD,mpierr) !DH

    if (.not. lprecursor) return

    if (lprecursor) Nsim = 2
    if (lprecursor) turid = 2
    if (.not. lstatref .and. Nsim == 2) statid = 2

    ! Which of these are really needed?
    ! Maybe save some memory here
    allocate(umsave(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(vmsave(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(wmsave(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thlmsave(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(e12msave(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(qtmsave(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(svmsave(2-ih:i1+ih,2-jh:j1+jh,k1,1:nsv))

    allocate(u0save(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(v0save(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(w0save(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thl0save(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(qt0save(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(ql0save(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(ql0hsave(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(e120save(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(dthvdzsave(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(ekmsave(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(tmp0save(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(eslsave(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(qvslsave(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(qvsisave(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(sv0save(2-ih:i1+ih,2-jh:j1+jh,k1,1:nsv))

    allocate(thv0hsave(2-ih:i1+ih,2-jh:j1+jh,k1))

    allocate(presfsave(k1))
    allocate(preshsave(k1))

    allocate(thvhsave(k1))

    allocate(u0avsave(k1))
    allocate(v0avsave(k1))
    allocate(thl0avsave(k1))
    allocate(qt0avsave(k1))
    allocate(sv0avsave(k1,nsv))

    umsave(:,:,:) = um(:,:,:) ! Copy
    vmsave(:,:,:) = vm(:,:,:)
    wmsave(:,:,:) = wm(:,:,:)
    e12msave(:,:,:) = e12m(:,:,:)
    thlmsave(:,:,:) = thlm(:,:,:)
    qtmsave(:,:,:) = qtm(:,:,:)

    u0save(:,:,:) = u0(:,:,:)
    v0save(:,:,:) = v0(:,:,:)
    w0save(:,:,:) = w0(:,:,:)
    thl0save(:,:,:) = thl0(:,:,:)
    qt0save(:,:,:) = qt0(:,:,:)
    ql0save(:,:,:) = ql0(:,:,:)
    ql0hsave(:,:,:) = ql0h(:,:,:)
    e120save(:,:,:) = e120(:,:,:)
    dthvdzsave(:,:,:) = dthvdz(:,:,:)
    ekmsave(:,:,:) = ekm(:,:,:)
    tmp0save(:,:,:) = tmp0(:,:,:)
    eslsave(:,:,:) = esl(:,:,:)
    qvslsave(:,:,:) = qvsl(:,:,:)
    qvsisave(:,:,:) = qvsi(:,:,:)

    thv0hsave(:,:,:) = thv0h(:,:,:)

    presfsave(:) = presf(:)
    preshsave(:) = presf(:)

    thvhsave(:) = thvh(:)

    u0avsave(:) = u0av(:)
    v0avsave(:) = v0av(:)
    thl0avsave(:) = thl0av(:)
    qt0avsave(:) = qt0av(:)

    ! Not using a loop here causes a segfault for some reason
    do n = 1, nsv
      sv0save(:,:,:,n) = sv0(:,:,:,n)
      svmsave(:,:,:,n) = svm(:,:,:,n)
      sv0avsave(:,n) = sv0av(:,n)
    end do

    allocate(fnudgeglob(1-ih:itot+ih,1-jh:jtot+jh,1:k1))
    allocate(fnudgeloc(2-ih:i1+ih,2-jh:j1+jh,k1))

    call calcfnudge
    
    if (lprecursor) then
      if (myid == 0 ) then
        write(*,*) 'nudgedepthgr = ', nudgedepthgr 
        write(*,*) 'Succesfully initialized modprecursor'
      end if
    end if
  end subroutine init_precursor

  subroutine calcfnudge
    use modglobal, only : pi, itot, jtot, ih, jh, k1, j1, i1, kmax
    use modfields, only : u0av, v0av
    use modmpi, only : myidx, myidy

    implicit none
    integer i,j,k
    real(field_r) fnudge

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

  subroutine swap_fields()
    use modfields, only : u0, v0, w0, e120, thl0, qt0, ql0, ql0h, tmp0, &
                          um, vm, wm, e12m, thlm, qtm, &
                          presf, presh, dthvdz, esl, qvsl, qvsi, thv0h, &
                          u0av, v0av, thl0av, qt0av, thvh
    use modsubgrid, only : ekm

    implicit none

    call swap(um, umsave)
    call swap(vm, vmsave)
    call swap(wm, wmsave)
    call swap(e12m, e12msave)
    call swap(thlm, thlmsave)
    call swap(qtm, qtmsave)
    call swap(u0, u0save)
    call swap(v0, v0save)
    call swap(w0, w0save)
    call swap(thl0, thl0save)
    call swap(qt0, qt0save)
    call swap(ql0, ql0save)
    call swap(ql0h, ql0hsave)
    call swap(e120, e120save)
    call swap(dthvdz, dthvdzsave)
    call swap(ekm, ekmsave)
    call swap(tmp0, tmp0save)
    call swap(esl, eslsave)
    call swap(qvsl, qvslsave)
    call swap(qvsi, qvsisave)
    call swap(thv0h, thv0hsave)
    call swap(presf, presfsave)
    call swap(presh, preshsave)
    call swap(thvh, thvhsave)
    call swap(u0av, u0avsave)
    call swap(v0av, v0avsave)
    call swap(thl0av, thl0avsave)
    call swap(qt0av, qt0avsave)
    call swap(sv0, sv0save)
    call swap(svm, svmsave)
    call swap(sv0av, sv0avsave)

  end subroutine swap_fields

  subroutine precursor_nudge_boundary 
    use modglobal, only : kmax, i1, j1, rdt, nsv
    use modfields, only : u0, v0, w0, thl0, e120, qt0, up, vp, wp, thlp, qtp, e12p

    implicit none
    integer i,j,k, s

    do k=1,kmax
      do j=2,j1
        do i=2,i1
          up(i,j,k) = (1-fnudgeloc(i,j,k))*up(i,j,k) + fnudgeloc(i,j,k)*(u0save(i,j,k)-u0(i,j,k))/rdt
          vp(i,j,k) = (1-fnudgeloc(i,j,k))*vp(i,j,k) + fnudgeloc(i,j,k)*(v0save(i,j,k)-v0(i,j,k))/rdt
          wp(i,j,k) = (1-fnudgeloc(i,j,k))*wp(i,j,k) + fnudgeloc(i,j,k)*(w0save(i,j,k)-w0(i,j,k))/rdt
          thlp(i,j,k) = (1-fnudgeloc(i,j,k))*thlp(i,j,k) + fnudgeloc(i,j,k)*(thl0save(i,j,k)-thl0(i,j,k))/rdt
          qtp(i,j,k) = (1-fnudgeloc(i,j,k))*qtp(i,j,k) + fnudgeloc(i,j,k)*(qt0save(i,j,k)-qt0(i,j,k))/rdt
          e12p(i,j,k) = (1-fnudgeloc(i,j,k))*e12p(i,j,k) + fnudgeloc(i,j,k)*(e120save(i,j,k)-e120(i,j,k))/rdt
        end do
      end do
    end do

    do s = 1, nsv
      do k = 1, kmax
        do j = 2, j1
          do i = 2, i1
            svp(i,j,k,s) = (1-fnudgeloc(i,j,k))*svp(i,j,k,s) + fnudgeloc(i,j,k)*(sv0save(i,j,k,s)-sv0(i,j,k,s))/rdt
          end do
        end do
      end do
    end do

  end subroutine precursor_nudge_boundary

  subroutine exit_precursor
    deallocate(fnudgeglob,fnudgeloc)
    deallocate(umsave, vmsave, wmsave, thlmsave, qtmsave, e12msave)
    deallocate(u0save, v0save, w0save, thl0save, qt0save, e120save)
    deallocate(ql0save, ql0hsave, dthvdzsave,ekmsave,tmp0save,eslsave,qvslsave,qvsisave,presfsave,preshsave, thv0hsave)
    deallocate(u0avsave,v0avsave,thl0avsave,qt0avsave, thvhsave)
    deallocate(svmsave, sv0save, sv0avsave)
  end subroutine exit_precursor

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

  subroutine swap_4d(a, b)

    real(field_r), pointer, intent(inout) :: a(:,:,:,:)
    real(field_r), pointer, intent(inout) :: b(:,:,:,:)

    real(field_r), pointer :: temp(:,:,:,:)

    temp => a
    a => b
    b => temp

  end subroutine swap_4d

  subroutine swap_3d(a, b)

    real(field_r), pointer, intent(inout) :: a(:,:,:)
    real(field_r), pointer, intent(inout) :: b(:,:,:)

    real(field_r), pointer :: temp(:,:,:)

    temp => a
    a => b
    b => temp

  end subroutine swap_3d

  subroutine swap_2d(a, b)

    real(field_r), pointer, intent(inout) :: a(:,:)
    real(field_r), pointer, intent(inout) :: b(:,:)

    real(field_r), pointer :: temp(:,:)

    temp => a
    a => b
    b => temp

  end subroutine swap_2d

  subroutine swap_1d(a, b)

    real(field_r), pointer, intent(inout) :: a(:)
    real(field_r), pointer, intent(inout) :: b(:)

    real(field_r), pointer :: temp(:)

    temp => a
    a => b
    b => temp

  end subroutine swap_1d

end module modprecursor