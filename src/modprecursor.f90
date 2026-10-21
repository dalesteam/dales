!> Concurrent precursor method.
!!
!! @author Pim van Dorp, TU Delft, 2015
!! @author Caspar Jungbacker, TU Delft, 2026
module modprecursor

  use fortran_support, only: nnml_output
  use modglobal,       only: longint, nsv, checknamelisterror, ifnamopt, pi, &
                             itot, i1, jtot, j1, kmax, k1, ih, jh, rdt
  use modfields,       only: u0, v0, w0, thl0, qt0, e120, ql0, ql0h, tmp0, um, &
                             vm, wm, e12m, thlm, qtm, presf, presh, dthvdz, &
                             esl, qvsl, qvsi, thv0h, u0av, v0av, thl0av, &
                             qt0av, thvh, sv0, svm, svp, sv0av, e12p, qtp, &
                             thlp, up, vp, wp
  use modmpi,          only: myid, comm3d, d_mpi_bcast, myidx, myidy, nprocx, &
                             nprocy
  use modsubgriddata,  only: ekm
  use modtimer,        only: timer_tic, timer_toc
  use modprecision,    only: field_r, real32, real64

  implicit none 

  private

  public :: lprecursor
  public :: statid
  public :: refid
  public :: turid
  public :: nsim
  public :: precursor_read_namelist
  public :: init_precursor
  public :: precursor_nudge_boundary
  public :: swap_fields
  public :: swap
  public :: exit_precursor

  character(len=*), parameter :: modname = 'modprecursor'

  !> Pointer swap
  interface swap
    module procedure :: swap_4d_r4
    module procedure :: swap_4d_r8
    module procedure :: swap_3d_r4
    module procedure :: swap_3d_r8
    module procedure :: swap_2d_r4
    module procedure :: swap_2d_r8
    module procedure :: swap_1d_r4
    module procedure :: swap_1d_r8
  end interface swap
  logical :: lprecursor = .false. !< Switch for enabling the precursor method
  logical :: lstatref = .false.   !< Output statistics of the reference simulation
  integer :: nudgedepthgr = 10    !< Depth of the nudging layer

  integer :: nsim = 1   !< Number of simulations
  integer :: statid = 1 !< ID of the simulation that should output statistics
  integer :: refid = 1  !< ID of the reference simulation (undisturbed)
  integer :: turid = 1  !< ID of the disturbed simulation
  
  real(field_r), allocatable :: umsave(:,:,:)
  real(field_r), allocatable :: vmsave(:,:,:)
  real(field_r), allocatable :: wmsave(:,:,:)
  real(field_r), allocatable :: thlmsave(:,:,:)
  real(field_r), allocatable :: e12msave(:,:,:)
  real(field_r), allocatable :: qtmsave(:,:,:)
  real(field_r), allocatable :: u0save(:,:,:)
  real(field_r), allocatable :: v0save(:,:,:)
  real(field_r), allocatable :: w0save(:,:,:)
  real(field_r), allocatable :: thl0save(:,:,:)
  real(field_r), allocatable :: qt0save(:,:,:)
  real(field_r), allocatable :: ql0save(:,:,:)
  real(field_r), allocatable :: ql0hsave(:,:,:)
  real(field_r), allocatable :: e120save(:,:,:)
  real(field_r), allocatable :: dthvdzsave(:,:,:)  
  real(field_r), allocatable :: ekmsave(:,:,:)  
  real(field_r), allocatable :: tmp0save(:,:,:)  
  real(field_r), allocatable :: eslsave(:,:,:)  
  real(field_r), allocatable :: qvslsave(:,:,:)  
  real(field_r), allocatable :: qvsisave(:,:,:)  
  real(field_r), allocatable :: thv0hsave(:,:,:)  
  real(field_r), allocatable :: presfsave(:)  
  real(field_r), allocatable :: preshsave(:)  
  real(field_r), allocatable :: thvhsave(:)  
  real(field_r), allocatable :: u0avsave(:)
  real(field_r), allocatable :: v0avsave(:)
  real(field_r), allocatable :: thl0avsave(:)
  real(field_r), allocatable :: qt0avsave(:)
  real(field_r), allocatable :: svmsave(:,:,:,:)
  real(field_r), allocatable :: sv0save(:,:,:,:)
  real(field_r), allocatable :: sv0avsave(:,:)

contains

  !> Compute the nudging factor.
  pure function nudge_fac(i, depth) result(fac)

    integer, intent(in) :: i     !< Index of grid point in the nudging layer (1 at the boundary, increasing towards the interior)
    integer, intent(in) :: depth !< Depth of the nudging layer

    real(field_r) :: fac !< Nudging factor

    fac = 0.5 + 0.5 * cos((pi / (depth - 1)) * (i - 1))

  end function nudge_fac

  !> Compute the updated tendency of a field after nudging at the boundary.
  pure function calc_nudged_tend(tend_in, nudge_target, field, dt, depth, i) &
    result(tend_out)

    real(field_r), intent(in) :: tend_in      !< Original tendency of the field.
    real(field_r), intent(in) :: nudge_target !< Target values to nudge towards.
    real(field_r), intent(in) :: field        !< Current values of the field.
    real(field_r), intent(in) :: dt           !< Time step size
    integer,       intent(in) :: depth        !< Depth of the nudging layer
    integer,       intent(in) :: i            !< Index of grid point in the nudging layer (1 at the boundary, increasing towards the interior)

    real(field_r) :: nudge_fac_
    real(field_r) :: tend_out

    nudge_fac_ = nudge_fac(i, depth)
    tend_out = (1 - nudge_fac_) * tend_in &
               + nudge_fac_ * (nudge_target - field) / dt

  end function calc_nudged_tend
  
  !> Read precursor namelist options.
  subroutine precursor_read_namelist(nml_filename)

    character(len=*), intent(in) :: nml_filename !<  Name of namelist file

    integer :: ierr

    namelist /precursor/ lprecursor, lstatref, nudgedepthgr

    if (myid == 0) then
      open(ifnamopt, file=nml_filename, status='old', iostat=ierr)
      read(ifnamopt,precursor, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'precursor')
      write(nnml_output, precursor)
      close(ifnamopt)
    end if

    call d_mpi_bcast(lprecursor, 1, 0, comm3d, ierr)
    call d_mpi_bcast(lstatref, 1, 0, comm3d, ierr)
    call d_mpi_bcast(nudgedepthgr, 1, 0, comm3d, ierr)

  end subroutine precursor_read_namelist

  !> Initialize the precursor method
  subroutine init_precursor

    integer :: i, j, k, s

    if (lprecursor) then

      Nsim = 2
      turid = 2
      if (.not. lstatref) statid = 2

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

      ! Copy intial values
      umsave(:,:,:) = um(:,:,:)
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
      do s = 1, nsv
        sv0save(:,:,:,s) = sv0(:,:,:,s)
        svmsave(:,:,:,s) = svm(:,:,:,s)
        sv0avsave(:,s) = sv0av(:,s)
      end do

      !$acc enter data copyin(umsave(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   vmsave(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   wmsave(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   e12msave(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   thlmsave(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   qtmsave(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   u0save(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   v0save(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   w0save(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   thl0save(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   qt0save(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   ql0save(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   ql0hsave(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   e120save(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   dthvdzsave(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   ekmsave(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   tmp0save(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   eslsave(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   qvslsave(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   qvsisave(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   thv0hsave(2-ih:i1+ih,2-jh:j1+jh,1:k1), &
      !$acc                   svmsave(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv), &
      !$acc                   sv0save(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv), &
      !$acc                   sv0avsave(1:k1,1:nsv), &
      !$acc                   presfsave(1:k1), preshsave(1:k1), thvhsave(1:k1), &
      !$acc                   u0avsave(1:k1), v0avsave(1:k1), thl0avsave(1:k1), &
      !$acc                   qt0avsave(1:k1))
!$omp target enter data map(to:umsave(2-ih:i1+ih,2-jh:j1+jh,1:k1),&
!$omp vmsave(2-ih:i1+ih,2-jh:j1+jh,1:k1),wmsave(2-ih:i1+ih,2-jh:j1+jh,&
!$omp 1:k1),e12msave(2-ih:i1+ih,2-jh:j1+jh,1:k1),thlmsave(2-ih:i1+ih,&
!$omp 2-jh:j1+jh,1:k1),qtmsave(2-ih:i1+ih,2-jh:j1+jh,1:k1),&
!$omp u0save(2-ih:i1+ih,2-jh:j1+jh,1:k1),v0save(2-ih:i1+ih,2-jh:j1+jh,&
!$omp 1:k1),w0save(2-ih:i1+ih,2-jh:j1+jh,1:k1),thl0save(2-ih:i1+ih,&
!$omp 2-jh:j1+jh,1:k1),qt0save(2-ih:i1+ih,2-jh:j1+jh,1:k1),&
!$omp ql0save(2-ih:i1+ih,2-jh:j1+jh,1:k1),ql0hsave(2-ih:i1+ih,&
!$omp 2-jh:j1+jh,1:k1),e120save(2-ih:i1+ih,2-jh:j1+jh,1:k1),&
!$omp dthvdzsave(2-ih:i1+ih,2-jh:j1+jh,1:k1),ekmsave(2-ih:i1+ih,&
!$omp 2-jh:j1+jh,1:k1),tmp0save(2-ih:i1+ih,2-jh:j1+jh,1:k1),&
!$omp eslsave(2-ih:i1+ih,2-jh:j1+jh,1:k1),qvslsave(2-ih:i1+ih,&
!$omp 2-jh:j1+jh,1:k1),qvsisave(2-ih:i1+ih,2-jh:j1+jh,1:k1),&
!$omp thv0hsave(2-ih:i1+ih,2-jh:j1+jh,1:k1),svmsave(2-ih:i1+ih,&
!$omp 2-jh:j1+jh,1:k1,1:nsv),sv0save(2-ih:i1+ih,2-jh:j1+jh,1:k1,1:nsv),&
!$omp sv0avsave(1:k1,1:nsv),presfsave(1:k1),preshsave(1:k1),&
!$omp thvhsave(1:k1),u0avsave(1:k1),v0avsave(1:k1),thl0av(1:k1),&
!$omp qt0avsave(1:k1))

    end if

  end subroutine init_precursor

  !> Swap the pointers of the current and saved fields.
  subroutine swap_fields()

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

  !> Nudge the fields at the lateral boundaries towards the saved values.
  subroutine precursor_nudge_boundary 

    character(len=*), parameter :: routine = modname//'/precursor_nudge_boundary'

    integer :: s

    call timer_tic(routine, 0)

    call nudge_field_at_boundary(u0, u0save, rdt, up)
    call nudge_field_at_boundary(v0, v0save, rdt, vp)
    call nudge_field_at_boundary(w0, w0save, rdt, wp)
    call nudge_field_at_boundary(thl0, thl0save, rdt, thlp)
    call nudge_field_at_boundary(qt0, qt0save, rdt, qtp)
    call nudge_field_at_boundary(e120, e120save, rdt, e12p)

    do s = 1, nsv
      call nudge_field_at_boundary(sv0(:,:,:,s), sv0save(:,:,:,s), rdt, &
                                   svp(:,:,:,s))
    end do

    !$acc wait

    call timer_toc(routine)

  end subroutine precursor_nudge_boundary

  !> Nudge a prognostic field at the lateral boundaries
  subroutine nudge_field_at_boundary(field, nudge_target, dt, tend)

    real(field_r), intent(in) :: field(2-ih:,2-jh:,:)        !< Field to nudge at the boundary.
    real(field_r), intent(in) :: nudge_target(2-ih:,2-jh:,:) !< Target values to nudge towards.
    real(field_r), intent(in) :: dt                          !< Time step size.

    real(field_r), intent(inout) :: tend(2-ih:,2-jh:,:) !< Tendency of field.

    character(len=*), parameter :: routine = modname//'/nudge_field_at_boundary'

    integer :: i, j, k

    call timer_tic(routine, 1)

    ! North
    if (myidy == 0) then
      !$acc parallel loop gang vector collapse(3) default(present) async
!!$omp target teams loop collapse(3) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
      do k = 1, kmax
        do j = 1, nudgedepthgr
          do i = 2, i1
            tend(i,j,k) = calc_nudged_tend(tend(i,j,k), nudge_target(i,j,k), &
                                           field(i,j,k), dt, nudgedepthgr, j)
          end do
        end do
      end do
    end if

    ! South
    if (myidy == nprocy - 1) then
      !$acc parallel loop gang vector collapse(3) default(present) async
!!$omp target teams loop collapse(3) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
      do k = 1, kmax
        do j = j1 - nudgedepthgr + 1, j1
          do i = 2, i1
            tend(i,j,k) = calc_nudged_tend(tend(i,j,k), nudge_target(i,j,k), &
                                           field(i,j,k), dt, nudgedepthgr, &
                                           j1 - j + 1)
          end do
        end do
      end do
    end if

    ! West
    if (myidx == 0) then
      !$acc parallel loop gang vector collapse(3) default(present) async
!!$omp target teams loop collapse(3) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
      do k = 1, kmax
        do j = 2, j1
          do i = 1, nudgedepthgr
            tend(i,j,k) = calc_nudged_tend(tend(i,j,k), nudge_target(i,j,k), &
                                           field(i,j,k), dt, nudgedepthgr, i)
          end do
        end do
      end do
    end if

    ! East
    if (myidx == nprocx - 1) then
      !$acc parallel loop gang vector collapse(3) default(present) async
!!$omp target teams loop collapse(3) defaultmap(present:aggregate)&
!!$omp defaultmap(present:allocatable)
      do k = 1, kmax
        do j = 2, j1
          do i = i1 - nudgedepthgr + 1, i1
            tend(i,j,k) = calc_nudged_tend(tend(i,j,k), nudge_target(i,j,k), &
                                           field(i,j,k), dt, nudgedepthgr, &
                                           i1 - i + 1)
          end do
        end do
      end do
    end if

    call timer_toc(routine)

  end subroutine nudge_field_at_boundary

  !> Deallocate the arrays used for the precursor method.
  subroutine exit_precursor

    !$acc exit data delete(umsave, vmsave, wmsave, thlmsave, qtmsave, &
    !$acc                  e12msave, u0save, v0save, w0save, thl0save, &
    !$acc                  qt0save, e120save, ql0save, ql0hsave, dthvdzsave, &
    !$acc                  ekmsave, tmp0save, eslsave, qvslsave, qvsisave, &
    !$acc                  presfsave, preshsave, thv0hsave, u0avsave, &
    !$acc                  v0avsave,thl0avsave,qt0avsave, thvhsave, &
    !$acc                  svmsave, sv0save, sv0avsave)
!$omp target exit data map(delete:umsave,vmsave,wmsave,thlmsave,&
!$omp qtmsave,e12msave,u0save,v0save,w0save,thl0save,qt0save,e120save,&
!$omp ql0save,ql0hsave,dthvdzsave,ekmsave,tmp0save,eslsave,qvslsave,&
!$omp qvsisave,presfsave,preshsave,thv0hsave,u0avsave,v0avsave,&
!$omp thl0avsave,qt0avsave,thvhsave,svmsave,sv0save,sv0avsave)

    deallocate(umsave, vmsave, wmsave, thlmsave, qtmsave, e12msave, &
               u0save, v0save, w0save, thl0save, qt0save, e120save, &
               ql0save, ql0hsave, dthvdzsave, ekmsave, tmp0save, eslsave, &
               qvslsave, qvsisave, presfsave, preshsave, thv0hsave, &
               u0avsave, v0avsave, thl0avsave, qt0avsave, thvhsave, &
               svmsave, sv0save, sv0avsave)

  end subroutine exit_precursor

  !> Swap the pointers of two 4D arrays.
  subroutine swap_4d_r4(a, b)

    real(real32), allocatable, intent(inout) :: a(:,:,:,:)
    real(real32), allocatable, intent(inout) :: b(:,:,:,:)

    real(real32), allocatable :: temp(:,:,:,:)
    
    call move_alloc(a, temp)
    call move_alloc(b, a)
    call move_alloc(temp, b)

  end subroutine swap_4d_r4
  !> Swap the pointers of two 4D arrays.
  subroutine swap_4d_r8(a, b)

    real(real64), allocatable, intent(inout) :: a(:,:,:,:)
    real(real64), allocatable, intent(inout) :: b(:,:,:,:)

    real(real64), allocatable :: temp(:,:,:,:)
    
    call move_alloc(a, temp)
    call move_alloc(b, a)
    call move_alloc(temp, b)

  end subroutine swap_4d_r8
  !> Swap the pointers of two 3D arrays.
  subroutine swap_3d_r4(a, b)

    real(real32), allocatable, intent(inout) :: a(:,:,:)
    real(real32), allocatable, intent(inout) :: b(:,:,:)

    real(real32), allocatable :: temp(:,:,:)

    call move_alloc(a, temp)
    call move_alloc(b, a)
    call move_alloc(temp, b)

  end subroutine swap_3d_r4
  !> Swap the pointers of two 3D arrays.
  subroutine swap_3d_r8(a, b)

    real(real64), allocatable, intent(inout) :: a(:,:,:)
    real(real64), allocatable, intent(inout) :: b(:,:,:)

    real(real64), allocatable :: temp(:,:,:)

    call move_alloc(a, temp)
    call move_alloc(b, a)
    call move_alloc(temp, b)

  end subroutine swap_3d_r8
  !> Swap the pointers of two 2D arrays.
  subroutine swap_2d_r4(a, b)

    real(real32), allocatable, intent(inout) :: a(:,:)
    real(real32), allocatable, intent(inout) :: b(:,:)

    real(real32), allocatable :: temp(:,:)

    call move_alloc(a, temp)
    call move_alloc(b, a)
    call move_alloc(temp, b)

  end subroutine swap_2d_r4
  !> Swap the pointers of two 2D arrays.
  subroutine swap_2d_r8(a, b)

    real(real64), allocatable, intent(inout) :: a(:,:)
    real(real64), allocatable, intent(inout) :: b(:,:)

    real(real64), allocatable :: temp(:,:)

    call move_alloc(a, temp)
    call move_alloc(b, a)
    call move_alloc(temp, b)

  end subroutine swap_2d_r8
  !> Swap the pointers of two 1D arrays.
  subroutine swap_1d_r4(a, b)

    real(real32), allocatable, intent(inout) :: a(:)
    real(real32), allocatable, intent(inout) :: b(:)

    real(real32), allocatable :: temp(:)

    call move_alloc(a, temp)
    call move_alloc(b, a)
    call move_alloc(temp, b)

  end subroutine swap_1d_r4
  !> Swap the pointers of two 1D arrays.
  subroutine swap_1d_r8(a, b)

    real(real64), allocatable, intent(inout) :: a(:)
    real(real64), allocatable, intent(inout) :: b(:)

    real(real64), allocatable :: temp(:)

    call move_alloc(a, temp)
    call move_alloc(b, a)
    call move_alloc(temp, b)

  end subroutine swap_1d_r8
end module modprecursor
