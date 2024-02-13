!> \file modmsebudg.f90
!! Compute moist static energy budget
!! as required by RCE-MIP experiments
!!
!! The statistics is done in two steps, msebudg1 is called before tstep_integrate,
!! msebudg2 is called after tstep_integrate.
!!
!! \author Stephan de Roode
!! \author Fredrik Jansson

module modmsebudg
  use modglobal, only : longint, output_prefix
  use modprecision, only : longint, field_r

implicit none
PUBLIC :: initmsebudg, msebudg1, msebudg2, exitmsebudg
save

integer :: Nmse_2D = 4
real(field_r), allocatable :: mse0(:,:,:)      !<   rcemip, moist static energy
real(field_r), allocatable :: msem(:,:,:)
real(field_r), allocatable :: field_mse_2D(:,:,:)

logical :: lmsebudg = .false. !< switch to enable the MSE budget (on/off)
real    :: dtav
integer(kind=longint) :: idtav,tnext

!NetCDF variables
integer,parameter :: nmsevar = 4
integer :: ncid2,ncid3,nrec2 = 0, nrec3 = 0
character(80) :: fname = 'mse_budget.xxx.xxx.xxx.nc'
character(80),dimension(nmsevar,4) :: ncmsename
character(80),dimension(1,4) :: tmsencname

contains

  subroutine initmsebudg
    use modmpi,   only :myid,comm3d,myidx,myidy,d_mpi_bcast
    use modglobal,only :i1,ih,j1,jh,k1,imax,jmax,cexpnr,ifnamopt,fname_options,dtmax,dtav_glob,ladaptive,dt_lim,btime,tres
    use modstat_nc,only : lnetcdf,open_nc, define_nc,ncinfo,writestat_dims_nc
    implicit none
    integer :: ierr

     namelist/NAMMSEBUDG/ &
    dtav,lmsebudg

    dtav=dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMMSEBUDG,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMMSEBUDG'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMMSEBUDG'
      endif
      write(6 ,NAMMSEBUDG)
      close(ifnamopt)
    end if
    call D_MPI_BCAST(dtav      ,1, 0, comm3d,ierr)
    call D_MPI_BCAST(lmsebudg  ,1, 0, comm3d,ierr)
    idtav = int(dtav/tres,longint)
    tnext      = idtav   +btime

    if(.not.(lmsebudg)) return

    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be an integer multiple of dtmax'
    end if

    allocate(field_mse_2D(2-ih:i1+ih,2-jh:j1+jh,Nmse_2D ))
    allocate(mse0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(msem(2-ih:i1+ih,2-jh:j1+jh,k1))

         write(fname,'(A,i3.3,A,i3.3,A)') 'mse_budget.', myidx, '.', myidy, '.xxx.nc'    !rce table 6 2D hourly averaged variables
      fname(20:22) = cexpnr
      call ncinfo(tmsencname(1,:),'time','Time','s','time')
      call ncinfo(ncmsename( 1,:),'fmse'     ,'mass-weighted vertical integral of frozen moist static energy','J/m2','tt0t')
      call ncinfo(ncmsename( 2,:),'hadvmse'  ,'mass-weighted vertical integral of horizontal advective tendency of frozen moist static energy','J/m2/','tt0t')
      call ncinfo(ncmsename( 3,:),'vadvmse'  ,'mass-weighted vertical integral of vertical advective tendency of frozen moist static energy','J/m2/s','tt0t')
      call ncinfo(ncmsename( 4,:),'tnfmse'   ,'total tendency of mass-weighted vertical integral of frozen moist static energy','J/m2/s','tt0t')
!      call ncinfo(ncmsename( 5,:),'tnfmsevar','total tendency of spatial variance of mass-weighted vertical integral of frozen moist static energy','J2/m4/s','tt0t')

      call open_nc(trim(output_prefix)//fname, ncid3,nrec3,n1=imax,n2=jmax)
      if (nrec3==0) then
        call define_nc( ncid3, 1, tmsencname)
        call writestat_dims_nc(ncid3)
      end if
     call define_nc( ncid3, nmsevar, ncmsename)

   end subroutine initmsebudg

  ! call one of the advection functions to advect a scalar
  subroutine advect_scalar(a_in, a_out, iadv)
    use modglobal, only: i1,ih,j1,jh,k1,iadv_cd2,iadv_5th,iadv_52,     &
                         iadv_cd6,iadv_62,iadv_kappa,   &
                         iadv_hybrid,iadv_hybrid_f

    use advec_2nd,      only : advecc_2nd
    use advec_52,       only : advecc_52
    use advec_5th,      only : advecc_5th
    use advec_62,       only : advecc_62
    use advec_6th,      only : advecc_6th
    use advec_hybrid,   only : advecc_hybrid
    use advec_hybrid_f, only : advecc_hybrid_f
    use advec_kappa,    only : advecc_kappa
    use advec_upw,      only : advecc_upw
    implicit none

    real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the cell centered field
    real(field_r), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency
    integer :: iadv

    select case(iadv)
      case(iadv_cd2)
        call advecc_2nd(a_in,a_out)
      case(iadv_5th)
        call advecc_5th(a_in,a_out)
      case(iadv_52)
        call advecc_52(a_in,a_out)
      case(iadv_cd6)
        call advecc_6th(a_in,a_out)
      case(iadv_62)
        call advecc_62(a_in,a_out)
      case(iadv_kappa)
        call advecc_kappa(a_in,a_out)
      case(iadv_hybrid)
        call advecc_hybrid(a_in,a_out)
      case(iadv_hybrid_f)
        call advecc_hybrid_f(a_in,a_out)
      case default
          stop "Unknown advection scheme "
    end select
  end subroutine advect_scalar

! Step 1 of MSE budget. Called before tstep_integrate.
    subroutine msebudg1
      use modglobal, only : ih,i1,jh,j1,k1,dzf,iadv_thl, &
                            rk3step,timee,cp,grav,rlv,zf
      use modfields, only : rhobf,tmp0,qt0,ql0,u0,v0,w0
      use modmpi,    only : excjs

      integer k

      real(field_r), allocatable :: msep(:,:,:), &
                                    u0_save(:,:,:), &    !container for u
                                    v0_save(:,:,:), &    !container for v
                                    w0_save(:,:,:)       !container for w

       if (.not. lmsebudg) return
       if(timee<tnext) return

       if (rk3step==1) then
          ! start of first substep.
          ! qt0 = qtm, etc
          ! save msem - value at start of this full time step

          field_mse_2D = 0 ! reset all

          do k=1,k1
             msem(:,:,k)  = cp * tmp0(:,:,k) + grav * zf(k) + rlv * (qt0(:,:,k) - ql0(:,:,k))
             field_mse_2D (2:i1,2:j1,4) = field_mse_2D (2:i1,2:j1,4) + rhobf(k) * dzf(k) * msem (2:i1,2:j1,k)     ! m-field for total tendency
          end do
          !write (*,*) '*** MSE1 step:', rk3step, 'time:', timee
       end if

       if (rk3step==3) then
          ! start of 3rd substep
          ! tendencies at this point are used to do the final step
          allocate( &
               msep       (2-ih:i1+ih,2-jh:j1+jh,k1), &
               u0_save    (2-ih:i1+ih,2-jh:j1+jh,k1), &
               v0_save    (2-ih:i1+ih,2-jh:j1+jh,k1), &
               w0_save    (2-ih:i1+ih,2-jh:j1+jh,k1)  )

          do k=1,k1
             ! don't need to calculate for ghost cells here - probably their tmp,ql is not right
             mse0(:,:,k)  = cp * tmp0(:,:,k) + grav * zf(k) + rlv * (qt0(:,:,k) - ql0(:,:,k))
          end do

          call excjs(mse0, 2,i1,2,j1,1,k1,ih,jh)   ! get mse halo for advection

          u0_save = u0
          v0_save = v0
          w0_save = w0

          w0 = 0
          msep = 0
          !advect uv
          call advect_scalar(mse0, msep, iadv_thl)

          do k=1,k1
             field_mse_2D (2:i1,2:j1,2) = field_mse_2D (2:i1,2:j1,2) + rhobf(k) * dzf(k) * msep(2:i1,2:j1,k) ! horizontal adv tend
          enddo

          w0 = w0_save
          u0 = 0
          v0 = 0
          msep = 0
          ! advect w
          call advect_scalar(mse0, msep, iadv_thl)
          do k=1,k1
             field_mse_2D (2:i1,2:j1,3) = field_mse_2D (2:i1,2:j1,3) + rhobf(k) * dzf(k) * msep(2:i1,2:j1,k)  ! vertical adv tend
          enddo

          u0 = u0_save
          v0 = v0_save

          deallocate (msep,u0_save,v0_save,w0_save)
          !write (*,*) '*** MSE1 step:', rk3step, 'time:', timee
       end if


     end subroutine msebudg1


! Step 2 of MSE budget. Called after tstep_integrate.
     subroutine msebudg2
       use modglobal, only : i1,j1,k1,imax,jmax,rdt,dzf,&
                             rk3step,timee,rtimee,dt_lim,cp,grav,rlv,zf
       use modfields, only : rhobf,tmp0,qt0,ql0
       use modstat_nc, only : lnetcdf, writestat_nc
       implicit none

       integer k
       real, allocatable :: vars(:,:,:)

       if (.not. lmsebudg) return
       if (rk3step/=3) return

       if(timee<tnext) then
          dt_lim = min(dt_lim,tnext-timee)
          return
       end if

       tnext = tnext+idtav
       dt_lim = min(dt_lim, tnext-timee)

       do k=1,k1
          mse0(:,:,k)  = cp * tmp0(:,:,k) + grav * zf(k) + rlv * (qt0(:,:,k) - ql0(:,:,k))
          field_mse_2D (2:i1,2:j1,1) = field_mse_2D (2:i1,2:j1,1) + rhobf(k) * dzf(k) * mse0(2:i1,2:j1,k)     ! current value from 0-field
       end do

       field_mse_2D (2:i1,2:j1,4) =  (field_mse_2D (2:i1,2:j1,1) - field_mse_2D (2:i1,2:j1,4)) / rdt  ! total tendency = (mse0-msem)/dt


       if (lnetcdf) then
          allocate(vars(imax,jmax,Nmse_2D))

          vars(1:imax,1:jmax,:) = field_mse_2D (2:i1,2:j1,:)

          call writestat_nc(ncid3,1,tmsencname,(/rtimee/),nrec3,.true.)
          call writestat_nc(ncid3,nmsevar,ncmsename,vars,nrec3,imax,jmax)
          deallocate(vars)
       end if

       !write (*,*) '*** MSE2 step:', rk3step, 'time:', timee, 'rdt:', rdt
     end subroutine msebudg2


  subroutine exitmsebudg
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none

    if(.not.(lmsebudg)) return

    deallocate(field_mse_2D)
    deallocate(mse0,msem)

    if(lmsebudg .and. lnetcdf) call exitstat_nc(ncid3)
  end subroutine exitmsebudg


end module modmsebudg
