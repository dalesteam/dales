#define imax 500
#define jmax 500
#define kmax 130
subroutine advecu_2nd
  !$tuner initialize
  real(8), allocatable, dimension(:,:,:) :: a_out, u0, v0, w0 
  real(8), allocatable, dimension(:) :: rhobf
  real(8) :: dxiq, dyiq, dziq
  integer :: i, j, k

  allocate(a_out(imax, jmax, kmax))
  allocate(u0(imax, jmax, kmax))
  allocate(v0(imax, jmax, kmax))
  allocate(w0(imax, jmax, kmax))
  allocate(rhobf(kmax))

  call random_number(a_out)
  call random_number(u0)
  call random_number(v0)
  call random_number(w0)
  call random_number(rhobf)

  dxiq = 0.05
  dyiq = 0.05
  dziq = 0.05

  !$acc enter data copyin(a_out(:imax,:jmax,:kmax))
  !$acc enter data copyin(u0(:imax,:jmax,:kmax))
  !$acc enter data copyin(v0(:imax,:jmax,:kmax))
  !$acc enter data copyin(w0(:imax,:jmax,:kmax))
  !$tuner stop

  !$tuner start advecu_2nd_horz
  !$acc parallel loop collapse(ncollapse) vector_length(vlength) default(present)
  do k = 1, kmax
    do j = 2, jmax - 1
      do i = 2, imax - 1
        a_out(i,j,k)  = a_out(i,j,k)- ( &
                ( &
                (u0(i,j,k)+u0(i+1,j,k))*(u0(i,j,k)+u0(i+1,j,k)) &
                -(u0(i,j,k)+u0(i-1,j,k))*(u0(i,j,k)+u0(i-1,j,k)) &
                )*dxiq &
                +(  &
                (u0(i,j,k)+u0(i,j+1,k))*(v0(i,j+1,k)+v0(i-1,j+1 ,k)) &
                -(u0(i,j,k)+u0(i,j-1,k))*(v0(i,j  ,k)+v0(i-1,j  ,k)) &
                )*dyiq )
      end do
    end do
  end do
  !$tuner stop

  !$tuner start advecu_2nd_vert_bot
  !$acc parallel loop collapse(ncollapse) vector_length(vlength) default(present)
  do j = 2, jmax - 1
    do i = 2, imax - 1
      a_out(i,j,1) = a_out(i,j,1)-(1/rhobf(1))*( &
                     ( rhobf(2) * a_in(i,j,2) + rhobf(1) * a_in(i,j,1))*( w0(i,j,2)+ w0(i-1,j,2) ) &
                     ) *dziq
    end do
  end do
  !$tuner stop

  !$tuner start advecu_2nd_vert
  !$acc parallel loop collapse(ncollapse) vector_length(vlength) default(present)
  do k = 2, kmax
    do j = 2, jmax - 1
      do i = 2, imax - 1
        a_out(i,j,k) = a_out(i,j,k)- (1/rhobf(k))*( &
                       (rhobf(k) * a_in(i,j,k) + rhobf(k+1) * a_in(i,j,k+1) )*(w0(i,j,k+1)+w0(i-1,j,k+1)) &
                       -(rhobf(k) * a_in(i,j,k) + rhobf(k-1) * a_in(i,j,k-1) )*(w0(i,j,k )+w0(i-1,j,k )) &
                           )*dziq
      end do
    end do
  end do
  !$tuner stop

  !$tuner deinitialize

  !$acc exit data delete(a_out(:imax,:jmax,:kmax))
  !$acc exit data delete(u0(:imax,:jmax,:kmax))
  !$acc exit data delete(v0(:imax,:jmax,:kmax))
  !$acc exit data delete(w0(:imax,:jmax,:kmax))
  !$tuner stop
end subroutine advecu_2nd

