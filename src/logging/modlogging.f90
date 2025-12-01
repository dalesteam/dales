
module modlogging
    use iso_fortran_env,   only: error_unit
    use modmpi, only: myid, D_MPI_ALLREDUCE, comm3d, MPI_SUM
    use mpi, only: mpi_abort
    use fortran_support, only: util_abort, open_nml_output, close_nml_output, find_next_free_unit,fg_green, fg_default
    use modloghelpers, only: init_logger, set_msg_timestamp, fs_message=>message, fs_finish=>finish, fs_warning=>warning
    integer :: profile_output

    save

    logical :: is_initializing

    contains

    subroutine abort_with_mpi
      use modmpi, only: mpicomm=>MPI_COMM_WORLD, mpierr
      use mpi
      implicit none
      integer, parameter :: return_code = 1

      call MPI_ABORT(mpicomm%MPI_VAL, return_code, mpierr)
      
    end subroutine abort_with_mpi

    subroutine initlogging
      implicit none
      logical :: write_output = .false.
      if (myid == 0) write_output = .true.

      call init_logger(proc_id=myid, l_write_output=write_output, nerr_unit=error_unit, callback_abort=abort_with_mpi)
      call open_nml_output("namelist_used.nml")
      call open_profile_input("profiles_used.txt")

    end subroutine initlogging

    subroutine exitlogging
      implicit none
      call close_nml_output
      call close_profile_output
    end subroutine exitlogging

    subroutine handle_exit
      implicit none
    end subroutine handle_exit

  !>
  !! close the ASCII output that contains all the profile information
  subroutine open_profile_input(file)
    implicit none
    character(len=*), intent(in) :: file
    integer :: istat

    profile_output = find_next_free_unit(10, 20)

    open (profile_output, FILE=TRIM(file), IOSTAT=istat)

    if (istat /= 0) THEN
      call finish('open_nml_output', 'Could not open '//TRIM(file))
    end if

    end subroutine open_profile_input

  !>
  !!  close the ASCII output that contains all the profile information
  !!
  subroutine close_profile_output
    implicit none
    integer :: istat

    close (profile_output, IOSTAT=istat)

    if (istat /= 0) THEN
      call finish('close_profile_output', 'Could not close the profile output')
    end if

  end subroutine close_profile_output


  subroutine finish(name, text1, text2, text3, text4, text5, text6, text7, text8, text9, text10)
    use modmpi, only: myid
    implicit none
    character(len=*), intent(in) :: name
    class(*), intent(in), optional :: text1, text2, text3, text4, text5
    class(*), intent(in), optional :: text6, text7, text8, text9, text10
    logical :: write_backtrace
    write_backtrace = (.not.is_initializing)
    call fs_finish(name=name, text=convertstring(text1, text2, text3, text4, text5, text6, text7, text8, text9, text10), write_backtrace=write_backtrace)
  end subroutine finish

  subroutine warning(name, text1, text2, text3, text4, text5, text6, text7, text8, text9, text10)
    implicit none
    character(len=*), intent(in) :: name
    class(*), intent(in), optional :: text1, text2, text3, text4, text5
    class(*), intent(in), optional :: text6, text7, text8, text9, text10

    call fs_warning(name, text=convertstring(text1, text2, text3, text4, text5, text6, text7, text8, text9, text10))
  end subroutine warning

  subroutine message(name, text1, text2, text3, text4, text5, text6, text7, text8, text9, text10, all_print)
    implicit none
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: all_print
    class(*), intent(in), optional :: text1, text2, text3, text4, text5
    class(*), intent(in), optional :: text6, text7, text8, text9, text10

    call fs_message(name, text=convertstring(text1, text2, text3, text4, text5, text6, text7, text8, text9, text10), all_print=all_print)
  end subroutine message

  subroutine enable_init_error_logging
    implicit none
    is_initializing = .true.
  end subroutine enable_init_error_logging

  subroutine disable_init_error_logging
    implicit none
    is_initializing = .false.
  end subroutine disable_init_error_logging
 
 
  ! based on public domain source code from https://fortranwiki.org/fortran/show/tostring
  function convertstring(text1, text2, text3, text4, text5, text6, text7, text8, text9, text10)
    implicit none
    class(*), intent(in), optional :: text1, text2, text3, text4, text5
    class(*), intent(in), optional :: text6, text7, text8, text9, text10
    character(len=:), allocatable :: convertstring
    character(len=4096)        :: line
    character(len=4096)        :: curpart
    integer                    :: istart
    integer i


    istart=1
    if (present(text1)) call print_part(text1)
    if (present(text2)) call print_part(text2)
    if (present(text3)) call print_part(text3)
    if (present(text4)) call print_part(text4)
    if (present(text5)) call print_part(text5)
    if (present(text6)) call print_part(text6)
    if (present(text7)) call print_part(text7)
    if (present(text8)) call print_part(text8)
    if (present(text9)) call print_part(text9)
    if (present(text10)) call print_part(text10)
    convertstring=trim(line)

    contains

    subroutine print_part(text)
      use, intrinsic :: iso_fortran_env, only : int8, int16, int32, int64, real32, real64, real128

      class(*), intent(in), optional :: text
          select type(text)
          type is (integer(kind=int8))
            write(line(istart:),'(i0)') text
          type is (integer(kind=int16))
            write(line(istart:),'(i0)') text
          type is (integer(kind=int32))
            write(line(istart:),'(i0)') text
          type is (integer(kind=int64))
            write(line(istart:),'(i0)') text
          type is (real(kind=real32))
            write(line(istart:),'(1pg0)') text
          type is (real(kind=real64))
            write(line(istart:),'(1pg0)') text
          type is (real(kind=real128))
            write(line(istart:),'(1pg0)') text
          type is (logical)
            write(line(istart:),'(1l)') text
          type is (character(len=*))
            write(line(istart:),'(a)') text
            istart = istart + len(text)
            return
          type is (complex)
            write(line(istart:),'("(",1pg0,",",1pg0,")")') text
          end select
          istart = len_trim(line) + 1
    end subroutine print_part

  end function convertstring

end module modlogging
