
module modlogging
    use iso_fortran_env,   only: error_unit
    use modmpi, only: mpierr, myid, commwrld, d_mpi_bcast
    use mpi, only: mpi_abort, MPI_FINALIZE
    use fortran_support, only: init_logger, open_nml_output, close_nml_output, nnml_output, find_next_free_unit, & 
                               fg_green, fg_default, fs_message=>message, fs_finish=>finish, fs_warning=>warning, filename_max


    character(len=*), parameter :: modname = 'modlogging'

    integer :: profile_output !< the unit for the ascii output of profiles/fields/tracers.

    character(filename_max) :: profile_output_file = "used_profiles.txt"    !< path of the profile output file
    character(filename_max) :: namelist_output_file = "used_namelist.txt"   !< path of the namelist output file


    contains

    subroutine abort_with_mpi
      implicit none
      integer, parameter :: return_code = 1

      call MPI_ABORT(commwrld%MPI_VAL, return_code, mpierr)
      
    end subroutine abort_with_mpi

    subroutine initlogging
      implicit none
      logical :: write_output = .false.
      if (myid == 0) write_output = .true.

      call init_logger(proc_id=myid, l_write_output=write_output, nerr_unit=error_unit, callback_abort=abort_with_mpi)

      if (myid == 0) then
        call open_nml_output(namelist_output_file)
      endif
      call open_profile_input(profile_output_file)

    end subroutine initlogging

    subroutine exitlogging
      implicit none
      if (myid == 0) then
        call close_nml_output
      end if
      call close_profile_output
    end subroutine exitlogging

  !>
  !! opens the ASCII output file that contains all the profile information
  subroutine open_profile_input(file)
    implicit none

    character(len=*), parameter :: routine = modname//'/open_profile_input'

    character(len=*), intent(in) :: file
    integer :: istat

    profile_output = find_next_free_unit(10, 20)

    open (profile_output, FILE=TRIM(file), IOSTAT=istat)

    if (istat /= 0) THEN
      call finish(routine, 'Could not open '//TRIM(file))
    end if

    end subroutine open_profile_input

  !>
  !!  close the ASCII output file that contains all the profile information
  subroutine close_profile_output
    implicit none

    character(len=*), parameter :: routine = modname//'/close_profile_input'

    integer :: istat

    close (profile_output, IOSTAT=istat)

    if (istat /= 0) THEN
      call finish(routine, 'Could not close the profile output')
    end if

  end subroutine close_profile_output

  !>
  !!  wrapper around the fortran-support finish function that enables writing also numbers without having to format and define an extra character array
  subroutine finish(name, text1, text2, text3, text4, text5, text6, text7, text8, text9, text10)
    implicit none
    character(len=*), intent(in) :: name !< the name of the routine which caused an error
    class(*), intent(in), optional :: text1, text2, text3, text4, text5
    class(*), intent(in), optional :: text6, text7, text8, text9, text10
  
    call fs_finish(name=name, text=convertstring(text1, text2, text3, text4, text5, text6, text7, text8, text9, text10))

  end subroutine finish

  !>
  !!  wrapper around the fortran-support warning function that enables writing also numbers without having to format and define an extra character array
  subroutine warning(name, text1, text2, text3, text4, text5, text6, text7, text8, text9, text10)
    implicit none
    character(len=*), intent(in) :: name !< the name of the routine which sends this warning
    class(*), intent(in), optional :: text1, text2, text3, text4, text5
    class(*), intent(in), optional :: text6, text7, text8, text9, text10

    call fs_warning(name, text=convertstring(text1, text2, text3, text4, text5, text6, text7, text8, text9, text10))

  end subroutine warning

  !>
  ! wrapper around the fortran-support message function that enables writing also numbers without having to format and define an extra character array
  subroutine message(name, text1, text2, text3, text4, text5, text6, text7, text8, text9, text10, all_print)
    implicit none
    character(len=*), intent(in) :: name !< the name of the routine which sends this message
    logical, intent(in), optional :: all_print
    class(*), intent(in), optional :: text1, text2, text3, text4, text5
    class(*), intent(in), optional :: text6, text7, text8, text9, text10

    call fs_message(name, text=convertstring(text1, text2, text3, text4, text5, text6, text7, text8, text9, text10), all_print=all_print)

  end subroutine message

  !> converts up to 10 optional arguments and returns a line with all of them concatenated
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
      use, intrinsic :: iso_fortran_env, only : int8, int16, int32, int64, real32, real64

      class(*), intent(in), optional :: text
          select type(text)
          type is (integer(kind=int8))
            write(line(istart:),'(i0) ') text
          type is (integer(kind=int16))
            write(line(istart:),'(i0) ') text
          type is (integer(kind=int32))
            write(line(istart:),'(i0) ') text
          type is (integer(kind=int64))
            write(line(istart:),'(i0) ') text
          type is (real(kind=real32))
            write(line(istart:),'(1pg0) ') text
          type is (real(kind=real64))
            write(line(istart:),'(1pg0) ') text
          type is (logical)
            write(line(istart:),'(1l) ') text
          type is (character(len=*))
            write(line(istart:),'(a) ') text
            istart = istart + len(text)
            return
          type is (complex)
            write(line(istart:),'("(",1pg0,",",1pg0,")") ') text
          end select
          istart = len_trim(line) + 1
    end subroutine print_part

  end function convertstring

end module modlogging
