!###############################################################################
!
! WARNING: do not confuse this file with go.F90 from LOTOS-EUROS!
! This is an interface between le_drydepos_gas_depac.f90 and the rest of DALES.
! The DEPAC routines use 'gol', 'goErr' and 'goPr' that were defined in the GO
! routines for LOTOS-EUROS. This file contains an implementation of these
! functions. This way, le_drydepos_depac.f90 did not have to be adapted to DALES
! and it is directly transferrable between LE and DALES.
!
! L. Geers, June 2022
!###############################################################################

module GO
    implicit none

    public :: gol, goPr, goErr

    ! buffer size:
    integer, parameter :: lengol = 1024

    ! buffer for standard output
    character(len=lengol) :: gol

contains
    !> Print ordinary message to stdout
    !!
    !! Subroutine prints `gol` buffer to stdout.
    subroutine goPr
        ! write buffer to standard output (stdout)
        write(*, '(a)') trim(gol)
        ! clear buffer
        gol = ''
    end subroutine goPr

    !> Print error message.
    !!
    !! Currently, errors are also sent to stdout. Need to adapt this to stderr.
    subroutine goErr
        call goPr
    end subroutine goErr

    ! *****************************************************


    subroutine goSplitString_s( line, n, values, status, sep )

      ! --- in/out --------------------------------

      character(len=*), intent(in)            ::  line
      integer, intent(out)                    ::  n
      character(len=*), intent(out)           ::  values(:)
      integer, intent(out)                    ::  status
      character(len=1), intent(in), optional  ::  sep

      ! --- const ----------------------------

      ! character(len=*), parameter   ::  rname = mname//'/goSplitString_s'

      ! --- local ---------------------------------

      integer                     ::  l
      character(len=1)            ::  the_sep
      character(len=len(line))    ::  line_curr
      character(len=len(line))    ::  val

      ! --- begin ---------------------------------

      ! length:
      l = len_trim(line)
      ! any?
      if ( l > 0 ) then

      ! seperation character:
      the_sep = ' '
      if ( present(sep) ) the_sep = sep

      ! copy input:
      line_curr = line

      ! no parts extracted yet:
      n = 0

      ! loop until all elements in line_curr are processed:
      do
          ! empty ? then finished:
          if ( len_trim(line_curr) == 0 ) exit
          ! next number:
          n = n + 1
          ! storage problem ?
          if ( n > size(values) ) then
          write (gol,'("output array is too small:")'); call goErr
          write (gol,'("  input line   : ",a )') trim(line); call goErr
          write (gol,'("  size(values) : ",i4)') size(values); call goErr
          ! TRACEBACK; status=1; return
          end if
          ! extract leading name:
          call ReadFromLine_s( line_curr, val, status, sep=the_sep )
          ! IF_NOTOK_RETURN(status=1)
          ! store value in output:
          values(n) = adjustl(val)
      end do

      ! adhoc: might end with seperation character ..
      if ( line(l:l) == the_sep ) then
          n = n + 1
          values(n) = ''
      end if

      else
      ! empty:
      n = 0
      end if ! l > 0

      ! ok
      status = 0

    end subroutine goSplitString_s

    ! ***


    subroutine ReadFromLine_s( s, ss, status, sep )

      ! --- in/out --------------------------

      character(len=*), intent(inout)         ::  s
      character(len=*), intent(out)           ::  ss
      integer, intent(out)                    ::  status

      character(len=1), intent(in), optional  ::  sep

      ! --- const ------------------------------

      ! character(len=*), parameter  ::  rname = mname//'/ReadFromLine_s'

      ! --- local ----------------------------

      character(len=len(s))     ::  s1, s2
      character(len=1)          ::  thesep
      integer                   ::  l, ll

      ! --- begin ----------------------------

      ! default seperation character provided as argument:
      thesep = ','
      if (present(sep)) thesep = sep

      ! split at seperation character:
      call goSplitLine( s, s1, thesep, s2, status )
      ! IF_ERROR_RETURN(status=1)

      ! check storage:
      l = len_trim(s1)
      ll = len(ss)
      if ( ll < l ) then
        write (gol,'("size of output string not sufficient:")'); call goErr
        write (gol,'("  first part of input : ",a )') trim(s1) ; call goErr
        write (gol,'("  output length       : ",i4)') ll       ; call goErr
        ! TRACEBACK; status=1; return
      end if
      ! store:
      ss = trim(s1)

      ! return remainder:
      s = s2

      ! ok
      status = 0

    end subroutine ReadFromLine_s


     !**********************************************************************


    subroutine goSplitLine( line, s1, c, s2, status )

      ! --- in/out ----------------------------

      character(len=*), intent(in)      ::  line
      character(len=*), intent(out)     ::  s1
      character(len=1), intent(in)      ::  c
      character(len=*), intent(out)     ::  s2
      integer, intent(out)              ::  status

      ! --- local -----------------------------

      integer                     ::  l, pos
      character(len=len(line))    ::  s

      ! --- begin -----------------------------

      s = line
      l = len_trim(s)

      pos = scan(s,c)
      if ( (pos<1) .or. (pos>l) ) then
        ! s='abcd'  -> s1='abcd', s2=''
        !call AdjustLeft( s1, s(1:l) )
        s1 = AdjustL( s(1:l) )
        s2 = ''
      else if (pos==1) then
        ! s=',' or s=',abcd'  ->  s1='', s2='' or 'abcd'
        s1 = ''
        if (l==1) then
          ! s=','
          s2 = ''
        else
          !call AdjustLeft( s2, s(pos+1:l) )
          s2 = AdjustL( s(pos+1:l) )
        end if
      else
        ! s='ab,' or s='ab,cd'
        !call AdjustLeft( s1, s(1:pos-1) )
        s1 = AdjustL( s(1:pos-1) )
        if (pos==l) then
          ! s='ab,'
          s2 = ''
        else
          ! s='ab,cd'
          !call AdjustLeft( s2, s(pos+1:l) )
          s2 = AdjustL( s(pos+1:l) )
        end if
      end if

      ! ok
      status = 0

    end subroutine goSplitLine

end module ! GO
