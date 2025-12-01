! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2025, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE modloghelpers

  ! Serial and parallel logger

  ! LOG-LEVEL                 EXAMPLE CALL
  !-------------------------------------------------------------------------
  !
  !                           CALL debug_on()
  ! debug:                    CALL message('graupel','20K')
  !                           CALL debug_off()
  !
  ! info:                     CALL message('graupel','20K')
  !
  ! info(print all procs):    CALL message('graupel','20K',all_print=.TRUE.)
  !
  ! info(print to own unit):  CALL message_to_own_unit('graupel','20K',99)
  !
  ! param:                    CALL print_value('graupel',variable)
  !
  ! warning:                  CALL warning('graupel','Wrong value')
  !
  ! error:                    CALL finish('graupel','To high qv')
  !
  !-------------------------------------------------------------------------
  !
  ! Logger is able to handle a callback:
  !
  !       abort: routine to end the code in case finish is called

  USE mo_iconlib_kind, ONLY: wp, dp, sp, i8
  USE mo_io_units, ONLY: filename_max
  USE mo_util_backtrace, ONLY: ftn_util_backtrace

  IMPLICIT NONE
  PRIVATE

  ! init
  PUBLIC :: init_logger
  ! getter/setters
  PUBLIC :: debug_on, debug_off, set_msg_timestamp, enable_logging
  ! logging
  PUBLIC :: message, warning, finish, message_to_own_unit
!   PUBLIC :: message, warning, finish, print_value, message_to_own_unit

  ! character to write message for logger, only used outside of mo_exception
  PUBLIC :: message_text

  ! output unit
  PRIVATE :: nerr
  ! helper routine for all levels
  PRIVATE :: write_to_unit

  ! logical/functions to define behaviour of logs
  PRIVATE :: lwrite, ldebug, lprint, lmsg_timestamp
  PRIVATE :: lvl_info_is_active, lvl_debug_is_active

  ! helpers for prefix
  PRIVATE :: extra_prefix
  PRIVATE :: debug_prefix, info_prefix_default, info_prefix_extra, param_prefix, error_prefix

!   INTERFACE print_value !< report on a parameter value
!     MODULE PROCEDURE print_lvalue !< logical
!     MODULE PROCEDURE print_ivalue !< integer
!     MODULE PROCEDURE print_i8value !< integer(i8)
!     MODULE PROCEDURE print_rvalue_dp !< real64
!     MODULE PROCEDURE print_rvalue_sp !< real32
!   END INTERFACE

  INTERFACE
    SUBROUTINE callback_function()
      IMPLICIT NONE
    END SUBROUTINE callback_function
  END INTERFACE

  INTERFACE
    FUNCTION prefix_function() RESULT(pre)
      IMPLICIT NONE
      CHARACTER(len=:), ALLOCATABLE :: pre
      CHARACTER(LEN=100) :: tmp
    END FUNCTION prefix_function
  END INTERFACE

  INTEGER :: id, nerr

  ! highest order flag -> passed at init
  ! defines if a PE writes output
  ! for lvl_info/lvl_param
  LOGICAL :: lwrite = .FALSE.

  ! flag to switch on/off log for lvl_info/lvl_param
  ! for all PE's
  LOGICAL :: ldebug = .FALSE.

  ! flag to switch on/off log for lvl_info/lvl_param
  ! only once
  LOGICAL :: lprint = .FALSE.

  !> Flag. If .TRUE., precede output messages by time stamp.
  !  This parameters is set by mo_run_nml
  LOGICAL :: lmsg_timestamp = .FALSE.

  ! character to write message for logger, only used outside of mo_exception
  CHARACTER(LEN=filename_max) :: message_text

  CHARACTER(LEN=filename_max) :: extra_prefix

  PROCEDURE(callback_function), POINTER :: model_abort => NULL()
  PROCEDURE(prefix_function), POINTER :: info_prefix => NULL()

  INTEGER, PARAMETER :: lvl_debug = 0 !< debug message
  INTEGER, PARAMETER :: lvl_info = 1 !< informational message
  INTEGER, PARAMETER :: lvl_param = 2 !< report parameter
  INTEGER, PARAMETER :: lvl_warn = 3 !< warning message
  INTEGER, PARAMETER :: lvl_error = 4 !< error message

CONTAINS

  SUBROUTINE init_logger(proc_id, l_write_output, nerr_unit, callback_abort, &
                         l_extra_output, extra_info_prefix)

    INTEGER, INTENT(IN) :: proc_id
    INTEGER, INTENT(IN) :: nerr_unit
    LOGICAL, INTENT(IN) :: l_write_output
    PROCEDURE(callback_function) :: callback_abort
    LOGICAL, OPTIONAL, INTENT(IN) :: l_extra_output
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: extra_info_prefix

    info_prefix => info_prefix_default

    id = proc_id
    nerr = nerr_unit

    IF (PRESENT(l_extra_output)) THEN
      lwrite = l_write_output .OR. l_extra_output
    ELSE
      lwrite = l_write_output
    END IF

    model_abort => callback_abort

    IF (PRESENT(l_extra_output)) THEN
      IF (l_extra_output) THEN
        ! use user-defined prefix
        IF (PRESENT(extra_info_prefix)) THEN
          extra_prefix = extra_info_prefix
          info_prefix => info_prefix_extra
        END IF
      END IF
    END IF

  END SUBROUTINE init_logger

  SUBROUTINE set_msg_timestamp(ltimestamp)
    LOGICAL, INTENT(IN) :: ltimestamp
    lmsg_timestamp = ltimestamp
  END SUBROUTINE set_msg_timestamp

  SUBROUTINE enable_logging(l_write_output)
    LOGICAL, INTENT(IN) :: l_write_output
    lwrite = l_write_output
  END SUBROUTINE enable_logging

  SUBROUTINE debug_on
    ldebug = .TRUE.
  END SUBROUTINE debug_on

  SUBROUTINE debug_off
    ldebug = .FALSE.
  END SUBROUTINE debug_off

  LOGICAL FUNCTION lvl_info_is_active()
    lvl_info_is_active = lwrite .OR. ldebug .OR. lprint
  END FUNCTION lvl_info_is_active

  LOGICAL FUNCTION lvl_debug_is_active()
    lvl_debug_is_active = ldebug
  END FUNCTION lvl_debug_is_active

  FUNCTION debug_prefix() RESULT(pre)
    CHARACTER(len=:), ALLOCATABLE :: pre
    CHARACTER(LEN=100) :: tmp
    WRITE (tmp, '(4x,a,i5)') 'DEBUG PE: ', id
    pre = TRIM(tmp)
  END FUNCTION debug_prefix

  FUNCTION info_prefix_default() RESULT(pre)
    CHARACTER(len=:), ALLOCATABLE :: pre
    CHARACTER(LEN=100) :: tmp
    WRITE (tmp, '(4x,a)') ''
    pre = TRIM(tmp)
  END FUNCTION info_prefix_default

  FUNCTION info_prefix_extra() RESULT(pre)
    CHARACTER(len=:), ALLOCATABLE :: pre
    CHARACTER(LEN=100) :: tmp
    WRITE (tmp, '(4x,a,a,i5)') TRIM(extra_prefix), ' PE: ', id
    pre = TRIM(tmp)
  END FUNCTION info_prefix_extra

  FUNCTION param_prefix() RESULT(pre)
    CHARACTER(len=:), ALLOCATABLE :: pre
    CHARACTER(LEN=100) :: tmp
    WRITE (tmp, '(1x,a)') '---     '
    pre = TRIM(tmp)
  END FUNCTION param_prefix

  FUNCTION warning_prefix() RESULT(pre)
    CHARACTER(len=:), ALLOCATABLE :: pre
    CHARACTER(LEN=100) :: tmp
    WRITE (tmp, '(4x,a,i5)') 'WARNING PE: ', id
    pre = TRIM(tmp)
  END FUNCTION warning_prefix

  FUNCTION error_prefix() RESULT(pre)
    CHARACTER(len=:), ALLOCATABLE :: pre
    CHARACTER(LEN=100) :: tmp
    WRITE (tmp, '(4x,a,i5)') 'FINISH PE: ', id
    pre = TRIM(tmp)
  END FUNCTION error_prefix

  SUBROUTINE finish(name, text, write_backtrace)
    implicit none
    CHARACTER(len=*), INTENT(IN)           :: name
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: text
    CHARACTER(LEN=filename_max) :: tmp
    logical, intent(in) :: write_backtrace

    tmp(:) = ' '
    IF (PRESENT(text)) THEN
      tmp = text
    END IF

    
    CALL write_to_unit(name, tmp, lvl_error)

    if (lwrite.or.write_backtrace) then
        CALL ftn_util_backtrace()
    end if

    ! callback passed at init
    CALL model_abort()

  END SUBROUTINE finish

  SUBROUTINE message(name, text, all_print)

    CHARACTER(len=*), INTENT(IN) :: name
    CHARACTER(len=*), INTENT(IN) :: text
    LOGICAL, INTENT(IN), OPTIONAL :: all_print

    IF (PRESENT(all_print)) THEN
      lprint = all_print
    ELSE
      lprint = .FALSE.
    END IF

    IF (lvl_info_is_active() .AND. lvl_debug_is_active()) THEN
      CALL write_to_unit(name, text, lvl_debug)
    ELSE IF (lvl_info_is_active()) THEN
      CALL write_to_unit(name, text, lvl_info)
    END IF

  END SUBROUTINE message

  SUBROUTINE message_to_own_unit(name, text, nerr_unit)

    CHARACTER(len=*), INTENT(IN) :: name
    CHARACTER(len=*), INTENT(IN) :: text
    INTEGER, INTENT(IN) :: nerr_unit

    IF (lvl_info_is_active() .AND. lvl_debug_is_active()) THEN
      CALL write_to_unit(name, text, lvl_debug, nerr_unit)
    ELSE IF (lvl_info_is_active()) THEN
      CALL write_to_unit(name, text, lvl_info, nerr_unit)
    END IF

  END SUBROUTINE message_to_own_unit

  SUBROUTINE param(name, text)

    CHARACTER(len=*), INTENT(IN) :: name
    CHARACTER(len=*), INTENT(IN) :: text

    IF (lvl_info_is_active()) THEN
      CALL write_to_unit(name, text, lvl_param)
    END IF

  END SUBROUTINE param

  SUBROUTINE write_to_unit(name, text, level, nerr_unit)
    CHARACTER(len=*), INTENT(IN) :: name
    CHARACTER(len=*), INTENT(IN) :: text
    INTEGER, INTENT(IN) :: level
    INTEGER, INTENT(IN), OPTIONAL :: nerr_unit

    CHARACTER(len=filename_max) :: prefix
    CHARACTER(len=10) :: ctime
    CHARACTER(len=8)  :: cdate
    INTEGER :: unit_for_logging

    IF (PRESENT(nerr_unit)) THEN
      unit_for_logging = nerr_unit
    ELSE
      unit_for_logging = nerr
    END IF

    prefix(:) = ' ' ! init clean
    SELECT CASE (level)
    CASE (lvl_debug); prefix = debug_prefix()
    CASE (lvl_info); prefix = info_prefix()
    CASE (lvl_param); prefix = param_prefix()
    CASE (lvl_warn); prefix = warning_prefix()
    CASE (lvl_error); prefix = error_prefix()
    END SELECT

    IF (lmsg_timestamp) THEN
      CALL DATE_AND_TIME(date=cdate, time=ctime)
      prefix = '['//cdate//' '//ctime//'] '//TRIM(prefix)
    END IF

    IF (LEN_TRIM(name) > 0) THEN
      WRITE (unit_for_logging, '(a)') TRIM(prefix)//' '//TRIM(name)//': '//TRIM(text)
    ELSE
      WRITE (unit_for_logging, '(a)') TRIM(prefix)//' '//TRIM(text)
    END IF

  END SUBROUTINE write_to_unit

  SUBROUTINE warning(name, text)

    CHARACTER(len=*), INTENT(IN) :: name
    CHARACTER(len=*), INTENT(IN) :: text

    CALL write_to_unit(name, text, lvl_warn)

  END SUBROUTINE warning

END MODULE modloghelpers
