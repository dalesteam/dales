! This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1993-2026 The DALES team.
!

!> Microphysics abstraction layer.
!!  \author Hans Cuijpers, IMAU
!!  \author Thijs Heus,MPI-M
!!  \author Steef B\"oing, TU Delft
module modmicrophysics

  use modaerosol,        only: init_aerosol, laerosol
  use modglobal,         only: ifnamopt, checknamelisterror
  use moddrizzle,        only: drizzle
  use modbulkmicro,      only: initbulkmicro, &
                               exitbulkmicro, bulkmicro
  use modbulkmicro_data, only: l_lognormal, l_mur_cst, l_sb, mur_cst, sig_gr, &
                               l_sedc
  use modbulkmicro3,     only: initbulkmicro3, &
                               exitbulkmicro3, bulkmicro3
  use modmicrodata3,     only: l_sb_classic, l_sb_dumpall, l_sb_all_or, l_sb_dbg,   &
                               l_setclouds, l_setccn, l_corr_neg_qt, l_sb_lim_aggr, &
                               l_sb_stickyice, l_sb_conv_par, l_c_ccn,              &
                               l_sb_sat_max, l_sb_nuc_sat, l_sb_nuc_expl,           &
                               l_sb_nuc_diff, l_sb_inuc_sat, l_sb_inuc_expl,        &
                               l_sb_reisner, N_inuc, n_i_max, tmp_inuc, x_inuc,     &
                               N_inuc_R, c_inuc_R, a1_inuc_R, a2_inuc_R, c_ccn,     &
                               n_clmax, kappa_ccn, x_cnuc,sat_max, xc0_min, Nccn0,  &
                               l_statistics, l_tendencies, l_sb_tlimhetfreeze,      &
                               tlimhetfreeze
  use modmicrodata,      only: imicro, lstat, l_rain, Nc_0, sig_g
  use modsimpleice,      only: initsimpleice, &
                               exitsimpleice, simpleice
  use modsimpleice_data, only: l_berry, l_graupel, l_warm, l_mp, evapfactor, &
                               courantp
  use modmpi,            only: myid, D_MPI_BCAST, comm3d, print_info_stderr
  use modtimer,          only: timer_tic, timer_toc
  use moduser,           only: micro_user
  use modlogging,        only: finish
#ifdef USE_LCM
  use modlcm_adapter,    only: prepare_lcm, init_lcm, lcm_microphysics
#endif

  implicit none

  private

  public :: microphysics_read_namelist
  public :: initmicrophysics
  public :: initmicrophysics_state
  public :: microphysics
  public :: exitmicrophysics

  character(len=*), parameter :: modname = 'modmicrophysics'

  integer, parameter :: &
    imicro_none = 0,    & !< No microphysics.
    imicro_drizzle = 1, & !< Drizzle microphyics.
    imicro_bulk = 2,    & !< Double-moment warm microphysics.
    imicro_sice = 5,    & !< Single-moment mixed-phase microphysics.
    imicro_user = 10,   & !< User-provided microphysics.
    imicro_bulk3 = 11,  & !< Double-moment mixed-phase microphysics.
    imicro_lcm = 12       !< Lagrangian cloud microphysics.

contains

  !> Read microphysics namelist entry and broadcast settings.
  subroutine microphysics_read_namelist(nml_filename)
    use fortran_support, only: nnml_output

    character(len=*), intent(in) :: nml_filename

    character(len=*), parameter :: routine = modname//'/microphysics_read_namelist'

    integer :: ierr

    namelist /nammicrophysics/ &
      ! Common options
      imicro, lstat, l_rain, Nc_0, sig_g,                                       &
      ! Bulkmicro
      l_sb, l_sedc, l_mur_cst, l_lognormal, mur_cst, sig_gr,                    &
      ! Bulkmicro3
      l_sb_classic, l_sb_dumpall, l_sb_all_or, l_sb_dbg, l_setclouds, l_setccn, &
      l_corr_neg_qt, l_sb_lim_aggr, l_sb_stickyice, l_sb_conv_par, l_c_ccn,     &
      l_sb_sat_max, l_sb_nuc_sat, l_sb_nuc_expl, l_sb_nuc_diff, l_sb_inuc_sat,  &
      l_sb_inuc_expl, l_sb_reisner, N_inuc, n_i_max, tmp_inuc, x_inuc,          &
      N_inuc_R, c_inuc_R, a1_inuc_R, a2_inuc_R, c_ccn, n_clmax, kappa_ccn,      &
      x_cnuc,sat_max, xc0_min, Nccn0, l_statistics, l_tendencies,               &
      l_sb_tlimhetfreeze, tlimhetfreeze,                                        &
      ! Simpleice
      l_berry, l_graupel, l_warm, l_mp, evapfactor, courantp

    if (myid == 0) then
      open(ifnamopt, file=nml_filename, status='old', iostat=ierr)
      read(ifnamopt, nammicrophysics, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'nammicrophysics')
      write(nnml_output, nammicrophysics)
      close(ifnamopt)
    end if

    ! Common
    call D_MPI_BCAST(imicro, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(lstat, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(Nc_0, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(sig_g, 1, 0, comm3d, ierr)
    ! Bulkmicro
    call D_MPI_BCAST(l_sb, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_sedc, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_mur_cst, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_lognormal, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(mur_cst, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(sig_gr, 1, 0, comm3d, ierr)
    ! Bulkmicro3
    call D_MPI_BCAST(l_sb_classic, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_sb_dumpall, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_sb_all_or, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_sb_dbg, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_corr_neg_qt, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_sb_lim_aggr, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_sb_stickyice, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_sb_conv_par, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_c_ccn, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_sb_sat_max, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_sb_nuc_sat, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_sb_nuc_expl, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_sb_nuc_diff, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_sb_inuc_sat, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_sb_inuc_expl, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_sb_reisner, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(N_inuc_R, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(c_inuc_R, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(a1_inuc_R, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(a2_inuc_R, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(n_i_max, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(N_inuc, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(tmp_inuc, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(x_inuc, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(c_ccn, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(n_clmax, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(kappa_ccn, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(sat_max, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(x_cnuc, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(xc0_min, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(Nccn0, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_sb_tlimhetfreeze, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(tlimhetfreeze, 1, 0, comm3d, ierr)
    ! Simpleice
    call D_MPI_BCAST(l_berry, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_graupel, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_warm, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(l_mp, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(evapfactor, 1, 0, comm3d, ierr)
    call D_MPI_BCAST(courantp, 1, 0, comm3d, ierr)

    ! Perform some checks
    if (Nc_0 < 1e4) then
      call finish(routine, &
        'Nc_0 is suspiciously small (unit should be number per m3).')
    end if

  end subroutine microphysics_read_namelist

  !> Call the initialization routine of the selected microphysical scheme.
  subroutine initmicrophysics()
#ifndef USE_LCM
    character(len=*), parameter :: routine = modname//'/initmicrophysics'
#endif

    select case(imicro)
      case(imicro_bulk)
        call initbulkmicro
      case(imicro_sice)
        call initsimpleice
      case(imicro_bulk3)
        call initbulkmicro3
      case(imicro_lcm)
#ifdef USE_LCM
        call prepare_lcm
#else
        call finish(routine, &
          'LCM microphysics selected, but DALES was built without USE_LCM.')
#endif
    end select

    if(laerosol) call init_aerosol()

  end subroutine initmicrophysics

  !> Initialize selected microphysics after atmospheric fields are available.
  subroutine initmicrophysics_state()
#ifndef USE_LCM
    character(len=*), parameter :: routine = modname//'/initmicrophysics_state'
#endif

    select case(imicro)
      case(imicro_lcm)
#ifdef USE_LCM
        call init_lcm
#else
        call finish(routine, &
          'LCM microphysics selected, but DALES was built without USE_LCM.')
#endif
    end select

  end subroutine initmicrophysics_state

  !> Do the microphysics.
  subroutine microphysics

    character(len=*), parameter :: routine = modname//'/microphysics'

    call timer_tic(routine, 0)

    select case (imicro)
      case(imicro_drizzle)
        call drizzle
      case(imicro_bulk)
        call bulkmicro
      case(imicro_sice)
         call simpleice
      case(imicro_bulk3)
        call bulkmicro3
      case(imicro_user)
        call micro_user
      case(imicro_lcm)
#ifdef USE_LCM
        call lcm_microphysics
#else
        call finish(routine, &
          'LCM microphysics selected, but DALES was built without USE_LCM.')
#endif
    end select

    call timer_toc(routine)

  end subroutine microphysics

  !> Calls the clean-up routine for the selected microphysical scheme.
  subroutine exitmicrophysics

    select case (imicro)
      case(imicro_bulk)
        call exitbulkmicro
      case(imicro_sice)
        call exitsimpleice
      case(imicro_bulk3)
        call exitbulkmicro3
    end select

  end subroutine exitmicrophysics

end module modmicrophysics
