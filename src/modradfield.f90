!  This file is part of DALES.
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
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!
!> Dumps 2D fields of several (radiation related) variables.
module modradfield

  use fortran_support,   only: finish, nnml_output
  use modfields,         only: rhof, qt0, ql0, tmp0, u0, v0, presf
  use modglobal,         only: itot, jtot, i1, j1, kmax, dzf, ifnamopt, &
                               dtav_glob, timeav_glob, cu, cv, tup, tdn, cp, &
                               rlv, checknamelisterror
  use modmpi,            only: D_MPI_BCAST, commwrld, myid
  use modnetcdf_file_t,  only: cross_section_file_t
  use modprecision,      only: field_r
  use modraddata,        only: lwd, lwu, swd, swu, lwdca, lwuca, swdca, swuca, &
                               swdir, swdif, sw_up_toa, sw_dn_toa, lw_up_toa, &
                               sw_up_ca_toa, lw_up_ca_toa
  use modstat_nc_files,  only: add_output_file, is_sampling_timestep
  use modsurfdata,       only: qtflux, thlflux
  use modthermodynamics, only: calc_qsat
  use modtimer,          only: timer_tic, timer_toc

  implicit none

  character(len=*), parameter :: modname = 'modradfield'

  public :: radfield_read_namelist
  public :: initradfield
  public :: radfield

  type(cross_section_file_t) :: ofile
  integer                    :: ofile_id

  real    :: dtav
  real    :: timeav
  real    :: nsamples
  logical :: lradfield

contains 

  !> Read radfield namelist.
  subroutine radfield_read_namelist(nml_filename)

    character(len=*), intent(in) :: nml_filename

    integer :: ierr

    namelist /NAMRADFIELD/ dtav, timeav, lradfield

    dtav = dtav_glob
    timeav = timeav_glob

    if (myid == 0) then
      open(ifnamopt, file=nml_filename, status='old', iostat=ierr)
      read(ifnamopt, NAMRADFIELD, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMRADFIELD')
      write(nnml_output, NAMRADFIELD)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(dtav       ,1,0,commwrld,ierr)
    call D_MPI_BCAST(timeav     ,1,0,commwrld,ierr)
    call D_MPI_BCAST(lradfield  ,1,0,commwrld,ierr)

  end subroutine radfield_read_namelist

  !> Open radfield NetCDF file and add variables.
  subroutine initradfield

    character(len=*), parameter :: routine = modname//'/initradfield'

    if (lradfield) then
      call timer_tic(routine, 1)

      nsamples = timeav / dtav

      ofile = cross_section_file_t('radfield', nx=itot, ny=jtot)
      call add_output_file(ofile, dtav, ofile_id, dt_write=timeav)

      call ofile%add_var('hfls','surface upward latent heat flux','W/m2','tt0t')
      call ofile%add_var('hfss','surface upward sensible heat flux','W/m2','tt0t')
      call ofile%add_var('rlds','surface downwellling longwave flux','W/m2','tt0t')
      call ofile%add_var('rlus','surface upwelling longwave flux','W/m2','tt0t')
      call ofile%add_var('rsds','surface downwellling shortwave flux','W/m2','tt0t')
      call ofile%add_var('rsus','surface upwellling shortwave flux','W/m2','tt0t')
      call ofile%add_var('rsdtm','TOM incoming shortwave flux','W/m2','tt0t')
      call ofile%add_var('rldtm','TOM incoming longwave flux','W/m2','tt0t')
      call ofile%add_var('rsutm','TOM outgoing shortwave flux','W/m2','tt0t')
      call ofile%add_var('rlutm','TOM outgoing longwave flux','W/m2','tt0t')
      call ofile%add_var('rsdscs','surface downwelling shortwave flux - clear sky','W/m2','tt0t')
      call ofile%add_var('rsuscs','surface upwelling shortwave flux - clear sky','W/m2','tt0t')
      call ofile%add_var('rldscs','surface downwelling longwave flux - clear sky','W/m2','tt0t')
      call ofile%add_var('rluscs','surface upwelling longwave flux - clear sky','W/m2','tt0t')
      call ofile%add_var('rsutmcs','TOM outgoing shortwave flux - clear sky','W/m2','tt0t')
      call ofile%add_var('rlutmcs','TOM outgoing longwave flux - clear sky','W/m2','tt0t')
      call ofile%add_var('rsds_dir','surface downwellling shortwave direct flux','W/m2','tt0t')
      call ofile%add_var('rsds_dif','surface downwellling shortwave diffuse flux','W/m2','tt0t')

      call ofile%add_var('prw','water vapor path','kg/m2','tt0t')
      call ofile%add_var('clwvi','condensed water path','kg/m2','tt0t')
      call ofile%add_var('clivi','ice water path','kg/m2','tt0t')
      call ofile%add_var('spwr','saturated water vapor path','kg/m2','tt0t')
      call ofile%add_var('uabot','eastward wind at lowest model level','m/s','mt0t')
      call ofile%add_var('vabot','northward wind at lowest model level','m/s','tm0t')
      call ofile%add_var('tabot','air temperature at lowest model level','K','tt0t')

      call ofile%add_var('rsdt','TOA incoming shortwave flux','W/m2','tt0t')
      call ofile%add_var('rsut','TOA outgoing shortwave flux','W/m2','tt0t')
      call ofile%add_var('rlut','TOA outgoing longwave flux','W/m2','tt0t')
      call ofile%add_var('rsutcs','TOA outgoing shortwave flux - clear sky','W/m2','tt0t')
      call ofile%add_var('rlutcs','TOA outgoing longwave flux - clear sky','W/m2','tt0t')

      call timer_toc(routine)
    end if

  end subroutine initradfield

  !> Sample radiation fields.
  subroutine radfield

    character(len=*), parameter :: routine = modname//'/radfield'

    integer       :: i, j, k
    real(field_r) :: ilratio

    real(field_r), pointer :: hfls(:,:)
    real(field_r), pointer :: hfss(:,:)
    real(field_r), pointer :: rlds(:,:)
    real(field_r), pointer :: rlus(:,:)
    real(field_r), pointer :: rsds(:,:)
    real(field_r), pointer :: rsus(:,:)
    real(field_r), pointer :: rsdtm(:,:)
    real(field_r), pointer :: rldtm(:,:)
    real(field_r), pointer :: rsutm(:,:)
    real(field_r), pointer :: rlutm(:,:)
    real(field_r), pointer :: rsdscs(:,:)
    real(field_r), pointer :: rsuscs(:,:)
    real(field_r), pointer :: rldscs(:,:)
    real(field_r), pointer :: rluscs(:,:)
    real(field_r), pointer :: rsutmcs(:,:)
    real(field_r), pointer :: rlutmcs(:,:)
    real(field_r), pointer :: rsds_dir(:,:)
    real(field_r), pointer :: rsds_dif(:,:)
    real(field_r), pointer :: prw(:,:)
    real(field_r), pointer :: clwvi(:,:)
    real(field_r), pointer :: clivi(:,:)
    real(field_r), pointer :: spwr(:,:)
    real(field_r), pointer :: uabot(:,:)
    real(field_r), pointer :: vabot(:,:)
    real(field_r), pointer :: tabot(:,:)
    real(field_r), pointer :: rsdt(:,:)
    real(field_r), pointer :: rsut(:,:)
    real(field_r), pointer :: rlut(:,:)
    real(field_r), pointer :: rsutcs(:,:)
    real(field_r), pointer :: rlutcs(:,:)

    if (lradfield) then
      if (is_sampling_timestep(ofile_id)) then
        call timer_tic(routine, 1) 

        call ofile%get_pointer('hfls', hfls)
        call ofile%get_pointer('hfss', hfss)
        call ofile%get_pointer('rlds', rlds)
        call ofile%get_pointer('rlus', rlus)
        call ofile%get_pointer('rsds', rsds)
        call ofile%get_pointer('rsus', rsus)
        call ofile%get_pointer('rsdtm', rsdtm)
        call ofile%get_pointer('rldtm', rldtm)
        call ofile%get_pointer('rsutm', rsutm)
        call ofile%get_pointer('rlutm', rlutm)
        call ofile%get_pointer('rsdscs', rsdscs)
        call ofile%get_pointer('rsuscs', rsuscs)
        call ofile%get_pointer('rldscs', rldscs)
        call ofile%get_pointer('rluscs', rluscs)
        call ofile%get_pointer('rsutmcs', rsutmcs)
        call ofile%get_pointer('rlutmcs', rlutmcs)
        call ofile%get_pointer('rsds_dir', rsds_dir)
        call ofile%get_pointer('rsds_dif', rsds_dif)
        call ofile%get_pointer('prw', prw)
        call ofile%get_pointer('clwvi', clwvi)
        call ofile%get_pointer('clivi', clivi)
        call ofile%get_pointer('spwr', spwr)
        call ofile%get_pointer('uabot', uabot)
        call ofile%get_pointer('vabot', vabot)
        call ofile%get_pointer('tabot', tabot)
        call ofile%get_pointer('rsdt', rsdt)
        call ofile%get_pointer('rsut', rsut)
        call ofile%get_pointer('rlut', rlut)
        call ofile%get_pointer('rsutcs', rsutcs)
        call ofile%get_pointer('rlutcs', rlutcs)

        do j = 2, j1
          do i = 2, i1
            hfls(i,j) = hfls(i,j) + rhof(1) * rlv * qtflux(i,j)
            hfss(i,j) = hfss(i,j) + rhof(1) * cp * thlflux(i,j)
          end do
        end do

        rlds(:,:) = rlds(:,:) + abs(lwd(2:i1,2:j1,1)) 
        rlus(:,:) = rlus(:,:) + abs(lwu(2:i1,2:j1,1))
        rsds(:,:) = rsds(:,:) + abs(swd(2:i1,2:j1,1))
        rsus(:,:) = rsus(:,:) + abs(swu(2:i1,2:j1,1))
        rsdtm(:,:) = rsdtm(:,:) + abs(swd(2:i1,2:j1,kmax))
        rldtm(:,:) = rldtm(:,:) + abs(lwd(2:i1,2:j1,kmax))
        rsutm(:,:) = rsutm(:,:) + abs(swu(2:i1,2:j1,kmax))
        rlutm(:,:) = rlutm(:,:) + abs(lwu(2:i1,2:j1,kmax))
        rsdscs(:,:) = rsdscs(:,:) + abs(swdca(2:i1,2:j1,1))
        rsuscs(:,:) = rsuscs(:,:) + abs(swuca(2:i1,2:j1,1))
        rldscs(:,:) = rldscs(:,:) + abs(lwdca(2:i1,2:j1,1))
        rluscs(:,:) = rluscs(:,:) + abs(lwuca(2:i1,2:j1,1))
        rsutmcs(:,:) = rsutmcs(:,:) + abs(swuca(2:i1,2:j1,kmax))
        rlutmcs(:,:) = rlutmcs(:,:) + abs(lwuca(2:i1,2:j1,kmax))
        rsds_dir(:,:) = rsds_dir(:,:) + abs(swdir(2:i1,2:j1,1))
        rsds_dif(:,:) = rsds_dif(:,:) + abs(swdif(2:i1,2:j1,1))

        do k = 1, kmax
          do j = 2, j1
            do i = 2, i1
              prw(i,j) = prw(i,j) + rhof(k) * (qt0(i,j,k) - ql0(i,j,k)) * dzf(k)
              clwvi(i,j) = clwvi(i,j) + rhof(k) *  ql0(i,j,k) * dzf(k)
              ilratio=max(0._field_r,min(1._field_r,(tmp0(i,j,k)-tdn)/(tup-tdn)))! cloud water vs cloud ice partitioning
              clivi(i,j) = clivi(i,j) + rhof(k) *  ql0(i,j,k) * dzf(k) *(1-ilratio)
              spwr(i,j) = spwr(i,j) + rhof(k) *  calc_qsat(tmp0(i,j,k), presf(k)) * dzf(k)
            end do
          end do
        end do

        uabot(2:i1,2:j1) = uabot(2:i1,2:j1) + u0(2:i1,2:j1,1) + cu
        vabot(2:i1,2:j1) = vabot(2:i1,2:j1) + v0(2:i1,2:j1,1) + cv
        tabot(2:i1,2:j1) = tabot(2:i1,2:j1) + tmp0(2:i1,2:j1,1)

        rsdt(:,:) = rsdt(:,:) + abs(SW_dn_TOA(2:i1,2:j1))
        rsut(:,:) = rsut(:,:) + abs(SW_up_TOA(2:i1,2:j1))
        rlut(:,:) = rlut(:,:) + abs(LW_up_TOA(2:i1,2:j1))
        rsutcs(:,:) = rsutcs(:,:) + abs(SW_up_ca_TOA(2:i1,2:j1))
        rlutcs(:,:) = rlutcs(:,:) + abs(LW_up_ca_TOA(2:i1,2:j1))

        call timer_toc(routine)
      end if
    end if

  end subroutine radfield

end module modradfield