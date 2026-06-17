!> \file modlsmcrosssection.f90
!!   Dumps an instantenous lsmcrosssection of the field
!! lsmcrosssections in the yz-plane and in the xy-plane            |
!        of u,v,w,thl,thv,qt,ql. Written to movv_*.expnr and movh_*.expnr
!! If netcdf is true, this module leads the cross.myid.expnr.nc output

!!  \par Revision list
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
module modlsmcrosssection

  use modglobal,          only : longint
  use modlogging,         only : finish, warning
  use modnetcdf_file_t,   only : cross_section_file_t
  use modprecision,       only : field_r
  use modstat_nc_files,   only : add_output_file
  use modsurfdata,        only : ksoilmax

  implicit none

  character(len=*), parameter :: modname = 'modlsmcrosssection'

  private

  public :: initlsmcrosssection, lsmcrosssection

  save

  real    :: dtav
  integer(kind=longint) :: idtav, tnext
  logical :: lcross = .false.     !< switch for doing the lsmcrosssection (on/off)
  logical :: lcrosssoil = .false. !< switch for doing vertical soil crosssection (on/off)
  integer :: crossplane = 2       !< Location of the xz lsmcrosssection
  integer :: crossplane_local = -1 !< Local j-index of xz crosssection on this rank
  integer :: crossheight = 2      !< Height of the xy lsmcrosssection
  character(4) :: cheight

  type(cross_section_file_t) :: soil_xz_file
  type(cross_section_file_t) :: soil_xy_file
  type(cross_section_file_t) :: surf_file

  integer :: soil_xz_file_id = 0
  integer :: soil_xy_file_id = 0
  integer :: surf_file_id = 0

  logical :: soil_xz_enabled = .false.
  logical :: soil_xy_enabled = .false.
  logical :: surf_enabled = .false.
  logical :: lenable_gradients = .true.
  logical :: lenable_temp_more = .true.

contains

  !> Initializing lsmcrosssection. Read out the namelist, initializing the variables
  subroutine initlsmcrosssection
    use modmpi,     only : myid, myidy, mpierr, comm3d, D_MPI_BCAST
    use modglobal,  only : ifnamopt, fname_options, dtmax, dtav_glob, ladaptive, &
      j1, jmax, dy, y0, dt_lim, tres, btime, checknamelisterror, itot, jtot
    use modstat_nc, only : lnetcdf
    use modsurfdata, only : isurf
    use modlsm,     only : lags
    use fortran_support, only : nnml_output

    implicit none

    character(len=*), parameter :: routine = modname//'/initlsmcrossection'

    integer :: ierr, crossplane_global
    real(field_r) :: loc
    character(len=4) :: cloc

    namelist/NAMLSMCROSSSECTION/ &
    lcross, lcrosssoil, dtav, crossheight, crossplane, &
    lenable_gradients, lenable_temp_more

    dtav = dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMLSMCROSSSECTION,iostat=ierr)
       call checknamelisterror(ierr, ifnamopt, 'NAMLSMCROSSSECTION')
      write(nnml_output ,NAMLSMCROSSSECTION)
      close(ifnamopt)
    end if

    if (lcross .and. .not. (isurf == 1 .or. isurf == 2 .or. isurf == 11)) then
       lcross = .FALSE.
       call warning(routine, "Ignoring lcross, lsmcrossection currently implemented only for isurf==1, 2, or 11.")
    endif

    if (lcrosssoil .and. .not. (isurf == 1 .or. isurf == 11)) then
       lcrosssoil = .FALSE.
       call warning (routine, "Ignoring lcrosssoil, lsm soil crossection currently implemented only for isurf==1 or 11.")
    endif

    if (lcross .and. .not. lnetcdf) then
       lcross = .FALSE.
       call warning (routine, "Ignoring lcross, lcross output implemented only for netcdf output.")
    endif

     if (lenable_temp_more .and. .not. (isurf == 1 .or. isurf == 11)) then
       lenable_temp_more = .FALSE.
       call warning(routine, "Ignoring lenable_temp_more, extra LSM diagnostics are only available for isurf==1 or 11.")
     endif

    call D_MPI_BCAST(dtav       ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lcross     ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lcrosssoil ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(crossheight,1,0,comm3d,mpierr)
    call D_MPI_BCAST(crossplane ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lenable_gradients,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lenable_temp_more,1,0,comm3d,mpierr)

    idtav = int(dtav / tres, kind=kind(idtav))
    tnext   = idtav+btime
    if(.not.(lcross .or. lcrosssoil)) return
    dt_lim = min(dt_lim,tnext)

    if (lcrosssoil) then
      crossplane_global = crossplane
      if (crossheight < 1 .or. crossheight > ksoilmax .or. crossplane_global < 1 .or. crossplane_global > jtot + 1) then
        call finish(routine, 'lsmcrosssection: lsmcrosssection out of range')
      end if
      crossplane_local = crossplane_global - myidy * jmax + 1
      if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
        call finish(routine, 'lsmcrosssection: dtav should be a integer multiple of dtmax')
      end if

      if (.not. lnetcdf) return

      if (crossplane_local >= 2 .and. crossplane_local <= j1) then
        write(cloc, '(i4.4)') crossplane_global
        loc = y0 + dy * (crossplane_global - 1) + 0.5_field_r * dy
        soil_xz_file = cross_section_file_t('lsmcrossxz.'//cloc, nx=itot, nzs=ksoilmax, loc=loc, lgpu=.false.)
        call soil_xz_file%add_var('tsoil', 'xz crosssection of the Soil temperature', 'K', 't0tts')
        call soil_xz_file%add_var('phiw', 'xz crosssection of the Soil moisture', 'm3/m3', 't0tts')
        call add_output_file(soil_xz_file, dtav, soil_xz_file_id)
        soil_xz_enabled = .true.
      end if

      write(cheight, '(i4.4)') crossheight
      soil_xy_file = cross_section_file_t('lsmcrossxy.'//cheight, nx=itot, ny=jtot, lgpu=.false.)
      call soil_xy_file%add_var('tsoil', 'xy crosssection of the Soil temperature', 'K', 'tt0t')
      call soil_xy_file%add_var('phiw', 'xy crosssection of the Soil moisture', 'm3/m3', 'tt0t')
      call add_output_file(soil_xy_file, dtav, soil_xy_file_id)
      soil_xy_enabled = .true.
    end if

    if (lcross.and.lnetcdf) then
      surf_file = cross_section_file_t('surfcross', nx=itot, ny=jtot, lgpu=.false.)
      surf_enabled = .true.
      call add_output_file(surf_file, dtav, surf_file_id)

      if (isurf == 1) then
        call surf_file%add_var('Qnet', 'Net radiation', 'W/m^2', 'tt0t')
        call surf_file%add_var('H', 'Sensible heat flux', 'W/m^2', 'tt0t')
        call surf_file%add_var('LE', 'Latent heat flux', 'W/m^2', 'tt0t')
        call surf_file%add_var('G0', 'Ground heat flux', 'W/m^2', 'tt0t')
        call surf_file%add_var('tskin', 'Skin temperature', 'K', 'tt0t')
        call surf_file%add_var('tendskin', 'Skin tendency', 'W/m^2', 'tt0t')
        call surf_file%add_var('rs', 'Surface resistance', 's/m', 'tt0t')
        call surf_file%add_var('ra', 'Aerodynamic resistance', 's/m', 'tt0t')
        call surf_file%add_var('cliq', 'Fraction of vegetated surface covered with liquid water', '-', 'tt0t')
        call surf_file%add_var('Wl', 'Liquid water reservoir', 'm', 'tt0t')
        call surf_file%add_var('rssoil', 'Soil evaporation resistance', 's/m', 'tt0t')
        call surf_file%add_var('rsveg', 'Vegetation resistance', 's/m', 'tt0t')
      else if (isurf == 2) then
        call surf_file%add_var('hfss', 'Surface upward sensible heat flux', 'W/m^2', 'tt0t')
        call surf_file%add_var('hfls', 'Surface upward latent heat flux', 'W/m^2', 'tt0t')
        call surf_file%add_var('obuk', 'Obukhov length', 'm', 'tt0t')
        call surf_file%add_var('ustar', 'Friction velocity', 'm/s^-1', 'tt0t')
        call surf_file%add_var('Cs', 'Drag coefficient for scalars', '-', 'tt0t')
        call surf_file%add_var('Cm', 'Drag coefficient for momentum', '-', 'tt0t')
        call surf_file%add_var('z0h', 'Surface roughness length for heat', 'm', 'tt0t')
        call surf_file%add_var('z0m', 'Surface roughness length for momentum', 'm', 'tt0t')
      else if (isurf == 11) then
        call surf_file%add_var('H', 'Sensible heat flux', 'W/m^2', 'tt0t')
        call surf_file%add_var('LE', 'Latent heat flux', 'W/m^2', 'tt0t')
        call surf_file%add_var('G0', 'Ground heat flux', 'W/m^2', 'tt0t')
        call surf_file%add_var('tskin', 'Skin temperature', 'K', 'tt0t')
        call surf_file%add_var('obuk', 'Obukhov length', 'm', 'tt0t')
        call surf_file%add_var('ustar', 'Friction velocity', 'm/s^-1', 'tt0t')
        call surf_file%add_var('cliq', 'Fraction of vegetated surface covered with liquid water', '-', 'tt0t')
        call surf_file%add_var('wl', 'Liquid water reservoir', 'm', 'tt0t')
        call surf_file%add_var('ra', 'Aerodynamic resistance', 's/m', 'tt0t')
        call surf_file%add_var('rssoil', 'Soil evaporation resistance', 's/m', 'tt0t')
        call surf_file%add_var('rsveg', 'Vegetation resistance', 's/m', 'tt0t')
        call surf_file%add_var('f1', 'f1(SWD) function vegetation resistance', 's/m', 'tt0t')
        call surf_file%add_var('f2_b', 'f2(theta) function soil resistance', 's/m', 'tt0t')
        call surf_file%add_var('Qnet', 'Net radiation', 'W/m^2', 'tt0t')
        if (lags) then
          call surf_file%add_var('an_co2', 'Net CO2 assimilation', 'ppm m s-1', 'tt0t')
          call surf_file%add_var('resp_co2', 'CO2 respiration soil + plant', 'ppm m s-1', 'tt0t')
        end if
      end if

      if (lenable_gradients) then
        call surf_file%add_var('dudz', 'U-wind gradient in surface layer', 's^-1', 'tt0t')
        call surf_file%add_var('dvdz', 'V-wind gradient in surface layer', 's^-1', 'tt0t')
        call surf_file%add_var('dqtdz', 'Specific humidity gradient in surface layer', 'kg kg^-1 m^-1', 'tt0t')
        call surf_file%add_var('dthldz', 'Liquid water potential temperature gradient in surface layer', 'K m^-1', 'tt0t')
      end if

      if (lenable_temp_more) then
        call surf_file%add_var('thlskin', 'Grid-cell mean skin liquid water potential temperature', 'K', 'tt0t')
        call surf_file%add_var('qtskin', 'Grid-cell mean skin specific humidity', 'kg/kg', 'tt0t')
        call surf_file%add_var('db', 'Grid-cell mean buoyancy difference surface-atmosphere', 'm s^-2', 'tt0t')
        call surf_file%add_var('thl0_1', 'Lowest model level liquid water potential temperature', 'K', 'tt0t')
        call surf_file%add_var('qt0_1', 'Lowest model level specific humidity', 'kg/kg', 'tt0t')
      end if
    end if

  end subroutine initlsmcrosssection
!>Run lsmcrosssection. Mainly timekeeping
  subroutine lsmcrosssection
    use modglobal, only : rk3step, timee, dt_lim

    implicit none


    if (.not. (lcross .or. lcrosssoil)) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    if (lcrosssoil) then
       call wrtvert
       call wrthorz
    end if
    if (lcross) then
       call wrtsurf
    end if
  end subroutine lsmcrosssection


!> Do the xz lsmcrosssections and dump them to file
  subroutine wrtvert
    use modglobal,   only : i1, j1, cexpnr, ifoutput
    use modstat_nc,  only : lnetcdf
    use modsurfdata, only : tsoil, phiw

    implicit none

    integer :: i, k
    real(field_r), pointer :: tsoil_ptr(:,:), phiw_ptr(:,:)

    if (crossplane_local < 2 .or. crossplane_local > j1) return

    open(ifoutput,file='movv_tsoil.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((tsoil(i,crossplane_local,k),i=2,i1),k=1,ksoilmax)
    close(ifoutput)

    open(ifoutput,file='movv_phiw.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((phiw(i,crossplane_local,k),i=2,i1),k=1,ksoilmax)
    close(ifoutput)

    if (.not. (lnetcdf .and. soil_xz_enabled)) return

    call soil_xz_file%get_pointer('tsoil', tsoil_ptr)
    call soil_xz_file%get_pointer('phiw', phiw_ptr)

    tsoil_ptr(:,:) = tsoil(2:i1, crossplane_local, 1:ksoilmax)
    phiw_ptr(:,:) = phiw(2:i1, crossplane_local, 1:ksoilmax)
  end subroutine wrtvert

!> Do the xy lsmcrosssections and dump them to file
  subroutine wrthorz
    use modglobal,   only : i1, j1, cexpnr, ifoutput
    use modstat_nc,  only : lnetcdf
    use modsurfdata, only : tsoil, phiw

    implicit none

    integer :: i, j
    real(field_r), pointer :: tsoil_ptr(:,:), phiw_ptr(:,:)

    write(cheight,'(i4.4)') crossheight
    open(ifoutput,file='movh_tsoil.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((tsoil(i,j,crossheight),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_phiw.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((phiw(i,j,crossheight),i=2,i1),j=2,j1)
    close(ifoutput)

    if (.not. (lnetcdf .and. soil_xy_enabled)) return

    call soil_xy_file%get_pointer('tsoil', tsoil_ptr)
    call soil_xy_file%get_pointer('phiw', phiw_ptr)

    tsoil_ptr(:,:) = tsoil(2:i1, 2:j1, crossheight)
    phiw_ptr(:,:) = phiw(2:i1, 2:j1, crossheight)
  end subroutine wrthorz

  !> Do the xy lsmcrosssections and dump them to file
  subroutine wrtsurf
    use modglobal,   only : i1, j1, cp, rlv
    use modfields,   only : rhof, thl0, qt0
    use modlsm,      only : f1, f2b, lags, an_co2, resp_co2
    use modlsmdata,  only : tile, nlu
    use modstat_nc,  only : lnetcdf
    use modsurfdata, only : Qnet, H, LE, G0, rs, ra, tskin, tendskin, &
                            cliq, rsveg, rssoil, Wl, isurf, obl, ustar, &
                            Cs, Cm, z0h, z0m, qtflux, thlflux, dudz, dvdz, dqtdz, dthldz, qskin

    implicit none

    integer :: i, j, ilu_idx
    real(field_r), pointer :: qnet_ptr(:,:), h_ptr(:,:), le_ptr(:,:), g0_ptr(:,:), tskin_ptr(:,:), tendskin_ptr(:,:), &
                              rs_ptr(:,:), ra_ptr(:,:), cliq_ptr(:,:), wl_ptr(:,:), rssoil_ptr(:,:), rsveg_ptr(:,:), &
                              hfss_ptr(:,:), hfls_ptr(:,:), obuk_ptr(:,:), ustar_ptr(:,:), cs_ptr(:,:), cm_ptr(:,:), &
                  z0h_ptr(:,:), z0m_ptr(:,:), f1_ptr(:,:), f2_b_ptr(:,:), an_co2_ptr(:,:), resp_co2_ptr(:,:), &
                  dudz_ptr(:,:), dvdz_ptr(:,:), dqtdz_ptr(:,:), dthldz_ptr(:,:), thlskin_ptr(:,:), qtskin_ptr(:,:), &
                  db_ptr(:,:), thl0_1_ptr(:,:), qt0_1_ptr(:,:)

    if (.not. (lnetcdf .and. surf_enabled)) return

    if (isurf == 1) then
      call surf_file%get_pointer('Qnet', qnet_ptr)
      call surf_file%get_pointer('H', h_ptr)
      call surf_file%get_pointer('LE', le_ptr)
      call surf_file%get_pointer('G0', g0_ptr)
      call surf_file%get_pointer('tskin', tskin_ptr)
      call surf_file%get_pointer('tendskin', tendskin_ptr)
      call surf_file%get_pointer('rs', rs_ptr)
      call surf_file%get_pointer('ra', ra_ptr)
      call surf_file%get_pointer('cliq', cliq_ptr)
      call surf_file%get_pointer('Wl', wl_ptr)
      call surf_file%get_pointer('rssoil', rssoil_ptr)
      call surf_file%get_pointer('rsveg', rsveg_ptr)

      qnet_ptr(:,:) = Qnet(2:i1,2:j1)
      h_ptr(:,:) = H(2:i1,2:j1)
      le_ptr(:,:) = LE(2:i1,2:j1)
      g0_ptr(:,:) = G0(2:i1,2:j1)
      tskin_ptr(:,:) = tskin(2:i1,2:j1)
      tendskin_ptr(:,:) = tendskin(2:i1,2:j1)
      rs_ptr(:,:) = rs(2:i1,2:j1)
      ra_ptr(:,:) = ra(2:i1,2:j1)
      cliq_ptr(:,:) = cliq(2:i1,2:j1)
      wl_ptr(:,:) = Wl(2:i1,2:j1)
      rssoil_ptr(:,:) = rssoil(2:i1,2:j1)
      rsveg_ptr(:,:) = rsveg(2:i1,2:j1)
    else if (isurf == 2) then
      call surf_file%get_pointer('hfss', hfss_ptr)
      call surf_file%get_pointer('hfls', hfls_ptr)
      call surf_file%get_pointer('obuk', obuk_ptr)
      call surf_file%get_pointer('ustar', ustar_ptr)
      call surf_file%get_pointer('Cs', cs_ptr)
      call surf_file%get_pointer('Cm', cm_ptr)
      call surf_file%get_pointer('z0h', z0h_ptr)
      call surf_file%get_pointer('z0m', z0m_ptr)

      hfss_ptr(:,:) = rhof(1) * cp * thlflux(2:i1,2:j1)
      hfls_ptr(:,:) = rhof(1) * rlv * qtflux(2:i1,2:j1)
      obuk_ptr(:,:) = obl(2:i1,2:j1)
      ustar_ptr(:,:) = ustar(2:i1,2:j1)
      cs_ptr(:,:) = Cs(2:i1,2:j1)
      cm_ptr(:,:) = Cm(2:i1,2:j1)
      z0h_ptr(:,:) = z0h(2:i1,2:j1)
      z0m_ptr(:,:) = z0m(2:i1,2:j1)
    else if (isurf == 11) then
      !$acc update host(H, LE, G0, tskin, qskin, obl, ustar, cliq, Wl, ra, rssoil, rsveg, f1, f2b, dudz, dvdz, dqtdz, dthldz, thl0, qt0)
      call surf_file%get_pointer('H', h_ptr)
      call surf_file%get_pointer('LE', le_ptr)
      call surf_file%get_pointer('G0', g0_ptr)
      call surf_file%get_pointer('tskin', tskin_ptr)
      call surf_file%get_pointer('obuk', obuk_ptr)
      call surf_file%get_pointer('ustar', ustar_ptr)
      call surf_file%get_pointer('cliq', cliq_ptr)
      call surf_file%get_pointer('wl', wl_ptr)
      call surf_file%get_pointer('ra', ra_ptr)
      call surf_file%get_pointer('rssoil', rssoil_ptr)
      call surf_file%get_pointer('rsveg', rsveg_ptr)
      call surf_file%get_pointer('f1', f1_ptr)
      call surf_file%get_pointer('f2_b', f2_b_ptr)
      call surf_file%get_pointer('Qnet', qnet_ptr)

      h_ptr(:,:) = H(2:i1,2:j1)
      le_ptr(:,:) = LE(2:i1,2:j1)
      g0_ptr(:,:) = G0(2:i1,2:j1)
      tskin_ptr(:,:) = tskin(2:i1,2:j1)
      obuk_ptr(:,:) = obl(2:i1,2:j1)
      ustar_ptr(:,:) = ustar(2:i1,2:j1)
      cliq_ptr(:,:) = cliq(2:i1,2:j1)
      wl_ptr(:,:) = Wl(2:i1,2:j1)
      ra_ptr(:,:) = ra(2:i1,2:j1)
      rssoil_ptr(:,:) = rssoil(2:i1,2:j1)
      rsveg_ptr(:,:) = rsveg(2:i1,2:j1)
      f1_ptr(:,:) = f1(2:i1,2:j1)
      f2_b_ptr(:,:) = f2b(2:i1,2:j1)
      qnet_ptr(:,:) = Qnet(2:i1,2:j1)
      if (lags) then
        call surf_file%get_pointer('an_co2', an_co2_ptr)
        call surf_file%get_pointer('resp_co2', resp_co2_ptr)

        an_co2_ptr(:,:) = an_co2(2:i1,2:j1)
        resp_co2_ptr(:,:) = resp_co2(2:i1,2:j1)
      end if
    end if

    if (lenable_gradients) then
      call surf_file%get_pointer('dudz', dudz_ptr)
      call surf_file%get_pointer('dvdz', dvdz_ptr)
      call surf_file%get_pointer('dqtdz', dqtdz_ptr)
      call surf_file%get_pointer('dthldz', dthldz_ptr)

      dudz_ptr(:,:) = dudz(2:i1,2:j1)
      dvdz_ptr(:,:) = dvdz(2:i1,2:j1)
      dqtdz_ptr(:,:) = dqtdz(2:i1,2:j1)
      dthldz_ptr(:,:) = dthldz(2:i1,2:j1)
    end if

    if (lenable_temp_more) then
      do ilu_idx = 1, nlu
        !$acc update host(tile(ilu_idx)%thlskin, tile(ilu_idx)%qtskin, tile(ilu_idx)%db, tile(ilu_idx)%frac)
      end do

      call surf_file%get_pointer('thlskin', thlskin_ptr)
      call surf_file%get_pointer('qtskin', qtskin_ptr)
      call surf_file%get_pointer('db', db_ptr)
      call surf_file%get_pointer('thl0_1', thl0_1_ptr)
      call surf_file%get_pointer('qt0_1', qt0_1_ptr)

      thlskin_ptr(:,:) = tskin(2:i1,2:j1)
      qtskin_ptr(:,:) = qskin(2:i1,2:j1)
      thl0_1_ptr(:,:) = thl0(2:i1,2:j1,1)
      qt0_1_ptr(:,:) = qt0(2:i1,2:j1,1)

      db_ptr(:,:) = 0._field_r
      do ilu_idx = 1, nlu
        do j = 2, j1
          do i = 2, i1
            db_ptr(i,j) = db_ptr(i,j) + tile(ilu_idx)%frac(i,j) * tile(ilu_idx)%db(i,j)
          end do
        end do
      end do
    end if


  end subroutine wrtsurf

end module modlsmcrosssection
