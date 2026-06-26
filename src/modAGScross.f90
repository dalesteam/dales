!> \file modAGScross.f90
!!   Dumps an instantenous AGScross of the field

!>
!! Dumps an instantenous AGScross of the field.
!>
!! AGScrosss in the yz-plane and in the xy-plane            |
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
module modAGScross


  use modglobal, only : longint, kmax, itot, jtot
  use modlogging, only: finish, warning
  use modnetcdf_file_t, only : cross_section_file_t
  use modprecision, only : field_r
  use modstat_nc_files, only : add_output_file

implicit none

character(len=*), parameter :: modname = 'modAGScross'

private
PUBLIC :: initAGScross, AGScross,exitAGScross
save
  type(cross_section_file_t) :: ags_file
  integer :: ags_file_id = 0
  logical :: ags_file_enabled = .false.

  real    :: dtav
  integer(kind=longint) :: idtav,tnext
  logical :: lAGScross = .false. !< switch for doing the AGScross (on/off)

contains
!> Initializing AGScross. Read out the namelist, initializing the variables
  subroutine initAGScross
    use modmpi,   only :myid,mpierr,comm3d, D_MPI_BCAST
    use modglobal,only :ifnamopt,fname_options,dtmax, dtav_glob,ladaptive,dt_lim,tres,btime,checknamelisterror,timee
    use modstat_nc,only : lnetcdf
    use modsurfdata, only : lrsAgs, ksoilmax,lsplitleaf
    use modraddata,only   : irad_par,irad_rrtmg,irad_rte_rrtmgp,iradiation
    use fortran_support, only: nnml_output
   implicit none

    integer :: ierr
    character(len=*), parameter :: routine = modname//'/initAGScross'

    namelist/NAMAGScross/ &
    lAGScross, dtav

    dtav = dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMAGScross,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMAGScross')
      write(nnml_output ,NAMAGScross)
      close(ifnamopt)
    end if

    if (.not. lrsAgs) lAGScross = .false.
    if (lAGScross .and. .not. lnetcdf) then
      lAGScross = .false.
      call warning(routine, 'Ignoring lAGScross, AGScross output implemented only for netcdf output.')
    end if

    call D_MPI_BCAST(dtav     ,1 ,0,comm3d,mpierr)
    call D_MPI_BCAST(lAGScross,1 ,0,comm3d,mpierr)

    idtav = int(dtav / tres, kind=kind(idtav))
    tnext   = idtav+btime
    if(.not.(lAGScross)) return
    dt_lim = min(dt_lim,tnext - timee)
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      call finish(routine, 'AGScross: dtav should be a integer multiple of dtmax')
    end if
    if (ksoilmax /= 4) call finish(routine, 'ksoilmax is not equal to 4... this can give problems with AGScross.f90... update this file as well')

    ags_file = cross_section_file_t('crossAGS', nx=itot, ny=jtot, lgpu=.false.)
    call add_output_file(ags_file, dtav, ags_file_id)
    ags_file_enabled = .true.

    call ags_file%add_var('An', 'xy AGScross of An', 'mg/m2/s', 'tt0t')
    call ags_file%add_var('Resp', 'xy AGScross of Resp', 'mg/m2/s', 'tt0t')
    call ags_file%add_var('wco2', 'xy AGScross of wco2', 'ppm m/s', 'tt0t')
    call ags_file%add_var('rs', 'xy AGScross of diagnosed rs', 's/m', 'tt0t')
    call ags_file%add_var('ra', 'xy AGScross of ra', 's/m', 'tt0t')
    call ags_file%add_var('rsCO2', 'xy AGScross of rsCO2', 's/m', 'tt0t')
    call ags_file%add_var('rsveg', 'xy AGScross of rsveg=rsAgs', 's/m', 'tt0t')
    call ags_file%add_var('rssoil', 'xy AGScross of rssoil', 's/m', 'tt0t')
    call ags_file%add_var('fstr', 'xy AGScross of stress fnct.', '-', 'tt0t')
    call ags_file%add_var('phiw1', 'xy AGScross of phiw top', '-', 'tt0t')
    call ags_file%add_var('phiw2', 'xy AGScross of phiw level 2', '-', 'tt0t')
    call ags_file%add_var('phiw3', 'xy AGScross of phiw level 3', '-', 'tt0t')
    call ags_file%add_var('phiw4', 'xy AGScross of phiw level 4', '-', 'tt0t')
    call ags_file%add_var('CO2', 'xy AGScross of CO2 (grid 1)', 'ppm', 'tt0t')
    call ags_file%add_var('tskin', 'xy AGScross of curr. tskin', 'K', 'tt0t')
    call ags_file%add_var('tskinm', 'xy AGScross of prev. tskin', 'K', 'tt0t')
    call ags_file%add_var('tsoil1', 'xy AGScross of tsoil top', 'K', 'tt0t')
    call ags_file%add_var('tsoil2', 'xy AGScross of tsoil lvl 2', 'K', 'tt0t')
    call ags_file%add_var('tsoil3', 'xy AGScross of tsoil lvl 3', 'K', 'tt0t')
    call ags_file%add_var('tsoil4', 'xy AGScross of tsoil lvl 4', 'K', 'tt0t')
    call ags_file%add_var('wtheta', 'xy AGScross of kin. heat fl', 'K m/s', 'tt0t')
    call ags_file%add_var('wq', 'xy AGScross of kin. wat. fl', '- m/s', 'tt0t')
    call ags_file%add_var('lwp', 'xy AGScross of liq. wat. p.', 'kg/m2', 'tt0t')
    call ags_file%add_var('tau', 'xy AGScross of opt. thickn.', '-', 'tt0t')
    call ags_file%add_var('swd', 'xy AGScross of SW down rad.', 'W/m2', 'tt0t')
    call ags_file%add_var('swu', 'xy AGScross of SW up rad.', 'W/m2', 'tt0t')
    call ags_file%add_var('lwd', 'xy AGScross of LW down rad.', 'W/m2', 'tt0t')
    call ags_file%add_var('lwu', 'xy AGScross of LW up rad.', 'W/m2', 'tt0t')
    call ags_file%add_var('ci', 'xy AGScross of int CO2 conc', 'mg/m3', 'tt0t')
    call ags_file%add_var('gc_CO2', 'xy AGScross of gc_CO2', 'mm/s?', 'tt0t')
    call ags_file%add_var('PAR', 'xy AGScross of PAR', 'W/m2', 'tt0t')
    call ags_file%add_var('Qnet', 'xy AGScross of Qnet', 'W/m2', 'tt0t')
    call ags_file%add_var('LE', 'xy AGScross of LE', 'W/m2', 'tt0t')
    call ags_file%add_var('H', 'xy AGScross of H', 'W/m2', 'tt0t')
    call ags_file%add_var('G0', 'xy AGScross of G0', 'W/m2', 'tt0t')
    if (iradiation == irad_par .or. iradiation == irad_rrtmg .or. iradiation == irad_rte_rrtmgp) then
      call ags_file%add_var('swdir', 'xy AGScross of SW dir rad.', 'W/m2', 'tt0t')
      call ags_file%add_var('swdif', 'xy AGScross of SW diff rad.', 'W/m2', 'tt0t')
      if (lsplitleaf) then
        call ags_file%add_var('PARdir', 'xy AGScross of direct PAR', 'W/m2', 'tt0t')
        call ags_file%add_var('PARdif', 'xy AGScross of diffuse PAR', 'W/m2', 'tt0t')
      end if
    end if

  end subroutine initAGScross
!>Run AGScross. Mainly timekeeping
  subroutine AGScross
    use modglobal, only : rk3step,timee,dt_lim
    implicit none


    if (.not. lAGScross) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    call AGShorz

  end subroutine AGScross


!> Do the xy AGScrosss and dump them to file
  subroutine AGShorz
    use modglobal, only : i1,j1,dzf
    use modsurfdata, only : AnField, RespField, wco2Field,phiw,fstrField, rs, ra, rsco2Field, rsveg, rssoil, &
                            indCO2, tskin, tskinm, tsoil, thlflux, qtflux, tauField, ciField, gcco2Field, &
                            PARField,Qnet,LE,H,G0,PARdirField,PARdifField,lsplitleaf
    use modfields, only   : svm, rhof, ql0
    use modraddata,only   : swd, swu, lwd, lwu,swdir,swdif,irad_par,iradiation,irad_rrtmg,irad_rte_rrtmgp,lwc
    implicit none


    ! LOCAL
    integer i,j
    real(field_r) :: lwp(2:i1,2:j1)
    real(field_r), pointer :: an_ptr(:,:), resp_ptr(:,:), wco2_ptr(:,:), rs_ptr(:,:), ra_ptr(:,:), rsco2_ptr(:,:), &
                  rsveg_ptr(:,:), rssoil_ptr(:,:), fstr_ptr(:,:), phiw1_ptr(:,:), phiw2_ptr(:,:), &
                  phiw3_ptr(:,:), phiw4_ptr(:,:), co2_ptr(:,:), tskin_ptr(:,:), tskinm_ptr(:,:), &
                  tsoil1_ptr(:,:), tsoil2_ptr(:,:), tsoil3_ptr(:,:), tsoil4_ptr(:,:), wtheta_ptr(:,:), &
                  wq_ptr(:,:), lwp_ptr(:,:), tau_ptr(:,:), swd_ptr(:,:), swu_ptr(:,:), lwd_ptr(:,:), &
                  lwu_ptr(:,:), ci_ptr(:,:), gcco2_ptr(:,:), par_ptr(:,:), qnet_ptr(:,:), le_ptr(:,:), &
                  h_ptr(:,:), g0_ptr(:,:), swdir_ptr(:,:), swdif_ptr(:,:), pardir_ptr(:,:), pardif_ptr(:,:)

    if (.not. ags_file_enabled) return

    do i = 2,i1
      do j = 2,j1
        if (iradiation == irad_rrtmg) then
          lwp(i,j) = sum(lwc(i,j,1:kmax))*1.e-3 ! we get the already calculated lwc from RRTMG
        else
          lwp(i,j) = sum(ql0(i,j,1:kmax)*rhof(1:kmax)*dzf(1:kmax))
        end if
      enddo
    enddo

    call ags_file%get_pointer('An', an_ptr)
    call ags_file%get_pointer('Resp', resp_ptr)
    call ags_file%get_pointer('wco2', wco2_ptr)
    call ags_file%get_pointer('rs', rs_ptr)
    call ags_file%get_pointer('ra', ra_ptr)
    call ags_file%get_pointer('rsCO2', rsco2_ptr)
    call ags_file%get_pointer('rsveg', rsveg_ptr)
    call ags_file%get_pointer('rssoil', rssoil_ptr)
    call ags_file%get_pointer('fstr', fstr_ptr)
    call ags_file%get_pointer('phiw1', phiw1_ptr)
    call ags_file%get_pointer('phiw2', phiw2_ptr)
    call ags_file%get_pointer('phiw3', phiw3_ptr)
    call ags_file%get_pointer('phiw4', phiw4_ptr)
    call ags_file%get_pointer('CO2', co2_ptr)
    call ags_file%get_pointer('tskin', tskin_ptr)
    call ags_file%get_pointer('tskinm', tskinm_ptr)
    call ags_file%get_pointer('tsoil1', tsoil1_ptr)
    call ags_file%get_pointer('tsoil2', tsoil2_ptr)
    call ags_file%get_pointer('tsoil3', tsoil3_ptr)
    call ags_file%get_pointer('tsoil4', tsoil4_ptr)
    call ags_file%get_pointer('wtheta', wtheta_ptr)
    call ags_file%get_pointer('wq', wq_ptr)
    call ags_file%get_pointer('lwp', lwp_ptr)
    call ags_file%get_pointer('tau', tau_ptr)
    call ags_file%get_pointer('swd', swd_ptr)
    call ags_file%get_pointer('swu', swu_ptr)
    call ags_file%get_pointer('lwd', lwd_ptr)
    call ags_file%get_pointer('lwu', lwu_ptr)
    call ags_file%get_pointer('ci', ci_ptr)
    call ags_file%get_pointer('gc_CO2', gcco2_ptr)
    call ags_file%get_pointer('PAR', par_ptr)
    call ags_file%get_pointer('Qnet', qnet_ptr)
    call ags_file%get_pointer('LE', le_ptr)
    call ags_file%get_pointer('H', h_ptr)
    call ags_file%get_pointer('G0', g0_ptr)

    an_ptr(:,:) = AnField(2:i1,2:j1)
    resp_ptr(:,:) = RespField(2:i1,2:j1)
    wco2_ptr(:,:) = wco2Field(2:i1,2:j1)
    rs_ptr(:,:) = rs(2:i1,2:j1)
    ra_ptr(:,:) = ra(2:i1,2:j1)
    rsco2_ptr(:,:) = rsco2Field(2:i1,2:j1)
    rsveg_ptr(:,:) = rsveg(2:i1,2:j1)
    rssoil_ptr(:,:) = rssoil(2:i1,2:j1)
    fstr_ptr(:,:) = fstrField(2:i1,2:j1)
    phiw1_ptr(:,:) = phiw(2:i1,2:j1,1)
    phiw2_ptr(:,:) = phiw(2:i1,2:j1,2)
    phiw3_ptr(:,:) = phiw(2:i1,2:j1,3)
    phiw4_ptr(:,:) = phiw(2:i1,2:j1,4)
    co2_ptr(:,:) = svm(2:i1,2:j1,1,indCO2) / 1000.0_field_r
    tskin_ptr(:,:) = tskin(2:i1,2:j1)
    tskinm_ptr(:,:) = tskinm(2:i1,2:j1)
    tsoil1_ptr(:,:) = tsoil(2:i1,2:j1,1)
    tsoil2_ptr(:,:) = tsoil(2:i1,2:j1,2)
    tsoil3_ptr(:,:) = tsoil(2:i1,2:j1,3)
    tsoil4_ptr(:,:) = tsoil(2:i1,2:j1,4)
    wtheta_ptr(:,:) = thlflux(2:i1,2:j1)
    wq_ptr(:,:) = qtflux(2:i1,2:j1)
    lwp_ptr(:,:) = lwp(2:i1,2:j1)
    tau_ptr(:,:) = tauField(2:i1,2:j1)
    swd_ptr(:,:) = swd(2:i1,2:j1,1)
    swu_ptr(:,:) = swu(2:i1,2:j1,1)
    lwd_ptr(:,:) = lwd(2:i1,2:j1,1)
    lwu_ptr(:,:) = lwu(2:i1,2:j1,1)
    ci_ptr(:,:) = ciField(2:i1,2:j1)
    gcco2_ptr(:,:) = gcco2Field(2:i1,2:j1)
    par_ptr(:,:) = PARField(2:i1,2:j1)
    qnet_ptr(:,:) = Qnet(2:i1,2:j1)
    le_ptr(:,:) = LE(2:i1,2:j1)
    h_ptr(:,:) = H(2:i1,2:j1)
    g0_ptr(:,:) = G0(2:i1,2:j1)

    if (iradiation == irad_par .or. iradiation == irad_rrtmg .or. iradiation == irad_rte_rrtmgp) then
      call ags_file%get_pointer('swdir', swdir_ptr)
      call ags_file%get_pointer('swdif', swdif_ptr)
      swdir_ptr(:,:) = swdir(2:i1,2:j1,1)
      swdif_ptr(:,:) = swdif(2:i1,2:j1,1)
      if (lsplitleaf) then
        call ags_file%get_pointer('PARdir', pardir_ptr)
        call ags_file%get_pointer('PARdif', pardif_ptr)
        pardir_ptr(:,:) = PARdirField(2:i1,2:j1)
        pardif_ptr(:,:) = PARdifField(2:i1,2:j1)
      end if
    end if

  end subroutine AGShorz


!> Clean up when leaving the run
  subroutine exitAGScross
    implicit none

  end subroutine exitAGScross

end module modAGScross
