!> \file moddrydeposition.f90
!! Calculates the (dry) deposition of gaseous components
!!
!! if `ldeposition` is .true., the deposition is calculated for the scalars
!! svskip+1 to nsv. The names of these scalars can be given in the namoptions
!! file in `tracernames`, and these will align with `emisnames` in case all emitted
!! species are also deposited.
!
! Copyright (c) 2020-2022 TNO
!
! This file is part of DALES
!
! DALES is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DALES.  If not, see <http://www.gnu.org/licenses/>.
!
module moddrydeposition
  use modlsm,    only : tile, nlu
  use modfields, only : svp   ! tracer tendency array
  use modglobal, only : nsv, i1, j1
  use modlsm, only : llsm
  use modtracers, only: tracer_prop

  implicit none

  save
  public  :: initdrydep, drydep, exitdrydep

  logical :: ldrydep           !< On/Off switch dry deposition
  real, allocatable :: depfield(:,:,:) !< deposition flux (i,j,sv) [ug * m / (g * s)]
  logical, dimension(100) :: ldeptracers = .false. !< List of switches determining which of the tracers to deposit
  integer  :: ndeptracers = 0  !< Number of tracers that deposits
  integer  :: iname

  private :: Rc, Rb, vd

  real, allocatable :: Rb(:, :)  ! Temporary storage of the quasilaminar sublayer resistance
  real, allocatable :: Rc(:, :)  ! Temporary storage of the total canopy resistance
  real, allocatable :: vd(:, :)  ! Temporary storage of the deposition velocity

  ! Diffusion coefficients (up to 'co', originating from DEPAC, rest from Carl L. Yaws, "Transport Properties of Chemicals
  ! and Hydrocarbons. Viscosity, Thermal Conductivity, and Diffusivity of C1 to C1000 Organics and Ac to Zr Inorganics", 2009)
  character(*), parameter :: species(13) = (/'o3    ', 'so2   ', 'no2   ', 'no    ', 'nh3   ', 'co    ', &
                                             'no3   ', 'hno3  ', 'n2o5  ', 'h2o2  ', 'co2   ', 'h2o   ', &
                                             'n2o   '/)
  real, parameter :: diffusivity(13) =     (/ 0.13e-4,  0.11e-4,  0.13e-4,  0.16e-4,  0.21e-4,  0.16e-4, &
                                              0.11e-4,  0.14e-4,  0.11e-4,  0.18e-4,  0.16e-4,  0.25e-4, &
                                              0.18e-4/)

! Unfortunately, there is no other way to get the parameters from DEPAC (lai_par, laitype)
#include "depac_lu.inc"
 
contains

!> Initialize deposition calculation.
!!
!! Read the namelist NAMDEPOSITION from the namoptions file, distribute
!! the parameters to all processes and allocate the deposition flux field.
subroutine initdrydep
  ! read namelist    
  ! read LU parameters from file
  ! init drydep fields

  use modglobal, only : i2, j2, nsv, ifnamopt, fname_options, &
                        checknamelisterror
  use modmpi,    only : myid, comm3d, mpi_logical

  implicit none

  ! Auxiliary variables
  integer  :: ierr, isv

  ! ---------------------------------------------------------------------!
  ! Main variables                                                       !
  ! ---------------------------------------------------------------------!

  ! --- Read & broadcast namelist DEPOSITION -----------------------------------
  namelist/NAMDEPOSITION/ ldrydep

  if (myid == 0) then
    open(ifnamopt,file=fname_options,status='old',iostat=ierr)
    read (ifnamopt,nml=NAMDEPOSITION,iostat=ierr)
    call checknamelisterror(ierr, ifnamopt, 'NAMDEPOSITION')
    write(6, NAMDEPOSITION)
    close(ifnamopt)
  endif

  call mpi_bcast(ldrydep,              1, mpi_logical,   0, comm3d, ierr)

  do isv = 1,nsv
    if (.not. tracer_prop(isv)%ldep) cycle
    ndeptracers = ndeptracers + 1
    write(*,*) 'tracer is deposited: ', tracer_prop(isv)%tracname
  enddo

  ! --- Local pre-calculations and settings
  if (ldrydep .and. ndeptracers == 0) then
    write (*,*) "initdrydep: WARNING .. drydeposition switched on, but no tracers to deposit. &
      Continuing without deposition model"
  end if
  if (.not. (ldrydep) .or. .not. (llsm) .or. ndeptracers == 0)  return

  allocate(depfield(i2, j2, ndeptracers))
  allocate(Rb(i2, j2))
  allocate(Rc(i2, j2))
  allocate(vd(i2, j2))

  Rb = 0.0
  Rc = 0.0
  vd = 0.0

end subroutine initdrydep

!> Run deposition calculation.
!!
!! This routine calls `calc_depfield` to calculate the deposition flux and
!! adds it to the scalar variables tendency
!!
!! @see calc_depfield
subroutine drydep  ! called in program.f90
  use modglobal, only : dzf

  implicit none
  integer :: isv, idt, i, j ! isv = sv index, idt = deptracer index

  ! if (ldrydep .and. ndeptracers == 0) then
  !   write (*, *) "drydep: Skipping deposition calculation, since no tracers selected (ndeptracers = 0)"
  ! end if
  if (.not. llsm .or. .not. ldrydep .or. ndeptracers == 0) return  ! Dry deposition cannot be run if LSM not activated
  
  call calc_depfield

  ! tracer tendency due to deposition. The calculated deposition flux `depfield` is
  ! negative. Only tracers that will deposit (ldeptracers(isv) == .true.) have an entry in depfield to
  ! limit memory usage. 
  idt = 1  ! deposition tracers have their own index
  do isv = 1, nsv
    if (.not. tracer_prop(isv)%ldep) cycle
    do j = 2, j1
      do i = 2, i1
        svp(i,j,1,tracer_prop(isv)%trac_idx) = svp(i,j,1,tracer_prop(isv)%trac_idx) + depfield(i,j,idt) / dzf(1)
      end do
    end do
    idt = idt + 1
  end do

end subroutine drydep    

!> Wrapper function around call to DEPAC to calculate total canopy resistance
!!
!! Rc is calculated with the DEPAC model (v3.21 with minor changes by TNO) using the
!! classical approach without compensation points. 
!!
!! The DEPAC specific calculations are performed by `DryDepos_Gas_DEPAC`.
!!
!! @param[in] ilu Tile index number (defined by land use class)
!!
!! @see M.C. van Zanten et al., "Description of the DEPAC module", RIVM report nr. 680180001/2010
!! @see DryDepos_Gas_DEPAC
subroutine depac_call(ilu, species)
  use modlsm, only : tile
  use modglobal, only : i1, j1, xday, xlat, xlon, xtime, rtimee
  use modfields, only : thl0, exnf, qt0, qsat
  use le_drydepos_gas_depac, only : DryDepos_Gas_DEPAC
  use modraddata, only : zenith, swd
  use go, only : to_upper
  implicit none

  integer, intent(in) :: ilu
  character(*), intent(in) :: species
  character(len=6) :: depac_species
  integer :: i, j, nwet = 0, status, depac_ilu
  real :: T, RH, ccomp_tot, sinphi, lai, sai

  ! Temporary values, until something better is available
  ! for now, assuming low NH3/SO2 ratios
  integer, parameter :: iratns = 1

  ! TEMPORARY HACK, Assume NOx=NO2, since DEPAC needs NO/NO2
  if (trim(species) == 'nox') then
    depac_species = 'NO2' 
  else 
    depac_species = to_upper(trim(species))
  end if

  depac_ilu = get_depac_luindex(tile(ilu)%lushort)

  ! Currently used conventions:
  ! - To calculate RH, qt0 and qsat are used. There may be a better way, in fact, I believe qt is the total 
  !   specific humidity, which includes liquid moisture. Should ql be subtracted?
  ! - Latitude fed to DEPAC is the latitude of the lower left corner of the domain
  ! - As temperature, the temperature in the lowest layer is used. It is disputable whether this is correct, might need
  !   to be changed to 2m temperature or leaf temperature (provided they are available)
  ! - The `nwet` parameter is only used to discern between dry and wet; snow is not covered yet.
  ! - The parameter `iratns` (ratio NH3/SO2) is set to 1 (low)
  ! - The parameter `ccomp_tot` needs to be provided, but is not used (default 0)
  ! - Parameters after `ccomp_tot` are: `hlaw` for Henry's law constant, which is not used for non-VBS
  !   species, `react` for reactivity (not used for non-VBS species) and `status` for registering errors,
  !   but these are already handled in the depac routine.
  call calc_lai_sai(tile(ilu)%lushort, xday, xlat, tile(ilu)%SAI_a(2, 2), &
                    tile(ilu)%SAI_b(2, 2), lai, sai)  ! running the LAI & SAI here save imax x jmax - 
  sinphi = zenith(xtime*3600 + rtimee, xday, xlat, xlon)
  call flush(6)
  do i = 2, i1
    do j = 2, j1
      T = thl0(i, j, 1) * exnf(1)
      RH = qt0(i, j, 1) / qsat(i, j, 1) * 100
      ! swd needs to be negated, since it is pointing downward.
      ! tsea is a temperature DEPAC needs in case of water LU classes
      call DryDepos_Gas_DEPAC(depac_species, int(xday), xlat, T, &
                              tile(ilu)%ustar(i, j), -swd(i, j, 1), sinphi, RH, lai, sai, nwet, &
                              depac_ilu, iratns, Rc(i, j), ccomp_tot, 0.0, 0.0, status, tsea=tile(ilu)%tskin(i, j))
      ! check for missing Rc values, i.e. -9999, and return huge resistance, so virtually no deposition takes place
      if (missing_real(Rc(i, j), -9999.)) then
        Rc(i,j) = 1.e5
      endif
    end do
  end do

  ! ! DEBUG feedback
  ! write (6, '("DEPAC: Land use class= ",a)') tile(ilu)%lushort
  ! write (6, '("DEPAC: species = ",a)') to_upper(trim(species))
  ! write (6, '("DEPAC: xday = ",i3)') int(xday)
  ! write (6, '("DEPAC: xlat= ",f10.4)') xlat
  ! write (6, '("DEPAC: ustar(2, 2) = ",e12.4)') tile(ilu)%ustar(2, 2)
  ! write (6, '("DEPAC: -swd(2, 2, 1)= ",f10.2)') -swd(2, 2, 1)
  ! write (6, '("DEPAC: sinphi= ",f10.4)') sinphi
  ! write (6, '("DEPAC: temperature= ",f10.2)') T
  ! write (6, '("DEPAC: qt0(2, 2, 1)= ",f10.4)') qt0(2, 2, 1)
  ! write (6, '("DEPAC: qsat(2, 2, 1)= ",f10.4)') qsat(2, 2, 1)
  ! write (6, '("DEPAC: RH= ",f10.3)') RH
  ! write (6, '("DEPAC: LAI= ",f10.1)') lai
  ! write (6, '("DEPAC: SAI= ",f10.1)') sai
  ! write (6, '("DEPAC: DEPAC ilu= ",i3)') depac_ilu
  ! write (6, '("DEPAC: Rc(2, 2)= ",f12.2)') Rc(2, 2)

end subroutine depac_call

!> Finalize deposition calculation.
!!
!! Clean up the deposition flux field.
subroutine exitdrydep
  implicit none

  ! if (.not. llsm .or. .not. ldrydep .or. count(ldeptracers) == 0) return
  if (.not. llsm .or. .not. ldrydep .or. ndeptracers == 0) return

  deallocate(depfield)
  deallocate(Rb)
  deallocate(Rc)
  deallocate(vd)

end subroutine exitdrydep

!> Calculate the deposition flux for species other than water.
!!
!! This is the core of the module. The deposition flux is
!! calculated from the aerodynamic resistance, the quasilaminar layer resistance 
!! and the canopy resistance.
subroutine calc_depfield
  use modglobal, only : i1, j1, fkar
  use modfields, only : sv0
  use modlsm, only : tile, nlu
  ! Necessary to retrieve/calculate all parameters necessary for the DEPAC routine
  implicit none
 
  integer :: ilu, isv, idt, i, j ! Indices for land use type (ilu), scalar variable (isv) and deposited tracer (idt)
  real :: Sc, ScPrfac
  real, parameter :: Pr_air = 0.71  !< Prandtl number of air at atmospheric conditions (+/- 1.5%)
  real, parameter :: nu_air = 15e-6 !< Kinematic viscosity of air [m2/s]
  ! The variation of values of nu_air in literature is quite substantial, so for now temperature dependence ignored.
  ! This will only have a minor effect, since Rb is usually much smaller than Ra+Rc

  idt = 1
  do isv = 1, nsv
    if (.not. tracer_prop(isv)%ldep) cycle
    depfield(:, :, idt) = 0.0
    Sc = nu_air / findval(tracer_prop(isv)%tracname, species, diffusivity, defltvalue=.22e-4)  ! Default is O2 in air
    ScPrfac = (Sc/Pr_air) ** (2/3.)
    ! HACK: now running only on non-wet tiles. Wet surfaces still to be covered.
    do ilu = 1, nlu - 1
      ! ilu = 5
      call depac_call(ilu, tracer_prop(isv)%tracname) ! Update Rc
      ! Quasilaminar sublayer resistance according to Hicks et al, Water Air Soil Pollut., v35, p311-330, 1987
      do j = 2, j1
        do i = 2, i1
          if (tile(ilu)%frac(i,j) > 0.) then
            Rb(i,j) = 1/(fkar*tile(ilu)%ra(i,j)) * ScPrfac
            vd(i,j) = (tile(ilu)%ra(i,j) + Rb(i,j) + Rc(i,j)) ** (-1)
            ! Deposition flux in ppb m s-1
            depfield(i,j,idt) = depfield(i,j,idt) - tile(ilu)%frac(i,j) * vd(i,j) * sv0(i,j,1,tracer_prop(isv)%trac_idx)
          endif
        end do
      end do
    end do  ! ilu = 1, nlu
    idt = idt + 1
  end do  ! isv = i, nsv

end subroutine calc_depfield

!> Find a value in an array of values based on a key in an array of keys
!!
!! Lookup a value in one array (`values`) at the index position of `key` in the array `keys`.
!! Provide the optional default value (`defltvalue`), that is returned when `key` isn't found.
!! When no default is given, zero is returned.
!!
!! @param[in] key The key to be looked up
!! @param[in] keys The array of keys
!! @param[in] values The array of values
!! @param[in] defltvalue The default value (optional, default 0.0)
!! @return The value corresponding to the key
pure real function findval(key, keys, values, defltvalue)
  implicit none
  character(*), intent(in) :: key 
  character(*), intent(in) :: keys(:)
  real, intent(in) :: values(:)
  real, intent(in), optional :: defltvalue
  character(6) :: fkey
  integer :: idx(1)

  fkey = trim(key)

  idx = findloc(keys, fkey)
  if (idx(1) /= 0 .and. idx(1) <= size(values)) then
    findval = values(idx(1))
  else
    if (present(defltvalue)) then
      findval = defltvalue 
    else
      findval = 0.0 
    end if
  end if
end function findval

!> Calculate the leaf area and surface area indices as a function of growing season
!!
!! This is (hopefully) a temporary function to be used until the land use classes of
!! the LSM model directly relate to the DEPAC classes. It uses data defined in
!! `depac_lu.inc` for now.
!!
!! @param[in] luclass The acronym of the LU class to calculate the LAI for
!! @param[in] doy Day of the year
!! @param[in] latitude The latitude of the current location
!! @param[in] SAI_a Parameter to calculate SAI from LAI (from tile)
!! @param[in] SAI_b Parameter to calculate SAI from LAI (from tile)
!! @param[out] lai The leaf area index for the tile
!! @param[out] sai The surface area index for the tile
subroutine calc_lai_sai(luclass, doy, latitude, SAI_a, SAI_b, lai, sai)
  implicit none
  character(len=3), intent(in) :: luclass
  type(laitype) :: tab_data
  real, intent(in) :: doy, latitude, SAI_a, SAI_b
  real, intent(out) :: lai, sai
  integer :: idx
  real :: sgs, mgs, egs

  idx = get_depac_luindex(luclass)  ! If the index is not found, returns zero
  if (idx == 0) then
    write (6, *) "moddrydeposition ERROR: Land use class not found in DEPAC"
    stop  ! ... so stop the calculations
  end if
  tab_data = lai_par(idx)

  ! In case the values are not defined, return zero for LAI and SAI
  if (missing_int(tab_data%sgs50, -999)) then
    lai = 0.0
    sai = 0.0
    return
  end if

  call SGS_MGS_EGS(tab_data, latitude, sgs, mgs, egs)

  ! Calculating LAI
  if (doy < sgs) then
    lai = 0.0
  else if (doy >= sgs .and. doy < mgs) then
    lai = tab_data%laimin + (tab_data%laimax - tab_data%laimin) / tab_data%s_lai_len * (doy - sgs)
  else if (doy >= mgs .and. doy < egs - tab_data%e_lai_len) then
    lai = tab_data%laimax
  else if (doy >= egs - tab_data%e_lai_len .and. doy < egs) then
    lai = tab_data%laimax + (tab_data%laimin - tab_data%laimax) / tab_data%e_lai_len * (doy - egs + tab_data%e_lai_len)
  else
    lai = 0.0
  end if
  
  ! Calculating SAI
  if (luclass == 'ara') then
    if (doy < sgs) then
      sai = lai
    else
      sai = max(5 / 3.5 * lai, lai + 1.5)
    end if
  else
    sai = SAI_a * lai + SAI_b
  end if

end subroutine calc_lai_sai

!> Calculate the latitude dependent start, maximum and end of growing season
!!
!! From DEPAC routines
!!
!! @param[in] tab_data List of LAI parameters from `depac_lu.inc`
!! @param[in] latitude The latitude of the location
!! @param[out] sgs Start of the growing season
!! @param[out] mgs Maximum of the growin season
!! @param[out] egs End of the growing season
subroutine SGS_MGS_EGS(tab_data, latitude, sgs, mgs, egs)
  implicit none
  type(laitype), intent(in) :: tab_data
  real, intent(in) :: latitude
  real, intent(out) :: sgs, mgs, egs

  ! If no data is defined, return sensible data
  if (missing_int(tab_data%sgs50, -999)) then
    sgs = 0.0
    mgs = 0.0
    egs = 365.0
    return
  end if

  sgs = tab_data%sgs50 + tab_data%dsgs * (latitude - 50.0)
  mgs = sgs + tab_data%s_lai_len
  egs = tab_data%egs50 + tab_data%degs * (latitude - 50.0)

end subroutine SGS_MGS_EGS

!> Retrieve the three letter LU class acronym from the tile and return the corresponding index in DEPAC arrays
!!
!! This function is an interface between the current LSM implementation in DALES and the one used in DEPAC.
!! Once both DALES and DEPAC incorporate an LSM model that is not hard linked to data, this function is no longer needed.
!!
!! @param[in] luclass The LU class acronym from the tile (`lushort`)
!! @returns The corresponding index in DEPAC arrays defined in `depac_lu.inc`  
pure integer function get_depac_luindex(luclass)
  implicit none
  character(len=3), intent(in) :: luclass
  integer :: idx(1)
  character(len=3) :: depac_lu

  select case (luclass)
    case ('fcd', 'fce')
      depac_lu = 'cnf'
    case ('fbd', 'fbe')
      depac_lu = 'dec'
    case ('aqu')
      depac_lu = 'wai'
    case ('sem')
      depac_lu = 'grs'
    case ('brn')
      depac_lu = 'dsr'
    case default
      depac_lu = luclass
  end select

  idx = findloc(lu_name_abbr, depac_lu)
  get_depac_luindex = idx(1)
end function get_depac_luindex

!> Auxiliary function (source: le_drydepos_gas_depac.f90)
!!
!! Apalling way to check for missing values...
!!
!! @param[in] x Check whether x is 'missing'
!! @param[in] val Value to check x against. e.g. -999
logical function missing_real(x, val)
  implicit none
  real, intent(in) :: x, val
  ! bandwidth for checking (in)equalities of floats
  real, parameter :: EPS = 1.0e-5
  missing_real = (abs(x - val) .le. EPS)
end function missing_real

logical function missing_int(x, val)
  implicit none
  integer, intent(in) :: x, val
  missing_int = missing_real(real(x), real(val))
end function missing_int

end module moddrydeposition
