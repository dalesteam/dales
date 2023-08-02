!> \file modtracdata.f90
!! Parameter definitions and auxilary routines for tracers

!>
!! This routine should have no dependency on any other routine, save
!perhaps modglobal or modfields.
!!  \author Ruud Janssen, TNO
!!  \todo Documentation
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
!  Copyright 1993-2023 Delft University of Technology, Wageningen
!  University, Utrecht University, KNMI, Vrije Universiteit Amsterdam, TNO

module modtracdata
  
  implicit none
  save

  integer, parameter :: nsv_all = 30       ! number of all available tracers

  ! tracer names:
  character(len=6), parameter  ::  tracname_short(nsv_all) = (/ &
                              'qr    ', &   ! qr
                              'nr    ', &   ! nr
                              'co2   ', &   ! co2
                              'co2veg', &   ! co2veg
                              'co2sum', &   ! co2sum
                              'nox   ', &   ! nox
                              'nh3   ', &   ! nh3
                              'so2   ', &   ! so2
                              'inert ', &   ! INERT
                              'o3    ', &   ! O3
                              'o1d   ', &   ! O1D
                              'no    ', &   ! NO
                              'no2   ', &   ! NO2
                              'ch4   ', &   ! CH4
                              'ch2o  ', &   ! CH2O
                              'ch3o2 ', &   ! CH3O2
                              'mvk   ', &   ! MVK
                              'iso   ', &   ! ISO
                              'ro2   ', &   ! RO2
                              'oh    ', &   ! OH
                              'ho2   ', &   ! HO2
                              'co    ', &   ! CO
                              'h2o   ', &   ! H2O
                              'produc', &   ! PRODUC
                              'o2    ', &   ! O2
                              'n2    ', &   ! N2
                              'hno3  ', &   ! HNO3
                              'h2o2  ', &   ! H2O2
                              'no3   ', &   ! NO3
                              'n2o5  ' /)   ! N2O5

  ! tracer long names:
  character(len=30), parameter  ::  tracname_long(nsv_all) = (/ &
                              'Rain water mixing ratio       ', &   ! qr
                              'Cloud droplet number          ', &   ! nr
                              'Carbon dioxide                ', &   ! co2
                              'Carbon dioxide - vegetation   ', &   ! co2veg
                              'Carbon dioxide - sum          ', &   ! co2sum
                              'Nitrogen oxides               ', &   ! nox
                              'Ammonia                       ', &   ! nh3
                              'Sulphur dioxide               ', &   ! so2
                              'INERT                         ', &   ! INERT
                              'Ozone                         ', &   ! O3
                              'Excited oxygen atom           ', &   ! O1D
                              'Nitrogen oxide                ', &   ! NO
                              'Nitrogen dioxide              ', &   ! NO2
                              'Methane                       ', &   ! CH4
                              'Formaldehyde                  ', &   ! CH2O
                              'Methyl peroxy radical         ', &   ! CH3O2
                              'Methyl vinyl ketone           ', &   ! MVK
                              'Isoprene                      ', &   ! ISO
                              'Peroxy radical                ', &   ! RO2
                              'Hydroxyl radical              ', &   ! OH
                              'Hydroperoxyl radical          ', &   ! HO2
                              'Carbon monoxide               ', &   ! CO
                              'Water                         ', &   ! H2O
                              'PRODUC                        ', &   ! PRODUC
                              'Oxygen                        ', &   ! O2
                              'Nitrogen                      ', &   ! N2
                              'Nitric acid                   ', &   ! HNO3
                              'Hydrogen peroxide             ', &   ! H2O2
                              'Nitrate radical               ', &   ! NO3
                              'Dinitrogen pentoxide          '/)    ! N2O5

  ! tracer units:
  character(len=8), parameter  ::  tracer_unit(nsv_all) = (/ &
                              'kg/kg   ', &   ! qr
                              '-       ', &   ! nr
                              'kg/kg   ', &   ! co2
                              'kg/kg   ', &   ! co2veg
                              'kg/kg   ', &   ! co2sum
                              'kg/kg   ', &   ! nox
                              'kg/kg   ', &   ! nh3
                              'kg/kg   ', &   ! so2
                              'kg/kg   ', &   ! INERT
                              'kg/kg   ', &   ! O3
                              'kg/kg   ', &   ! O1D
                              'kg/kg   ', &   ! NO
                              'kg/kg   ', &   ! NO2
                              'kg/kg   ', &   ! CH4
                              'kg/kg   ', &   ! CH2O
                              'kg/kg   ', &   ! CH3O2
                              'kg/kg   ', &   ! MVK
                              'kg/kg   ', &   ! ISO
                              'kg/kg   ', &   ! RO2
                              'kg/kg   ', &   ! OH
                              'kg/kg   ', &   ! HO2
                              'kg/kg   ', &   ! CO
                              'kg/kg   ', &   ! H2O
                              'kg/kg   ', &   ! PRODUC
                              'kg/kg   ', &   ! O2
                              'kg/kg   ', &   ! N2
                              'kg/kg   ', &   ! HNO3
                              'kg/kg   ', &   ! H2O2
                              'kg/kg   ', &   ! NO3
                              'kg/kg   '/)    ! N2O5

  ! flag to check if tracer has "emission" property:
  logical, parameter  ::  tracer_is_emitted(nsv_all) = (/ &
                            .false.  , &   ! qr
                            .false.  , &   ! nr
                            .true.   , &   ! co2
                            .false.  , &   ! co2veg
                            .false.  , &   ! co2sum
                            .true.   , &   ! nox
                            .true.   , &   ! nh3
                            .true.   , &   ! so2
                            .false.  , &   ! INERT
                            .false.  , &   ! O3
                            .false.  , &   ! O1D
                            .true.   , &   ! NO
                            .true.   , &   ! NO2
                            .true.   , &   ! CH4
                            .false.  , &   ! CH2O
                            .false.  , &   ! CH3O2
                            .false.  , &   ! MVK
                            .false.  , &   ! ISO
                            .false.  , &   ! RO2
                            .false.  , &   ! OH
                            .false.  , &   ! HO2
                            .true.   , &   ! CO
                            .false.  , &   ! H2O
                            .false.  , &   ! PRODUC
                            .false.  , &   ! O2
                            .false.  , &   ! N2
                            .false.  , &   ! HNO3
                            .false.  , &   ! H2O2
                            .false.  , &   ! NO3
                            .false.  /)    ! N2O5

  ! flag to check if tracer has "reactivity" property (i.e. it is part of the chemical mechanism):
  logical, parameter  ::  tracer_is_reactive(nsv_all) = (/ &
                            .false. , &   ! qr
                            .false. , &   ! nr
                            .false. , &   ! co2
                            .false. , &   ! co2veg
                            .false. , &   ! co2sum
                            .false. , &   ! nox
                            .false. , &   ! nh3
                            .false. , &   ! so2
                            .true.  , &   ! INERT
                            .true.  , &   ! O3
                            .true.  , &   ! O1D
                            .true.  , &   ! NO
                            .true.  , &   ! NO2
                            .true.  , &   ! CH4
                            .true.  , &   ! CH2O
                            .true.  , &   ! CH3O2
                            .true.  , &   ! MVK
                            .true.  , &   ! ISO
                            .true.  , &   ! RO2
                            .true.  , &   ! OH
                            .true.  , &   ! HO2
                            .true.  , &   ! CO
                            .true.  , &   ! H2O
                            .true.  , &   ! PRODUC
                            .true.  , &   ! O2
                            .true.  , &   ! N2
                            .true.  , &   ! HNO3
                            .true.  , &   ! H2O2
                            .true.  , &   ! NO3
                            .true.  /)    ! N2O5

  ! flag to check if tracer has "deposition" property:
  logical, parameter  ::  tracer_is_deposited(nsv_all) = (/ &
                            .false. , &   ! qr
                            .false. , &   ! nr
                            .false. , &   ! co2
                            .false. , &   ! co2veg
                            .false. , &   ! co2sum
                            .true.  , &   ! nox
                            .true.  , &   ! nh3
                            .true.  , &   ! so2
                            .false. , &   ! INERT
                            .true.  , &   ! O3
                            .false. , &   ! O1D
                            .true.  , &   ! NO
                            .true.  , &   ! NO2
                            .false. , &   ! CH4
                            .false.  , &   ! CH2O ! todo: not available in DEPAC yet
                            .false.  , &   ! CH3O2 ! todo: not available in DEPAC yet
                            .false.  , &   ! MVK ! todo: not available in DEPAC yet
                            .false. , &   ! ISO
                            .false. , &   ! RO2
                            .false. , &   ! OH
                            .false. , &   ! HO2
                            .true.  , &   ! CO
                            .false. , &   ! H2O
                            .false. , &   ! PRODUC
                            .false. , &   ! O2
                            .false. , &   ! N2
                            .true.  , &   ! HNO3
                            .true.  , &   ! H2O2
                            .false. , &   ! NO3
                            .true.  /)    ! N2O5

  ! flag to check if tracer has "photosynthesis" property, i.e. it is included in A-gs:
  logical, parameter  ::  tracer_is_photosynthesized(nsv_all) = (/ &
                            .false. , &   ! qr
                            .false. , &   ! nr
                            .true.  , &   ! co2
                            .true.  , &   ! co2veg
                            .false. , &   ! co2sum
                            .false. , &   ! nox
                            .false. , &   ! nh3
                            .false. , &   ! so2
                            .false. , &   ! INERT
                            .false. , &   ! O3
                            .false. , &   ! O1D
                            .false. , &   ! NO
                            .false. , &   ! NO2
                            .false. , &   ! CH4
                            .false. , &   ! CH2O
                            .false. , &   ! CH3O2
                            .false. , &   ! MVK
                            .false. , &   ! ISO
                            .false. , &   ! RO2
                            .false. , &   ! OH
                            .false. , &   ! HO2
                            .false. , &   ! CO
                            .false. , &   ! H2O
                            .false. , &   ! PRODUC
                            .false. , &   ! O2
                            .false. , &   ! N2
                            .false. , &   ! HNO3
                            .false. , &   ! H2O2
                            .false. , &   ! NO3
                            .false. /)    ! N2O5

  ! flag to check if tracer has "micro physics" property:
  logical, parameter  ::  tracer_is_microphys(nsv_all) = (/ &
                            .true.  , &   ! qr
                            .true.  , &   ! nr
                            .false. , &   ! co2
                            .false. , &   ! co2veg
                            .false. , &   ! co2sum
                            .false. , &   ! nox
                            .false. , &   ! INERT
                            .false. , &   ! nh3
                            .false. , &   ! so2
                            .false. , &   ! O3
                            .false. , &   ! O1D
                            .false. , &   ! NO
                            .false. , &   ! NO2
                            .false. , &   ! CH4
                            .false. , &   ! CH2O
                            .false. , &   ! CH3O2
                            .false. , &   ! MVK
                            .false. , &   ! ISO
                            .false. , &   ! RO2
                            .false. , &   ! OH
                            .false. , &   ! HO2
                            .false. , &   ! CO
                            .false. , &   ! H2O
                            .false. , &   ! PRODUC
                            .false. , &   ! O2
                            .false. , &   ! N2
                            .false. , &   ! HNO3
                            .false. , &   ! H2O2
                            .false. , &   ! NO3
                            .false. /)    ! N2O5

end module modtracdata