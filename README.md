# DALES 4.1 with explicit aerosol

The following document describes the installation of the DALES (Dutch Atmospheric Large Eddy Simulation) model
version 4.1 as used in de Bruine et al. (TBA), doi:TBA. 

The main repository can be found at github.com/dalesteam/dales

## Files specific for this branch:

### Source code
- [modbulkmicro.f90](src/modbulkmicro.f90)
- [modbulkmicrostat.f90](src/modbulkmicrostat.f90)
- [modcrosssection.f90](src/modcrosssection.f90)
- [modfielddump.f90](src/modfielddump.f90)
- [modfields.f90](src/modfields.f90)
- [modmicrodata.f90](src/modmicrodata.f90)
- [modmicrophysics.f90](src/modmicrophysics.f90)
- [modsamptend.f90](src/modsamptend.f90)
- [tstep.f90](src/tstep.f90)

### Aerosol initialisation
- [scalar.inp.AERO](scalar.inp.AERO)
- [CDNC1_cloudscavenging_005.dat](CDNC1_cloudscavenging_005.dat)
- [CDNC1_cloudscavenging_010.dat](CDNC1_cloudscavenging_010.dat)
- [CDNC1_cloudscavenging_015.dat](CDNC1_cloudscavenging_015.dat)
- [CDNC1_cloudscavenging_020.dat](CDNC1_cloudscavenging_020.dat)
- [CDNC1_cloudscavenging_025.dat](CDNC1_cloudscavenging_025.dat)
- [CDNC1_cloudscavenging_030.dat](CDNC1_cloudscavenging_030.dat)
- [CDNC1_cloudscavenging_035.dat](CDNC1_cloudscavenging_035.dat)
- [CDNC1_cloudscavenging_040.dat](CDNC1_cloudscavenging_040.dat)
- [CDNC1_cloudscavenging_045.dat](CDNC1_cloudscavenging_045.dat)
- [CDNC1_cloudscavenging_050.dat](CDNC1_cloudscavenging_050.dat)

- [CDNC1_ncloudscavenging_005.dat](CDNC1_ncloudscavenging_005.dat)
- [CDNC1_ncloudscavenging_010.dat](CDNC1_ncloudscavenging_010.dat)
- [CDNC1_ncloudscavenging_015.dat](CDNC1_ncloudscavenging_015.dat)
- [CDNC1_ncloudscavenging_020.dat](CDNC1_ncloudscavenging_020.dat)
- [CDNC1_ncloudscavenging_025.dat](CDNC1_ncloudscavenging_025.dat)
- [CDNC1_ncloudscavenging_030.dat](CDNC1_ncloudscavenging_030.dat)
- [CDNC1_ncloudscavenging_035.dat](CDNC1_ncloudscavenging_035.dat)
- [CDNC1_ncloudscavenging_040.dat](CDNC1_ncloudscavenging_040.dat)
- [CDNC1_ncloudscavenging_045.dat](CDNC1_ncloudscavenging_045.dat)
- [CDNC1_ncloudscavenging_050.dat](CDNC1_ncloudscavenging_050.dat)

- [belowcloud_m.dat](belowcloud_m.dat)
- [belowcloud_n.dat](belowcloud_n.dat)
- [rbar_aerosol](rbar_aerosol)
- [rbar_aerosol_belowcloud.dat](rbar_aerosol_belowcloud.dat)

## Installation

How to INSTALL/RUN using git (following the DALES installation manual)

1. git clone https://github.com/dalesteam/dales.git
2. git checkout origin 4.1

Then follow the following instructions:
3. cd dales
4. mkdir build
5. cd build
6. cmake .. -DCASE="rico"
7. make install

## Simulation instructions

1. Make new run directory
2. Copy the following files from the rico case directory (dales/cases/rico)
- lscale.inp.001
- prof.inp.001
- scalar.inp.001

3. Copy executable dales4 from dales/build and run the executable.

## Additional simulation specific instructions
Two simulations use the base version of DALES4.1 and only requires an
updated version of the job script and namelist files, specific setting:

BASE   : Set Nc_0 to 70e6
BASE30 : Set Nc_0 to 30e6

The remaining 4 simulations require other job scripts and namelists,
as well as new source files and files for initialisation of aerosol
specified above. Specific settings:

KAPPA  : Set l_kohler = True and Ssat = 0.2
SAT02  : Set l_kohler = True and Ssat = 0.4
SAT10  : Set l_kohler = True and Ssat = 1.0
PN     : Set l_kohler = False

## New namelist options
- l_kohler (default: True)
Switch for activation parameterization:
True   Kappa-Kohler based activation described in de Bruine et al.
False  Updraft based activation described in Pousse-Nottelmann, 2015

- Ssat (default = 0.4)
Value for supersaturation used in Kappa-Kohler based activation, value in %.
