# DALES 4.1 with explicit aerosol

The following document describes the installation of the DALES model 
version 4.1 as used in de Bruine et al. (TBA), doi:TBA.

## Included files in directory:

### Runfiles
- [job.BASE](job.BASE)
- [job.AERO](job.AERO)
- [namoptions.BASE](namoptions.BASE)
- [namoptions.AERO](namoptions.AERO)
- [scalar.inp.AERO](scalar.inp.AERO)

### Source code
- [modbulkmicro.f90](modbulkmicro.f90)
- [modbulkmicrostat.f90](modbulkmicrostat.f90)
- [modcrosssection.f90](modcrosssection.f90)
- [modfielddump.f90](modfielddump.f90)
- [modfields.f90](modfields.f90)
- [modmicrodata.f90](modmicrodata.f90)
- [modmicrophysics.f90](modmicrophysics.f90)
- [modsamptend.f90](modsamptend.f90)
- [tstep.f90](tstep.f90)

### Aerosol initialisation
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
Two simulations use the base version of DALES4.1 and only require an 
updated version of the job script and namelist files, with suffixes 
replaced by the number specified in the job script.
- [job.BASE](job.BASE)
- [namoptions.BASE](namoptions.BASE)

BASE   : No additional instructions  
BASE30 : Set Nc_0 to 30e6  

The remaining 4 simulations require other job scripts and namelists,
as well as new source files and files for initialisation of aerosol 
specified above. Suffixes **AERO** should again be replaced by a number 
specified in the job script. 
- [job.AERO](jon.AERO)
- [namoptions.AERO](namoptions.AERO)
- [scalar.inp.AERO](scalar.inp.AERO)

KAPPA  : No additional instructions  
SAT02  : Set Ssat = 0.2  
SAT10  : Set Ssat = 1.0  
PN     : Set l_kohler = False  

## New namelist options
- l_kohler (default: True)  
Switch for activation parameterization:  
True   Kappa-Kohler based activation described in de Bruine et al.  
False  Updraft based activation described in Pousse-Nottelmann, 2015  

- Ssat (default = 0.4)  
Value for supersaturation used in Kappa-Kohler based activation, value in %.  

