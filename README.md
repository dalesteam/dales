# dales_nh3
 The addition of ammonia to DALES
 
 This README file notes the changes made to the code
 
 - _v00_dales-master: 
	Original DALES code. 
	* BUGFIX: modradiation.f90. 
		line 253:	use modsurfdata,  only : tauField. 
					--> use modsurfdata,  only : tauField, lrsAgs. 
		line 342: 	tauField(i,j) = tauc. 
					--> if (lrAgs) tauField(i,j) = tauc. 

- _v01_dales-master_nh3_ruben:
	Update date: 31-07-2020. 
	Continuation based on "_v00_dales-master". 
	Includes NH3 and additional outputs. 

- _v02_dales-master_hetero-sfc:
	update date: 08-10-2020.
	Continuation based on "_v01_dales-master_nh3_ruben".
	Adds option for surface heterogeneity. 
	* BUGFIX modsurface.f90 & modsurfdata.f90. 
		when lhetero = true: rssoilmin = rsmin (line 1580) was incorrect and is fixed. 
		Now separate rssoilmin input is added to surface.interactive.inp. 
	* BUGFIX: modsurface.f90. 
		Line 432: 	thls         = 300. 
					---> !thls         = 300. 
	* The bugfix in dales-master will show up here in Github. 
	- UPDATE ON 26-02-2021! 
	  Added lpatchoutput which makes it possible to NOT write patch specific outputs
	  files for lheteor = .true., like "tmser1patchiiixjjj.[cexpnr]", which is written 
	  for each patch and can significantly slow down the model when many patches are
	  used.
		* modsurfdata.f90: 	Added lpatchoutput = .true.
		* modsurface.f90:	Added lpatchoutput as an input in NAMOPTIONS for the 
							NAMSURFACE section. 
		* modtimestat.f90: 	Added lpatchoutput as an additional condition for writing
							the patch-specific output files. 
	
- _v03_dales-master_non-periodic BCs:
	Update data: 07-01-2021. 
	Continuation based on "_v02_dales-master_hetero-sfc". 
	- Makes the boundary conditions (BCs) at x = 0 and x = itot non-periodic. 
		Changed subroutine "cyclich" in modboundary.f90. 
		Force the ghost-cells to be zero. 
	- Add two new input files to namoptions, under &DYNAMICS:
		* lnonperiodbc_sv: 	Array (of length 1000) with logicals for the non-periodic 
							BC scalars. 
		*** Changes are made to modglobal.f90, modstartup.f90 & modboundary.f90.
			Search for (ctrl+F) "Ruben" to find the changes. 
	- Increased the maximum of patches for lhetero. 
		Changed mpatches from 16 to 1000 in modsurfdata. 
		
- _v04_dales-master_split-flux: 
	Update data: 26-01-2021. 
	Continuation based on "_v03_dales-master_non-periodic BCs". 
	Add the option for "split-flux functionality", where the prescribed constant surface
	fluxes are split between multiple scalars. 
	- modsurfdata.f90: Add new (NAMOPTIONS input) variables to DALES:
		* Added: lpartialflux, sf_dim1, sf_dim2, sf_scalars. 
	- modsurface.f90: subroutine Surface. 
		* Call variables svm, lpartialflux, sf_dim1, sf_dim2, sf_scalars. 
		  Start at line 709. 
		* Define new variables: sf_counter, sf_i, sf_j, sf_svm, sf_flux. 
		  Start at line 733. 
		* Implement new script to allow for calculation of splitted fluxes. 
			+ Documentation in "DALES split-flux instructions.docx". 
			+ Start at line 870. 
