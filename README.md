# dales_nh3
 The addition of ammonia to DALES
 
 This README file notes the changes made to the code
 
 - _v00_dales-master:
	Original DALES code
	* BUGFIX: modradiation.f90 
		line 253:	use modsurfdata,  only : tauField 
					--> use modsurfdata,  only : tauField, lrsAgs
		line 342: 	tauField(i,j) = tauc 
					--> if (lrAgs) tauField(i,j) = tauc

- _v01_dales-master_nh3_ruben:
	Update date: 31-07-2020
	Continuation based on "_v00_dales-master"
	Includes NH3 and additional outputs

- _v02_dales-master_hetero-sfc
	update date: 08-10-2020
	Continuation based on "_v01_dales-master_nh3_ruben"
	Adds option for surface heterogeneity 
	* BUGFIX modsurface.f90 & modsurfdata.f90 
		when lhetero = true: rssoilmin = rsmin (line 1580) was incorrect and is fixed. 
		Now separate rssoilmin input is added to surface.interactive.inp
	* BUGFIX: modsurface.f90 
		Line 432: 	thls         = 300
					---> !thls         = 300
	* The bugfix in dales-master will show up here in Github
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
