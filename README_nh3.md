# dales_nh3
 The addition of ammonia to DALES
 
 This README file notes the changes made to the code
 
- _v08_final: 
	Update data on 23-05-2025. 
	Continuation based on "_v07c_dales-master_Snellius-CMakeList". 
	This version inclused the updates from _v07c_dales-master_Snellius-CMakeList, which 
	are not listed in this document. These updates make sure that DALES will run on 
	the new Snellius supercomputer. 
	The only other updates are:
	- Added matlab scripts to process and analyze the data
	- Added example input files for a meteo spin-up run and an analysis run
 
 - _v07b_dales-master_reduce-output: 
	Update data on 23-04-2021. 
	Continuation based on "_v07_dales-master_back-to-v4-2". 
	For larger runs, the output resulting from namcrosssection & namfielddump becomes 
	very large. This version removes several meteorological variables from these outputs,
	keeping the scalar concentrations, w, thl, qt & e120.
	- modcrosssection: remove multiple variables
		* Line 43:		Decreased nvar: 	nvar = 12 --> nvar = 4
		* Line 155:  	remove variables from xz cross-section output initialization
		* Line 186:  	remove variables from xy cross-section output initialization
		* Line 215:  	remove variables from yz cross-section output initialization
		* Line 334:  	remove variables from xz cross-section output writing
		* Line 437:  	remove variables from xy cross-section output writing
		* Line 537:  	remove variables from yz cross-section output writing
	- modfielddump: remove multiple variables
		* Line 40:		Decreased nvar: 	nvar = 7 --> nvar = 3
		* Line 114:  	remove variables from the fielddump output initialization
		* line 171 to 300	Remove several variables from fielddump output writing 
	
- _v07_dales-master_back-to-v4-2:
	Update data on 30-03-2021.
	Continuation based on "_v05_dales-master_perc_chem". 
	Going back from DALES v4.3 to DALES v4.2. 
	!!! IMPORTANT !!! 
		Remove the "rssoilmin" variable from surface.interactive.inp
	!!! IMPORTANT !!! 
	- Upgrading from DALES v4.2 to DALES v4.3 comes with multiple unexpected issues which 
	  appear not to be easily fixed. The issues listed below are the reasons we revert 
	  back to v4.2: 
		* lcoriol = .false. but the coriolis force still acts upon the domain.
		* Meteorological parameters, mainly the CBL height (based on potent. temp.), are 
		  different when using non-equidistant vertical grids, when compared to equidistant 
		  vertical grids. This was tested using testruns aerosolrad 410, 411, 412 & 413. No 
		  reason for the differences were found. 
	- Updates on bugfixes: 
		* "thls bug" (fixed in version _v02): 
			!thls = 300 --> line is removed
		* "rssoilmin(i,j) = rsmin_patch bug" (fixed in version _v02):
			Newly added input rssoilmin to surface.interactive.inp is removed again.
			New fix: rssoilmin(i,j) = rssoilminav. The new input variable made this 
			DALES version incompatible with other versions. This way, the 
			surface.interactive.inp input file has the same input variables again.
			Changes per line:
				* modsurface.f90
				line 275: remove "allocate(rssoilmin_patch(xpatches,ypatches))"
				line 304: remove "rssoilmin_patch  = -1"
				line 363: remove ", rssoilmin_land(i)" from line
				line 259: remove "rssoilminav  = 0"
				line 487: remove "rssoilmin_patch(i,j)  = rssoilmin_land(landindex)"
				line 508: remove  line "rssoilminav   = rssoilminav  +  ..."
				line 1670: Add use modsurfdata, only : rssoilminav	
				line 1758: change "rssoilmin  (i,j) = rssoilmin_patch (tempx,tempy)" into 
					"rssoilmin  (i,j) = rssoilminav" 
				* modsurfdata
				line 287: remove "real, allocatable :: rssoilmin_patch(:,:)..."

- _v06_dales-master_dales-v4-3:
	Update data on 25-03-2021. 
	Continuation based on "_v05_dales-master_perc_chem". 
	This is an update from DALES v4.2 to DALES v4.3. 
	All of the above changes are applied to the source code of DALES v4.3. 
	- The main reason for the update is that non-equidistant vertical grids are not supported 
	  for the kappa advection scheme (nr. 7) in v4.2, but it is supported in v4.3 This 
	  allows us to save costs on the research runs. 
	- modsurface.f90: Changed the bugfix for thls, mentioned in v02 of this branch, to a 
	  cleaner fix (line 475). 

- _v05_dales-master_perc_chem: 
	Update data on 10-02-2021. 
	Continuation based on "_v04_dales-master_split-flux". 
	Add the option for "percentage-chemistry functionality", where specific scalars 
	have an additional loss/source term in a percentage of the concentration per hour,
	representing simplified chemistry. 
	- modglobal.f90: Add new (NAMOPTIONS input) variables to DALES. 
		* Added: lprec_chem, pc_chemrate. 
	- modfields.f90: allocate pc_chemrate to determine the size of the variable of length "nsv". 
		* allocate pc_chemrate to set the size of the variable to be "nsv". 
		* Set default value of pc_chemrate to  be 0. 
		* deallocate pc_chemrate. 
	- modstartup.f90: Read the new variables under the &DYNAMICS list. 
		* Added: lprec_chem, pc_chemrate. 
	- tstep.f90: Add the percentage-chemistry functionallity to scalar concentration. 
		* Added equation svp(:,:,:,n) = svp(:,:,:,n) + svm(:,:,:,n) * pc_chemrate(n)/3600. 
		
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

- _v01_dales-master_nh3_ruben:
	Update date: 31-07-2020. 
	Continuation based on "_v00_dales-master". 
	Includes NH3 and additional outputs. 
	
- _v00_dales-master: 
	Original DALES code. 
	* BUGFIX: modradiation.f90. 
		line 253:	use modsurfdata,  only : tauField. 
					--> use modsurfdata,  only : tauField, lrsAgs. 
		line 342: 	tauField(i,j) = tauc. 
					--> if (lrAgs) tauField(i,j) = tauc. 
	
	