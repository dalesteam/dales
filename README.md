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
