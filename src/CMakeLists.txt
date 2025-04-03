# DALES source files
set( dales_srcs
	advec_2nd.f90
	advec_5th.f90
	advec_6th.f90
	advec_hybrid.f90
	advec_hybrid_f.f90
	advec_kappa.f90
	advec_upw.f90
	bulkmicro_kk.f90
	bulkmicro_sb.f90
	daleslib.f90
	fftnew.f90
	go.f90
	le_drydepos_gas_depac.f90
	modAGScross.f90
	modadvection.f90
	modboundary.f90
	modbudget.f90
	modbulkmicro.f90
	modbulkmicro3.f90
	modbulkmicro3_column.f90
	modbulkmicro3_point.f90
	modbulkmicrostat.f90
	modbulkmicrostat3.f90
	modcanopy.f90
	modcape.f90
	modchecksim.f90
	modchem.f90
	modcloudfield.f90
	modcrosssection.f90
	modcufft.f90
	moddatetime.f90
	moddepcrosssection.f90
	moddrydeposition.f90
	modemisdata.f90
	modemission.f90
	modfft2d.f90
	modfftw.f90
	modfielddump.f90
	modfields.f90
	modforces.f90
	modgenstat.f90
	modglobal.f90
	modgpu.f90
	modgpumpiinterface.f90
	modheterostats.f90
	modhypre.f90
	modlsm.f90
	modlsmcrosssection.f90
	modlsmdata.f90
	modlsmstat.f90
	modmicrodata.f90
	modmicrodata3.f90
	modmicrophysics.f90
	modmpi.f90
	modmpiinterface.f90
	modmsebudg.f90
	modnetcdf.f90
	modnudge.f90
	modnudgeboundary.f90
	modnvtx.f90
	modopenboundary.f90
	modpois.f90
	modprecision.f90
	modquadrant.f90
	modraddata.f90
	modradfield.f90
	modradfull.f90
	modradiation.f90
	modradrrtmg.f90
	modradrte_rrtmgp.f90
	modradstat.f90
	modsampdata.f90
	modsampling.f90
	modsamptend.f90
	modsimpleice.f90
	modsimpleice2.f90
	modstartup.f90
	modstat_nc.f90
	modstattend.f90
	modsubgrid.f90
	modsubgriddata.f90
	modsurface.f90
	modsurfdata.f90
	modsynturb.f90
	modtestbed.f90
	modthermodynamics.f90
	modtimedep.f90
	modtimedepsv.f90
	modtimer.f90
	modtimestat.f90
	modtracers.f90
	moduser.f90
	modvarbudget.f90
	program.f90
	rad_rndnmb.f90
	shr_orb_mod.f90
	tstep.f90
	utils.f90
)

# RRTMG - longwave
set( rrtmg_lw_srcs
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/parkind.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/parrrtm.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_cld.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_con.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_kg01.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_kg02.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_kg03.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_kg04.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_kg05.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_kg06.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_kg07.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_kg08.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_kg09.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_kg10.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_kg11.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_kg12.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_kg13.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_kg14.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_kg15.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_kg16.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_ncpar.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_ref.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_tbl.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_vsn.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/modules/rrlw_wvn.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/src/rrtmg_lw_cldprop.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/src/rrtmg_lw_init.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/src/rrtmg_lw_rad.nomcica.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/src/rrtmg_lw_read_nc.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/src/rrtmg_lw_rtrn.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/src/rrtmg_lw_rtrnmc.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/src/rrtmg_lw_rtrnmr.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/src/rrtmg_lw_setcoef.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/src/rrtmg_lw_taumol.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/src/src_link
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_LW/src/util_link
)

# RRTMG - shortwave
set( rrtmg_sw_srcs
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/parrrsw.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_aer.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_cld.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_con.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_kg16.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_kg17.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_kg18.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_kg19.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_kg20.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_kg21.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_kg22.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_kg23.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_kg24.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_kg25.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_kg26.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_kg27.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_kg28.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_kg29.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_ncpar.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_ref.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_tbl.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_vsn.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/modules/rrsw_wvn.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/src/rrtmg_sw_cldprop.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/src/rrtmg_sw_init.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/src/rrtmg_sw_rad.nomcica.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/src/rrtmg_sw_read_nc.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/src/rrtmg_sw_reftra.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/src/rrtmg_sw_setcoef.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/src/rrtmg_sw_spcvrt.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/src/rrtmg_sw_taumol.f90
	${PROJECT_SOURCE_DIR}/external/RRTMG/RRTMG_SW/src/rrtmg_sw_vrtqdr.f90
)

# RTE-RRTMGP
set( rrtmgp_srcs
	# RTE
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/rte-kernels/accel/mo_optical_props_kernels.F90
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/rte-kernels/accel/mo_rte_solver_kernels.F90
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/rte-kernels/mo_rte_util_array.F90
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/rte-kernels/mo_fluxes_broadband_kernels.F90
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/rte-kernels/mo_rte_kind.F90
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/rte-frontend/mo_fluxes.F90
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/rte-frontend/mo_optical_props.F90
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/rte-frontend/mo_rte_config.F90
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/rte-frontend/mo_rte_lw.F90
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/rte-frontend/mo_rte_sw.F90
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/rte-frontend/mo_source_functions.F90
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/rte-frontend/mo_rte_util_array_validation.F90
	# RRTMGP
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/rrtmgp-kernels/accel/mo_gas_optics_rrtmgp_kernels.F90
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/rrtmgp-kernels/mo_cloud_optics_rrtmgp_kernels.F90
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/rrtmgp-frontend/mo_gas_optics_rrtmgp.F90
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/rrtmgp-frontend/mo_cloud_optics_rrtmgp.F90
	# Gas optics
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/gas-optics/mo_gas_concentrations.F90
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/gas-optics/mo_gas_optics.F90
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/gas-optics/mo_gas_optics_constants.F90
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/gas-optics/mo_gas_optics_util_string.F90
	# Examples (why are these needed?)
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/examples/all-sky/mo_load_cloud_coefficients.F90
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/examples/mo_load_coefficients.F90
	${PROJECT_SOURCE_DIR}/external/RTE-RRTMGP/examples/mo_simple_netcdf.F90
)

if( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" )
	set_source_files_properties( ${rrtmg_lw_srcs} PROPERTIES COMPILE_FLAGS -Wno-error=implicit-interface )
	set_source_files_properties( ${rrtmg_sw_srcs} PROPERTIES COMPILE_FLAGS -Wno-error=implicit-interface )
endif()

# Use git-version.cmake to create modversion.f90, containing a version string from git, e.g. "4.2-34-g62b85a-dirty"
add_custom_target( tag_git_version ALL
  COMMAND ${CMAKE_COMMAND} -D TARGET_DIR=${CMAKE_CURRENT_BINARY_DIR} -P ${CMAKE_CURRENT_SOURCE_DIR}/git-version.cmake
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/modversion.f90
)

# Temporary, because we use CUDA-Fortran in these places
if( ENABLE_ACC )
  set_source_files_properties( modmpi.f90 modgpumpiinterface.f90 PROPERTIES COMPILE_OPTIONS -cuda )
endif()

# Stand-alone DALES program
ecbuild_add_executable( TARGET dales
                        SOURCES
                          ${dales_srcs}
                          $<$<BOOL:${ENABLE_RRTMG}>:${rrtmg_lw_srcs}>
                          $<$<BOOL:${ENABLE_RRTMG}>:${rrtmg_sw_srcs}>
                          ${rrtmgp_srcs}
                          ${CMAKE_CURRENT_BINARY_DIR}/modversion.f90
                        DEPENDS tag_git_version
                        INCLUDES
                          ${NetCDF_Fortran_INCLUDE_DIR}
                          $<$<BOOL:${ENABLE_FFTW}>:${FFTW_INCLUDE_DIRS}>
                        LIBS
                          NetCDF::NetCDF_Fortran
                          NetCDF::NetCDF_C
                          $<$<BOOL:${ENABLE_ACC}>:OpenACC::OpenACC_Fortran>
                          $<$<BOOL:${ENABLE_FFTW}>:FFTW::fftw3>
                          $<$<BOOL:${ENABLE_FFTW}>:FFTW::fftw3f>
                          $<$<BOOL:${ENABLE_HYPRE}>:${HYPRE_LIB}> 
                          $<$<BOOL:${ENABLE_ACC}>:CUDA::cufft> 
                          $<$<BOOL:${ENABLE_NVTX}>:CUDA::nvToolsExt>
)
# Temporary
if( ENABLE_ACC )
  target_link_options( dales PUBLIC -cuda )
endif()

# DALES library, e.g. for use with OMUSE. Only built if -DENABLE_DALESLIB is turned on.
ecbuild_add_library( TARGET libdales
                     SOURCES
                       ${dales_srcs}
                       ${rrtmg_lw_srcs}
                       ${rrtmg_sw_srcs}
                       ${rrtmgp_srcs}
                       ${CMAKE_CURRENT_BINARY_DIR}/modversion.f90
                     DEPENDS tag_git_version
                     PUBLIC_INCLUDES
                       ${NetCDF_Fortran_INCLUDE_DIR}
                       $<$<BOOL:${ENABLE_FFTW}>:${FFTW_INCLUDE_DIRS}>
                     PUBLIC_LIBS
                       NetCDF::NetCDF_Fortran
                       $<$<BOOL:${ENABLE_FFTW}>:FFTW::fftw3>
                       $<$<BOOL:${ENABLE_FFTW}>:FFTW::fftw3f>
                       $<$<BOOL:${ENABLE_HYPRE}>:${HYPRE_LIB}>
                     CONDITION ENABLE_DALESLIB
)
