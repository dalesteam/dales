!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_read_nc.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.1 $
!     created:   $Date: 2009/05/22 21:02:13 $
!

!===============================================================================
! rrtmg_lw_read_nc.f90
!
! Description: This program reads all of the RRTM longwave data from a NetCDF
!              file in band by band subroutines, as a replacement for the
!              rrtmg_lw_k_g.f90 data statements.
!
! Written By: Patrick Hofmann
! Last Update: 1/23/2009
!===============================================================================

module rrtmg_lw_read_nc
contains
!*******************************************************************************
subroutine lw_kgb01
    use rrlw_kg01, only : fracrefao, fracrefbo, kao, kbo, kao_mn2, kbo_mn2, &
                        selfrefo, forrefo, no1
    use rrlw_ncpar
    use netcdf

    implicit none
    save

    integer(kind=im) :: ab
    integer(kind=im), parameter :: bandNumber = 1, numGPoints = no1
    integer(kind=im), parameter :: gPointSetNumber = 1
    integer(kind=im) :: ncid, varID,sts(17)

    sts(:)   = nf90_NoErr
    sts(1)   = nf90_open('rrtmg_lw.nc',nf90_nowrite,ncid)

    sts(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
    sts(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))

    sts(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
    sts(5)   = nf90_get_var(ncid, varID, fracrefbo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))

    sts(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
    sts(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,plower,numGPoints,1,1/))

    sts(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
    sts(9)   = nf90_get_var(ncid, varID, kbo, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,pupper,numGPoints,1,1/))

    sts(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
    sts(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))

    sts(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
    sts(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))

    !Get absorber index for N2
    call getAbsorberIndex('N2',ab)
    sts(14)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
    sts(15)  = nf90_get_var(ncid, varID, kao_mn2, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))

    sts(16)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsUpperAtmos",varID)
    sts(17)  = nf90_get_var(ncid, varID, kbo_mn2, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))

    if(any(sts(:) /= nf90_NoErr)) stop  "Error reading variables from file"

    sts(1) = nf90_close(ncid)

end subroutine lw_kgb01
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb02
        use rrlw_kg02, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no2
        use rrlw_ncpar
    use netcdf

        implicit none
        save

       integer(kind=im), parameter :: bandNumber = 2
    integer(kind=im), parameter :: numGPoints = no2
    integer(kind=im), parameter :: gPointSetNumber = 1
    integer(kind=im) :: ncid, varID,sts(13)

    sts(:)   = nf90_NoErr
    sts(1)   = nf90_open('rrtmg_lw.nc',nf90_nowrite,ncid)

    sts(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
    sts(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))

    sts(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
    sts(5)   = nf90_get_var(ncid, varID, fracrefbo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))

    sts(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
    sts(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,plower,numGPoints,1,1/))

    sts(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
    sts(9)   = nf90_get_var(ncid, varID, kbo, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,pupper,numGPoints,1,1/))

    sts(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
    sts(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))

    sts(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
    sts(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))

    if(any(sts(:) /= nf90_NoErr)) stop  "Error reading variables from file"

    sts(1) = nf90_close(ncid)

end subroutine lw_kgb02
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb03
        use rrlw_kg03, only : fracrefao, fracrefbo, kao, kbo, kao_mn2o, &
                          kbo_mn2o, selfrefo, forrefo, no3
    use rrlw_ncpar
    use netcdf

        implicit none
        save

    integer(kind=im) :: ab
        integer(kind=im), parameter :: bandNumber = 3
    integer(kind=im), parameter :: numGPoints = no3
    integer(kind=im), parameter :: gPointSetNumber = 1
    integer(kind=im) :: ncid, varID,sts(17)

    sts(:)   = nf90_NoErr
    sts(1)   = nf90_open('rrtmg_lw.nc',nf90_nowrite,ncid)

    sts(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
    sts(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,keylower,1,1/))

    sts(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
    sts(5)   = nf90_get_var(ncid, varID, fracrefbo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,keyupper,1,1/))

    sts(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
    sts(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/keylower,Tdiff,plower,numGPoints,1,1/))

    sts(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
    sts(9)   = nf90_get_var(ncid, varID, kbo, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/keyupper,Tdiff,pupper,numGPoints,1,1/))

    sts(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
    sts(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))

    sts(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
    sts(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))

    !Get absorber index for N2
    call getAbsorberIndex('N2O',ab)
    sts(14)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
    sts(15)  = nf90_get_var(ncid, varID, kao_mn2o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/keylower,T,numGPoints,1,1,1/))

    sts(16)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsUpperAtmos",varID)
    sts(17)  = nf90_get_var(ncid, varID, kbo_mn2o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/keyupper,T,numGPoints,1,1,1/))

    if(any(sts(:) /= nf90_NoErr)) stop  "Error reading variables from file"

    sts(1) = nf90_close(ncid)

end subroutine lw_kgb03
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb04
        use rrlw_kg04, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no4
    use rrlw_ncpar
    use netcdf

        implicit none
        save

        integer(kind=im), parameter :: bandNumber = 4
    integer(kind=im), parameter :: numGPoints = no4
    integer(kind=im), parameter :: gPointSetNumber = 1
    integer(kind=im) :: ncid, varID,sts(13)

    sts(:)   = nf90_NoErr
    sts(1)   = nf90_open('rrtmg_lw.nc',nf90_nowrite,ncid)

    sts(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
    sts(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,keylower,1,1/))

    sts(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
    sts(5)   = nf90_get_var(ncid, varID, fracrefbo(:,1:5), &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,keyupper,1,1/))

    sts(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
    sts(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/keylower,Tdiff,plower,numGPoints,1,1/))

    sts(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
    sts(9)   = nf90_get_var(ncid, varID, kbo, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/keyupper,Tdiff,pupper,numGPoints,1,1/))

    sts(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
    sts(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))

    sts(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
    sts(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))

    if(any(sts(:) /= nf90_NoErr)) stop  "Error reading variables from file"

    sts(1) = nf90_close(ncid)

end subroutine lw_kgb04
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb05
        use rrlw_kg05, only : fracrefao, fracrefbo, kao, kbo, kao_mo3, &
                          selfrefo, forrefo, ccl4o, no5
    use rrlw_ncpar
    use netcdf

        implicit none
        save

        integer(kind=im) :: ab
        integer(kind=im), parameter :: bandNumber = 5
    integer(kind=im), parameter :: numGPoints = no5
    integer(kind=im), parameter :: gPointSetNumber = 1
    integer(kind=im) :: ncid, varID,sts(17)

    sts(:)   = nf90_NoErr
    sts(1)   = nf90_open('rrtmg_lw.nc',nf90_nowrite,ncid)

    sts(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
    sts(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,keylower,1,1/))

    sts(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
    sts(5)   = nf90_get_var(ncid, varID, fracrefbo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,keyupper,1,1/))

    sts(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
    sts(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/keylower,Tdiff,plower,numGPoints,1,1/))

    sts(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
    sts(9)   = nf90_get_var(ncid, varID, kbo, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/keyupper,Tdiff,pupper,numGPoints,1,1/))

    sts(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
    sts(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))

    sts(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
    sts(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))

    !Get absorber index for O3
    call getAbsorberIndex('O3',ab)
    sts(14)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
    sts(15)  = nf90_get_var(ncid, varID, kao_mo3, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/keylower,T,numGPoints,1,1,1/))

    !Get absorber index for CCL4
    call getAbsorberIndex('CCL4',ab)
    sts(16)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
    sts(17)  = nf90_get_var(ncid, varID, ccl4o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,1,numGPoints,1,1,1/))

    if(any(sts(:) /= nf90_NoErr)) stop  "Error reading variables from file"

    sts(1) = nf90_close(ncid)

end subroutine lw_kgb05
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb06
        use rrlw_kg06, only : fracrefao, kao, kao_mco2, selfrefo, forrefo, &
                          cfc11adjo, cfc12o, no6
    use rrlw_ncpar
    use netcdf

        implicit none
        save

    integer(kind=im) :: ab
        integer(kind=im), parameter :: bandNumber = 6
    integer(kind=im), parameter :: numGPoints = no6
    integer(kind=im), parameter :: gPointSetNumber = 1
    integer(kind=im) :: ncid, varID,sts(17)

    sts(:)   = nf90_NoErr
    sts(1)   = nf90_open('rrtmg_lw.nc',nf90_nowrite,ncid)

    sts(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
    sts(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))

    sts(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
    sts(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,plower,numGPoints,1,1/))

    sts(8)   = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
    sts(9)   = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))

    sts(10)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
    sts(11)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))

    !Get absorber index for CO2
    call getAbsorberIndex('CO2',ab)
    sts(12)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
    sts(13)  = nf90_get_var(ncid, varID, kao_mco2, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))

    !Get absorber index for CFC11
    call getAbsorberIndex('CFC11',ab)
    sts(14)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
    sts(15)  = nf90_get_var(ncid, varID, cfc11adjo, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,1,numGPoints,1,1,1/))

    !Get absorber index for CFC12
    call getAbsorberIndex('CFC12',ab)
    sts(16)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
    sts(17)  = nf90_get_var(ncid, varID, cfc12o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,1,numGPoints,1,1,1/))

    if(any(sts(:) /= nf90_NoErr)) stop  "Error reading variables from file"

    sts(1) = nf90_close(ncid)

end subroutine lw_kgb06
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb07
        use rrlw_kg07, only : fracrefao, fracrefbo, kao, kbo, kao_mco2, &
                          kbo_mco2, selfrefo, forrefo, no7
    use rrlw_ncpar
    use netcdf

        implicit none
        save

        integer(kind=im) :: ab
       integer(kind=im), parameter :: bandNumber = 7
    integer(kind=im), parameter :: numGPoints = no7
    integer(kind=im), parameter :: gPointSetNumber = 1
    integer(kind=im) :: ncid, varID,sts(17)

    sts(:)   = nf90_NoErr
    sts(1)   = nf90_open('rrtmg_lw.nc',nf90_nowrite,ncid)

    sts(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
    sts(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,keylower,1,1/))

    sts(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
    sts(5)   = nf90_get_var(ncid, varID, fracrefbo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))

    sts(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
    sts(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/keylower,Tdiff,plower,numGPoints,1,1/))

    sts(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
    sts(9)   = nf90_get_var(ncid, varID, kbo, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,pupper,numGPoints,1,1/))

    sts(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
    sts(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))

    sts(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
    sts(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))

    !Get absorber index for CO2
    call getAbsorberIndex('CO2',ab)
    sts(14)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
    sts(15)  = nf90_get_var(ncid, varID, kao_mco2, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/keylower,T,numGPoints,1,1,1/))

    sts(16)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsUpperAtmos",varID)
    sts(17)  = nf90_get_var(ncid, varID, kbo_mco2, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))

    if(any(sts(:) /= nf90_NoErr)) stop  "Error reading variables from file"

    sts(1) = nf90_close(ncid)

end subroutine lw_kgb07
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb08
        use rrlw_kg08, only : fracrefao, fracrefbo, kao, kao_mco2, kao_mn2o, &
                          kao_mo3, kbo, kbo_mco2, kbo_mn2o, selfrefo, forrefo, &
                          cfc12o, cfc22adjo, no8
    use rrlw_ncpar
    use netcdf

        implicit none
        save

        integer(kind=im) :: ab
        integer(kind=im), parameter :: bandNumber = 8
    integer(kind=im), parameter :: numGPoints = no8
        integer(kind=im), parameter :: gPointSetNumber = 1
    integer(kind=im) :: ncid, varID,sts(27)

        sts(:)   = nf90_NoErr
    sts(1)   = nf90_open('rrtmg_lw.nc',nf90_nowrite,ncid)

    sts(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
    sts(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))

    sts(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
    sts(5)   = nf90_get_var(ncid, varID, fracrefbo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))

    sts(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
    sts(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,plower,numGPoints,1,1/))

    sts(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
    sts(9)   = nf90_get_var(ncid, varID, kbo, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,pupper,numGPoints,1,1/))

    sts(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
    sts(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))

    sts(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
    sts(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))

    !Get absorber index for O3
    call getAbsorberIndex('O3',ab)
    sts(14)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
    sts(15)  = nf90_get_var(ncid, varID, kao_mo3, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))

    !Get absorber index for CO2
    call getAbsorberIndex('CO2',ab)
    sts(16)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
    sts(17)  = nf90_get_var(ncid, varID, kao_mco2, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))

    sts(18)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsUpperAtmos",varID)
    sts(19)  = nf90_get_var(ncid, varID, kbo_mco2, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))

    !Get absorber index for N2O
    call getAbsorberIndex('N2O',ab)
    sts(20)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
    sts(21)  = nf90_get_var(ncid, varID, kao_mn2o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))

    sts(22)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsUpperAtmos",varID)
    sts(23)  = nf90_get_var(ncid, varID, kbo_mn2o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))

    !Get absorber index for CFC12
    call getAbsorberIndex('CFC12',ab)
    sts(24)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
    sts(25)  = nf90_get_var(ncid, varID, cfc12o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,1,numGPoints,1,1,1/))

    !Get absorber index for CFC22
    call getAbsorberIndex('CFC22',ab)
    sts(26)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
    sts(27)  = nf90_get_var(ncid, varID, cfc22adjo, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,1,numGPoints,1,1,1/))

    if(any(sts(:) /= nf90_NoErr)) stop  "Error reading variables from file"

    sts(1) = nf90_close(ncid)

end subroutine lw_kgb08
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb09
        use rrlw_kg09, only : fracrefao, fracrefbo, kao, kbo, kao_mn2o, &
                            kbo_mn2o, selfrefo, forrefo, no9
    use rrlw_ncpar
    use netcdf

        implicit none
        save

        integer(kind=im) :: ab
        integer(kind=im), parameter :: bandNumber = 9
    integer(kind=im), parameter :: numGPoints = no9
    integer(kind=im), parameter :: gPointSetNumber = 1
    integer(kind=im) :: ncid, varID,sts(17)

    sts(:)   = nf90_NoErr
    sts(1)   = nf90_open('rrtmg_lw.nc',nf90_nowrite,ncid)

    sts(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
    sts(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,keylower,1,1/))

    sts(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
    sts(5)   = nf90_get_var(ncid, varID, fracrefbo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))

    sts(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
    sts(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/keylower,Tdiff,plower,numGPoints,1,1/))

    sts(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
    sts(9)   = nf90_get_var(ncid, varID, kbo, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,pupper,numGPoints,1,1/))

    sts(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
    sts(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))

    sts(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
    sts(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))

    !Get absorber index for N2O
    call getAbsorberIndex('N2O',ab)
    sts(14)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
    sts(15)  = nf90_get_var(ncid, varID, kao_mn2o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/keylower,T,numGPoints,1,1,1/))

    sts(16)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsUpperAtmos",varID)
    sts(17)  = nf90_get_var(ncid, varID, kbo_mn2o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))

    if(any(sts(:) /= nf90_NoErr)) stop  "Error reading variables from file"

    sts(1) = nf90_close(ncid)

end subroutine lw_kgb09
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb10
        use rrlw_kg10, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no10
    use rrlw_ncpar
    use netcdf

        implicit none
        save

       integer(kind=im), parameter :: bandNumber = 10
    integer(kind=im), parameter :: numGPoints = no10
    integer(kind=im), parameter :: gPointSetNumber = 1
    integer(kind=im) :: ncid, varID,sts(13)

    sts(:)   = nf90_NoErr
    sts(1)   = nf90_open('rrtmg_lw.nc',nf90_nowrite,ncid)

    sts(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
    sts(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))

    sts(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
    sts(5)   = nf90_get_var(ncid, varID, fracrefbo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))

    sts(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
    sts(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,plower,numGPoints,1,1/))

    sts(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
    sts(9)   = nf90_get_var(ncid, varID, kbo, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,pupper,numGPoints,1,1/))

    sts(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
    sts(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))

    sts(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
    sts(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))

    if(any(sts(:) /= nf90_NoErr)) stop  "Error reading variables from file"

    sts(1) = nf90_close(ncid)

end subroutine lw_kgb10
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb11
        use rrlw_kg11, only : fracrefao, fracrefbo, kao, kbo, kao_mo2, &
                          kbo_mo2, selfrefo, forrefo, no11
    use rrlw_ncpar
    use netcdf

        implicit none
        save

    integer(kind=im) :: ab
    integer(kind=im), parameter :: bandNumber = 11
    integer(kind=im), parameter :: numGPoints = no11
    integer(kind=im), parameter :: gPointSetNumber = 1
    integer(kind=im) :: ncid, varID,sts(17)

    sts(:)   = nf90_NoErr
    sts(1)   = nf90_open('rrtmg_lw.nc',nf90_nowrite,ncid)

    sts(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
    sts(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))

    sts(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
    sts(5)   = nf90_get_var(ncid, varID, fracrefbo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))

    sts(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
    sts(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,plower,numGPoints,1,1/))

    sts(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
    sts(9)   = nf90_get_var(ncid, varID, kbo, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,pupper,numGPoints,1,1/))

    sts(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
    sts(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))

    sts(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
    sts(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))

    !Get absorber index for O2
    call getAbsorberIndex('O2',ab)
    sts(14)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
    sts(15)  = nf90_get_var(ncid, varID, kao_mo2, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))

    sts(16)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsUpperAtmos",varID)
    sts(17)  = nf90_get_var(ncid, varID, kbo_mo2, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))

    if(any(sts(:) /= nf90_NoErr)) stop  "Error reading variables from file"

    sts(1) = nf90_close(ncid)

end subroutine lw_kgb11
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb12
        use rrlw_kg12, only : fracrefao, kao, selfrefo, forrefo, no12
    use rrlw_ncpar
    use netcdf

        implicit none
        save

    integer(kind=im), parameter :: bandNumber = 12
    integer(kind=im), parameter :: numGPoints = no12
    integer(kind=im), parameter :: gPointSetNumber = 1
    integer(kind=im) :: ncid, varID,sts(9)

    sts(:) = nf90_NoErr
    sts(1) = nf90_open('rrtmg_lw.nc',nf90_nowrite,ncid)

    sts(2) = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
    sts(3) = nf90_get_var(ncid, varID, fracrefao, &
                    start = (/1,1,bandNumber,gPointSetNumber/), &
                    count = (/numGPoints,keylower,1,1/))

    sts(4) = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
    sts(5) = nf90_get_var(ncid, varID, kao, &
                    start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                    count = (/keylower,Tdiff,plower,numGPoints,1,1/))

    sts(6) = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
    sts(7) = nf90_get_var(ncid, varID, selfrefo, &
                    start = (/1,1,bandNumber,gPointSetNumber/), &
                    count = (/Tself,numGPoints,1,1/))

    sts(8) = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
    sts(9) = nf90_get_var(ncid, varID, forrefo, &
                    start = (/1,1,bandNumber,gPointSetNumber/), &
                    count = (/Tforeign,numGPoints,1,1/))

    if(any(sts(:) /= nf90_NoErr)) stop  "Error reading variables from file"

    sts(1) = nf90_close(ncid)

end subroutine lw_kgb12
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb13
        use rrlw_kg13, only : fracrefao, fracrefbo, kao, kao_mco2, kao_mco, &
                          kbo_mo3, selfrefo, forrefo, no13
    use rrlw_ncpar
    use netcdf

        implicit none
        save

        integer(kind=im) :: ab
    integer(kind=im), parameter :: bandNumber = 13
    integer(kind=im), parameter :: numGPoints = no13
    integer(kind=im), parameter :: gPointSetNumber = 1
    integer(kind=im) :: ncid, varID,sts(17)

    sts(:)  = nf90_NoErr
    sts(1)  = nf90_open('rrtmg_lw.nc',nf90_nowrite,ncid)

    sts(2)  = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
    sts(3)  = nf90_get_var(ncid, varID, fracrefao, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,keylower,1,1/))

    sts(4)  = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
    sts(5)  = nf90_get_var(ncid, varID, fracrefbo, &
                     start = (/1,1,bandNumber,gPointSetNumber/),  &
                     count = (/numGPoints,1,1,1/))

    sts(6)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
    sts(7)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/keylower,Tdiff,plower,numGPoints,1,1/))

    sts(8)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
    sts(9)  = nf90_get_var(ncid, varID, selfrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tself,numGPoints,1,1/))

    sts(10) = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
    sts(11) = nf90_get_var(ncid, varID, forrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tforeign,numGPoints,1,1/))

    !Get absorber index for O3
    call getAbsorberIndex('O3',ab)
    sts(12) = nf90_inq_varid(ncid,"AbsorptionCoefficientsUpperAtmos",varID)
    sts(13) = nf90_get_var(ncid, varID, kbo_mo3, &
                     start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                     count = (/1,T,numGPoints,1,1,1/))

    !Get absorber index for CO2
    call getAbsorberIndex('CO2',ab)
    sts(14) = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
    sts(15) = nf90_get_var(ncid, varID, kao_mco2, &
                     start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                     count = (/keylower,T,numGPoints,1,1,1/))

    !Get absorber index for CO
    call getAbsorberIndex('CO',ab)
    sts(16) = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
        sts(17) = nf90_get_var(ncid, varID, kao_mco, &
                     start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                     count = (/keylower,T,numGPoints,1,1,1/))

    if(any(sts(:) /= nf90_NoErr)) stop  "Error reading variables from file"

    sts(1) = nf90_close(ncid)

end subroutine lw_kgb13
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb14
        use rrlw_kg14, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no14
    use rrlw_ncpar
    use netcdf

        implicit none
        save

    integer(kind=im), parameter :: bandNumber = 14
    integer(kind=im), parameter :: numGPoints = no14
    integer(kind=im), parameter :: gPointSetNumber = 1
    integer(kind=im) :: ncid, varID,sts(13)

    sts(:)  = nf90_NoErr
    sts(1)  = nf90_open('rrtmg_lw.nc',nf90_nowrite,ncid)

    sts(2)  = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
    sts(3)  = nf90_get_var(ncid, varID, fracrefao, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,1,1,1/))

    sts(4)  = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
    sts(5)  = nf90_get_var(ncid, varID, fracrefbo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,1,1,1/))

    sts(6)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
    sts(7)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,Tdiff,plower,numGPoints,1,1/))

    sts(8)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
    sts(9)  = nf90_get_var(ncid, varID, kbo, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,Tdiff,pupper,numGPoints,1,1/))

    sts(10) = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
    sts(11) = nf90_get_var(ncid, varID, selfrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tself,numGPoints,1,1/))

    sts(12) = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
    sts(13) = nf90_get_var(ncid, varID, forrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tforeign,numGPoints,1,1/))

    if(any(sts(:) /= nf90_NoErr)) stop  "Error reading variables from file"

    sts(1) = nf90_close(ncid)

end subroutine lw_kgb14
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb15
        use rrlw_kg15, only : fracrefao, kao, kao_mn2, selfrefo, forrefo, no15
    use rrlw_ncpar
    use netcdf

        implicit none
        save

        integer(kind=im) :: ab
     integer(kind=im), parameter :: bandNumber = 15
    integer(kind=im), parameter :: numGPoints = no15
    integer(kind=im), parameter :: gPointSetNumber = 1
    integer(kind=im) :: ncid, varID,sts(11)

    sts(:)  = nf90_NoErr
    sts(1)  = nf90_open('rrtmg_lw.nc',nf90_nowrite,ncid)

    sts(2)  = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
    sts(3)  = nf90_get_var(ncid, varID, fracrefao, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,keylower,1,1/))

    sts(4)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
    sts(5)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/keylower,Tdiff,plower,numGPoints,1,1/))

    sts(6)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
    sts(7)  = nf90_get_var(ncid, varID, selfrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tself,numGPoints,1,1/))

    sts(8)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
    sts(9)  = nf90_get_var(ncid, varID, forrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tforeign,numGPoints,1,1/))

    !Get absorber index for N2
    call getAbsorberIndex('N2',ab)
    sts(10) = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
    sts(11) = nf90_get_var(ncid, varID, kao_mn2, &
                     start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                     count = (/keylower,T,numGPoints,1,1,1/))

    if(any(sts(:) /= nf90_NoErr)) stop  "Error reading variables from file"

    sts(1) = nf90_close(ncid)

end subroutine lw_kgb15
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb16
        use rrlw_kg16, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no16
    use rrlw_ncpar
    use netcdf

        implicit none
        save

        integer(kind=im), parameter :: bandNumber = 16
     integer(kind=im), parameter :: numGPoints = no16
        integer(kind=im), parameter :: gPointSetNumber = 1
    integer(kind=im) :: ncid, varID,sts(13)

        sts(:)  = nf90_NoErr
    sts(1)  = nf90_open('rrtmg_lw.nc',nf90_nowrite,ncid)

    sts(2)  = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
    sts(3)  = nf90_get_var(ncid, varID, fracrefao, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,keylower,1,1/))

    sts(4)  = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
    sts(5)  = nf90_get_var(ncid, varID, fracrefbo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,1,1,1/))

    sts(6)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
    sts(7)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/keylower,Tdiff,plower,numGPoints,1,1/))

    sts(8)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
    sts(9)  = nf90_get_var(ncid, varID, kbo, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,Tdiff,pupper,numGPoints,1,1/))

    sts(10) = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
    sts(11) = nf90_get_var(ncid, varID, selfrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tself,numGPoints,1,1/))

    sts(12) = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
    sts(13) = nf90_get_var(ncid, varID, forrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tforeign,numGPoints,1,1/))

    if(any(sts(:) /= nf90_NoErr)) stop  "Error reading variables from file"

    sts(1) = nf90_close(ncid)

end subroutine lw_kgb16
!*******************************************************************************
end module rrtmg_lw_read_nc
