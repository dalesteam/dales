#if defined(_OPENACC)
module modgpu
  use modprecision, only: pois_r
  implicit none

save
  real(pois_r), allocatable, target :: workspace(:)

contains
  !< Initialize fields on the GPU post-startup  
  subroutine initgpu
    use modfields, only: um, vm, wm, thlm, e12m, qtm, &
                         u0, v0, w0, thl0, thl0h, qt0h, e120, qt0, &
                         up, vp, wp, thlp, e12p, qtp, &
                         svm, sv0, svp, &
                         rhobf, rhobh, &
                         ql0, ql0h, tmp0, thv0h, dthvdz, &
                         whls, ug, vg, thvf, thvh, &
                         presf, presh,exnf, exnh, &
                         rhof, &
                         qt0av, ql0av, thl0av, u0av, v0av, &
                         dpdxl, dpdyl, &
                         dthldxls, dthldyls, dthldtls, &
                         dqtdxls, dqtdyls, dqtdtls, &
                         dudxls, dudyls, dudtls, &
                         dvdxls, dvdyls, dvdtls, &
                         thlpcar, sv0av, &
                         qvsl, qvsi, esl, qsat

    use modglobal, only: dzf, dzh, zh, zf, delta, deltai, &
                         rd, rv, esatmtab, esatitab, esatltab

    use modsurfdata, only: ustar, dudz, dvdz, &
                           thlflux, qtflux, dqtdz, dthldz, &
                           svflux, svs

    use modsubgriddata, only: ekm, ekh, zlt, csz, anis_fac, &
                              sbdiss, sbshr, sbbuo

    use modradiation, only: thlprad

    use modthermodynamics, only: th0av, thv0

    use modboundary, only: tsc

    implicit none

    ! Prognostic variables
    !$acc enter data copyin(um, vm, wm, thlm, e12m, qtm)
    !$acc enter data copyin(u0, v0, w0, thl0, thl0h, qt0h, e120, qt0)
    !$acc enter data copyin(up, vp, wp, thlp, e12p, qtp)
    !$acc enter data copyin(svm, sv0, svp)

    ! Base state variables
    !$acc enter data copyin(rhobf, rhobh)

    ! Diagnostic variables
    !$acc enter data copyin(ql0, ql0h, tmp0, thv0h, dthvdz)
    !$acc enter data copyin(whls, ug, vg)
    !$acc enter data copyin(thvf, thvh)
    !$acc enter data copyin(qt0av, ql0av, thl0av, u0av, v0av, sv0av)
    !$acc enter data copyin(dpdxl, dpdyl)
    !$acc enter data copyin(dthldxls, dthldyls, dthldtls)
    !$acc enter data copyin(dqtdxls, dqtdyls, dqtdtls)
    !$acc enter data copyin(dudxls, dudyls, dudtls)
    !$acc enter data copyin(dvdxls, dvdyls, dvdtls)

    ! Global
    !$acc enter data copyin(dzf, dzh, zh, zf, delta, deltai)

    ! Surface
    !$acc enter data copyin(ustar, dudz, dvdz)
    !$acc enter data copyin(thlflux, qtflux, dqtdz, dthldz)
    !$acc enter data copyin(svflux, svs)

    ! Boundary
    !$acc enter data copyin(tsc)

    ! Subgrid
    !$acc enter data copyin(ekm, ekh, zlt, csz, anis_fac)
    !$acc enter data copyin(sbdiss, sbshr, sbbuo)
    
    ! Radiation
    !$acc enter data copyin(thlpcar, thlprad)

    ! Thermodynamics
    !$acc enter data copyin(th0av, thv0)
    !$acc enter data copyin(presf, presh, exnf, exnh)
    !$acc enter data copyin(rhof) 
    !$acc enter data copyin(qvsl, qvsi, esl, qsat)
    !$acc enter data copyin(esatmtab, esatitab, esatltab)

  end subroutine initgpu

  !> Allocate reusable workspace for transposes and FFT
  subroutine allocate_workspace(nx, ny, nz)
    implicit none
    
    integer, intent(in) :: nx, ny, nz

    allocate(workspace(nx*ny*nz))

    !$acc enter data create(workspace)
  
  end subroutine allocate_workspace
end module modgpu
#endif
