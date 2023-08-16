module modgpu
  use modprecision, only: pois_r
  implicit none

save
  real(pois_r), allocatable, target :: workspace(:)

#if defined(_OPENACC)
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

    use modchecksim, only: courxl, couryl, courzl, courtotl, peclettotl

    implicit none

    ! Prognostic variables
    !$acc enter data copyin(um, vm, wm, thlm, e12m, qtm) async
    !$acc enter data copyin(u0, v0, w0, thl0, thl0h, qt0h, e120, qt0) async
    !$acc enter data copyin(up, vp, wp, thlp, e12p, qtp) async
    !$acc enter data copyin(svm, sv0, svp) async

    ! Base state variables
    !$acc enter data copyin(rhobf, rhobh) async

    ! Diagnostic variables
    !$acc enter data copyin(ql0, ql0h, tmp0, thv0h, dthvdz) async
    !$acc enter data copyin(whls, ug, vg) async
    !$acc enter data copyin(thvf, thvh) async
    !$acc enter data copyin(qt0av, ql0av, thl0av, u0av, v0av, sv0av) async
    !$acc enter data copyin(dpdxl, dpdyl) async
    !$acc enter data copyin(dthldxls, dthldyls, dthldtls) async
    !$acc enter data copyin(dqtdxls, dqtdyls, dqtdtls) async
    !$acc enter data copyin(dudxls, dudyls, dudtls) async
    !$acc enter data copyin(dvdxls, dvdyls, dvdtls) async

    ! Global
    !$acc enter data copyin(dzf, dzh, zh, zf, delta, deltai) async

    ! Surface
    !$acc enter data copyin(ustar, dudz, dvdz) async
    !$acc enter data copyin(thlflux, qtflux, dqtdz, dthldz) async
    !$acc enter data copyin(svflux, svs) async

    ! Boundary
    !$acc enter data copyin(tsc) async

    ! Subgrid
    !$acc enter data copyin(ekm, ekh, zlt, csz, anis_fac) async
    !$acc enter data copyin(sbdiss, sbshr, sbbuo) async
    
    ! Radiation
    !$acc enter data copyin(thlpcar, thlprad) async

    ! Thermodynamics
    !$acc enter data copyin(th0av, thv0) async
    !$acc enter data copyin(presf, presh, exnf, exnh) async
    !$acc enter data copyin(rhof) async
    !$acc enter data copyin(qvsl, qvsi, esl, qsat) async
    !$acc enter data copyin(esatmtab, esatitab, esatltab) async

    !$acc wait

    ! Statistics
    !$acc enter data copyin(courxl, couryl, courzl, courtotl, &
    !$acc&                  peclettotl)

  end subroutine initgpu

<<<<<<< HEAD
end module modgpu
=======
  !> Allocate reusable workspace for transposes and FFT
  subroutine allocate_workspace(nx, ny, nz)
    implicit none
    
    integer, intent(in) :: nx, ny, nz

    allocate(workspace(nx*ny*nz))

    !$acc enter data create(workspace)
  
  end subroutine allocate_workspace
>>>>>>> OpenACC-data
#endif
end module modgpu
