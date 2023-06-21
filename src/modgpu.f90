#if defined(_OPENACC)
module modgpu

public :: initgpu

contains
  !< Initialize fields on the GPU post-startup  
  subroutine initgpu
    use modfields, only: um, vm, wm, thlm, e12m, qtm, &
                         u0, v0, w0, thl0, thl0h, qt0h, e120, qt0, &
                         up, vp, wp, thlp, e12p, qtp, &
                         svm, sv0, svp, &
                         rhobf, rhobh, &
                         ql0, ql0h, tmp0, thv0h, dthvdz, &
                         whls, ug, vg, &
                         qt0av, ql0av, thl0av, u0av, v0av, &
                         dpdxl, dpdyl, &
                         dthldxls, dthldyls, dthldtls, &
                         dqtdxls, dqtdyls, dqtdtls, &
                         dudxls, dudyls, dudtls, &
                         dvdxls, dvdyls, dvdtls

    use modglobal, only: dzf, dzh, zh, zf, delta, deltai

    use modsurfdata, only: ustar, dudz, dvdz, &
                           thlflux, qtflux, dqtdz, dthldz, &
                           svflux, svs

    use modsubgriddata, only: ekm, ekh, zlt, csz, anis_fac, &
                              sbdiss, sbshr, sbbuo
    implicit none

    ! Prognostic variables
    !$acc update device(um, vm, wm, thlm, e12m, qtm)
    !$acc update device(u0, v0, w0, thl0, thl0h, qt0h, e120, qt0)
    !$acc update device(up, vp, wp, thlp, e12p, qtp)
    !$acc update device(svm, sv0, svp)

    ! Base state variables
    !$acc update device(rhobf, rhobh)

    ! Diagnostic variables
    !$acc update device(ql0, ql0h, tmp0, thv0h, dthvdz)
    !$acc update device(whls, ug, vg)
    !$acc update device(qt0av, ql0av, thl0av, u0av, v0av)
    !$acc update device(dpdxl, dpdyl)
    !$acc update device(dthldxls, dthldyls, dthldtls)
    !$acc update device(dqtdxls, dqtdyls, dqtdtls)
    !$acc update device(dudxls, dudyls, dudtls)
    !$acc update device(dvdxls, dvdyls, dvdtls)

    ! Global
    !$acc update device(dzf, dzh, zh, zf, delta, deltai)

    ! Surface
    !$acc update device(ustar, dudz, dvdz)
    !$acc update device(thlflux, qtflux, dqtdz, dthldz)
    !$acc update device(svflux, svs)

    ! Subgrid
    !$acc update device(ekm, ekh, zlt, csz, anis_fac)
    !$acc update device(sbdiss, sbshr, sbbuo)

    
  end subroutine initgpu
end module modgpu
#endif
