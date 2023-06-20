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
                         dpdxl, dpdyl, &
                         dthldxls, dthldyls, dthldtls, &
                         dqtdxls, dqtdyls, dqtdtls, &
                         dudxls, dudyls, dudtls, &
                         dvdxls, dvdyls, dvdtls
    use modglobal, only: dzf, dzh, zh, zf, delta, deltai
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
    !$acc update device(dpdxl, dpdyl)
    !$acc update device(dthldxls, dthldyls, dthldtls, &
    !$acc update device(dqtdxls, dqtdyls, dqtdtls)
    !$acc update device(dudxls, dudyls, dudtls)
    !$acc update device(dvdxls, dvdyls, dvdtls)

    ! Global
    !$acc update device(dzf, dzh, zh, zf, delta, deltai)


    write(*,*) "GPU fields updated"
    
  end subroutine initgpu
end module modgpu
#endif
