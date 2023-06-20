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
                         rhobf, rhobh
    use modglobal, only: dzf, dzh, zh, zf, delta, deltai
    implicit none

    ! FIELDS
    !$acc update device(um, vm, wm, thlm, e12m, qtm)
    !$acc update device(u0, v0, w0, thl0, thl0h, qt0h, e120, qt0)
    !$acc update device(up, vp, wp, thlp, e12p, qtp)
    !$acc update device(svm, sv0, svp)
    !$acc update device(rhobf, rhobh)

    ! GLOBAL
    !$acc update device(dzf, dzh, zh, zf, delta, deltai)

    write(*,*) "GPU fields updated"
    
  end subroutine initgpu
end module modgpu
#endif
