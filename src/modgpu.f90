module modgpu
  use modprecision, only: pois_r
  implicit none

save
  real(pois_r), allocatable, target :: workspace(:)

#if defined(_OPENACC)
contains
  !> @brief Copies fields and arrays to GPU  
  subroutine update_gpu
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

    !$acc update device(um, vm, wm, thlm, e12m, qtm, &
    !$acc&              u0, v0, w0, thl0, thl0h, qt0h, e120, qt0, &
    !$acc&              up, vp, wp, thlp, e12p, qtp, &
    !$acc&              svm, sv0, svp, &
    !$acc&              rhobf, rhobh, &
    !$acc&              ql0, ql0h, tmp0, thv0h, dthvdz, &
    !$acc&              whls, ug, vg, &
    !$acc&              thvf, thvh, &
    !$acc&              qt0av, ql0av, thl0av, u0av, v0av, sv0av, &
    !$acc&              dpdxl, dpdyl, &
    !$acc&              dthldxls, dthldyls, dthldtls, &
    !$acc&              dqtdxls, dqtdyls, dqtdtls, &
    !$acc&              dudxls, dudyls, dudtls, &
    !$acc&              dvdxls, dvdyls, dvdtls, &
    !$acc&              dzf, dzh, zh, zf, delta, deltai, &
    !$acc&              ustar, dudz, dvdz, &
    !$acc&              thlflux, qtflux, dqtdz, dthldz, &
    !$acc&              tsc, &
    !$acc&              ekm, ekh, zlt, sbdiss, sbshr, sbbuo, csz, anis_fac, &
    !$acc&              svflux, &
    !$acc&              thlpcar, thlprad, &
    !$acc&              presf, presh, exnf, exnh, &
    !$acc&              rhof, &
    !$acc&              qvsl, qvsi, esl, qsat, &
    !$acc&              esatmtab, esatitab, esatltab)

  end subroutine update_gpu
  
  !> @brief Copies data from GPU to host, mostly for debugging
  subroutine update_host
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

    !$acc update self(um, vm, wm, thlm, e12m, qtm, &
    !$acc&            u0, v0, w0, thl0, thl0h, qt0h, e120, qt0, &
    !$acc&            up, vp, wp, thlp, e12p, qtp, &
    !$acc&            svm, sv0, svp, &
    !$acc&            rhobf, rhobh, &
    !$acc&            ql0, ql0h, tmp0, thv0h, dthvdz, &
    !$acc&            whls, ug, vg, &
    !$acc&            thvf, thvh, &
    !$acc&            qt0av, ql0av, thl0av, u0av, v0av, sv0av, &
    !$acc&            dpdxl, dpdyl, &
    !$acc&            dthldxls, dthldyls, dthldtls, &
    !$acc&            dqtdxls, dqtdyls, dqtdtls, &
    !$acc&            dudxls, dudyls, dudtls, &
    !$acc&            dvdxls, dvdyls, dvdtls, &
    !$acc&            dzf, dzh, zh, zf, delta, deltai, &
    !$acc&            ustar, dudz, dvdz, &
    !$acc&            thlflux, qtflux, dqtdz, dthldz, &
    !$acc&            tsc, &
    !$acc&            ekm, ekh, zlt, sbdiss, sbshr, sbbuo, csz, anis_fac, &
    !$acc&            svflux, &
    !$acc&            thlpcar, thlprad, &
    !$acc&            presf, presh, exnf, exnh, &
    !$acc&            rhof, &
    !$acc&            qvsl, qvsi, esl, qsat, &
    !$acc&            esatmtab, esatitab, esatltab)

  end subroutine update_host

  !> @brief Allocate reusable workspace for transposes and FFT
  subroutine allocate_workspace(nx, ny, nz)
    implicit none
    
    integer, intent(in) :: nx, ny, nz

    allocate(workspace(nx*ny*nz))

    !$acc enter data create(workspace)
  
  end subroutine allocate_workspace
#endif
end module modgpu
