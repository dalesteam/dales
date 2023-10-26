module modgpu
  use modprecision, only: pois_r
  implicit none

save
  real(pois_r), allocatable, target :: workspace(:)
  logical :: host_is_updated = .false.

#if defined(_OPENACC)
contains
  !> @brief Copies fields and arrays to GPU  
  subroutine update_gpu
    use modfields, only: um, u0, up, vm, v0, vp, wm, w0, wp, &
                         thlm, thl0, thlp, qtm, qt0, qtp, &
                         e12m, e120, e12p, svm, sv0, svp, &
                         rhobf, rhobh, ql0, tmp0, ql0h, thv0h, &
                         thl0h, qt0h, presf, presh, exnf, exnh, &
                         thvh, thvf, rhof, qt0av, ql0av, thl0av, &
                         u0av, v0av, sv0av, ug, vg, dpdxl, dpdyl, &
                         wfls, whls, thlpcar, dthldxls, dthldyls, &
                         dthldtls, dqtdxls, dqtdyls, dqtdtls, &
                         dudxls, dudyls, dudtls, dvdxls, dvdyls, &
                         dvdtls, dthvdz, qvsl, qvsi, esl, qsat
    use modglobal, only: dzf, dzh, zh, zf, delta, deltai, &
                         rd, rv, esatmtab, esatitab, esatltab
    use modsurfdata, only: z0m, z0h, obl, tskin, qskin, Cm, Cs, &
                           ustar, dudz, dvdz, thlflux, qtflux, &
                           dqtdz, dthldz, svflux, svs, horv
    use modsubgriddata, only: ekm, ekh, zlt, csz, anis_fac, &
                              sbdiss, sbshr, sbbuo
    use modradiation, only: thlprad
    use modthermodynamics, only: th0av, thv0, thetah, qth, qlh
    use modboundary, only: tsc
    use modchecksim, only: courxl, couryl, courzl, courtotl, peclettotl

    implicit none

    !$acc update device(um, u0, up, vm, v0, vp, wm, w0, wp, &
    !$acc&              thlm, thl0, thlp, qtm, qt0, qtp, &
    !$acc&              e12m, e120, e12p, svm, sv0, svp, &
    !$acc&              rhobf, rhobh, ql0, tmp0, ql0h, thv0h, &
    !$acc&              thl0h, qt0h, presf, presh, exnf, exnh, &
    !$acc&              thvh, thvf, rhof, qt0av, ql0av, thl0av, &
    !$acc&              u0av, v0av, sv0av, ug, vg, dpdxl, dpdyl, &
    !$acc&              wfls, whls, thlpcar, dthldxls, dthldyls, &
    !$acc&              dthldtls, dqtdxls, dqtdyls, dqtdtls, &
    !$acc&              dudxls, dudyls, dudtls, dvdxls, dvdyls, &
    !$acc&              dvdtls, dthvdz, qvsl, qvsi, esl, qsat, &
    !$acc&              dzf, dzh, zh, zf, delta, deltai, &
    !$acc&              z0m, z0h, obl, tskin, qskin, Cm, Cs, &
    !$acc&              ustar, dudz, dvdz, thlflux, qtflux, &
    !$acc&              dqtdz, dthldz, svflux, svs, horv, &
    !$acc&              ekm, ekh, zlt, sbdiss, sbshr, sbbuo, csz, &
    !$acc&              anis_fac, tsc, thlpcar, thlprad, presf, & 
    !$acc&              presh, exnf, exnh, rhof, thetah, &
    !$acc&              qvsl, qvsi, esl, qsat, qth, qlh, &
    !$acc&              esatmtab, esatitab, esatltab)

  end subroutine update_gpu
  
  !> @brief Copies data from GPU to host, mostly for debugging
  subroutine update_host
    use modfields, only: um, u0, up, vm, v0, vp, wm, w0, wp, &
                         thlm, thl0, thlp, qtm, qt0, qtp, &
                         e12m, e120, e12p, svm, sv0, svp, &
                         rhobf, rhobh, ql0, tmp0, ql0h, thv0h, &
                         thl0h, qt0h, presf, presh, exnf, exnh, &
                         thvh, thvf, rhof, qt0av, ql0av, thl0av, &
                         u0av, v0av, sv0av, ug, vg, dpdxl, dpdyl, &
                         wfls, whls, thlpcar, dthldxls, dthldyls, &
                         dthldtls, dqtdxls, dqtdyls, dqtdtls, &
                         dudxls, dudyls, dudtls, dvdxls, dvdyls, &
                         dvdtls, dthvdz, qvsl, qvsi, esl, qsat
    use modglobal, only: dzf, dzh, zh, zf, delta, deltai, &
                         rd, rv, esatmtab, esatitab, esatltab
    use modsurfdata, only: z0m, z0h, obl, tskin, qskin, Cm, Cs, &
                           ustar, dudz, dvdz, thlflux, qtflux, &
                           dqtdz, dthldz, svflux, svs, horv
    use modsubgriddata, only: ekm, ekh, zlt, csz, anis_fac, &
                              sbdiss, sbshr, sbbuo
    use modradiation, only: thlprad
    use modthermodynamics, only: th0av, thv0, thetah, qth, qlh
    use modboundary, only: tsc
    use modchecksim, only: courxl, couryl, courzl, courtotl, peclettotl

    implicit none

    if (host_is_updated) return

    !$acc update self(um, u0, up, vm, v0, vp, wm, w0, wp, &
    !$acc&            thlm, thl0, thlp, qtm, qt0, qtp, &
    !$acc&            e12m, e120, e12p, svm, sv0, svp, &
    !$acc&            rhobf, rhobh, ql0, tmp0, ql0h, thv0h, &
    !$acc&            thl0h, qt0h, presf, presh, exnf, exnh, &
    !$acc&            thvh, thvf, rhof, qt0av, ql0av, thl0av, &
    !$acc&            u0av, v0av, sv0av, ug, vg, dpdxl, dpdyl, &
    !$acc&            wfls, whls, thlpcar, dthldxls, dthldyls, &
    !$acc&            dthldtls, dqtdxls, dqtdyls, dqtdtls, &
    !$acc&            dudxls, dudyls, dudtls, dvdxls, dvdyls, &
    !$acc&            dvdtls, dthvdz, qvsl, qvsi, esl, qsat, &
    !$acc&            dzf, dzh, zh, zf, delta, deltai, &
    !$acc&            z0m, z0h, obl, tskin, qskin, Cm, Cs, &
    !$acc&            ustar, dudz, dvdz, thlflux, qtflux, &
    !$acc&            dqtdz, dthldz, svflux, svs, horv, &
    !$acc&            ekm, ekh, zlt, sbdiss, sbshr, sbbuo, csz, &
    !$acc&            anis_fac, tsc, thlpcar, thlprad, presf, & 
    !$acc&            presh, exnf, exnh, rhof, thetah, &
    !$acc&            qvsl, qvsi, esl, qsat, qth, qlh, &
    !$acc&            esatmtab, esatitab, esatltab)

    host_is_updated = .true.

  end subroutine update_host

  !> @brief Allocate reusable workspace for transposes and FFT
  subroutine allocate_workspace(n)
    implicit none
    
    integer, intent(in) :: n

    allocate(workspace(n))

    workspace = 0

    !$acc enter data copyin(workspace)
  
  end subroutine allocate_workspace

  !> @brief Deallocate GPU workspaces
  subroutine deallocate_workspace
    implicit none

    !$acc exit data delete(workspace)
    deallocate(workspace)

  end subroutine deallocate_workspace
#endif
end module modgpu
