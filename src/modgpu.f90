module modgpu
  use modprecision, only: pois_r
  implicit none

save
  real(pois_r), allocatable, target :: workspace_0(:), workspace_1(:)
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
                         dvdtls, dthvdz, qvsl, qvsi, esl, qsat, w0av
    use modglobal, only: dzf, dzh, zh, zf, delta, deltai, &
                         rd, rv
    use modsurfdata, only: z0m, z0h, obl, tskin, qskin, Cm, Cs, &
                           ustar, dudz, dvdz, thlflux, qtflux, &
                           dqtdz, dthldz, svflux, svs, horv, ra, rs, wsvsurf
    use modsubgriddata, only: ekm, ekh, zlt, csz, anis_fac, &
                              sbdiss, sbshr, sbbuo
    use modraddata, only: thlprad, lwd, lwu, swd, swu, lwc, swdir, swdif, &
                          lwdca, lwuca, swdca, swuca, &
                          LW_dn_TOA, LW_up_TOA, SW_dn_TOA, SW_up_TOA, &
                          LW_dn_ca_TOA, LW_up_ca_TOA, SW_dn_ca_TOA, SW_up_ca_TOA
    use modthermodynamics, only: th0av, thv0, thetah, qth, qlh
    use modboundary, only: tsc
    use modmicrodata, only: precep, thlpmcr, qtpmcr
    use modibm,      only: fluid_mask, iobst, ixw_p, ixw_m, iyw_p, iyw_m, izw_p

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
    !$acc&              dqtdz, dthldz, svflux, svs, horv, ra, rs, wsvsurf, &
    !$acc&              ekm, ekh, zlt, sbdiss, sbshr, sbbuo, csz, &
    !$acc&              anis_fac, tsc, thlpcar, presf, &
    !$acc&              presh, exnf, exnh, thetah, &
    !$acc&              qvsl, qvsi, esl, qsat, qth, qlh, &
    !$acc&              th0av, thv0, thetah, qth, qlh, &
    !$acc&              precep, thlpmcr, qtpmcr, &
    !$acc&              thlprad, lwd, lwu, swd, swu, lwc, swdir, swdif, &
    !$acc&              lwdca, lwuca, swdca, swuca, &
    !$acc&              LW_dn_TOA, LW_up_TOA, SW_dn_TOA, SW_up_TOA, &
    !$acc&              LW_dn_ca_TOA, LW_up_ca_TOA, SW_dn_ca_TOA, SW_up_ca_TOA, &
    !$acc&              fluid_mask, iobst, ixw_p, ixw_m, iyw_p, iyw_m, izw_p, &
    !$acc&              w0av)

  end subroutine update_gpu

  !> @brief Copies surface data to GPU containers
  subroutine update_gpu_surface

    use modsurfdata, only: obl, tskin, qskin, ra, rs

    implicit none

    !$acc update device(tskin, qskin, ra, rs, obl)

  end subroutine update_gpu_surface
  
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
                         dvdtls, dthvdz, qvsl, qvsi, esl, qsat, w0av
    use modglobal, only: dzf, dzh, zh, zf, delta, deltai, &
                         rd, rv
    use modsurfdata, only: z0m, z0h, obl, tskin, qskin, Cm, Cs, &
                           ustar, dudz, dvdz, thlflux, qtflux, &
                           dqtdz, dthldz, svflux, svs, horv, ra, rs, wsvsurf
    use modsubgriddata, only: ekm, ekh, zlt, csz, anis_fac, &
                              sbdiss, sbshr, sbbuo
    use modraddata, only: thlprad, lwd, lwu, swd, swu, lwc, swdir, swdif, &
                          lwdca, lwuca, swdca, swuca, &
                          LW_dn_TOA, LW_up_TOA, SW_dn_TOA, SW_up_TOA, &
                          LW_dn_ca_TOA, LW_up_ca_TOA, SW_dn_ca_TOA, SW_up_ca_TOA
    use modthermodynamics, only: th0av, thv0, thetah, qth, qlh
    use modboundary, only: tsc
    use modmicrodata, only: precep, thlpmcr, qtpmcr
    use modibm,      only: fluid_mask, iobst, ixw_p, ixw_m, iyw_p, iyw_m, izw_p

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
    !$acc&            dqtdz, dthldz, svflux, svs, horv, ra, rs, wsvsurf, &
    !$acc&            ekm, ekh, zlt, sbdiss, sbshr, sbbuo, csz, &
    !$acc&            anis_fac, tsc, thlpcar, presf, &
    !$acc&            presh, exnf, exnh, thetah, &
    !$acc&            qvsl, qvsi, esl, qsat, qth, qlh, &
    !$acc&            th0av, thv0, thetah, qth, qlh, &
    !$acc&            precep, thlpmcr, qtpmcr, &
    !$acc&            thlprad, lwd, lwu, swd, swu, lwc, swdir, swdif, &
    !$acc&            lwdca, lwuca, swdca, swuca, &
    !$acc&            LW_dn_TOA, LW_up_TOA, SW_dn_TOA, SW_up_TOA, &
    !$acc&            LW_dn_ca_TOA, LW_up_ca_TOA, SW_dn_ca_TOA, SW_up_ca_TOA, &
    !$acc&            fluid_mask, iobst, ixw_p, ixw_m, iyw_p, iyw_m, izw_p, &
    !$acc&            w0av)

    host_is_updated = .true.

  end subroutine update_host

  !> @brief Copies surface data to host containers
  subroutine update_host_surface

    use modsurfdata, only: obl, tskin, qskin

    implicit none

    !$acc update self(tskin, qskin, obl)

  end subroutine update_host_surface

  !> @brief Allocate reusable workspace for transposes and FFT
  subroutine allocate_workspace(n)
    use modmpi, only: nprocs
    implicit none
    
    integer, intent(in) :: n

    allocate(workspace_0(n))
    workspace_0 = 0
    !$acc enter data copyin(workspace_0)

    ! Allocate another workspace for the all-to-all operations
    if (nprocs > 1) then
      allocate(workspace_1(n))
      workspace_1 = 0
      !$acc enter data copyin(workspace_1)
    end if

  end subroutine allocate_workspace

  !> @brief Deallocate GPU workspaces
  subroutine deallocate_workspace
    use modmpi, only: nprocs
    implicit none

    !$acc exit data delete(workspace_0)
    deallocate(workspace_0)

    if (nprocs > 1) then
      !$acc exit data delete (workspace_1)
      deallocate(workspace_1)
    end if

  end subroutine deallocate_workspace
#else
contains
  ! Dummies because CI/CD complains
  subroutine update_host
    implicit none
  end subroutine update_host

  subroutine update_gpu
    implicit none
  end subroutine update_gpu
#endif
end module modgpu
