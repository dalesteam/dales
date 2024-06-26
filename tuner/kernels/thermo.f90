#define imax 1000
#define jmax 500
#define kmax 130
module thermo
contains
  subroutine icethermo
    implicit none
    
    !$tuner initialize
    integer :: i, j, k, m
    real(8) :: qsat_, Tl, b, ql, qt, Tl_min
    real(8) :: ttab, esatltab, esatitab, ilratio, tdn, tup
    real(8) :: tlo, thi, es, T
    integer :: tlonr
    real(8), parameter :: rv = 461.5
    real(8), parameter :: rd = 287.04
    real(8), parameter :: cp = 1004.
    real(8), parameter :: rlv = 2.53e6
    real(8), allocatable, dimension(:, :, :) :: thl0, qt0, ql0
    real(8), allocatable, dimension(:) :: presf, exnf
    real(8), dimension(1:2000) :: esatmtab

    allocate(thl0(imax,jmax,kmax), qt0(imax,jmax,kmax), ql0(imax,jmax,kmax))
    allocate(presf(kmax), exnf(kmax))
  
    thl0(:,:,:) = 300
    qt0(:,:,:) = 0.0165
    ql0(:,:,:) = 0.
    presf(:) = 80000.
    exnf(:) = (presf(:) / 100000)**(rd/cp)
    Tl_min = 300 * exnf(1)

    do m=1,2000
      ttab=150.+0.2*m
      esatltab=exp(54.842763-6763.22/ttab-4.21*log(ttab)+0.000367*ttab+&
           tanh(0.0415*(ttab-218.8))*(53.878-1331.22/ttab-9.44523*log(ttab)+ 0.014025*ttab))

      esatitab=exp(9.550426-5723.265/ttab+3.53068*log(ttab)-0.00728332*ttab)
      ilratio = max(0.,min(1.,(ttab-tdn)/(tup-tdn)))
      esatmtab(m) = ilratio*esatltab + (1-ilratio)*esatitab
    end do
  
    !$acc enter data copyin(thl0(:imax,:jmax,:kmax))
    !$acc enter data copyin(qt0(:imax,:jmax,:kmax))
    !$acc enter data copyin(ql0(:imax,:jmax,:kmax))
    !$acc enter data copyin(presf(:kmax))
    !$acc enter data copyin(exnf(:kmax))
    !$acc enter data copyin(esatmtab(:2000))
    !$tuner stop
     
    !$tuner start icethermo
    !$acc parallel loop collapse(ncoll) private(Tl, qsat_, qt, ql, b, T, tlo, thi, tlonr, es) &
    !$acc& default(present) vector_length(nthreads)
    do k = 1, kmax
      do j = 2, jmax - 1
        do i = 2, imax - 1
          Tl = exnf(k)*thl0(i,j,k)
          qt = qt0(i,j,k)
  
          ! first step
          T = Tl
          tlonr=int((T-150.)*5.)
          tlo = 150. + 0.2*tlonr
          thi = tlo + 0.2
          es = (thi-T)*5.*esatmtab(tlonr)+(T-tlo)*5.*esatmtab(tlonr+1)

          qsat_ = (rd/rv) * es / (presf(k) - (1.-rd/rv)*es)
          b = rlv**2 / (rv * cp * Tl**2)
          qsat_ = qsat_ * (1 + b * qt) / (1 + b * qsat_)
  
          ql = max(qt0(i,j,k) - qsat_, 0.)
  
          ! update the starting point
          Tl = Tl + (rlv/cp) * ql
          qt = qt - ql
  
          ! second step
          T = Tl
          tlonr=int((T-150.)*5.)
          tlo = 150. + 0.2*tlonr
          thi = tlo + 0.2
          es = (thi-T)*5.*esatmtab(tlonr)+(T-tlo)*5.*esatmtab(tlonr+1)

          qsat_ = (rd/rv) * es / (presf(k) - (1.-rd/rv)*es)
          b = rlv**2 / (rv * cp * Tl**2)
          qsat_ = qsat_ * (1 + b * qt) / (1 + b * qsat_)
  
          ! save results
          ql = max(qt0(i,j,k) - qsat_, 0.)
          ql0(i,j,k) = ql
        end do
      end do
    end do
    !$tuner stop
  
    !$tuner deinitialize
    !$acc exit data delete(thl0(:imax,:jmax,:kmax))
    !$acc exit data delete(qt0(:imax,:jmax,:kmax))
    !$acc exit data delete(ql0(:imax,:jmax,:kmax))
    !$acc exit data delete(presf(:kmax))
    !$acc exit data delete(exnf(:kmax))
    !$tuner stop
  end subroutine icethermo
end module thermo
