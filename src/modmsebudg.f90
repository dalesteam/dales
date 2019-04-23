!> \file modmsebudg.f90
!! Compute moist static energy budget

module modmsebudg

implicit none
PRIVATE :: icethermom, u_advecc_52, w_advecc_52
PUBLIC :: msebudg
save

real, allocatable :: msem (:,:,:)


contains

  subroutine u_advecc_52(putin, putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi,dyi,dzf,dzh,leq
  use modfields, only : u0, v0, w0,rhobf
  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the cell centered field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency
  real                                       :: inv2dzfk, rhobf_p, rhobf_m

  integer :: i,j,k

  if (leq) then
   
  k = 1
  inv2dzfk = 1./(2. * dzf(k))
  rhobf_p = rhobf(k+1)/rhobf(k)
  
  do j=2,j1
     do i=2,i1

        putout(i,j,k)  = putout(i,j,k)- ( &
             ( &
             u0(i+1,j,k)/60.&
             *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
             -abs(u0(i+1,j,k))/60.&
             *(10.*(putin(i+1,j,k)-putin(i,j,k))-5.*(putin(i+2,j,k)-putin(i-1,j,k))+(putin(i+3,j,k)-putin(i-2,j,k)))&
             -u0(i,j,k)/60.&
             *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
             +abs(u0(i,j,k))/60.&
             *(10.*(putin(i,j,k)-putin(i-1,j,k))-5.*(putin(i+1,j,k)-putin(i-2,j,k))+(putin(i+2,j,k)-putin(i-3,j,k)))&
             )*dxi&
             +(&
             v0(i,j+1,k)/60.&
             *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
             -abs(v0(i,j+1,k))/60.&
             *(10.*(putin(i,j+1,k)-putin(i,j,k))-5.*(putin(i,j+2,k)-putin(i,j-1,k))+(putin(i,j+3,k)-putin(i,j-2,k)))&
             -v0(i,j,k)/60.&
             *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
             +abs(v0(i,j,k))/60.&
             *(10.*(putin(i,j,k)-putin(i,j-1,k))-5.*(putin(i,j+1,k)-putin(i,j-2,k))+(putin(i,j+2,k)-putin(i,j-3,k))) &
             )* dyi &
         )    
     end do
  end do

  do k=2,kmax
     inv2dzfk = 1./(2. * dzf(k))
     rhobf_p = rhobf(k+1)/rhobf(k)
     rhobf_m = rhobf(k-1)/rhobf(k)
     do j=2,j1
        do i=2,i1

            putout(i,j,k)  = putout(i,j,k)- (  &
                  ( &
                      u0(i+1,j,k)/60.&
                      *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                      -abs(u0(i+1,j,k))/60.&
                      *(10.*(putin(i+1,j,k)-putin(i,j,k))-5.*(putin(i+2,j,k)-putin(i-1,j,k))+(putin(i+3,j,k)-putin(i-2,j,k)))&
                      -u0(i,j,k)/60.&
                      *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                      +abs(u0(i,j,k))/60.&
                      *(10.*(putin(i,j,k)-putin(i-1,j,k))-5.*(putin(i+1,j,k)-putin(i-2,j,k))+(putin(i+2,j,k)-putin(i-3,j,k)))&
                  )*dxi&
                +(&
                      v0(i,j+1,k)/60.&
                      *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                      -abs(v0(i,j+1,k))/60.&
                      *(10.*(putin(i,j+1,k)-putin(i,j,k))-5.*(putin(i,j+2,k)-putin(i,j-1,k))+(putin(i,j+3,k)-putin(i,j-2,k)))&
                      -v0(i,j,k)/60.&
                      *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                      +abs(v0(i,j,k))/60.&
                      *(10.*(putin(i,j,k)-putin(i,j-1,k))-5.*(putin(i,j+1,k)-putin(i,j-2,k))+(putin(i,j+2,k)-putin(i,j-3,k)))&
                  )* dyi &
              )
      end do
    end do
 end do

 else ! non-equidistant grid
  k = 1
  inv2dzfk = 1./(2. * dzf(k))
  rhobf_p = rhobf(k+1)/rhobf(k)
  
  do j=2,j1
     do i=2,i1

        putout(i,j,k)  = putout(i,j,k)- ( &
             ( &
             u0(i+1,j,k)/60.&
             *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
             -abs(u0(i+1,j,k))/60.&
             *(10.*(putin(i+1,j,k)-putin(i,j,k))-5.*(putin(i+2,j,k)-putin(i-1,j,k))+(putin(i+3,j,k)-putin(i-2,j,k)))&
             -u0(i,j,k)/60.&
             *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
             +abs(u0(i,j,k))/60.&
             *(10.*(putin(i,j,k)-putin(i-1,j,k))-5.*(putin(i+1,j,k)-putin(i-2,j,k))+(putin(i+2,j,k)-putin(i-3,j,k)))&
             )*dxi&
             +(&
             v0(i,j+1,k)/60.&
             *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
             -abs(v0(i,j+1,k))/60.&
             *(10.*(putin(i,j+1,k)-putin(i,j,k))-5.*(putin(i,j+2,k)-putin(i,j-1,k))+(putin(i,j+3,k)-putin(i,j-2,k)))&
             -v0(i,j,k)/60.&
             *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
             +abs(v0(i,j,k))/60.&
             *(10.*(putin(i,j,k)-putin(i,j-1,k))-5.*(putin(i,j+1,k)-putin(i,j-2,k))+(putin(i,j+2,k)-putin(i,j-3,k))) &
             )* dyi &
             )       
     end do
  end do

  do k=2,kmax
     inv2dzfk = 1./(2. * dzf(k))
     rhobf_p = rhobf(k+1)/rhobf(k)
     rhobf_m = rhobf(k-1)/rhobf(k)
     do j=2,j1
        do i=2,i1

            putout(i,j,k)  = putout(i,j,k)- (  &
                  ( &
                      u0(i+1,j,k)/60.&
                      *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                      -abs(u0(i+1,j,k))/60.&
                      *(10.*(putin(i+1,j,k)-putin(i,j,k))-5.*(putin(i+2,j,k)-putin(i-1,j,k))+(putin(i+3,j,k)-putin(i-2,j,k)))&
                      -u0(i,j,k)/60.&
                      *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                      +abs(u0(i,j,k))/60.&
                      *(10.*(putin(i,j,k)-putin(i-1,j,k))-5.*(putin(i+1,j,k)-putin(i-2,j,k))+(putin(i+2,j,k)-putin(i-3,j,k)))&
                  )*dxi&
                +(&
                      v0(i,j+1,k)/60.&
                      *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                      -abs(v0(i,j+1,k))/60.&
                      *(10.*(putin(i,j+1,k)-putin(i,j,k))-5.*(putin(i,j+2,k)-putin(i,j-1,k))+(putin(i,j+3,k)-putin(i,j-2,k)))&
                      -v0(i,j,k)/60.&
                      *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                      +abs(v0(i,j,k))/60.&
                      *(10.*(putin(i,j,k)-putin(i,j-1,k))-5.*(putin(i,j+1,k)-putin(i,j-2,k))+(putin(i,j+2,k)-putin(i,j-3,k)))&
                  )* dyi &
               )
      end do
    end do
 end do

 end if

end subroutine u_advecc_52


 

subroutine w_advecc_52(putin, putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi,dyi,dzf,dzh,leq
  use modfields, only : u0, v0, w0,rhobf
  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the cell centered field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency
  real                                       :: inv2dzfk, rhobf_p, rhobf_m

  integer :: i,j,k

  if (leq) then
   
  k = 1
  inv2dzfk = 1./(2. * dzf(k))
  rhobf_p = rhobf(k+1)/rhobf(k)
  
  do j=2,j1
     do i=2,i1

        putout(i,j,k)  = putout(i,j,k)- ( &
                         ( &
             w0(i,j,k+1) * (rhobf_p * putin(i,j,k+1) + putin(i,j,k)) &
             ) * inv2dzfk  &
             )
     end do
  end do

  do k=2,kmax
     inv2dzfk = 1./(2. * dzf(k))
     rhobf_p = rhobf(k+1)/rhobf(k)
     rhobf_m = rhobf(k-1)/rhobf(k)
     do j=2,j1
        do i=2,i1

            putout(i,j,k)  = putout(i,j,k)- (  &
                  ( &
                  w0(i,j,k+1) * (rhobf_p * putin(i,j,k+1) +  putin(i,j,k)) &
                  -w0(i,j,k)  * (rhobf_m * putin(i,j,k-1) +  putin(i,j,k)) &
                  ) * inv2dzfk &
                  )


      end do
    end do
 end do

 else ! non-equidistant grid
  k = 1
  inv2dzfk = 1./(2. * dzf(k))
  rhobf_p = rhobf(k+1)/rhobf(k)
  
  do j=2,j1
     do i=2,i1

        putout(i,j,k)  = putout(i,j,k)- ( &
              ( &
             w0(i,j,k+1) * (rhobf_p * putin(i,j,k+1) * dzf(k) +  putin(i,j,k) * dzf(k+1) ) / dzh(k+1) &
             ) * inv2dzfk  &
             )       
     end do
  end do

  do k=2,kmax
     inv2dzfk = 1./(2. * dzf(k))
     rhobf_p = rhobf(k+1)/rhobf(k)
     rhobf_m = rhobf(k-1)/rhobf(k)
     do j=2,j1
        do i=2,i1

            putout(i,j,k)  = putout(i,j,k)- (  &
                  ( &
                w0(i,j,k+1) * (rhobf_p * putin(i,j,k+1) * dzf(k) +  putin(i,j,k) * dzf(k+1) ) / dzh(k+1) &
                -w0(i,j,k ) * (rhobf_m * putin(i,j,k-1) * dzf(k) +  putin(i,j,k) * dzf(k-1) ) / dzh(k) &
                  ) * inv2dzfk &
                  )
      end do
    end do
 end do

 end if

end subroutine w_advecc_52




  subroutine icethermom
!> Calculates liquid water content.and temperature
!! \author Steef B\"oing

  use modglobal, only : i1,j1,k1,ih,jh,rd,rv,rlv,tup,tdn,cp,ttab,esatltab,esatitab,grav,rlv,zf
  use modfields, only : qtm,thlm,exnf,presf,esl,thl0,qt0
  implicit none

  integer i, j, k
  real :: ilratio, esl1,esi1,qvsl1,qvsi1,qsatur, thlguess, thlguessmin,tlo,thi,ttry
  real :: Tnr,Tnr_old
  integer :: niter,nitert,tlonr,thinr
  real :: tmpm,qlm
  real :: qvsl,qvsi

!     calculation of T with Newton-Raphson method
!     first guess is Tnr=tl
      nitert = 0
      niter = 0
      do k=1,k1
      do j=2-jh,j1+jh
      do i=2-ih,i1+ih
         if (thlm(i,j,k).eq.0) then 
           write(6,*) 'thlm =0',i,j,k,thlm(i,j,k),thl0(i,j,k)
         endif
         if (qtm(i,j,k).eq.0) then 
            write(6,*) 'qtm =0',i,j,k,qtm(i,j,k),qt0(i,j,k)
               
         endif
        
            ! first guess for temperature
            Tnr=exnf(k)*thlm(i,j,k)
            ilratio = max(0.,min(1.,(Tnr-tdn)/(tup-tdn)))
            tlonr=int((Tnr-150.)*5.)
            thinr=tlonr+1
            tlo=ttab(tlonr)
            thi=ttab(thinr)
            esl1=(thi-Tnr)*5.*esatltab(tlonr)+(Tnr-tlo)*5.*esatltab(thinr)
            esi1=(thi-Tnr)*5.*esatitab(tlonr)+(Tnr-tlo)*5.*esatitab(thinr)
            qvsl1=(rd/rv)*esl1/(presf(k)-(1.-rd/rv)*esl1)
            qvsi1=(rd/rv)*esi1/(presf(k)-(1.-rd/rv)*esi1)
            qsatur = ilratio*qvsl1+(1.-ilratio)*qvsi1
            if(qtm(i,j,k)>qsatur) then
              Tnr_old=0.
              niter = 0
              thlguess = Tnr/exnf(k)-(rlv/(cp*exnf(k)))*max(qtm(i,j,k)-qsatur,0.)
              ttry=Tnr-0.002
              ilratio = max(0.,min(1.,(ttry-tdn)/(tup-tdn)))
              tlonr=int((ttry-150.)*5.)
              thinr=tlonr+1
              tlo=ttab(tlonr)
              thi=ttab(thinr)
              esl1=(thi-ttry)*5.*esatltab(tlonr)+(ttry-tlo)*5.*esatltab(thinr)
              esi1=(thi-ttry)*5.*esatitab(tlonr)+(ttry-tlo)*5.*esatitab(thinr)
              qsatur = ilratio*(rd/rv)*esl1/(presf(k)-(1.-rd/rv)*esl1)+(1.-ilratio)*(rd/rv)*esi1/(presf(k)-(1.-rd/rv)*esi1)
              thlguessmin = ttry/exnf(k)-(rlv/(cp*exnf(k)))*max(qtm(i,j,k)-qsatur,0.)

              Tnr = Tnr - (thlguess-thlm(i,j,k))/((thlguess-thlguessmin)*500.)
              do while ((abs(Tnr-Tnr_old) > 0.002).and.(niter<100))
                niter = niter+1
                Tnr_old=Tnr
                ilratio = max(0.,min(1.,(Tnr-tdn)/(tup-tdn)))
                tlonr=int((Tnr-150.)*5.)
                if(tlonr<1 .or.tlonr>1999) then
                  write(*,*) 'thermo crash in modmse: i,j,k,niter,thlm(i,j,k),qtm(i,j,k)'
                  write(*,*) i,j,k,niter,thlm(i,j,k),qtm(i,j,k)
                  write(*,*) thl0(i,j,k), qt0(i,j,k),presf(k),Tnr,exnf(k)*thlm(i,j,k)
                endif
                thinr=tlonr+1
                tlo=ttab(tlonr)
                thi=ttab(thinr)
                esl1=(thi-Tnr)*5.*esatltab(tlonr)+(Tnr-tlo)*5.*esatltab(thinr)
                esi1=(thi-Tnr)*5.*esatitab(tlonr)+(Tnr-tlo)*5.*esatitab(thinr)
                qsatur = ilratio*(rd/rv)*esl1/(presf(k)-(1.-rd/rv)*esl1)+(1.-ilratio)*(rd/rv)*esi1/(presf(k)-(1.-rd/rv)*esi1)
                thlguess = Tnr/exnf(k)-(rlv/(cp*exnf(k)))*max(qtm(i,j,k)-qsatur,0.)

                ttry=Tnr-0.002
                ilratio = max(0.,min(1.,(ttry-tdn)/(tup-tdn)))
                tlonr=int((ttry-150.)*5.)
                thinr=tlonr+1
                tlo=ttab(tlonr)
                thi=ttab(thinr)
                esl1=(thi-ttry)*5.*esatltab(tlonr)+(ttry-tlo)*5.*esatltab(thinr)
                esi1=(thi-ttry)*5.*esatitab(tlonr)+(ttry-tlo)*5.*esatitab(thinr)
                qsatur = ilratio*(rd/rv)*esl1/(presf(k)-(1.-rd/rv)*esl1)+(1.-ilratio)*(rd/rv)*esi1/(presf(k)-(1.-rd/rv)*esi1)
                thlguessmin = ttry/exnf(k)-(rlv/(cp*exnf(k)))*max(qtm(i,j,k)-qsatur,0.)

                Tnr = Tnr - (thlguess-thlm(i,j,k))/((thlguess-thlguessmin)*500.)
              enddo
              nitert =max(nitert,niter)
              tmpm= Tnr
              ilratio = max(0.,min(1.,(Tnr-tdn)/(tup-tdn)))
              tlonr=int((Tnr-150.)*5.)
              thinr=tlonr+1
              tlo=ttab(tlonr)
              thi=ttab(thinr)
              esl(i,j,k)=(thi-Tnr)*5.*esatltab(tlonr)+(Tnr-tlo)*5.*esatltab(thinr)
              esi1=(thi-Tnr)*5.*esatitab(tlonr)+(Tnr-tlo)*5.*esatitab(thinr)
              qvsl=rd/rv*esl(i,j,k)/(presf(k)-(1.-rd/rv)*esl(i,j,k))
              qvsi=rd/rv*esi1/(presf(k)-(1.-rd/rv)*esi1)
              qsatur = ilratio*qvsl+(1.-ilratio)*qvsi
            else
              tmpm= Tnr
              esl(i,j,k)=esl1
              esi1=esi1
            endif
            qlm = max(qtm(i,j,k)-qsatur,0.)
            msem(i,j,k)  = cp * tmpm + grav * zf(k) + rlv * (qtm(i,j,k) - qlm)
            !rce write (6,*) i,j,k,tmpm,qtm(i,j,k),thlm(i,j,k),qlm,msem(i,j,k)
      end do
      end do
      end do
      if(nitert>99) then
      write(*,*) 'thermowarning'
      endif

  end subroutine icethermom

  subroutine msebudg
      use modglobal, only : ih,i1,jh,j1,k1,kmax,rdt,dzf,ijtot,&
                            iadv_thl,iadv_cd2,iadv_5th,iadv_52,iadv_cd6,iadv_62,iadv_kappa,&
                            iadv_upw,iadv_hybrid,iadv_hybrid_f,iadv_null,leq
      use modfields, only : mse0, field_mse_2D,rhof, thlm,qtm
      use modmpi,    only : my_real,mpi_sum,mpi_max,mpi_min,comm3d,mpierr,myid, excjs
      integer i,j,k
      real msemintavl,  mse0intavl,  msemintav,  mse0intav
      real msemint2avl, mse0int2avl, msemint2av, mse0int2av
      !rce see above real, allocatable :: msem(:,:,:)
      real, allocatable :: msep(:,:,:)
      allocate(msem(2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(msep(2-ih:i1+ih,2-jh:j1+jh,k1))
      
      msem=0.
      msep=0.

      field_mse_2D = 0
      msemintavl = 0
      mse0intavl = 0
      msemintav = 0
      mse0intav = 0

      msemint2avl = 0
      mse0int2avl = 0
      msemint2av = 0
      mse0int2av = 0

      call excjs( mse0           , 2,i1,2,j1,1,k1,ih,jh)
      call excjs( thlm           , 2,i1,2,j1,1,k1,ih,jh)   
      call excjs( qtm            , 2,i1,2,j1,1,k1,ih,jh)
   
      write (6,*) 'do icethermo'
      call icethermom

!rce      write (6,*) 'dt ',rdt
      do k=1,kmax
      do j=2-ih,j1+ih
      do i=2-jh,i1+jh
         !write(6,*) i,j,k,mse0(i,j,k),msem(i,j,k)
     enddo
     enddo
     enddo
      write (6,*) 'done icethermo'


      if (iadv_thl.eq.iadv_52) then 
        call u_advecc_52(mse0,msep)

        do k=1,k1
        do j=2,j1
        do i=2,i1
          !rce write (6,*) i,j,k,msep(i,j,k) 
          field_mse_2D (i,j,2) = field_mse_2D (i,j,2) + rhof(k) * dzf(k) * msep(i,j,k) 
        enddo
        enddo
        enddo
      
        call w_advecc_52(mse0,msep)

        do k=1,k1
        do j=2,j1
        do i=2,i1
          field_mse_2D (i,j,1) = field_mse_2D (i,j,1) + rhof(k) * dzf(k) * mse0 (i,j,k) 
          field_mse_2D (i,j,3) = field_mse_2D (i,j,3) + rhof(k) * dzf(k) * msep(i,j,k)   
          field_mse_2D (i,j,4) = field_mse_2D (i,j,4) + rhof(k) * dzf(k) * msem (i,j,k) 
        enddo
        enddo
        enddo

        else
          ! null advection scheme 
          stop "budget not prepared for this advection option"
      end if


      !compute slab mean values, needed for variance
      !
 
      do j=2,j1
      do i=2,i1
         mse0intavl = msemintavl  + field_mse_2D (i,j,1)
         msemintavl = msemintavl  + field_mse_2D (i,j,4)
      enddo
      enddo
     
      call MPI_ALLREDUCE(mse0intavl, mse0intav, 1,    MY_REAL, &
                          MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(msemintavl, msemintav, 1,    MY_REAL, &
                          MPI_SUM, comm3d,mpierr)

      mse0intav = mse0intav / ijtot
      msemintav = msemintav / ijtot
  
      mse0int2avl = mse0int2avl + (field_mse_2D (i,j,1) - mse0intav )**2
      msemint2avl = msemint2avl + (field_mse_2D (i,j,4) - msemintav )**2
      
      call MPI_ALLREDUCE(mse0int2avl, mse0int2av, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(msemint2avl, msemint2av, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)

      mse0int2av = mse0int2av / ijtot
      msemint2av = msemint2av / ijtot

      do j=2,j1
      do i=2,i1
         field_mse_2D (i,j,4) =  (field_mse_2D (i,j,1) - field_mse_2D (i,j,4)) / rdt
          write (6,*) 'check' ,i,j,field_mse_2D (i,j,4), field_mse_2D (i,j,1)
      enddo
      enddo

      deallocate( msem,msep)
  end subroutine msebudg




end module modmsebudg






