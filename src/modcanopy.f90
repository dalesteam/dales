!> \file modcanopy.f90
!!  Canopy parameterization
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2014 Delft University of Technology, Wageningen University, Utrecht University, KNMI, MPIC
!
module modcanopy
  implicit none
  save

  ! Namoptions
  logical :: lcanopy   = .false.       !< Switch to enable canopy representation
  logical :: lhetcanopy= .false.       !< Switch to enable heterogeneous canopy representation

  integer :: ncanopy   = 10            !< Amount of layers to represent the canopy
  real    :: cd        = 0.15          !< Drag coefficient in the canopy
  real    :: lai       = 2             !< Leaf Area Index (or actually plant area index) of the canopy
  logical :: lpaddistr = .false.       !< Switch to customize the general plant area density distribution (at half levels)
  integer :: npaddistr = 11            !< (if lpaddistr): number of half levels for prescribed general plant area density distribution
  logical :: lsrcdistr = .false.       !< Switch to customized src distribution (for scalars only at the moment, only including the canopy source)
  integer :: maxsrcdistr = 50          !< (if lsrcdistr): number of half levels for prescribed general plant area density distribution
  integer :: nsrcdistr(100)  !< length of srcfactor for scalar i

  logical :: wth_total = .false.       !< Switch: prescribed SH flux is added to surface flux if .false., it contains (the effect of) the surface flux if .true.
  logical :: wqt_total = .false.       !< Switch: prescribed LE flux is added to surface flux if .false., it contains (the effect of) the surface flux if .true.
  logical :: wsv_total(100) = .false.  !< Switch: prescribed sv flux is added to surface flux if .false., it contains (the effect of) the surface flux if .true.

  real    :: wth_can  = 0.0            !< prescribed SH canopy flux
  real    :: wqt_can  = 0.0            !< prescribed LE canopy flux
  real    :: wsv_can(100) = 0.0        !< prescribed scalar canopy flux

  real    :: wth_alph = 0.6            !< Decay constant for flux as function of the vertically integrated PAI (from canopy top)
  real    :: wqt_alph = 0.6            !< Decay constant for flux as function of the vertically integrated PAI (from canopy top)
  real    :: wsv_alph(100) = 0.6       !< Decay constant for flux as function of the vertically integrated PAI (from canopy top)

  ! Fields
  real, allocatable :: padfactor(:)    !< prescribed weighing factor for plant area density
  real, allocatable :: ppad(:)         !< resulting prescribed plant area density
  real, allocatable :: zpad(:)         !< heights of prescribed plant area density
  real, allocatable :: padtemp(:)      !< temporary plant area density used for calculations
  real, allocatable :: padf(:)         !< plant area density field full level
  real, allocatable :: padh(:)         !< plant area density field full level
  real, allocatable :: pai(:)          !< plant area index of the column starting in this grid cell up to the canopy top

  real, allocatable :: srcfactor(:,:)  !< prescribed weighing factor for scalar source density

  real              :: f_lai_h         !< average plant area density [m2/m2 / m]

  ! For heterogeneous runs: canopy properties per patch
  integer, parameter:: max_canopy = 10             !< Number of canopy types that can be defined
  integer           :: ncanopies                   !< Number of canopy definitions (<= max_canopy)
  integer, parameter:: max_npad   = 50             !< Number of levels for which PAD can be prescribed (later rescaled to actual number of levels of canopy)
  integer           :: npaddistr_c(max_canopy)=0   !< (if lpaddistr_c): number of half levels for prescribed general plant area density distribution
  integer, parameter:: max_ncan   = 100            !< Number of levels for which PAD can be prescribed in model levels
  integer           :: act_max_ncan                !< actual maximum level of the canopy (should be <= max_ncan)
  integer           :: ncanopy_c(max_canopy) = 0   !< Number of layers to represent the canopy
  integer           :: canopytype(max_canopy) = -1 !< Type nr of the land in surface.inp.xxx
  character(len=10),dimension(max_canopy) :: canopyname = "none" !< Name of the land type in canopy.inp.xxx
  real              :: cd_c(max_canopy)            !< Drag coefficient in the canopy
  real              :: lai_c(max_canopy)           !< Leaf Area Index (or actually plant area index) of the canopy
  logical           :: lpaddistr_c(max_canopy)     !< Switch to customize the general plant area density distribution (at half levels)

  logical :: wth_total_c(max_canopy) = .false.    !< Switch: prescribed SH flux is added to surface flux if .false., it contains (the effect of) the surface flux if .true.
  logical :: wqt_total_c(max_canopy) = .false.    !< Switch: prescribed LE flux is added to surface flux if .false., it contains (the effect of) the surface flux if .true.
  logical :: wsv_total_c(100,max_canopy) = .false.!< Switch: prescribed sv flux is added to surface flux if .false., it contains (the effect of) the surface flux if .true.

  real :: wth_can_c(max_canopy) = 0.0    !< prescribed SH canopy flux
  real :: wqt_can_c(max_canopy) = 0.0   !< prescribed LE canopy flux
  real :: wsv_can_c(100,max_canopy) = 0.0  !< prescribed scalar canopy flux

  real :: wth_alph_c(max_canopy) = 0.6     !< Decay constant for flux as function of the vertically integrated PAI (from canopy top)
  real :: wqt_alph_c(max_canopy) = 0.6     !< Decay constant for flux as function of the vertically integrated PAI (from canopy top)
  real :: wsv_alph_c(100,max_canopy) = 0.6 !< Decay constant for flux as function of the vertically integrated PAI (from canopy top)

  ! Patch definitions
  integer, parameter   :: mpatch_c = 16   !< Maximum number of patches that can be defined 
  integer              :: patchtype(mpatch_c)
  real                 :: minx_p(mpatch_c) !< minimum x position of edge of patch (fraction of total domain)
  real                 :: maxx_p(mpatch_c) !< maximum x position of edge of patch (fraction of total domain)
  real                 :: miny_p(mpatch_c) !< minimum y position of edge of patch (fraction of total domain)
  real                 :: maxy_p(mpatch_c) !< maximum y position of edge of patch (fraction of total domain)
  integer, allocatable :: map_c(:,:)      !< Canopy type for each grid point in the horizontal direction

  ! These are read from, and constructed from, the pad files
  real, allocatable :: padfactor_c(:,:)!< prescribed weighing factor for plant area density
  real, allocatable :: ppad_c(:,:)     !< resulting prescribed plant area density
  real, allocatable :: zpad_c(:,:)     !< heights of prescribed plant area density
  real, allocatable :: padtemp_c(:,:)  !< temporary plant area density used for calculations
  real, allocatable :: padf_c(:,:)     !< plant area density field full level
  real, allocatable :: padh_c(:,:)     !< plant area density field full level
  real, allocatable :: pai_c(:,:)      !< plant area index of the column starting in this grid cell up to the canopy top


contains

!-----------------------------------------------------------------------------------------
! Return relative position along x-axis given local i index
function relxpos(ii)
    use modmpi,     only : myidx
    use modglobal,  only : imax,itot
    implicit none
    integer, intent(in) :: ii
    real                :: relxpos

    integer             :: positionx
  
    ! Converting the i position to the real i position by taking the processor number into account
    ! First grid point lies at i = 2. Make position = 0 for first grid point

    positionx = ii + (myidx * imax) - 2

    ! Convert position to relative position 
    relxpos  = (1.0*positionx)/itot

    return
end function

! Return relative position along y-axis given local j index
function relypos(jj)
    use modmpi,     only : myidy
    use modglobal,  only : jmax,jtot
    implicit none
    integer, intent(in) :: jj
    real                :: relypos

    integer             :: positiony

    ! Converting the j position to the real j position by taking the processor number into account
    ! First grid point lies at j = 2. Make position = 0 for first gridpoint
    positiony = jj + (myidy * jmax) - 2

    ! Convert position to relative position 
    relypos  = (1.0*positiony)/jtot

    return
end function

!-----------------------------------------------------------------------------------------
! Note that the source distribution is defined on full levels
SUBROUTINE read_srcdistr(maxnsrc, srcfact, nsrc, srcname)
    use modglobal, only : ifinput, cexpnr
    character(80) readstring
    integer       k, maxnsrc, nsrc, ierr
    real          srcfact(maxnsrc)
    character*(*)      srcname

    open (ifinput,file=srcname//'.inp.'//cexpnr)

    k = 0
    ierr = 0
    ! Initialize srcfact
    srcfact = 0
    do while (ierr == 0)
       read(ifinput, '(A)', iostat=ierr) readstring
       do while (readstring(1:1)=='#')     ! Skip the lines that are commented (like headers)
         read (ifinput,'(A)', iostat=ierr) readstring
       end do
       if (ierr == 0) then
         k = k + 1
         if (k > maxnsrc) then
           STOP 'Number of levels in srcdistr file is too large'
         else
           read(readstring, *) srcfact(k)
         endif
       endif
    end do
    nsrc = k

    close(ifinput)

    !And now, weigh it such that the array of srcfactor values is on average 1 (just in case the user's array does not fulfill that criterion)
    !Note that the weighing is different here than for the padfactor array, as
    !the srcfactors are defined on full levels, rather than half levels
    if (sum(srcfact(1:(nsrc-1)) + srcfact(2:nsrc)) > 0) then 
        srcfact = srcfact * (nsrc) / sum(srcfact(1:nsrc))
    else
        WRITE(*,*) 'The source has no non-zero values in file ', srcname//'.inp.'//cexpnr
        STOP 'Fix your source file'
    endif
    write(*,*) 'DEBUG ', srcfact(:nsrc)

end subroutine
!-----------------------------------------------------------------------------------------

SUBROUTINE read_paddistr(maxnpad, padfact, npad, padname)
    use modglobal, only : ifinput, cexpnr
    character(80) readstring
    integer       k, maxnpad, npad, ierr
    real          padfact(maxnpad)
    character*(*)      padname

    open (ifinput,file=padname//'.inp.'//cexpnr)

    k = 0
    ierr = 0
    ! Initialize padfact
    padfact = 0
    do while (ierr == 0)
       read(ifinput, '(A)', iostat=ierr) readstring
       do while (readstring(1:1)=='#')     ! Skip the lines that are commented (like headers)
         read (ifinput,'(A)', iostat=ierr) readstring
       end do
       if (ierr == 0) then
         k = k + 1
         if (k > maxnpad) then
           STOP 'Number of levels in paddistr file is too large'
         else
           read(readstring, *) padfact(k)
         endif
       endif
    end do
    npad = k

    close(ifinput)

    !And now, weigh it such that the array of averages of 2 adjacent padfactor values is on average 1 (just in case the user's array does not fulfill that criterion)
    padfact = padfact * (npad - 1) * 2 / sum(padfact(1:(npad-1))+padfact(2:npad))

    ! write(*,*) 'Prescribed weighing for plant area density from surface to canopy top (equidistant); normalized if necessary'
    ! do k=1,npad
    !   write (*,*) padfact(k)
    ! end do
end subroutine
!-----------------------------------------------------------------------------------------
  SUBROUTINE initcanopy
    use modmpi,    only : myid, mpi_logical, mpi_integer, my_real, comm3d, mpierr
    use modglobal, only : kmax,ifnamopt, fname_options, ifinput, cexpnr, zh, dzh, dzf, nsv, i2, j2

    implicit none

    integer       ::  ierr, i,j,k, kp, defined_cantypes = 0, canopytype_0 = -1, defined_canpatches = 0
    character(len=1500) :: readbuffer
    logical       ::       found_type
    character(len=5) :: srcname


    namelist/NAMCANOPY/ lcanopy, ncanopy, cd, lai, lpaddistr, npaddistr, &
                        wth_total, wqt_total, wsv_total, wth_can, wqt_can, wsv_can, &
                        wth_alph, wqt_alph, wsv_alph, &
                        lhetcanopy, &
                        lsrcdistr
  
    if(myid==0) then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMCANOPY,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMCANOPY'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMCANOPY'
      endif
      write(6 ,NAMCANOPY)
      close(ifnamopt)
  
      ncanopy = min(ncanopy,kmax)
    endif
  
  
    call MPI_BCAST(lcanopy   ,   1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(lhetcanopy,   1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(ncanopy   ,   1, mpi_integer , 0, comm3d, mpierr)
    call MPI_BCAST(cd        ,   1, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(lai       ,   1, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(lpaddistr ,   1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(npaddistr ,   1, mpi_integer , 0, comm3d, mpierr)
    call MPI_BCAST(wth_total ,   1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(wqt_total ,   1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(wsv_total , 100, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(wth_can   ,   1, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(wqt_can   ,   1, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(wsv_can   , 100, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(wth_alph  ,   1, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(wqt_alph  ,   1, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(wsv_alph  , 100, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(lsrcdistr ,   1, mpi_logical , 0, comm3d, mpierr)

    if (.not. (lcanopy)) return
    
    if (.not. lpaddistr) npaddistr = 11

    allocate(padfactor (npaddistr))
    allocate(ppad      (npaddistr))
    allocate(zpad      (npaddistr))
    allocate(padtemp   (npaddistr))
    allocate(padf      (ncanopy  ))
    allocate(padh      (ncanopy+1))
    allocate(pai       (ncanopy+1))

    if (lhetcanopy) then
      allocate(padfactor_c (max_npad,max_canopy))
      allocate(ppad_c      (max_npad,max_canopy))
      allocate(zpad_c      (max_npad,max_canopy))
      allocate(padtemp_c   (max_npad,max_canopy))
      allocate(padf_c      (max_ncan,max_canopy))
      allocate(padh_c      (max_ncan+1,max_canopy))
      allocate(pai_c       (max_ncan+1,max_canopy))
      allocate(map_c       (i2,j2))
    endif
  
    if (lsrcdistr) then
      allocate(srcfactor (maxsrcdistr,100))
    endif

    ! If heterogeneous canopy, read file that describes land types
    if (lhetcanopy) then
      if (myid == 0) then 
        open (ifinput,file='canopy.inp.'//cexpnr)
        ierr = 0
        do while (ierr == 0)
          read(ifinput, '(A)', iostat=ierr) readbuffer
          if (ierr == 0) then                               !So no end of file is encountered
            if (readbuffer(1:1)=='#') then
              print *,trim(readbuffer)
            else
              print *,trim(readbuffer)
              defined_cantypes = defined_cantypes + 1
              i = defined_cantypes

              read(readbuffer, *, iostat=ierr) canopytype(i), canopyname(i), ncanopy_c(i), &
                   cd_c(i), lai_c(i), lpaddistr_c(i), &
                   wth_can_c(i), wqt_can_c(i), wsv_can_c(1:nsv,i), &
                   wth_alph_c(i), wqt_alph_c(i), wsv_alph_c(1:nsv,i), &
                   wth_total_c(i), wqt_total_c(i), wsv_total_c(1:nsv,i)

              if (cd_c(i) .lt. 0) then
                stop "NAMCANOPY: A CD value in the canopy input file is negative"
              endif
              if (lai_c(i) .lt. 0) then
                stop "NAMCANOPY: A LAI value in the canopy input file is negative"
              endif
              if (wth_alph_c(i) .lt. 0) then
                stop "NAMCANOPY: The wth extinction factor is negative"
              endif
              if (wqt_alph_c(i) .lt. 0) then
                stop "NAMCANOPY: The wqt extinction factor is negative"
              endif
              do k = 1,nsv
                if (wsv_alph_c(k,i) .lt. 0) then
                  stop "NAMCANOPY: The wsc extinction factor is negative for a scalar"
                endif
              enddo

              if (canopytype(i) .eq. 0) canopytype_0 = i
              do j = 1, (i-1)
                if (canopytype(i) .eq. canopytype(j)) stop "NAMCANOPY: Two canopy types have the same type number"
              enddo
            endif
          endif
        enddo
        close(ifinput)
        if (canopytype_0 .eq. -1) then
          stop "NAMCANOPY: no standard canopy type (0) is defined"
        else
           print "(a,i2,a,i2)","There are ",defined_cantypes,&
           " canopy types defined in the canopy input file. The standard land type is defined by line ",canopytype_0
        endif
        act_max_ncan = maxval(ncanopy_c)
        ncanopies = defined_cantypes ! ncanopies is a global variable
      endif
      call MPI_BCAST(ncanopies   ,              1, mpi_integer , 0, comm3d, mpierr)
      call MPI_BCAST(act_max_ncan ,              1, mpi_integer , 0, comm3d, mpierr)
      call MPI_BCAST(defined_cantypes ,         1, mpi_integer , 0, comm3d, mpierr)
      call MPI_BCAST(ncanopy_c   ,     max_canopy, mpi_integer , 0, comm3d, mpierr)
      call MPI_BCAST(canopytype  ,     max_canopy, mpi_integer , 0, comm3d, mpierr)
      call MPI_BCAST(cd_c        ,     max_canopy, my_real     , 0, comm3d, mpierr)
      call MPI_BCAST(lai_c       ,     max_canopy, my_real     , 0, comm3d, mpierr)
      call MPI_BCAST(lpaddistr_c ,     max_canopy, mpi_logical , 0, comm3d, mpierr)
      call MPI_BCAST(npaddistr_c ,     max_canopy, mpi_integer , 0, comm3d, mpierr)
      call MPI_BCAST(wth_total_c ,     max_canopy, mpi_logical , 0, comm3d, mpierr)
      call MPI_BCAST(wqt_total_c ,     max_canopy, mpi_logical , 0, comm3d, mpierr)
      call MPI_BCAST(wsv_total_c , 100*max_canopy, mpi_logical , 0, comm3d, mpierr)
      call MPI_BCAST(wth_can_c   ,     max_canopy, my_real     , 0, comm3d, mpierr)
      call MPI_BCAST(wqt_can_c   ,     max_canopy, my_real     , 0, comm3d, mpierr)
      call MPI_BCAST(wsv_can_c   , 100*max_canopy, my_real     , 0, comm3d, mpierr)
      call MPI_BCAST(wth_alph_c  ,     max_canopy, my_real     , 0, comm3d, mpierr)
      call MPI_BCAST(wqt_alph_c  ,     max_canopy, my_real     , 0, comm3d, mpierr)
      call MPI_BCAST(wsv_alph_c  , 100*max_canopy, my_real     , 0, comm3d, mpierr)
    endif



    ! Determination of padfactor: relative weighing of plant area distribution inside canopy; equidistant from surface to canopy top
    if (lhetcanopy) then
      if (myid == 0) then 
        do i = 1,defined_cantypes
          if (lpaddistr_c(i)) then
             call read_paddistr(max_npad, padfactor_c(:,i), npaddistr_c(i), trim(canopyname(i)))
          else
             npaddistr_c(i) = 11
             padfactor_c(i,1:11) = (/ 0.4666666666666667, &
                       0.5307086614173228, &
                       0.6792650918635170, &
                       0.9548556430446193, &
                       1.3154855643044620, &
                       1.5490813648293960, &
                       1.5916010498687660, &
                       1.5275590551181100, &
                       1.2944881889763780, &
                       0.3236220472440945, &
                       0.0000000000000000  /)
          endif
        enddo
      endif
      call MPI_BCAST(padfactor_c, max_npad*max_canopy, my_real    , 0, comm3d, mpierr)
      call MPI_BCAST(npaddistr_c,          max_canopy, mpi_integer, 0, comm3d, mpierr)
    else
      if (lpaddistr) then  !< Profile prescribed by user in the file paddistr.inp.<expnr>
        if (myid==0) then
          call read_paddistr(npaddistr, padfactor, npaddistr, 'paddistr')
        endif

        call MPI_BCAST(padfactor, npaddistr, my_real , 0, comm3d, mpierr)

      else                 !< Standard profile fron Ned Patton
        padfactor = (/ 0.4666666666666667, &
                       0.5307086614173228, &
                       0.6792650918635170, &
                       0.9548556430446193, &
                       1.3154855643044620, &
                       1.5490813648293960, &
                       1.5916010498687660, &
                       1.5275590551181100, &
                       1.2944881889763780, &
                       0.3236220472440945, &
                       0.0000000000000000  /)
      endif
      if (lsrcdistr) then
      ! First only implement reading of source distribution for scalars
          do i=1,nsv 
              write(srcname,'(a,i3.3)') 'sv', i
              call read_srcdistr(maxsrcdistr, srcfactor(:,i), nsrcdistr(i), srcname)
          enddo
          call MPI_BCAST(srcfactor, 100*maxsrcdistr, my_real , 0, comm3d, mpierr)
          call MPI_BCAST(nsrcdistr, maxsrcdistr, my_real , 0, comm3d, mpierr)
      endif
    endif

    ! Make PAD profiles on grid levels 
    ! This is done locally on each processor
    if (lhetcanopy) then
      do i = 1,defined_cantypes
         f_lai_h = lai_c(i) / zh(1+ncanopy_c(i)) ! LAI of canopy divided by height of the top of the canopy
         ppad_c(:,i)  = f_lai_h * padfactor_c(:,i) ! prescribed PAD-values

         ! Define heights for pad
         zpad_c(:,i) = 0
         do k=1,npaddistr_c(i)
           zpad_c(k,i) = zh(1+ncanopy_c(i)) * real(k-1)/real(npaddistr_c(i)-1)
         end do
         call spline(zpad_c(:,i),ppad_c(:,i),npaddistr_c(i),padtemp_c(:,i))

         ! Interpolate from PAD definition levels to grid levels in the model
         padh_c(:,i) = 0
         do k=1,(1+ncanopy_c(i))
           call splint(zpad_c(:,i),ppad_c(:,i),padtemp_c(:,i),npaddistr_c(i),zh(k),padh_c(k,i))
         end do

         ! Interpolate plant area (index) density to full levels
         padf_c(:,i) = 0
         do k=1,ncanopy_c(i)
            kp      = k+1
            padf_c(k,i) = ( dzh(kp) * padh_c(k,i) + dzh(k) * padh_c(kp,i) ) / ( dzh(k) + dzh(kp) )
         enddo

         ! Vertically integrate the plant area density to arrive at plant area index 
         pai_c(:,i) = 0.0
         do k=ncanopy_c(i),1,-1
           pai_c(k,i) = pai_c(k+1,i) + dzf(k) * padf_c(k,i)
         end do
      enddo
    else
      f_lai_h = lai / zh(1+ncanopy) ! LAI of canopy divided by height of the top of the canopy
  
      ppad    = f_lai_h * padfactor ! prescribed PAD-values
      do k=1,npaddistr
        zpad(k) = zh(1+ncanopy) * real(k-1)/real(npaddistr-1)
      end do

      !interpolate the PAD values to the LES grid
      call spline(zpad,ppad,npaddistr,padtemp)
      do k=1,(1+ncanopy)
        call splint(zpad,ppad,padtemp,npaddistr,zh(k),padh(k))
      end do

      ! Interpolate plant area (index) density to full levels
      do k=1,ncanopy
        kp      = k+1
        padf(k) = ( dzh(kp) * padh(k) + dzh(k) * padh(kp) ) / ( dzh(k) + dzh(kp) )
      end do

      ! Vertically integrate the plant area density to arrive at plant area index 
      pai = 0.0
      do k=ncanopy,1,-1
        pai(k) = pai(k+1) + dzf(k) * padf(k)
      end do
    endif


    ! If heterogeneous, now read the patch definitions
    if (lhetcanopy) then
      if (myid == 0) then
        open (ifinput,file='canopypatch.inp.'//cexpnr)
        ierr = 0
        do while (ierr == 0)
          read(ifinput, '(A)', iostat=ierr) readbuffer
          if (ierr == 0) then                               !So no end of file is encountered
            if (readbuffer(1:1)=='#') then
              print *,trim(readbuffer)
            else
              print *,trim(readbuffer)
              defined_canpatches = defined_canpatches + 1
              i = defined_canpatches

              read(readbuffer, *, iostat=ierr) patchtype(i), &
                   minx_p(i), maxx_p(i), miny_p(i), maxy_p(i)

              ! Check that patchtype is known
              found_type = .false.
              do j = 1,defined_cantypes
                  if (patchtype(i) == canopytype(j)) found_type = .true.
              enddo
              if (.not. found_type) then
                stop "NAMCANOPY: unkown canopy type in patch definition"
              endif

              if (minx_p(i) .lt. 0) then
                stop "NAMCANOPY: A minx for patch is negative"
              endif
              if (maxx_p(i) .lt. 0) then
                stop "NAMCANOPY: A maxx for patch is negative"
              endif
              if (miny_p(i) .lt. 0) then
                stop "NAMCANOPY: A miny for patch is negative"
              endif
              if (maxy_p(i) .lt. 0) then
                stop "NAMCANOPY: A maxy for patch is negative"
              endif
              if (minx_p(i) .gt. 1) then
                stop "NAMCANOPY: A minx for patch is > 1 "
              endif
              if (maxx_p(i) .gt. 1) then
                stop "NAMCANOPY: A maxx for patch is > 1 "
              endif
              if (miny_p(i) .gt. 1) then
                stop "NAMCANOPY: A miny for patch is > 1 "
              endif
              if (maxy_p(i) .gt. 1) then
                stop "NAMCANOPY: A maxy for patch is > 1 "
              endif
              if (minx_p(i) > maxx_p(i)) then
                stop "NAMCANOPY: A minx > maxx for patch"
              endif
              if (miny_p(i) > maxy_p(i)) then
                stop "NAMCANOPY: A miny > maxy for patch"
              endif

            endif
          endif
        enddo
        close(ifinput)
      endif
      call MPI_BCAST(defined_canpatches, 1, mpi_integer, 0, comm3d, mpierr)
      call MPI_BCAST(patchtype,   mpatch_c, mpi_integer, 0, comm3d, mpierr)
      call MPI_BCAST(minx_p,      mpatch_c, my_real,     0, comm3d, mpierr)
      call MPI_BCAST(maxx_p,      mpatch_c, my_real,     0, comm3d, mpierr)
      call MPI_BCAST(miny_p,      mpatch_c, my_real,     0, comm3d, mpierr)
      call MPI_BCAST(maxy_p,      mpatch_c, my_real,     0, comm3d, mpierr)
    endif

    if (lhetcanopy) then
        ! Now fill map with canopy types
        ! First the standard canopy type
        map_c = canopytype(canopytype_0)
        do k=1,defined_canpatches
           do j=1,j2
               if ((relypos(j) > miny_p(k)) .and. (relypos(j) <= maxy_p(k))) then
                  do i=1,i2
                     if ((relxpos(i) > minx_p(k)) .and. (relxpos(i) <= maxx_p(k))) then
                        map_c(i,j) = patchtype(k)
                     endif
                  enddo
               endif
           enddo
        enddo
    endif
    return
  end subroutine initcanopy
  
  subroutine canopy
    use modfields,   only : up,vp,wp,e12p,thlp,qtp,svp
    use modsurfdata, only : thlflux, qtflux, svflux
    use modglobal,   only : nsv,i2,j2

    implicit none

    integer n, i
    real :: zeroar(i2,j2)
    
    zeroar = 0.0
  
    if (.not. (lcanopy)) return

!!  Momentum affected by trees
    call canopyu(up)
    call canopyv(vp)
    call canopyw(wp)
!!  TKE affected by trees
    call canopye(e12p)
!!  Emissions of heat, moisture and scalars by trees (effect on center of the grid)
 
    if (lhetcanopy) then
      do i=1,ncanopies
          if (wth_total_c(i)) then
            call canopyc_het(thlp,wth_can_c(i),thlflux,wth_alph_c(i),pai_c(:,i),i)
          else
            call canopyc_het(thlp,wth_can_c(i), zeroar,wth_alph_c(i),pai_c(:,i),i)
          endif
          if (wqt_total_c(i)) then
            call canopyc_het( qtp,wqt_can_c(i),qtflux,wqt_alph_c(i),pai_c(:,i),i)
          else
            call canopyc_het( qtp,wqt_can_c(i),zeroar,wqt_alph_c(i),pai_c(:,i),i)
          endif
          do n=1,nsv
            if (wsv_total_c(n,i)) then
              call canopyc_het(svp(:,:,:,n),wsv_can_c(n,i),svflux(:,:,n),wsv_alph_c(n,i),pai_c(:,i),i)
            else
              call canopyc_het(svp(:,:,:,n),wsv_can_c(n,i),       zeroar,wsv_alph_c(n,i),pai_c(:,i),i)
            endif
          end do
       enddo
    else
      
      if (wth_total) then
          call canopyc(thlp,wth_can,thlflux,wth_alph,pai)
      else
          call canopyc(thlp,wth_can, zeroar,wth_alph,pai)
      endif
      if (wqt_total) then
          call canopyc( qtp,wqt_can,qtflux,wqt_alph,pai)
      else
          call canopyc( qtp,wqt_can,zeroar,wqt_alph,pai)
      endif
      ! Here I have to implement a new canopyc that can use a prescribed source
      ! distibution, having srcdistr as an input, rather than the alpha and pai
      if (lsrcdistr) then
         do n=1,nsv
            if (wsv_total(n)) then
              call canopycd(svp(:,:,:,n),wsv_can(n),svflux(:,:,n),nsrcdistr(n), srcfactor(:,n))
            else
              call canopycd(svp(:,:,:,n),wsv_can(n),zeroar       ,nsrcdistr(n), srcfactor(:,n))
            endif
         end do
      else
         do n=1,nsv
            if (wsv_total(n)) then
              call canopyc(svp(:,:,:,n),wsv_can(n),svflux(:,:,n),wsv_alph(n),pai)
            else
              call canopyc(svp(:,:,:,n),wsv_can(n),       zeroar,wsv_alph(n),pai)
            endif
         end do
      endif
    endif
    return
  end subroutine canopy
  
  subroutine exitcanopy
    implicit none

    if (.not. (lcanopy)) return

    deallocate(padfactor)
    deallocate(ppad     )
    deallocate(zpad     )
    deallocate(padtemp  )
    deallocate(padf     )
    deallocate(padh     )
    deallocate(pai      )
    return
  end subroutine exitcanopy
  
  subroutine canopyu (putout)
    use modglobal, only  : i1, ih, j1, j2, jh, k1, cu, cv, dzh, imax, jmax
    use modfields, only  : u0, v0, w0
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ucor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: vcor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ftau  (imax,jmax)
    integer             :: k, kp, n

    ucor = u0 + cu
    vcor = v0 + cv

    if (lhetcanopy) then
      do k=1,act_max_ncan
         kp   = k+1
         ftau = 0
         do n=1,ncanopies 
            where (map_c(2:i1,2:j1) == canopytype(n)) 
                ftau = cd_c(n) * padf_c(k,n) * sqrt(ucor(2:i1,2:j1,k)**2 +  &
                ((vcor(1:imax,2:j1,k)+vcor(2:i1,2:j1,k)+vcor(1:imax,3:j2,k)+vcor(2:i1,3:j2,k))/4)**2 + &
                ((dzh(kp)*(w0(1:imax,2:j1,k)+w0(2:i1,2:j1,k))+dzh(k)*(w0(1:imax,2:j1,kp)+w0(2:i1,2:j1,kp)))/(2*(dzh(k)+dzh(kp))))**2 &
                )
            end where
         enddo
         putout(2:i1,2:j1,k) = putout(2:i1,2:j1,k) - ftau * ucor(2:i1,2:j1,k)
      end do
    else
      do k=1,ncanopy
        kp   = k+1
        ftau = cd * padf(k) * sqrt(ucor(2:i1,2:j1,k)**2 +  &
                ((vcor(1:imax,2:j1,k)+vcor(2:i1,2:j1,k)+vcor(1:imax,3:j2,k)+vcor(2:i1,3:j2,k))/4)**2 + &
                ((dzh(kp)*(w0(1:imax,2:j1,k)+w0(2:i1,2:j1,k))+dzh(k)*(w0(1:imax,2:j1,kp)+w0(2:i1,2:j1,kp)))/(2*(dzh(k)+dzh(kp))))**2 &
                )
      
        putout(2:i1,2:j1,k) = putout(2:i1,2:j1,k) - ftau * ucor(2:i1,2:j1,k)
      end do
    endif

    return
  end subroutine canopyu
    
  subroutine canopyv (putout)
    use modglobal, only  : i1, i2, ih, j1, jh, k1, cu, cv, dzh, imax, jmax
    use modfields, only  : u0, v0, w0
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ucor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: vcor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ftau  (imax,jmax)
    integer             :: k, kp, n

    ucor = u0 + cu
    vcor = v0 + cv

    if (lhetcanopy) then
      do k=1,act_max_ncan
         kp   = k+1
         ftau = 0
         do n=1,ncanopies 
            where (map_c(2:i1,2:j1) == canopytype(n)) 
                  ftau = cd_c(n) * padf_c(k,n) * sqrt(vcor(2:i1,2:j1,k)**2 +  &
                  ((ucor(3:i2,2:j1,k)+ucor(2:i1,2:j1,k)+ucor(3:i2,1:jmax,k)+ucor(2:i1,1:jmax,k))/4)**2 + &
                  ((dzh(kp)*(w0(2:i1,1:jmax,k)+w0(2:i1,2:j1,k))+dzh(k)*(w0(2:i1,1:jmax,kp)+w0(2:i1,2:j1,kp)))/(2*(dzh(k)+dzh(kp))))**2 &
                  )
            end where
         enddo
         putout(2:i1,2:j1,k) = putout(2:i1,2:j1,k) - ftau * vcor(2:i1,2:j1,k)
      end do
    else
      do k=1,ncanopy
        kp   = k+1
        ftau = cd * padf(k) * sqrt(vcor(2:i1,2:j1,k)**2 +  &
                  ((ucor(3:i2,2:j1,k)+ucor(2:i1,2:j1,k)+ucor(3:i2,1:jmax,k)+ucor(2:i1,1:jmax,k))/4)**2 + &
                  ((dzh(kp)*(w0(2:i1,1:jmax,k)+w0(2:i1,2:j1,k))+dzh(k)*(w0(2:i1,1:jmax,kp)+w0(2:i1,2:j1,kp)))/(2*(dzh(k)+dzh(kp))))**2 &
                  )
      
        putout(2:i1,2:j1,k) = putout(2:i1,2:j1,k) - ftau * vcor(2:i1,2:j1,k)
      end do
    endif
    return
  end subroutine canopyv
    
  subroutine canopyw (putout)
    use modglobal, only  : i1, i2, ih, j1, j2, jh, k1, cu, cv, dzh, dzf, imax, jmax
    use modfields, only  : u0, v0, w0
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ucor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: vcor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ftau  (imax,jmax)
    integer             :: k, km, n

    ucor = u0 + cu
    vcor = v0 + cv

    if (lhetcanopy) then
      do k=1,act_max_ncan
         km   = k-1
         ftau = 0
         do n=1,ncanopies 
            where (map_c(2:i1,2:j1) == canopytype(n)) 
               ftau = cd_c(n) * padh_c(k,n) * sqrt(w0(2:i1,2:j1,k)**2 +  &
                  ((dzf(km)*(ucor(2:i1,2:j1,k)+ucor(3:i2,2:j1,k))+dzf(k)*(ucor(2:i1,2:j1,km)+ucor(3:i2,2:j1,km)))/(4*dzh(k)))**2 + &
                  ((dzf(km)*(vcor(2:i1,2:j1,k)+vcor(2:i1,3:j2,k))+dzf(k)*(vcor(2:i1,2:j1,km)+vcor(2:i1,3:j2,km)))/(4*dzh(k)))**2 &
                  )
            end where
         enddo
         putout(2:i1,2:j1,k) = putout(2:i1,2:j1,k) - ftau * w0(2:i1,2:j1,k)
      end do
    else
      do k=2,(ncanopy+1)
        km   = k-1
        ftau = cd * padh(k) * sqrt(w0(2:i1,2:j1,k)**2 +  &
                  ((dzf(km)*(ucor(2:i1,2:j1,k)+ucor(3:i2,2:j1,k))+dzf(k)*(ucor(2:i1,2:j1,km)+ucor(3:i2,2:j1,km)))/(4*dzh(k)))**2 + &
                  ((dzf(km)*(vcor(2:i1,2:j1,k)+vcor(2:i1,3:j2,k))+dzf(k)*(vcor(2:i1,2:j1,km)+vcor(2:i1,3:j2,km)))/(4*dzh(k)))**2 &
                  )
        
        putout(2:i1,2:j1,k) = putout(2:i1,2:j1,k) - ftau * w0(2:i1,2:j1,k)
      end do
    endif
    return
  end subroutine canopyw
  
  subroutine canopye (putout)
    use modglobal, only  : i1, i2, ih, j1, j2, jh, k1, cu, cv, dzh, imax, jmax
    use modfields, only  : u0, v0, w0, e120
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ucor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: vcor  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: ftau  (imax,jmax)
    integer             :: k, kp,  n

    ucor = u0 + cu
    vcor = v0 + cv

    if (lhetcanopy) then
      do k=1,act_max_ncan
         kp   = k+1
         ftau = 0
         do n=1,ncanopies 
            where (map_c(2:i1,2:j1) == canopytype(n)) 
               ftau = cd_c(n) * padf_c(k,n) * sqrt(((ucor(3:i2,2:j1,k)+ucor(2:i1,2:j1,k))/2)**2 + &
                                                   ((vcor(2:i1,3:j2,k)+vcor(2:i1,2:j1,k))/2)**2 + &
                  ((dzh(kp)*w0(2:i1,2:j1,k)+dzh(k)*w0(2:i1,2:j1,kp))/(dzh(k)+dzh(kp)))**2 &
                  )
            end where
         enddo
         putout(2:i1,2:j1,k) = putout(2:i1,2:j1,k) - e120(2:i1,2:j1,k) * ftau
      end do
    else
      do k=1,ncanopy
        kp   = k+1
        ftau = cd * padf(k) * sqrt(((ucor(3:i2,2:j1,k)+ucor(2:i1,2:j1,k))/2)**2 + &
                                   ((vcor(2:i1,3:j2,k)+vcor(2:i1,2:j1,k))/2)**2 + &
                  ((dzh(kp)*w0(2:i1,2:j1,k)+dzh(k)*w0(2:i1,2:j1,kp))/(dzh(k)+dzh(kp)))**2 &
                  )
        putout(2:i1,2:j1,k) = putout(2:i1,2:j1,k) - e120(2:i1,2:j1,k) * ftau
      end do
    endif
    return
  end subroutine canopye
 
  subroutine canopyc_het (putout, flux_top, flux_surf, alpha, pai, ican)
    use modglobal, only  : i1, i2, ih, j1, j2, jh, k1, dzf, imax, jmax
    use modfields, only  : rhobh, rhobf
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real, intent(in   ) :: flux_top
    real, intent(in   ) :: flux_surf(i2,j2)
    real, intent(in   ) :: alpha
    real, intent(in   ) :: pai(max_ncan+1)
    integer, intent(in) :: ican
    real                :: flux_net (i2,j2)
    integer             :: k
    real                :: integratedcontribution(imax,jmax,max_ncan+1), tendency(imax,jmax,max_ncan)

    flux_net                        = flux_top * rhobh(ncanopy_c(ican)+1) - flux_surf * rhobh(1)
    integratedcontribution(:,:,:)   = 0.0
    tendency(:,:,:) = 0.0

    do k=2,(ncanopy_c(ican)+1)
      integratedcontribution(:,:,k) = flux_net(2:i1,2:j1) * exp(- alpha * pai(k))
    end do
    do k=1,ncanopy_c(ican)
      where (map_c(2:i1,2:j1) == canopytype(ican))
          tendency(:,:,k) = ( integratedcontribution(:,:,(k+1)) - integratedcontribution(:,:,k) ) / ( rhobf(k) * dzf(k) )
      elsewhere
          tendency(:,:,k) = 0
      end where
    end do

    putout(2:i1,2:j1,1:ncanopy_c(ican)) = putout(2:i1,2:j1,1:ncanopy_c(ican)) + tendency(:,:,1:ncanopy_c(ican))

    return
  end subroutine canopyc_het

  subroutine canopyc (putout, flux_top, flux_surf, alpha, pai)
    use modglobal, only  : i1, i2, ih, j1, j2, jh, k1, dzf, imax, jmax
    use modfields, only  : rhobh, rhobf
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real, intent(in   ) :: flux_top
    real, intent(in   ) :: flux_surf(i2,j2)
    real, intent(in   ) :: alpha
    real, intent(in   ) :: pai(ncanopy+1)
    real                :: flux_net (i2,j2)
    integer             :: k
    real                :: integratedcontribution(imax,jmax,ncanopy+1), tendency(imax,jmax,ncanopy)
    
    flux_net                        = flux_top * rhobh(ncanopy+1) - flux_surf * rhobh(1)
    integratedcontribution(:,:,1)   = 0.0

    do k=2,(ncanopy+1)
      integratedcontribution(:,:,k) = flux_net(2:i1,2:j1) * exp(- alpha * pai(k))
    end do
    do k=1,ncanopy
      tendency(:,:,k) = ( integratedcontribution(:,:,(k+1)) - integratedcontribution(:,:,k) ) / ( rhobf(k) * dzf(k) )
    end do

    putout(2:i1,2:j1,1:ncanopy) = putout(2:i1,2:j1,1:ncanopy) + tendency

    return
  end subroutine canopyc

! Scalar soure with prescribed source distribution
  subroutine canopycd (putout, flux_top, flux_surf, npoints, srcdistr)
    use modglobal, only  : i1, i2, ih, j1, j2, jh, k1, dzf, imax, jmax, zh
    use modfields, only  : rhobh, rhobf
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real, intent(in   ) :: flux_top
    real, intent(in   ) :: flux_surf(i2,j2)
    real, intent(in   ) :: srcdistr(:)
    integer, intent(in) :: npoints
    real                :: flux_net (i2,j2)
    integer             :: k
    real                :: integratedcontribution(imax,jmax,ncanopy+1), tendency(imax,jmax,ncanopy)
    
    flux_net                        = flux_top * rhobh(ncanopy+1) - flux_surf * rhobh(1)

    do k=1,npoints
      tendency(:,:,k) = srcdistr(k)*flux_net(2:i1,2:j1) / zh(1+ncanopy)
    end do

    putout(2:i1,2:j1,1:npoints) = putout(2:i1,2:j1,1:ncanopy) + tendency

    return
  end subroutine canopycd

! ======================================================================
! subroutine from Ned Patton ~ adapted where obvious
! ======================================================================
      subroutine spline(x,y,n,y2)
      implicit none

      integer n
      real :: yp1 = 2.0e30, ypn = 2.0e30  !< Values needed for the fit
      real x(n), y(n), y2(n)
      integer i, k
      real p, qn, sig, un
      real, allocatable :: u(:)

      allocate(u(max(n-1,1)))

      if(yp1 > .99e30) then
        y2(1) = 0.0
        u(1)  = 0.0
      else
        y2(1) = -0.5
        u(1)  = (3./(x(2) - x(1)))*((y(2) - y(1))/(x(2) - x(1)) - yp1)
      endif
      if(ypn > .99e+30) then
        qn = 0.0
        un = 0.0
      else
        qn = 0.5
        un = (3.0/(x(n) - x(n-1)))* (ypn - (y(n) - y(n-1))/(x(n) - x(n-1)))
      endif

      do i=2,n-1
        sig   = (x(i) - x(i-1))/(x(i+1) - x(i-1))
        p     = sig*y2(i-1) + 2.0
        y2(i) = (sig - 1.0)/p
        u(i)  = (6.0*((y(i+1) - y(i))/(x(i+1) - x(i)) - (y(i) - y(i-1)) &
              / (x(i) - x(i-1)))/(x(i+1) - x(i-1)) - sig*u(i-1))/p
      enddo

      y2(n) = (un - qn*u(n-1))/(qn*y2(n-1) + 1.0)

      do k=n-1,1,-1
        y2(k) = y2(k)*y2(k+1) + u(k)
      enddo

      deallocate(u)

      return
      end subroutine spline

! ======================================================================
! subroutine from Ned Patton ~ adapted where obvious
! ======================================================================

      subroutine splint(xa,ya,y2a,n,x,y)
      implicit none

      integer n
      real x,y, xa(n), y2a(n), ya(n)
      integer k,khi,klo
      real a,b,h
      klo = 1
      khi = n
    1 continue
      if(khi - klo > 1) then
         k = (khi + klo)/2
         if(xa(k) > x) then
           khi = k
         else
           klo = k
         endif
         go to 1
      endif
      h = xa(khi) - xa(klo)
      a = (xa(khi) - x)/h
      b = (x - xa(klo))/h
      y = a*ya(klo) + b*ya(khi) + ((a**3 - a)*y2a(klo) + (b**3 - b)*y2a(khi))*(h**2)/6.0

      return
      end subroutine splint

end module modcanopy
