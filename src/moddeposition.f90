!
! Copyright (c) 2020-2022 TNO
!
! This file is part of DALES
!
! DALES is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DALES.  If not, see <http://www.gnu.org/licenses/>.
!

module moddeposition
    use modlsm,    only : tile, nlu, llsm
    use modfields, only : svp   ! tracer tendency array
    use modglobal, only : nsv, i1, j1

    implicit none

    public  :: initdrydep, drydep, exitdrydep, ldrydep

    logical :: ldrydep         ! On/Off switch dry deposition 
    real, allocatable    :: depfield(:,:,:) ! deposition tendencies (i,j,sv)
    integer  :: svskip   =  0        ! no. scalars to exclude for deposition


contains

subroutine drydep
  implicit none
  integer :: ilu

  if (.not. llsm .or. .not. ldrydep) return ! Dry deposition cannot be run if LSM not activated

  !! for testing only
  do ilu=1, nlu
    write(*,*) 'luname ', trim(tile(ilu)%luname)
    write(*,*) 'frac   ', sum(tile(ilu)%base_frac)/(max(1,size(tile(ilu)%base_frac)))
    write(*,*) 'z0h    ', sum(tile(ilu)%z0h)/(max(1,size(tile(ilu)%z0h)))
    write(*,*) 'R_inc_b', sum(tile(ilu)%R_inc_b)/(max(1,size(tile(ilu)%R_inc_b)))
  end do

  depfield = 0
  ! call to drydep subroutines here
  
  ! tracer tendency due to deposition
  ! todo: convert units of depfield to ug g-1 s-1
  svp(2:i1, 2:j1, 1, svskip+1:nsv) = svp(2:i1, 2:j1, 1, svskip+1:nsv) + &
    depfield(2:i1, 2:j1, svskip+1:nsv) 
   
end subroutine drydep    

subroutine readdrydep
    implicit none
  ! read netCDF tables with
    ! species properties
      ! diffusivity
      ! Rsoilwet
      ! Rsoilfrozen
    ! land use properties
      ! IncPAR parameters
      ! Gs parameters
      ! SAI parameters
      ! fmin parameter
      ! alpha parameters
      ! temperature parameters
      ! Gsmax
      ! VPD parameters
      ! Rsoil parameters (per LU type AND species)
      ! gamma_stom
      ! gamma_soil
    ! radiation parameters
    ! leaf surface resistance parameters


end subroutine readdrydep

subroutine initdrydep
  ! read namelist    
  ! read LU parameters from file
  ! init drydep fields

  use modglobal, only : i2, j2, kmax, nsv, ifnamopt, fname_options, checknamelisterror
  use modmpi,    only : myid, comm3d, mpi_logical, mpi_integer, mpi_character

  implicit none

  ! Auxiliary variables
  integer  :: ierr
  integer  :: iname

  ! ---------------------------------------------------------------------!
  ! Namelist variables                                                   !
  ! ---------------------------------------------------------------------!
   
  logical  :: l_deposition = .false. ! logical deposition switch

  character(len = 6), dimension(100) :: & 
              depnames = (/ ('      ', iname=1, 100) /) ! list with scalar names,
                          ! each name must(!) be 6 characters long for now  

  ! --- Read & broadcast namelist DEPOSITION -----------------------------------
  namelist/NAMDEPOSITION/ ldrydep, svskip, depnames 

  if (myid == 0) then

    open(ifnamopt,file=fname_options,status='old',iostat=ierr)
    read (ifnamopt,NAMDEPOSITION,iostat=ierr)
    call checknamelisterror(ierr, ifnamopt, 'NAMDEPOSITION')
    write(6, NAMDEPOSITION)
    close(ifnamopt)

  endif

  call mpi_bcast(ldrydep,           1, mpi_logical,   0, comm3d, ierr)
  call mpi_bcast(svskip,            1, mpi_integer,   0, comm3d, ierr)
  call mpi_bcast(depnames(1:100), 100, mpi_character, 0, comm3d, ierr)

  ! --- Local pre-calculations and settings
  if (.not. (ldrydep)) return

  allocate(depfield(i2, j2, svskip+1:nsv))

  !call deposition calculation

end subroutine initdrydep

subroutine exitdrydep
  implicit none

  if (.not. llsm .or. .not. ldrydep) return

  deallocate(depfield)

end subroutine exitdrydep

subroutine drydep_Rc
    
  ! Calculate canopy resistance.

  ! In case there is no snow:

  ! .. math:: R_c = \left(\frac{1}{R_w} + \frac{1}{R_{soil,eff}} + \frac{1}{R_s}\right)^{-1}

  ! for :math:`R_s`, you take :math:`R_{stom} + R_{mes}`, but DEPAC does currently not provide values for mesophyll resistance
  ! since these are usually very low (and hence set to 0).

  ! And if there is snow:

  ! .. math:: 

  !     R_c=
  !     \begin{cases}
  !         500 &,\textrm{if}\ T < -1^\circ C \\
  !         70(2-T) &,\textrm{if}\ -1 \leq T \leq 1^\circ C \\
  !         70 &,\textrm{if}\ T > 1^\circ C \\
  !     \end{cases}

  ! Note
  ! -----
  ! When passing parameters as arrays, their lengths should be equal.

  ! Parameters
  ! ----------
  ! snow_present : bool
  !     Flag to represent snow presence
  ! tskin : float 
  !     The temperature of the canopy (:math:`^\circ C`)
  ! R_w : float 
  !     External leaf surface or water layer resistance (a.k.a. cuticular resistance) (:math:`s\ m^{-1}`)
  ! R_soileff : float 
  !     Effective soil resistance (:math:`s\ m^{-1}`)
  ! R_s : float 
  !     Stomatal resistance (:math:`s\ m^{-1}`)

  ! Returns
  ! -------
  ! R_c : float 
  !     Canopy resistance (:math:`s\ m^{-1}`)
  
  use modglobal, only : i1, j1, i2, j2
  use modfields, only : thl0, exnf
  implicit none

  !type(T_lsm_tile), intent(inout) :: tile
  integer :: i, j
  logical :: snow_present 
  real    :: R_w, R_soileff, R_s, tskin
  real    :: R_c(i2,j2)

  snow_present = .false.
  do j=2,j1
    do i=2,i1
      tskin = thl0(i,j,1) * exnf(1)
      !todo: get tksin from tile
      if (snow_present) then
        if (tskin < -1.) then
          R_c(i,j) = 500
        elseif ( tskin >= -1. .and. tskin <= 1.) then
          R_c(i,j) = 70 * (2-tskin)  
        elseif (tskin > 1.) then
          R_c(i,j) = 70
        endif
      else  
        ! no snow, normal calculation
        R_c = (1 / R_w + 1 / R_soileff + 1 / R_s) ** (-1)
      endif
    enddo
  enddo
end subroutine drydep_Rc

subroutine drydep_Rb
    implicit none

end subroutine drydep_Rb

subroutine drydep_Ra
    implicit none

end subroutine drydep_Ra

pure function f1(Sw_in) result(res)
    ! Dummy pure function for testing purposes
    implicit none
    real, intent(in) :: Sw_in
    real             :: res  
    ! Constants f1 calculation:
    real, parameter :: a_f1 = 0.81
    real, parameter :: b_f1 = 0.004
    real, parameter :: c_f1 = 0.05

    res = min(1., (b_f1*Sw_in + c_f1) / (a_f1 * (b_f1*Sw_in + 1.)))

end function f1

end module moddeposition
