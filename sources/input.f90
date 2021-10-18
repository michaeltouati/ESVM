!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!             ElectroStatic Vlasov-Maxwell (ESVM) code              !!
!!                                                                   !!
!! Copyright © 2015 Michaël J TOUATI                                 !!
!!                                                                   !!
!! This file is part of ESVM.                                        !!
!!                                                                   !!
!! ESVM is free software: you can redistribute it and/or modify      !!
!! it under the terms of the GNU General Public License as published !!
!! by the Free Software Foundation, either version 3 of the License, !!
!! or (at your option) any later version.                            !!
!!                                                                   !!
!! ESVM is distributed in the hope that it will be useful,           !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with ESVM. If not, see <https://www.gnu.org/licenses/>.     !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Initial commit written by Michaël J TOUATI - Dec. 2015
module input

use acuracy
use constants
use physics
use omp_lib

implicit none
public :: read_init_parameters

!======================================================================
!              Public variables known by the whole code 
!           (see the file input-deck for their definitions)
!======================================================================

character(len=60), public :: simu
integer,           public :: N_th
real(PR),          public :: T
real(PR),          public :: Z
real(PR),          public :: ni
real(PR),          public :: x_min
real(PR),          public :: x_max
real(PR),          public :: d_x
real(PR),          public :: vx_min
real(PR),          public :: vx_max
real(PR),          public :: d_vx
real(PR),          public :: cfl
real(PR),          public :: L_t
real(PR),          public :: dt_diag
integer,           public :: maxwell
integer,           public :: b_cond
integer,           public :: scheme
real(PR),          public :: b_VL
integer,           public :: perturb
real(PR),          public :: A_pert
real(PR),          public :: k_pert
real(PR),          public :: omega_pert
real(PR),          public :: v_d
integer,           public :: N_x
integer,           public :: N_vx
real(PR),          public :: n0
real(PR),          public :: Te

contains

subroutine read_init_parameters
  implicit none
  Character(len=80)         :: str
  Integer                   :: i, istr
  integer                   :: N_th_max

  ! read the file 'input-deck'
  open (unit=1, file='input-deck',&
  & position='rewind', access='sequential',&
  & form='formatted', action='read',status='old')

  read: do
    
    call get_str(str)
    
    if (str == 'end') exit read
    
    if (str(1:1) /= '#') cycle read
    
    istr = 80
    do i = 1, len(str)
      if (str(i:i) == ' ') then
        istr = i
        exit
      end if
    end do
    
    select case (str(1:istr))
      case ('#simu')
        simu=get_char(str(istr+1:))
      case ('#N_th')
        N_th = get_integer(str(istr+1:))
      case ('#T')  
        T = get_real(str(istr+1:))
      case ('#Z')  
        Z = get_real(str(istr+1:))
      case ('#ni')  
        ni = get_real(str(istr+1:))
      case ('#x_min')  
        x_min = get_real(str(istr+1:))
      case ('#x_max')  
        x_max = get_real(str(istr+1:))
      case ('#d_x')  
        d_x = get_real(str(istr+1:))
      case ('#vx_min')  
        vx_min = get_real(str(istr+1:))
      case ('#vx_max')  
        vx_max = get_real(str(istr+1:))
      case ('#d_vx')  
        d_vx = get_real(str(istr+1:))
      case ('#cfl')
        cfl = get_real(str(istr+1:))  
      case ('#L_t')
        L_t = get_real(str(istr+1:))
      case ('#dt_diag')
        dt_diag = get_real(str(istr+1:))
      case ('#maxwell')
        maxwell = get_integer(str(istr+1:))
      case ('#b_cond')
        b_cond = get_integer(str(istr+1:))
      case ('#scheme')
        scheme = get_integer(str(istr+1:))
      case ('#b_VL')
        b_VL = get_real(str(istr+1:))
      case ('#perturb')
        perturb = get_integer(str(istr+1:))
      case ('#A_pert')
        A_pert = get_real(str(istr+1:))
      case ('#k_pert')
        k_pert = get_real(str(istr+1:))
      case ('#omega_pert')
        omega_pert= get_real(str(istr+1:))
      case ('#v_d')
        v_d = get_real(str(istr+1:))
    end select

  end do read
  close(1)

  N_th_max = omp_get_max_threads()
  if (N_th == 0) N_th = N_th_max

  n0      = Z * ni
  Te      = T * eV

  N_x     = 1 + floor( (x_max-x_min)/d_x)
  N_vx    = 1 + floor( (vx_max-vx_min)/d_vx)

  if ((cfl < 0.) .or. (cfl > 0.999999)) then
    write(*,*)'The CFL parameter must be 0 < cfl < 1!'
    write(*,*)'here cfl = ',cfl
    stop
  end if

  if ((maxwell .ne. 1 ) .and. (maxwell .ne. 2)) then
    write(*,*)'maxwell must be 1 or 2 for Maxwell-Ampere or Maxwell-Poisson solver respectively!'
    write(*,*)'here maxwell = ',maxwell
    stop
  end if

  if ((b_cond .ne. 1 ) .and. (b_cond .ne. 2)) then
    write(*,*)'b_cond must be 1 or 2 for absorbing or periodic boundary conditions respectively!'
    write(*,*)'here b_cond = ',b_cond
    stop
  end if

  if ((b_VL < 1.) .or. (b_VL > 2.)) then
    write(*,*)'the Van Leer schemes parameter b_VL must be 1 < b < 2!'
    write(*,*)'here b_VL = ',b_VL
    stop
  end if

  write(*,*)'------------------------------------------'
  write(*,*)' Recapitulation of simulation parameters :'
  write(*,*)'------------------------------------------'
  write(*,*)'* Simulation :'
  write(*,*)'  simu    = ',trim(simu)
  write(*,*)'* Number of OpenMP threads :'
  write(*,'(A,1I4)')'   N_th    = ',N_th
  write(*,*)'------------------------------------------'
  write(*,*)'* Plasma properties (immobile ions)'
  write(*,'(A,1E21.14,A)')'   T          = ',T,' eV'
  write(*,'(A,1E21.14,A)')'              = ',Te,' K'
  write(*,'(A,1E21.14)')'   Z          = ',Z
  write(*,'(A,1E21.14,A)')'   ni         = ',ni,' /cm3'
  write(*,*)'* Deduced ESVM units'
  write(*,'(A,1E21.14,A)')'   Debye      = ',Debye(n0,Te),' cm'
  write(*,'(A,1E21.14,A)')'   me         = ',me,' g'
  write(*,'(A,1E21.14,A)')'   omega_p    = ',omega_pe(n0),' /s'
  write(*,'(A,1E21.14,A)')'   e          = ',e,' statC'
  write(*,'(A,1E21.14,A)')'   vTe0       = ',v_T(Te),' cm/s'
  write(*,'(A,1E21.14,A)')'   n0         = ',n0,' /cm3'
  write(*,*)'------------------------------------------'
  write(*,*)'* 1D-1V phase-space : '
  write(*,'(A,1E21.14)')'   x_min      = ',x_min
  write(*,'(A,1E21.14)')'   x_max      = ',x_max
  write(*,'(A,1E21.14)')'   d_x        = ',d_x
  write(*,'(A,1E21.14)')'   vx_min     = ',vx_min
  write(*,'(A,1E21.14)')'   vx_max     = ',vx_max
  write(*,'(A,1E21.14)')'   d_vx       = ',d_vx
  write(*,*)'------------------------------------------'
  write(*,*)'* Simulation properties : '
  write(*,'(A,1E21.14)')'   cfl        = ',cfl
  write(*,'(A,1E21.14)')'   L_t        = ',L_t
  write(*,'(A,1E21.14)')'   dt_diag    = ',dt_diag
  write(*,'(A,1I3)')'   maxwell    = ',maxwell
  write(*,'(A,1I3)')'   b_cond     = ',b_cond
  write(*,'(A,1I3)')'   scheme     = ',scheme
  write(*,'(A,1E21.14)')'   b_VL       = ',b_VL
  write(*,*)'------------------------------------------'
  write(*,*)'* Perturbation properties : '
  write(*,'(A,1I3)')'   perturb    = ',perturb
  if (perturb .ne. 0) then
    write(*,'(A,1E21.14)')'   A_pert     = ',A_pert
    write(*,'(A,1E21.14)')'   k_pert     = ',k_pert
    write(*,'(A,1E21.14)')'   omega_pert = ',omega_pert
    write(*,'(A,1E21.14)')'   v_d        = ',v_d
  else 
    write(*,*)'  no perturbation'
  end if
  write(*,*)'------------------------------------------'
  write(*,*)'* Deduced bins number : '
  write(*,'(A,1I7)')'   N_x        = ',N_x
  write(*,'(A,1I7)')'   N_vx       = ',N_vx
  write(*,*)'------------------------------------------' 
  write(*,*) ' '
end subroutine read_init_parameters

subroutine get_str(str)
  character(len=*), intent(inout) :: str
  read (1,'(A)',end=10) str
  str = adjustl(str)
  return
  continue ; 10 str = 'end'
end subroutine get_str

function get_logical(str)
  character(len=*), intent(in) :: str
  logical :: get_logical
  read (str, *) get_logical
  return
end function get_logical

function get_integer(str)
  character(len=*), intent(in) :: str
  integer :: get_integer
  read (str, *) get_integer
  return
end function get_integer

function get_real(str)
  character(len=*), intent(in) :: str
  real(PR) :: get_real
  read (str, *) get_real
  return
end function get_real

function get_char(str)
  character(len=*), intent(in) :: str
  Character(len=60)            :: get_char
  read (str, *) get_char
  get_char = adjustl(get_char)
  return
end function get_char

end module input
