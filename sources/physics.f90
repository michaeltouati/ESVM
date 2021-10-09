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
module physics

use acuracy
use constants
  
implicit none
  
public :: omega_pe, v_T, Debye
  
contains

function omega_pe(n0)
  implicit none
  real(PR), intent(in) :: n0       ! (/cm^3)
  real(PR)             :: omega_pe ! (/s)
  omega_pe = sqrt(4 * pi * n0 * (e**2._PR) / me)
end function omega_pe
  
function v_T(Te)
  implicit none
  real(PR), intent(in) :: Te  ! Te (K)
  real(PR)             :: v_T ! v_T (cm/s)
  v_T = (3._PR* kb * Te / me )**0.5_PR
end function v_T
  
function Debye(n0, Te)
  implicit none
  real(PR), intent(in) :: n0, Te ! ne (/cm^3) and T (K)
  real(PR)             :: Debye  ! Debye (cm)
  Debye = ( kB * Te / (4 * pi * n0 * (e**2._PR)) )**(0.5_PR) 
end function Debye

end module physics
