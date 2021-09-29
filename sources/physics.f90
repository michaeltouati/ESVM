!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!          1D-1V ElectroStatic Vlasov-Maxwell (ESV) code            !!
!!                                                                   !!
!!  Written by Dr Michaël J TOUATI - CLPU - 2020 - mtouati@clpu.es   !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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