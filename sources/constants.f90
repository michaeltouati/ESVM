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
module constants

  use acuracy
  
  implicit none
  
  !=======================================================================================
  !                              Mathematical constants
  !=======================================================================================

  real(PR), parameter, public :: pi   = 3.141592653589793_PR
  real(PR), parameter, public :: zero = 1.e-15_PR

  !=======================================================================================
  !                                Physical constants
  !=======================================================================================
  
  real(PR), parameter, public :: mu    = 1.66053907E-24_PR     ! Atomic mass unit [g]
  real(PR), parameter, public :: mp    = 1.67262192E-24_PR     ! Proton mass unit [g]
  real(PR), parameter, public :: me    = 9.10938370E-28_PR     ! Electron mass [g]
  real(PR), parameter, public :: e     = 4.803204E-10_PR       ! Elementary charge [statC]
  real(PR), parameter, public :: c     = 2.99792458E10_PR      ! Speed of light [cm/s]  
  real(PR), parameter, public :: kb    = 1.380649E-16_PR       ! Boltzmann constant [erg/K]
  real(PR), parameter, public :: h     = 6.62607015E-27_PR     ! Planck constant [erg.s]
  real(PR), parameter, public :: hbar  = h / (2._PR*pi)
  real(PR), parameter, public :: alpha = 1._PR / 137.035999_PR ! Fine Structure constant []

  !=======================================================================================
  !                                Conversion factors
  !=======================================================================================
  
  real(PR), parameter, public :: Joules  = 1.E7_PR             ! 1 J in erg
  real(PR), parameter, public :: qe      = 1.602176634E-19_PR  ! Elementary charge [C]
  real(PR), parameter, public :: eV      = qe * Joules / kB    ! 1 eV = 11604 K

  !=======================================================================================
  !                                Advection schemes
  !=======================================================================================

  integer , parameter, public :: L_donor_cell   = 1
  integer , parameter, public :: L_Lax_Wendroff = 2
  integer , parameter, public :: L_Beam_Warming = 3
  integer , parameter, public :: L_Fromm        = 4
  integer , parameter, public :: NL_minmod      = 5
  integer , parameter, public :: NL_superbee    = 6
  integer , parameter, public :: NL_Van_Leer    = 7
  integer , parameter, public :: NL_MUSCL1      = 8 
  integer , parameter, public :: NL_MUSCL2      = 9

end module constants
