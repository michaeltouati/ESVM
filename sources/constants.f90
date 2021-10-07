!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!          1D-1V ElectroStatic Vlasov-Maxwell (ESV) code            !!
!!                                                                   !!
!!  Written by Dr Michaël J TOUATI - CLPU - 2020 - mtouati@clpu.es   !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module constants

  use acuracy
  
  implicit none
  
  !=======================================================================================
  !                                         Mathematical constants
  !=======================================================================================

  real(PR),parameter,public :: pi = 3.14159_PR

  !=======================================================================================
  !                                         Physical_constants constants
  !=======================================================================================
  
  real(PR),parameter,public :: mu   = 1.6605E-24_PR           ! Atomic mass unit [g]
  real(PR),parameter,public :: mp   = 1.6726E-24_PR           ! Proton mass unit [g]
  real(PR),parameter,public :: me   = 9.1094E-28_PR           ! Electron mass [g] 
  real(PR),parameter,public :: e    = 4.8032E-10_PR           ! Elementary charge [statcoulomb]
  real(PR),parameter,public :: c    = 2.9979E10_PR            ! Speed of light in vacuum [cm/s]  
  real(PR),parameter,public :: kb   = 1.3807E-16_PR           ! Boltzmann constant [erg/K]
  real(PR),parameter,public :: h    = 6.6261E-27_PR           ! Placnk constant [erg.s]
  real(PR),parameter,public :: hbar = h / (2*pi)
  real(PR),parameter,public :: alpha = 1._PR / 137.035999_PR  ! Fine Structure constant []

  !=======================================================================================
  !                                        Conversion factors
  !=======================================================================================

  real(PR),parameter,public   :: eV = 1.6022E-19_PR / 1.3807E-23_PR       ! 1 eV = 11604 K
  real(PR), parameter, public :: fs = 1.E-15_PR
  real(PR), parameter, public :: microns = 1.E-4_PR
  real(PR), parameter, public :: Joules = 1.E7_PR                         ! 1 J in erg

  !=======================================================================================
  !                                      Advection schemes
  !=======================================================================================

  real(PR), parameter, public           :: zero = 1.e-14_PR
  integer, parameter, public            :: L_donor_cell   = 1
  integer, parameter, public            :: L_Lax_Wendroff = 2
  integer, parameter, public            :: L_Beam_Warming = 3
  integer, parameter, public            :: L_Fromm        = 4
  integer, parameter, public            :: NL_minmod      = 5
  integer, parameter, public            :: NL_superbee    = 6
  integer, parameter, public            :: NL_Van_Leer    = 7
  integer, parameter, public            :: NL_MUSCL1      = 8 
  integer, parameter, public            :: NL_MUSCL2      = 9

end module constants
