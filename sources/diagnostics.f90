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
module diagnostics

use acuracy
use constants
use input
use physics

implicit none

public  :: INIT_DIAG, DIAG_ENERGY
public  :: DIAG, CLOSE_DIAG

private :: DUMP_FIELD

integer, parameter, private :: Ndiag_timing = 1

integer, parameter, private :: Ndiag_fe  = 20
integer, parameter, private :: Ndiag_ne  = 21
integer, parameter, private :: Ndiag_je  = 22
integer, parameter, private :: Ndiag_ve  = 23
integer, parameter, private :: Ndiag_vTe = 24
integer, parameter, private :: Ndiag_Phi = 25
integer, parameter, private :: Ndiag_Ex  = 26

integer, parameter, private :: Ndiag_UK = 40
integer, parameter, private :: Ndiag_UT = 41
integer, parameter, private :: Ndiag_UE = 42

contains

! Subroutines

subroutine INIT_DIAG()
  implicit none
  open (unit=Ndiag_timing, file='results/'//trim(simu)//'/timing.dat', form='formatted', status='unknown')
  open (unit=Ndiag_fe    , file='results/'//trim(simu)//'/fe.dat'    , form='formatted', status='unknown')
  open (unit=Ndiag_ne    , file='results/'//trim(simu)//'/ne.dat'    , form='formatted', status='unknown')
  open (unit=Ndiag_je    , file='results/'//trim(simu)//'/je.dat'    , form='formatted', status='unknown')
  open (unit=Ndiag_ve    , file='results/'//trim(simu)//'/ve.dat'    , form='formatted', status='unknown')
  open (unit=Ndiag_vTe   , file='results/'//trim(simu)//'/vTe.dat'   , form='formatted', status='unknown')
  open (unit=Ndiag_Phi   , file='results/'//trim(simu)//'/Phi.dat'   , form='formatted', status='unknown')
  open (unit=Ndiag_Ex    , file='results/'//trim(simu)//'/Ex.dat'    , form='formatted', status='unknown')
  open (unit=Ndiag_UK    , file='results/'//trim(simu)//'/UK.dat'    , form='formatted', status='unknown')
  open (unit=Ndiag_UT    , file='results/'//trim(simu)//'/UT.dat'    , form='formatted', status='unknown')
  open (unit=Ndiag_UE    , file='results/'//trim(simu)//'/UE.dat'    , form='formatted', status='unknown')
end subroutine INIT_DIAG

subroutine DIAG_ENERGY(t0, UK, UT, UE)
  implicit none
  real(PR)                    , intent(in) :: t0
  real(PR)                    , intent(in) :: UK, UT, UE
  integer                                  :: k
  !
  !$omp PARALLEL DO DEFAULT(SHARED) PRIVATE(k) COLLAPSE(1)
  do k=Ndiag_UK,Ndiag_UE,1
    select case (k)
      case(Ndiag_UK)
        write(Ndiag_UK,'(2E23.15)') t0, UK
      case(Ndiag_UT)
        write(Ndiag_UT,'(2E23.15)') t0, UT
      case(Ndiag_UE)
        write(Ndiag_UE,'(2E23.15)') t0, UE
    end select
  end do
  !$omp END PARALLEL DO
end subroutine DIAG_ENERGY

subroutine DUMP_FIELD(Nval, t0, x0, val)
  implicit none
  integer,  intent(in) :: Nval
  real(PR), intent(in) :: t0, x0, val
  !
  if (abs(val).lt.zero) then
    write(Nval,'(3E23.15)') t0, x0, zero
  else
    write(Nval,'(3E23.15)') t0, x0, val
  end if
end subroutine DUMP_FIELD

subroutine DIAG(Nt, t0, Nx, x0, Nvx, vx0, &
              & positivity, UK, UT, UE, &
              & fn, ne, Exn, je, ve, vTe, phin)
  implicit none
  integer                               , intent(in) :: Nt, Nx, Nvx
  real(PR), dimension(-1:Nx+2)          , intent(in) :: x0
  real(PR), dimension(-1:Nvx+2)         , intent(in) :: vx0
  real(PR)                              , intent(in) :: t0
  logical                               , intent(in) :: positivity
  real(PR), dimension(-1:Nx+2,-1:N_vx+2), intent(in) :: fn
  real(PR), dimension(-1:Nx+2)          , intent(in) :: ne, je, ve, vTe
  real(PR), dimension(-1:Nx+2)          , intent(in) :: Exn, phin
  real(PR)                              , intent(in) :: UK, UT, UE
  real(PR)                                           :: Utot
  integer                                            :: l, i, k
  !
  write(*,*)'================================'
  write(*,'(A,1E14.7)')' time (/omega_p)  =',t0 
  write(*,'(A,1I10)')' Number of iterations :',Nt
  write(*,*)'================================'
  write(*,*)' '
  if (positivity.eqv..true.) write(*,*)'the distribution function became negative'
  write(*,'(A,1E14.7)')' Kinetic energy  (n0 Debye^3 me vTe0^2) =', UK
  write(*,'(A,1E14.7)')' Thermal energy  (n0 Debye^3 me vTe0^2) =', UT
  write(*,'(A,1E14.7)')' Electric energy (n0 Debye^3 me vTe0^2) =', UE
  write(*,*)'------------------------------------------------------'
  Utot = UK + UT + UE 
  write(*,'(A,1E14.7)')' Total energy    (n0 Debye^3 me vTe0^2) =', Utot
  write(*,*)' '
  !
  !$omp PARALLEL DO DEFAULT(SHARED) PRIVATE(k,i,l) COLLAPSE(1)
  do k=Ndiag_fe,Ndiag_Ex,1
    select case (k)
      case(Ndiag_fe)
        do l = 1,Nvx,1
          do i = 1,Nx,1
            if (abs(fn(i,l)).lt.zero) then
              write(Ndiag_fe,'(4E23.15)') t0, vx0(l), x0(i), zero
            else
              write(Ndiag_fe,'(4E23.15)') t0, vx0(l), x0(i), fn(i,l)
            end if
          end do
        end do
      case(Ndiag_ne)
        do i = 1,Nx,1
          call DUMP_FIELD(Ndiag_ne , t0, x0(i), ne(i)  )
        end do
      case(Ndiag_je)
        do i = 1,Nx,1
          call DUMP_FIELD(Ndiag_je , t0, x0(i), je(i)  )
        end do
      case(Ndiag_ve)
        do i = 1,Nx,1
          call DUMP_FIELD(Ndiag_ve , t0, x0(i), ve(i)  )
        end do
      case(Ndiag_vTe)
        do i = 1,Nx,1
          call DUMP_FIELD(Ndiag_vTe, t0, x0(i), vTe(i) )
        end do
      case(Ndiag_Phi)
        do i = 1,Nx,1
          call DUMP_FIELD(Ndiag_Phi, t0, x0(i), phin(i))
        end do        
      case(Ndiag_Ex)
        do i = 1,Nx,1
          call DUMP_FIELD(Ndiag_Ex , t0, x0(i), Exn(i) )
        end do
    end select
  end do
  !$omp END PARALLEL DO
end subroutine DIAG

subroutine CLOSE_DIAG(CPUtime0, clock0)
  implicit none
  real(PR), intent(in) :: CPUtime0,clock0
  !
  write(Ndiag_timing,'(A,1E14.7,A)')" Simulation CPUxtime        = ",CPUtime0," hours"
  write(Ndiag_timing,'(A,1E14.7,A)')" Simulation ellapsed time   = ",clock0," hours"
  !
  close(Ndiag_timing)
  close(Ndiag_fe)
  close(Ndiag_ne)
  close(Ndiag_je)
  close(Ndiag_ve)
  close(Ndiag_vTe)
  close(Ndiag_Phi)
  close(Ndiag_Ex)
  close(Ndiag_UK)
  close(Ndiag_UT)
  close(Ndiag_UE)
  !
  write (*,*)'================================================='
  write (*,'(A,1E14.7,A)')" Simulation CPUxtime        = ",CPUtime0," hours"
  write (*,'(A,1E14.7,A)')" Simulation ellapsed time   = ",clock0," hours"
end subroutine CLOSE_DIAG

end module diagnostics
