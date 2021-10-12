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
program ESVM
!
use acuracy
use constants
use input
use physics
use diagnostics
use library
use omp_lib
!
implicit none
real(PR)                              :: time, d_t
real(PR), dimension(:), allocatable   :: x, vx
real(PR), dimension(:,:), allocatable :: f_n, f_np1
real(PR), dimension(:), allocatable   :: n_e, j_e, v_e, vT_e
real(PR), dimension(:), allocatable   :: E_x_n, E_x_np1, phi_n
integer                               :: N_t, i, l 
real(PR)                              :: f_max, flux_im1, flux_ip1
real(PR), dimension(:), allocatable   :: dU_K, dU_T, dU_E 
real(PR)                              :: U_K, U_T, U_E 
logical                               :: test_positivity, save_results
real(PR)                              :: timer_start, timer_finish
!
call cpu_time(timer_start)
!
call read_init_parameters()
!
call system('mkdir -p results')
call system('mkdir -p results/'//trim(simu))
!
allocate(x(-1:N_x+2),vx(-1:N_vx+2))
allocate(f_n(-1:N_x+2,-1:N_vx+2),f_np1(-1:N_x+2,-1:N_vx+2))
allocate(n_e(-1:N_x+2),j_e(-1:N_x+2),v_e(-1:N_x+2),vT_e(-1:N_x+2))
allocate(phi_n(-1:N_x+2),E_x_n(-1:N_x+2),E_x_np1(-1:N_x+2))
allocate(dU_K(1:N_x),dU_T(1:N_x),dU_E(1:N_x))
!
call GRID(N_x,N_vx,d_x,d_vx,x_min,vx_min,x(-1:N_x+2),vx(-1:N_vx+2))
!
call INIT_VAR(f_n, f_np1, n_e, j_e, v_e, vT_e, E_x_np1, E_x_n , phi_n, &
            & dU_K, dU_T, dU_E, U_K, U_T, U_E, time, N_t, & 
            & test_positivity, save_results)
!
call INIT_SIMU(x(1:N_x), vx(1:N_vx), f_n)
!
call INIT_DIAG()
!
do while (time.lt.L_t)
  !
  f_max = maxval(maxval(f_n(:,:),dim=2))
  !
  call DENSITIES(vx(1:N_vx), f_n(-1:N_x+2,-1:N_vx+2),&
               & n_e(-1:N_x+2), j_e(-1:N_x+2),&
               & v_e(-1:N_x+2), vT_e(-1:N_x+2))
  ! Landau damping test-case
  if (perturb == 2) then
    if ( time < (6.* pi / omega_0) ) then
      call DRIVE(d_t, time, x, E_x_n, E_x_np1, phi_n)
    else 
      call MAXWELL_SOLVER(maxwell, N_t, d_t, d_x, j_e, n_e, E_x_n, E_x_np1, phi_n)
    end if
  ! All other cases
  else
    call MAXWELL_SOLVER(maxwell, N_t, d_t, d_x, j_e, n_e, E_x_n, E_x_np1, phi_n)
  end if
  !   
  d_t   = cfl*(0.5_PR/((maxval(vx(1:N_vx))/d_x)+(maxval(abs(E_x_n(1:N_x)))/d_vx)))
  !
  !$omp PARALLEL DO DEFAULT(SHARED) PRIVATE(l,i,flux_im1,flux_ip1) COLLAPSE(2)
  do l=1,N_vx,1
    do i = 1,N_x,1
      call fluxes(scheme, vx(l), f_max, d_t, d_x, &
                & f_n(i-2,l), f_n(i-1,l), f_n(i,l), f_n(i+1,l), f_n(i+2,l), &
                & flux_im1, flux_ip1)
      f_np1(i,l) = f_n(i,l) - (d_t * (flux_ip1 - flux_im1) / d_x)  
      !
      call fluxes(scheme,-E_x_n(i), f_max, d_t, d_vx,&
                & f_n(i,l-2), f_n(i,l-1), f_n(i,l), f_n(i,l+1), f_n(i,l+2), &
                & flux_im1, flux_ip1)
      f_np1(i,l) = f_np1(i,l) - (d_t * (flux_ip1 - flux_im1) / d_vx)
      if ((f_np1(i,l).lt.0._PR).and.(test_positivity.eqv..false.)) test_positivity = .true.
    end do
  end do
  !$omp END PARALLEL DO
  !
  call BOUNDARIES(f_np1(-1:N_x+2,-1:N_vx+2))
  !
  call DIAG_ENERGY(time, N_x, d_x, n_e, v_e, vT_e, E_x_n, &
                 & dU_K, dU_T, dU_E, U_K, U_T, U_E)
  !
  save_results = (mod(time,dt_diag).lt.d_t).and.(time.ge.d_t)
  save_results = save_results.or.(N_t.eq.1)
  save_results = save_results.or.((L_t-time).le.d_t)
  !
  if (save_results.eqv..true.) then
    call DIAG(N_t, time, N_x, x, N_vx, vx, test_positivity, U_K, U_T, U_E, &
           & f_n, n_e, E_x_n, j_e, v_e, vT_e, phi_n)
  end if
 !
 call INIT_NEXT_STEP(f_n, f_np1, n_e, j_e, v_e, vT_e, &
                   & phi_n, E_x_n, E_x_np1, N_t, d_t, time, U_K, U_T, U_E)
end do
!
deallocate(x,vx,f_n,f_np1,n_e,j_e,v_e,vT_e,E_x_n,E_x_np1)
!
call CLOSE_DIAG()
!
call cpu_time(timer_finish)
!
write (*,*)'=============================================='
write (*,'(A,1E11.3,A)')" Total simulation time duration = ",timer_finish-timer_start," s"
!
end program ESVM
