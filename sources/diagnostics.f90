!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!          1D-1V ElectroStatic Vlasov-Maxwell (ESV) code            !!
!!                                                                   !!
!!  Written by Dr Michaël J TOUATI - CLPU - 2020 - mtouati@clpu.es   !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module diagnostics

use acuracy
use constants
use input
use physics

implicit none
public  :: INIT_DIAG

contains

! Subroutines

subroutine INIT_DIAG()
  implicit none
  open (unit=1 ,file='results/fe.dat' ,form='formatted',status='unknown')
  open (unit=2 ,file='results/ne.dat' ,form='formatted',status='unknown')
  open (unit=3 ,file='results/Ex.dat' ,form='formatted',status='unknown')
  open (unit=4 ,file='results/je.dat' ,form='formatted',status='unknown')
  open (unit=5 ,file='results/ve.dat' ,form='formatted',status='unknown')
  open (unit=7 ,file='results/vTe.dat',form='formatted',status='unknown')
  open (unit=11,file='results/phi.dat',form='formatted',status='unknown')
  open (unit=8 ,file='results/UK.dat' ,form='formatted',status='unknown')
  open (unit=9 ,file='results/UT.dat' ,form='formatted',status='unknown')
  open (unit=10,file='results/UE.dat' ,form='formatted',status='unknown')
end subroutine INIT_DIAG

subroutine DIAG_ENERGY(time, N_x, d_x, n_e, v_e, vT_e, E_x_n, &
                     & dU_K, dU_T, dU_E, U_K, U_T, U_E)
  implicit none
  integer                      , intent(in)    :: N_x   
  real(PR)                     , intent(in)    :: d_x, time
  real(PR), dimension(-1:N_x+2), intent(in)    :: n_e, v_e, vT_e
  real(PR), dimension(-1:N_x+2), intent(in)    :: E_x_n
  real(PR), dimension(1:N_x)   , intent(inout) :: dU_K, dU_T, dU_E
  real(PR)                     , intent(inout) :: U_K, U_T, U_E
  integer                                      :: i
  !
  !$omp PARALLEL DO DEFAULT(SHARED) PRIVATE(i) COLLAPSE(1)
  do i = 1,N_x,1
    dU_K(i) = n_e(i) * (v_e(i)**2._PR)
    dU_T(i) = n_e(i) * (vT_e(i)**2._PR)
    dU_E(i) = (E_x_n(i)**2._PR)
  end do
  !$omp END PARALLEL DO
  !
  do i = 1,N_x,1
    U_K = U_K + (dU_K(i) * d_x)
    U_T = U_T + (dU_T(i) * d_x)
    U_E = U_E + (dU_E(i) * d_x)
  end do
  !
  if (abs(U_K).lt.zero) then
    write (8,'(4E22.14)') time, 0.
  else
    write(8,'(4E22.14)') time, U_K
  end if

  if (abs(U_T).lt.zero) then
    write (9,'(4E22.14)') time, 0.
  else
    write(9,'(4E22.14)') time, U_T
  end if

  if (abs(U_E).lt.zero) then
    write (10,'(4E22.14)') time, 0.
  else
    write(10,'(4E22.14)') time, U_E
  end if
end subroutine DIAG_ENERGY

subroutine DIAG(N_t, time, N_x, x, N_vx, vx, test_positivity, U_K, U_T, U_E, &
              & f_n, n_e, E_x_n, j_e, v_e, vT_e, phi_n)
  implicit none
  integer                                , intent(in) :: N_t, N_x, N_vx
  real(PR), dimension(-1:N_x+2)          , intent(in) :: x
  real(PR), dimension(-1:N_vx+2)         , intent(in) :: vx
  real(PR)                               , intent(in) :: time
  logical                                , intent(in) :: test_positivity
  real(PR), dimension(-1:N_x+2,-1:N_vx+2), intent(in) :: f_n
  real(PR), dimension(-1:N_x+2)          , intent(in) :: n_e, j_e, v_e, vT_e
  real(PR), dimension(-1:N_x+2)          , intent(in) :: E_x_n , phi_n
  real(PR)                               , intent(in) :: U_K, U_T, U_E
  real(PR)                                            :: U_tot
  integer                                             :: l, i
  !
  write (*,*)'==========================='
  write (*,'(A,1E22.14)')'time (/omega_p) =', time 
  write (*,'(A,1I7)')'Number of iteration :',N_t
  write (*,*)'==========================='
  write (*,*) ' '
  if (test_positivity.eqv..true.) print*, 'the distribution function became negative'
  if (abs(U_K).lt.zero) then
    write (*,'(A,1E22.2)')'Kinetic energy  (n0 Debye^3 me v_T^2 / 2) = ', U_K
  else
    write (*,'(A,1E22.2)')'Kinetic energy  (n0 Debye^3 me v_T^2 / 2) = ', 0.
  end if
  if (abs(U_T).lt.zero) then
    write (*,'(A,1E22.2)')'Thermal energy  (n0 Debye^3 me v_T^2 / 2) = ', U_T
  else
    write (*,'(A,1E22.2)')'Thermal energy  (n0 Debye^3 me v_T^2 / 2) = ', 0.
  end if
  if (abs(U_E).lt.zero) then
    write (*,'(A,1E22.2)')'Electric energy (n0 Debye^3 me v_T^2 / 2) = ', U_E
  else
    write (*,'(A,1E22.2)')'Electric energy (n0 Debye^3 me v_T^2 / 2) = ', 0.
  end if
  write (*,*)'--------------------------------------------'
  U_tot = U_K + U_T + U_E 
  if (abs(U_tot).lt.zero) then
    write (*,'(A,1E22.2)')'Total energy    (n0 Debye^3 me v_T^2 / 2) = ', U_tot
  else
    write (*,'(A,1E22.2)')'Total energy    (n0 Debye^3 me v_T^2 / 2) = ', 0.
  end if
  write (*,*) ' '
  !
  do l=1,N_vx,1
    do i = 1,N_x,1
      if (abs(f_n(i,l)).lt.zero) then
        write (1,'(4E22.14)') time, vx(l), x(i), 0.
      else
        write(1,'(4E22.14)') time, vx(l), x(i), f_n(i,l)
      end if
    end do
  end do
  do i = 1,N_x,1
    if (abs(n_e(i)).lt.zero) then
      write (2,'(4E22.14)') time, x(i), 0.
    else
      write(2,'(4E22.14)') time, x(i), n_e(i)
    end if
  end do
  do i = 1,N_x,1
    if (abs(E_x_n(i)).lt.zero) then
      write (3,'(4E22.14)') time, x(i), 0.
    else
      write(3,'(4E22.14)') time, x(i), E_x_n(i)
    end if
  end do
  do i = 1,N_x,1
    if (abs(j_e(i)).lt.zero) then
      write (4,'(4E22.14)') time, x(i), 0.
    else
      write(4,'(4E22.14)') time, x(i), j_e(i)
    end if
  end do
  do i = 1,N_x,1
    if (abs(v_e(i)).lt.zero) then
      write (5,'(4E22.14)') time, x(i), 0.
    else
      write(5,'(4E22.14)') time, x(i), v_e(i)
    end if
  end do
  do i = 1,N_x,1
    if (abs(vT_e(i)).lt.zero) then
      write (7,'(4E22.14)') time, x(i), 0.
    else
      write(7,'(4E22.14)') time, x(i), vT_e(i)
    end if
  end do
  do i = 1,N_x,1
    !
    if (abs(phi_n(i)).lt.zero) then
      write (11,'(4E22.14)') time, x(i), 0.
    else
      write(11,'(4E22.14)') time, x(i), phi_n(i)
    end if
  end do
end subroutine DIAG

subroutine CLOSE_DIAG()
  implicit none
  close(1)
  close(2)
  close(3)
  close(4)
  close(5)
  close(7)
  close(8)
  close(9)
  close(10)
  close(11)
end subroutine CLOSE_DIAG

end module diagnostics