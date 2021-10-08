!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!             ElectroStatic Vlasov-Maxwell (ESVM) code              !!
!!                                                                   !!
!!                  Written by Dr Michaël J TOUATI                   !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module library

use acuracy
use constants
use input
use physics
  
implicit none
public  :: GRID, INIT_VAR, INIT_SIMU
public  :: DENSITIES, POISSON, AMPERE
public  :: FLUXES, BOUNDARIES
private :: slope, minmod, minmod_3, maxmod, theta

contains

! Subroutines
  
subroutine GRID(Nx,Nv,dx,dv,xmin,vmin,x,vx)
  implicit none
  integer, intent(in)                         :: Nx, Nv
  real(PR), intent(in)                        :: dx, dv, xmin, vmin
  real(PR), dimension(-1:Nx+2), intent(out)   :: x
  real(PR), dimension(-1:Nv+2), intent(out)   :: vx
  integer                                     :: i, l
  do i=1,Nx+2,1
    x(i)  = xmin + real(i-1,PR)*dx
  end do
  x(0)  = x(1) - dx
  x(-1) = x(0) - dx 
  do l=1,Nv+2,1
    vx(l) = vmin + real(l-1,PR)*dv
  end do
    vx(0)  = vx(1) - dv
    vx(-1) = vx(0) - dx
  end subroutine GRID

subroutine INIT_VAR(f_n, f_np1, n_e, j_e, v_e, vT_e, E_x_np1, E_x_n , phi_n, &
            & dU_K, dU_T, dU_E, U_K, U_T, U_E, time, N_t, & 
            & test_positivity, save_results)
  implicit none
  real(PR), dimension(-1:N_x+2,-1:N_vx+2), intent(out)   :: f_n, f_np1
  real(PR), dimension(-1:N_x+2)          , intent(out)   :: n_e, j_e, v_e, vT_e
  real(PR), dimension(-1:N_x+2)          , intent(out)   :: E_x_np1, E_x_n , phi_n
  real(PR), dimension(1:N_x)             , intent(out)   :: dU_K, dU_T, dU_E
  real(PR)                               , intent(out)   :: U_K, U_T, U_E, time
  integer                                , intent(out)   :: N_t
  logical                                , intent(out)   :: test_positivity, save_results
  !
  f_n     = 0._PR
  f_np1   = 0._PR
  n_e     = 0._PR
  j_e     = 0._PR
  v_e     = 0._PR
  vT_e    = 0._PR
  E_x_np1 = 0._PR
  E_x_n   = 0._PR
  phi_n   = 0._PR
  !
  dU_K = 0._PR
  dU_T = 0._PR
  dU_E = 0._PR
  U_K  = 0._PR
  U_T  = 0._PR
  U_E  = 0._PR
  !
  time    = 0._PR
  N_t     = 1
  !
  test_positivity = .false.
  save_results    = .false.
end subroutine INIT_VAR

subroutine INIT_SIMU(x, vx, f_n)
  implicit none
  real(PR), dimension(1:N_x), intent(in)                 :: x
  real(PR), dimension(1:N_vx), intent(in)                :: vx
  real(PR), dimension(-1:N_x+2,-1:N_vx+2), intent(inout) :: f_n
  integer                                                :: l, i
  real(PR) :: xs, dx, dvx, X2
  if (perturb == 1) then
  ! Electrostatic wakefield test case :
    dx  = 0.25_PR  ! "particle size"
    dvx = 0.025_PR ! "particle size"
    xs  = x_min + ( (x_max - x_min) / 8. ) 
    do l=1,N_vx,1
      do i=1,N_x,1
        f_n(i,l) = (1.0_PR/sqrt(2._PR*pi))*&
        & exp(-(vx(l)**2._PR)/2._PR) 
        X2 = -0.5_PR * ( ( ( (x(i)-xs) / dx )**2._PR) + ( ( (vx(l)-vd) / dvx )**2._PR) )
        f_n(i,l) = f_n(i,l) + ( A * exp(X2)/ (2._PR * pi * dx * dvx ) )
      end do
    end do
  ! Landau damping test case :
  else if (perturb == 2) then
    do l=1,N_vx,1
      do i=1,N_x,1
        f_n(i,l) = (1._PR/sqrt(2._PR*pi))*&
                 & exp(-(vx(l)**2._PR)/2._PR) 
       end do
    end do
  else if (perturb == 3) then
  ! Two-stream instability test case :
    do l=1,N_vx,1
      do i=1,N_x,1
        f_n(i,l) = (0.5_PR/sqrt(2._PR*pi))*&
        & ((1._PR+A*sin(k*x(i)))*exp(-((vx(l)-vd)**2._PR)/2._PR) + &
        &  (1._PR-A*sin(k*x(i)))*exp(-((vx(l)+vd)**2._PR)/2._PR))
       end do
    end do
  else
  ! Maxwellian :
    do l=1,N_vx,1
      do i=1,N_x,1
        f_n(i,l) = (1._PR/sqrt(2._PR*pi))*&
                 & exp(-((vx(l)-vd)**2._PR)/2._PR) 
       end do
    end do
  end if
  ! Boundary conditions
  select case (b_cond)
    ! absorbing
    case (1)
      do l=-1,N_vx+2,1
        f_n(N_x+1,l) = 0._PR
        f_n(N_x+2,l) = 0._PR
        f_n(0,l)     = 0._PR
        f_n(-1,l)    = 0._PR
      end do
    ! periodic
    case (2)
      do l=-1,N_vx+2,1
        f_n(N_x+1,l) = f_n(1,l)
        f_n(N_x+2,l) = f_n(2,l)
        f_n(0,l)     = f_n(N_x,l)
        f_n(-1,l)    = f_n(N_x-1,l)
      end do
  end select
  do i=-1,N_x+2,1
    f_n(i,N_vx+1) = 0._PR
    f_n(i,N_vx+2) = 0._PR
    f_n(i,0)      = 0._PR
    f_n(i,-1)     = 0._PR
  end do
end subroutine INIT_SIMU

subroutine INIT_NEXT_STEP(f_n, f_np1, n_e, j_e, v_e, vT_e, &
                  & phi_n, E_x_n, E_x_np1, N_t, d_t, time, U_K, U_T, U_E)
  implicit none
  real(PR), dimension(-1:N_x+2,-1:N_vx+2), intent(inout) :: f_n, f_np1
  real(PR), dimension(-1:N_x+2)          , intent(inout) :: n_e, j_e, v_e, vT_e
  real(PR), dimension(-1:N_x+2)          , intent(inout) :: E_x_np1, E_x_n , phi_n
  real(PR)                               , intent(inout) :: U_K, U_T, U_E, d_t, time
  integer                                , intent(inout) :: N_t
  !
  f_n     = f_np1
  f_np1   = 0._PR
  n_e     = 0._PR
  j_e     = 0._PR
  v_e     = 0._PR
  vT_e     = 0._PR
  phi_n   = 0._PR
  E_x_n   = E_x_np1 
  E_x_np1 = 0._PR 
  N_t     = N_t + 1
  time    = time + d_t
  U_K = 0._PR
  U_T = 0._PR
  U_E = 0._PR
end subroutine INIT_NEXT_STEP

subroutine DENSITIES(vx, f_n, n_e, j_e, v_e, vT_e)
  implicit none
  real(PR), dimension(1:N_vx), intent(in)               :: vx
  real(PR), dimension(-1:N_x+2,-1:N_vx+2), intent(in)   :: f_n
  real(PR), dimension(-1:N_x+2), intent(out)            :: n_e, j_e, v_e, vT_e
  integer                                               :: i, l
  real(PR), dimension(:), allocatable                   :: F1,F2,F3
  !$omp PARALLEL DO DEFAULT(SHARED) PRIVATE(i,l,F1,F2,F3) COLLAPSE(1)
  do i=-1,N_x+2,1
    allocate(F1(1:N_vx),F2(1:N_vx),F3(1:N_vx))
    F1(1:N_vx) = f_n(i,1:N_vx)*d_vx
    do l=1,N_vx,1
      F2(l) = - f_n(i,l) * vx(l) * d_vx
      F3(l) = f_n(i,l) * (vx(l)**2._PR) * d_vx
    end do
    n_e(i)  = sum(F1(1:N_vx))
    j_e(i)  = sum(F2(1:N_vx))
    if (abs(n_e(i)).lt.zero) then
      v_e(i) = 0.
    else
      v_e(i)  = - j_e(i) / n_e(i)
    end if
    if (abs(n_e(i)).ne.(v_e(i)**2._PR)) then
      vT_e(i) = ((sum(F3(1:N_vx))/n_e(i))-(v_e(i)**2._PR))**0.5_PR
    else
      vT_e(i) = 0.
    end if
   deallocate(F1,F2,F3)
  end do
  !$omp END PARALLEL DO
end subroutine DENSITIES

subroutine MAXWELL_SOLVER(solver, N_t, d_t, d_x, j_e, n_e, E_x_n, E_x_np1, phi_n)
  implicit none
  integer, intent(in)                          :: solver, N_t
  real(PR), intent(in)                         :: d_t, d_x
  real(PR), dimension(-1:N_x+2), intent(in)    :: j_e, n_e
  real(PR), dimension(-1:N_x+2), intent(inout) :: E_x_n
  real(PR), dimension(-1:N_x+2), intent(out)   :: E_x_np1, phi_n
  if (solver == 1) then
    call AMPERE(N_t, d_t, d_x, j_e, n_e, E_x_n, E_x_np1, phi_n)
  else 
    call POISSON(d_x, n_e, E_x_n, E_x_np1, phi_n)
  end if
end subroutine MAXWELL_SOLVER

subroutine POISSON(d_x, n_e, E_x_n, E_x_np1, phi_n)
  implicit none
  real(PR), intent(in)                         :: d_x
  real(PR), dimension(-1:N_x+2), intent(in)    :: n_e
  real(PR), dimension(-1:N_x+2), intent(inout) :: E_x_n
  real(PR), dimension(-1:N_x+2), intent(out)   :: E_x_np1, phi_n
  integer                                      :: i
  real(PR), dimension(1:N_x+1)                 :: a, b, c, d, e, phi_temp 
  !omp PARALLEL DO DEFAULT(SHARED) PRIVATE(i) COLLAPSE(1)
  do i=1,N_x,1
    a(i)        = -1._PR
    b(i)        =  2._PR
    c(i)        = -1._PR
    d(i)        = (1._PR - n_e(i))*(d_x**2._PR)
    if (b_cond.eq.2) then
      e(i)        = 0._PR
      phi_temp(i) = 0._PR
    end if
  end do
  !omp END PARALLEL DO
  select case (b_cond)
    ! absorbing boundary conditions
    case (1)
      call SOLVE_TRIDIAG(a(1:N_x),b(1:N_x),c(1:N_x),d(1:N_x),phi_n(1:N_x)   ,N_x)
      phi_n(N_x+1) = 0._PR
      phi_n(N_x+2) = phi_n(N_x) - phi_n(N_x-1)
      phi_n(0)     = 0._PR
      phi_n(-1)    = phi_n(1) - phi_n(2)
    ! periodic boundary conditions
    case (2)
      call SOLVE_TRIDIAG(a(1:N_x-1),b(1:N_x-1),c(1:N_x-1),d(1:N_x-1),phi_n(1:N_x-1)   ,N_x-1)
      e(1)   = a(1)
      e(N_x-1) = c(N_x-1)
      call SOLVE_TRIDIAG(a(1:N_x-1),b(1:N_x-1),c(1:N_x-1),e(1:N_x-1),phi_temp(1:N_x-1),N_x-1)
      phi_n(N_x) = (d(N_x)-(c(N_x)*phi_n(1))   -(a(N_x)*phi_n(N_x-1))) &
                &/ (b(N_x)+(c(N_x)*phi_temp(1))+(a(N_x)*phi_temp(N_x-1))) 
      phi_n(N_x+1) = phi_n(1)
      phi_n(N_x+2) = phi_n(2)
      phi_n(0)  = phi_n(N_x)
      phi_n(-1) = phi_n(N_x-1)
  end select
  !omp PARALLEL DO DEFAULT(SHARED) PRIVATE(i) COLLAPSE(1)
  do i=1,N_x,1
    E_x_n(i) = (phi_n(i-1) - phi_n(i+1))/(2._PR*d_x)
  end do
  !omp END PARALLEL DO
  select case (b_cond)
    ! absorbing
    case (1)
      E_x_n(0)     = E_x_n(1)
      E_x_n(-1)    = E_x_n(0)
      E_x_n(N_x+1) = E_x_np1(N_x)
      E_x_n(N_x+2) = E_x_n(N_x+1)
      E_x_np1(0)     = E_x_np1(1)
      E_x_np1(-1)    = E_x_np1(0)
      E_x_np1(N_x+1) = E_x_np1(N_x)
      E_x_np1(N_x+2) = E_x_np1(N_x+1)
    ! periodic
    case (2)
      E_x_np1(0)     = E_x_np1(N_x)
      E_x_np1(-1)    = E_x_np1(N_x-1)
      E_x_np1(N_x+1) = E_x_np1(1)
      E_x_np1(N_x+2) = E_x_np1(2)
      E_x_n(0)     = E_x_n(N_x)
      E_x_n(-1)    = E_x_n(N_x-1)
      E_x_n(N_x+1) = E_x_n(1)
      E_x_n(N_x+2) = E_x_n(2)
  end select
end subroutine POISSON

subroutine AMPERE(N_t, d_t, d_x, j_e, n_e, E_x_n, E_x_np1, phi_n)
  implicit none
  integer, intent(in)                          :: N_t
  real(PR), intent(in)                         :: d_t, d_x
  real(PR), dimension(-1:N_x+2), intent(in)    :: j_e, n_e
  real(PR), dimension(-1:N_x+2), intent(inout) :: E_x_n
  real(PR), dimension(-1:N_x+2), intent(out)   :: E_x_np1, phi_n
  integer                                      :: i
  if (N_t.eq.1) then
    call POISSON(d_x, n_e, E_x_n, E_x_np1, phi_n)
  else
    !omp PARALLEL DO DEFAULT(SHARED) PRIVATE(i) COLLAPSE(1)
    do i=1,N_x,1
      E_x_np1(i) = E_x_n(i) - (d_t*j_e(i))
    end do
    !omp END PARALLEL DO
    select case (b_cond)
      ! absorbing
      case (1)
        E_x_n(0)     = E_x_n(1)
        E_x_n(-1)    = E_x_n(0)
        E_x_n(N_x+1) = E_x_np1(N_x)
        E_x_n(N_x+2) = E_x_n(N_x+1)
        E_x_np1(0)     = E_x_np1(1)
        E_x_np1(-1)    = E_x_np1(0)
        E_x_np1(N_x+1) = E_x_np1(N_x)
        E_x_np1(N_x+2) = E_x_np1(N_x+1)
      ! periodic
      case (2)
        E_x_np1(0)     = E_x_np1(N_x)
        E_x_np1(-1)    = E_x_np1(N_x-1)
        E_x_np1(N_x+1) = E_x_np1(1)
        E_x_np1(N_x+2) = E_x_np1(2)
        E_x_n(0)     = E_x_n(N_x)
        E_x_n(-1)    = E_x_n(N_x-1)
        E_x_n(N_x+1) = E_x_n(1)
        E_x_n(N_x+2) = E_x_n(2)
    end select
  end if
end subroutine AMPERE

subroutine DRIVE(d_t, time, x, E_x_n, E_x_np1)
  implicit none
  real(PR), intent(in)                         :: d_t, time
  real(PR), dimension(-1:N_x+2), intent(in)    :: x
  real(PR), dimension(-1:N_x+2), intent(inout) :: E_x_n
  real(PR), dimension(-1:N_x+2), intent(inout) :: E_x_np1
  integer                                      :: i
 !omp PARALLEL DO DEFAULT(SHARED) PRIVATE(i) COLLAPSE(1)
  do i=1,N_x
    E_x_n(i)   = A * sin( (omega_0 *  time     ) - (k*x(i)) )
    E_x_np1(i) = A * sin( (omega_0 * (time+d_t)) - (k*x(i)) )
  end do
  !omp END PARALLEL DO
  E_x_n(N_x+1)   = E_x_n(1)
  E_x_n(N_x+2)   = E_x_n(2)
  E_x_n(0)       = E_x_n(N_x-1)
  E_x_n(-1)      = E_x_n(N_x-2)
  E_x_np1(N_x+1) = E_x_np1(1)
  E_x_np1(N_x+2) = E_x_np1(2)
  E_x_np1(0)     = E_x_np1(N_x-1)
  E_x_np1(-1)    = E_x_np1(N_x-2) 
end subroutine DRIVE

subroutine FLUXES(scheme, vx, u_max, d_t, d_mu,&
                     & u_im2, u_im1, u_i, u_ip1, u_ip2, flux_l, flux_r) 
  implicit none
  integer, intent(in)   :: scheme
  real(PR), intent(in)  :: vx, u_max, d_t, d_mu, u_im2, u_im1, u_i, u_ip1, u_ip2
  real(PR), intent(out) :: flux_l, flux_r
  real(PR)              :: eps_l, eps_r
  real(PR)              :: t, sigma_im1, sigma_i, sigma_ip1
  if (scheme.eq.NL_MUSCL1) then 
    if (vx.ge.0._PR) then
      if ((u_ip1-u_i).gt.0._PR) then
        eps_r = min(1._PR,2._PR*u_i/(u_ip1-u_i))
      else 
        eps_r = 0._PR
      end if
      if ((u_i-u_im1).gt.0._PR) then
        eps_l = min(1._PR,2._PR*u_im1/(u_i-u_im1))
      else 
        eps_l = 0._PR
      end if
      flux_r = vx*(u_i  +(0.5_PR*eps_r*(u_ip1-u_i)))
      flux_l = vx*(u_im1+(0.5_PR*eps_l*(u_i  -u_im1)))
    else
      if ((u_ip1-u_i).lt.0._PR) then
        eps_r = min(1._PR,-2._PR*u_ip1/(u_ip1-u_i))
      else 
        eps_r = 0._PR
      end if
      if ((u_i-u_im1).lt.0._PR) then
        eps_l = min(1._PR,-2._PR*u_i/(u_i-u_im1))
      else 
        eps_l = 0._PR
      end if
      flux_r = vx*(u_ip1-(0.5_PR*eps_r*(u_ip1-u_i)))
      flux_l = vx*(u_i  -(0.5_PR*eps_l*(u_i-u_im1)))
    end if
  else if (scheme.eq.NL_MUSCL2) then
    if (vx.ge.0._PR) then
      if ((u_ip1-u_i)*(u_i-u_im1).le.0._PR) then
        eps_r = 0._PR
      else if ((u_ip1-u_i).lt.0._PR) then
        eps_r = min(1._PR,-2._PR*(u_max-u_i)/(u_ip1-u_i))
      else
        eps_r = min(1._PR,2._PR*u_i/(u_ip1-u_i))
      end if
      if ((u_ip1-u_i)*(u_i-u_im1).le.0._PR) then
        eps_l = 0._PR
      else if ((u_i-u_im1).lt.0._PR) then
        eps_l = min(1._PR,-2._PR*(u_max-u_im1)/(u_i-u_im1))
      else 
        eps_l = min(1._PR,2._PR*u_im1/(u_i-u_im1))
      end if
      flux_r = vx*(u_i  +(0.5_PR*eps_r*(u_ip1-u_i  )))
      flux_l = vx*(u_im1+(0.5_PR*eps_l*(u_i  -u_im1)))
    else
      if ((u_ip1-u_i)*(u_i-u_im1).le.0._PR) then
        eps_r = 0._PR
      else if ((u_ip1-u_i).gt.0._PR) then
        eps_r = min(1._PR,2._PR*(u_max-u_ip1)/(u_ip1-u_i))
      else
        eps_r = min(1._PR,-2._PR*u_ip1/(u_ip1-u_i))
      end if
      if ((u_ip1-u_i)*(u_i-u_im1).le.0._PR) then
        eps_l = 0._PR
      else if ((u_i-u_im1).gt.0._PR) then
        eps_l = min(1._PR,2._PR*(u_max-u_i)/(u_i-u_im1))
      else 
        eps_l = min(1._PR,-2._PR*u_i/(u_i-u_im1))
      end if
      flux_r = vx*(u_ip1-(0.5_PR*eps_r*(u_ip1-u_i  )))
      flux_l = vx*(u_i  -(0.5_PR*eps_l*(u_i  -u_im1)))
    end if
  else
    sigma_im1 = slope(scheme,u_im2,u_im1,u_i)
    sigma_i   = slope(scheme,u_im1,u_i,u_ip1)
    sigma_ip1 = slope(scheme,u_i,u_ip1,u_ip2)
    t = theta(vx)
    flux_l = ( vx * ( ( (1._PR+t) * u_im1 ) + ( (1._PR-t) * u_i ) ) / 2._PR )    &
         & + ( t * vx * ( 1._PR - ( t * vx * d_t / d_mu ) ) *                    &
         &   ( ( (1._PR+t) * sigma_im1 ) + ( (1._PR-t) * sigma_i   ) ) / 4._PR )
    flux_r = ( vx * ( ( (1._PR+t) * u_i ) + ( (1._PR-t) * u_ip1 ) ) / 2._PR )    &
         & + ( t * vx * ( 1._PR - ( t * vx * d_t / d_mu ) ) *                    &
         &   ( ( (1._PR+t) * sigma_i )   + ( (1._PR-t) * sigma_ip1 ) ) / 4._PR )
  end if
end subroutine FLUXES

subroutine BOUNDARIES(f_np1)  
   implicit none
  real(PR), dimension(-1:N_x+2,-1:N_vx+2),intent(inout) :: f_np1
  integer                                               :: l,i
  select case (b_cond)
    ! absorbing
    case (1)
      do l=-1,N_vx+2,1
        f_np1(N_x+1,l) = 0._PR
        f_np1(N_x+2,l) = 0._PR
        f_np1(0,l)     = 0._PR
        f_np1(-1,l)    = 0._PR
      end do
    ! periodic
    case (2)
      do l=-1,N_vx+2,1
        f_np1(N_x+1,l) = f_np1(1,l)
        f_np1(N_x+2,l) = f_np1(2,l)
        f_np1(0,l)     = f_np1(N_x,l)
        f_np1(-1,l)    = f_np1(N_x-1,l)
      end do
  end select
  do i=-1,N_x+2,1
    f_np1(i,N_vx+1) = 0._PR
    f_np1(i,N_vx+2) = 0._PR
    f_np1(i,0)      = 0._PR
    f_np1(i,-1)     = 0._PR
  end do
end subroutine BOUNDARIES

  subroutine SOLVE_TRIDIAG(a,b,c,d,x,n)
    implicit none
    integer,intent(in)                :: n
    real(PR),dimension(n),intent(in)  :: a,b,c,d
    real(PR),dimension(n),intent(out) :: x
    real(PR),dimension(n)             :: cp,dp
    real(PR)                          :: m
    integer                           :: i
    cp(1) = c(1)/b(1)
    dp(1) = d(1)/b(1)
    do i = 2,n,1
       m = b(i)-(cp(i-1)*a(i))
       cp(i) = c(i)/m
       dp(i) = (d(i)-(dp(i-1)*a(i)))/m
    enddo
    x(n) = dp(n)
    do i = n-1, 1, -1
       x(i) = dp(i)-cp(i)*x(i+1)
    end do 
  end subroutine SOLVE_TRIDIAG

! Functions
  
function slope(scheme,u_im1, u_i, u_ip1)
  implicit none
  integer, intent(in)  :: scheme
  real(PR), intent(in) :: u_im1, u_i, u_ip1
  real(PR)             :: slope
  select case (scheme)
    case DEFAULT
      slope = 0._PR
    case (L_donor_cell)
      slope = 0._PR
    case (L_Lax_Wendroff)
      slope = u_ip1 - u_i
    case (L_Beam_Warming)
      slope = u_i - u_im1
    case (L_Fromm)
      slope = (u_ip1 - u_im1) / 2._PR
    case (NL_minmod)
      slope = minmod(u_i - u_im1,u_ip1 - u_i)
    case (NL_superbee)
      slope = maxmod(minmod(u_ip1 - u_i,2._PR*(u_i - u_im1)),minmod(2._PR*(u_ip1 - u_i),u_i - u_im1))
    case (NL_Van_Leer)
      slope = minmod_3(b*(u_ip1 - u_i),0.5_PR*(u_ip1 - u_im1),b*(u_i - u_im1))
  end select
end function slope
  
function minmod(a, b)
   implicit none
  real(PR), intent(in) :: a,b
  real(PR)             :: minmod
  if ((a*b).le.0._PR) then
    minmod = 0._PR
  else
    if (abs(a).gt.abs(b)) then
      minmod = b
    else
      minmod = a
    end if
  end if
end function minmod
  
function minmod_3(a, b, c)
   implicit none
  real(PR), intent(in) :: a, b, c
  real(PR)             :: minmod_3
  minmod_3 = max(0._PR,min(a,b,c)) + min(0._PR,max(a,b,c))
end function minmod_3
  
function maxmod(a, b)
  implicit none
  real(PR), intent(in) :: a,b
  real(PR)             :: maxmod
  if ((a*b).le.0._PR) then
    maxmod = 0._PR
  else
    if (abs(a).gt.abs(b)) then
      maxmod = a
    else
      maxmod = b
    end if
  end if
end function maxmod
  
function theta(vx)
  implicit none
  real(PR), intent(in) :: vx
  real(PR)             :: theta
  if (vx.ge.0) then
    theta = 1._PR
  else
    theta = - 1._PR
  end if
end function theta
    
end module library
