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
module library

use acuracy
use constants
  
implicit none
public  :: GRID, INIT_VAR, INIT_SIMU
public  :: DENSITIES, ENERGIES
public  :: POISSON, AMPERE
public  :: INIT_NEXT_STEP, MAXWELL_SOLVER
public  :: DRIVE, FLUXES, FE_BOUNDARIES
private :: FIELD_BOUNDARIES, slope, minmod
private :: minmod_3, maxmod, theta

contains

! Subroutines
  
subroutine GRID(Nx, Nvx, dx, dvx, &
              & xmin, vmin, x0, vx0)
  implicit none
  integer,  intent(in)                        :: Nx, Nvx
  real(PR), intent(in)                        :: dx, dvx, xmin, vmin
  real(PR), dimension(-1:Nx+2),  intent(out)  :: x0
  real(PR), dimension(-1:Nvx+2), intent(out)  :: vx0
  integer                                     :: i, l
  !
  !$omp  PARALLEL DO DEFAULT(NONE) &
  !$omp& SHARED(Nx, xmin, dx, x0) &
  !$omp& PRIVATE(i) COLLAPSE(1)
  do i=-1,Nx+2,1
    x0(i)  = xmin + real(i-1,PR)*dx
  end do
  !$omp END PARALLEL DO
  !
  !$omp  PARALLEL DO DEFAULT(NONE) &
  !$omp& SHARED(Nvx, vmin, dvx, vx0) &
  !$omp& PRIVATE(l) COLLAPSE(1)
  do l=-1,Nvx+2,1
    vx0(l) = vmin + real(l-1,PR)*dvx
  end do
  !$omp END PARALLEL DO
end subroutine GRID

subroutine INIT_VAR(Nx, Nvx, fn, fnp1, &
                  & ne, je, ve, vTe, &
                  & Exnp1, Exn, phin, &
                  & dUK, dUT, dUE, &
                  & UK, UT, UE, t0, Nt, & 
                  & positivity, results)
  implicit none
  integer,                               intent(in)    :: Nx, Nvx
  real(PR), dimension(-1:Nx+2,-1:Nvx+2), intent(out)   :: fn, fnp1
  real(PR), dimension(-1:Nx+2)         , intent(out)   :: ne, je, ve, vTe
  real(PR), dimension(-1:Nx+2)         , intent(out)   :: Exnp1, Exn , phin
  real(PR), dimension(1:Nx)            , intent(out)   :: dUK, dUT, dUE
  real(PR)                             , intent(out)   :: UK, UT, UE, t0
  integer                              , intent(out)   :: Nt
  logical                              , intent(out)   :: positivity, results
  integer                                              :: i, l
  !
  !$omp  PARALLEL DO DEFAULT(NONE) &
  !$omp& SHARED(Nx, Nvx, fn, fnp1) &
  !$omp& PRIVATE(i,l) COLLAPSE(2)
  do l=-1,Nvx+2,1
    do i=-1,Nx+2,1
      fn(i,l)    = zero
      fnp1(i,l)  = zero
    end do
  end do
  !$omp END PARALLEL DO
  !
  !$omp  PARALLEL DO DEFAULT(NONE) &
  !$omp& SHARED(Nx, ne, je, ve, vTe, Exn, Exnp1, phin) &
  !$omp& PRIVATE(i) COLLAPSE(1)
  do i=-1,Nx+2,1
    ne(i)    = zero
    je(i)    = zero
    ve(i)    = zero
    vTe(i)   = zero
    Exnp1(i) = zero
    Exn(i)   = zero
    phin(i)  = zero
  end do
  !$omp END PARALLEL DO
  !
  !$omp  PARALLEL DO DEFAULT(NONE) &
  !$omp& SHARED(Nx, dUK, dUT, dUE) &
  !$omp& PRIVATE(i) COLLAPSE(1)
  do i=1,Nx,1
    dUK(i) = zero
    dUT(i) = zero
    dUE(i) = zero
  end do
  !$omp END PARALLEL DO
  !
  UK  = zero
  UT  = zero
  UE  = zero
  !
  t0    = 0._PR
  Nt    = 1
  !
  positivity = .false.
  results    = .false.
end subroutine INIT_VAR

subroutine FE_BOUNDARIES(bcond, Nx, Nvx, f0)  
  implicit none
  integer,                              intent(in)    :: bcond
  integer,                              intent(in)    :: Nx, Nvx
  real(PR), dimension(-1:Nx+2,-1:Nvx+2),intent(inout) :: f0
  integer                                             :: l,i
  !
  select case (bcond)
    ! absorbing
    case (1)
      !$omp  PARALLEL DO DEFAULT(NONE) &
      !$omp& SHARED(Nvx, Nx, f0) &
      !$omp& PRIVATE(l) COLLAPSE(1)
      do l=-1,Nvx+2,1
        f0(Nx+1,l) = zero
        f0(Nx+2,l) = zero
        f0(0,l)    = zero
        f0(-1,l)   = zero
      end do
      !$omp END PARALLEL DO
    ! periodic
    case (2)
      !$omp  PARALLEL DO DEFAULT(NONE) &
      !$omp& SHARED(Nvx, Nx, f0) &
      !$omp& PRIVATE(l) COLLAPSE(1)
      do l=-1,Nvx+2,1
        f0(Nx+1,l) = f0(1,l)
        f0(Nx+2,l) = f0(2,l)
        f0(0,l)    = f0(Nx,l)
        f0(-1,l)   = f0(Nx-1,l)
      end do
      !$omp END PARALLEL DO
  end select
  !$omp  PARALLEL DO DEFAULT(NONE) &
  !$omp& SHARED(Nvx, Nx, f0) &
  !$omp& PRIVATE(i) COLLAPSE(1)
  do i=-1,Nx+2,1
    f0(i,Nvx+1) = zero
    f0(i,Nvx+2) = zero
    f0(i,0)     = zero
    f0(i,-1)    = zero
  end do
  !$omp END PARALLEL DO
end subroutine FE_BOUNDARIES

subroutine INIT_SIMU(bcond, academic_case, &
                   & Ap, kp ,vd, &
                   & Nx, Nvx, x0, vx0, f0)
  implicit none
  integer,                               intent(in)    :: bcond
  integer,                               intent(in)    :: academic_case
  real(PR),                              intent(in)    :: Ap, kp ,vd
  integer,                               intent(in)    :: Nx, Nvx
  real(PR), dimension(-1:Nx+2),          intent(in)    :: x0
  real(PR), dimension(-1:Nvx+2),         intent(in)    :: vx0
  real(PR), dimension(-1:Nx+2,-1:Nvx+2), intent(inout) :: f0
  integer                                              :: l, i
  real(PR)                                             :: xs, dx
  real(PR)                                             :: dvx, X2
  real(PR)                                             :: norm
  real(PR)                                             :: norm1, norm2 
  !
  norm = 1.0_PR/sqrt(2._PR*pi)
  !
  select case (academic_case)
    case (1)
    ! Electrostatic wakefield test case :
      dx    = 0.25_PR  ! "particle size"
      dvx   = 0.025_PR ! "particle size"
      xs    = x0(1) + ( (x0(Nx) - x0(1)) / 8._PR ) 
      norm1 = Ap / ( 2._PR * pi * dx * dvx )
      !$omp  PARALLEL DO DEFAULT(NONE) &
      !$omp& SHARED(Nvx, Nx, f0, norm, norm1, vx0, vd, x0, xs, dx, dvx) &
      !$omp& PRIVATE(i,l,X2) COLLAPSE(2)
      do l=1,Nvx,1
        do i=1,Nx,1
          f0(i,l) = norm * exp(-(vx0(l)**2._PR)/2._PR) 
          X2 = -0.5_PR * &
             & (  ( ((x0(i) -xs)/dx )**2._PR) &
             &  + ( ((vx0(l)-vd)/dvx)**2._PR) )
          f0(i,l) = f0(i,l) + ( norm1 * exp(X2) )
        end do
      end do
      !$omp END PARALLEL DO
    ! Landau damping test case :
    case (2)
      !$omp  PARALLEL DO DEFAULT(NONE) &
      !$omp& SHARED(Nvx, Nx, f0, norm, vx0) &
      !$omp& PRIVATE(i,l) COLLAPSE(2)
      do l=1,Nvx,1
        do i=1,Nx,1
          f0(i,l) = norm * exp(-(vx0(l)**2._PR)/2._PR)
         end do
      end do
      !$omp END PARALLEL DO
    case (3)
    ! Two-stream instability test case :
      norm = 0.5_PR * norm
      !$omp  PARALLEL DO DEFAULT(NONE) &
      !$omp& SHARED(Nvx, Nx, f0, norm, Ap, kp, vx0, vd, x0) &
      !$omp& PRIVATE(i,l,norm1,norm2) COLLAPSE(2)
      do l=1,Nvx,1
        do i=1,Nx,1
          norm1   = Ap * sin(kp*x0(i))
          norm2   = 1._PR-norm1
          norm1   = 1._PR+norm1
          f0(i,l) = norm * &
          & ( (norm1*exp(-((vx0(l)-vd)**2._PR)/2._PR)) + &
          &   (norm2*exp(-((vx0(l)+vd)**2._PR)/2._PR)) )
         end do
      end do
      !$omp END PARALLEL DO
    case default
    ! Maxwellian :
      !$omp  PARALLEL DO DEFAULT(NONE) &
      !$omp& SHARED(Nvx, Nx, f0, norm, vx0, vd) &
      !$omp& PRIVATE(i,l) COLLAPSE(2)
      do l=1,Nvx,1
        do i=1,Nx,1
          f0(i,l) = norm * exp(-((vx0(l)-vd)**2._PR)/2._PR)
         end do
      end do
      !$omp END PARALLEL DO
  end select
  ! Boundary conditions
  call FE_BOUNDARIES(bcond, Nx, Nvx, f0)
  !
end subroutine INIT_SIMU

subroutine INIT_NEXT_STEP(Nx, Nvx, fn, fnp1, &
                        & Exn, Exnp1, &
                        & Nt, dt, t0)
  implicit none
  integer,                               intent(in)    :: Nx, Nvx
  real(PR), dimension(-1:Nx+2,-1:Nvx+2), intent(inout) :: fn, fnp1
  real(PR), dimension(-1:Nx+2)         , intent(inout) :: Exnp1, Exn
  real(PR)                             , intent(inout) :: dt, t0
  integer                              , intent(inout) :: Nt
  integer                                              :: l, i
  !
  !$omp  PARALLEL DO DEFAULT(NONE) &
  !$omp& SHARED(Nvx, Nx, fn, fnp1) &
  !$omp& PRIVATE(i,l) COLLAPSE(2)
  do l=-1,Nvx+2,1
    do i=-1,Nx+2,1
      fn(i,l)   = fnp1(i,l)
    end do 
  end do
  !$omp END PARALLEL DO
  !
  !$omp  PARALLEL DO DEFAULT(NONE) &
  !$omp& SHARED(Nx, Exn, Exnp1) &
  !$omp& PRIVATE(i) COLLAPSE(1)
  do i=-1,Nx+2,1
    Exn(i)   = Exnp1(i) 
  end do
  !$omp END PARALLEL DO
  !
  Nt = Nt + 1
  t0 = t0 + dt
end subroutine INIT_NEXT_STEP

subroutine DENSITIES(Nx, Nvx, dvx, vx0, fn, ne, je, ve, vTe)
  implicit none
  integer,                               intent(in)   :: Nx, Nvx
  real(PR),                              intent(in)   :: dvx
  real(PR), dimension(-1:Nvx+2),         intent(in)   :: vx0
  real(PR), dimension(-1:Nx+2,-1:Nvx+2), intent(in)   :: fn
  real(PR), dimension(-1:Nx+2),          intent(out)  :: ne, je
  real(PR), dimension(-1:Nx+2),          intent(out)  :: ve, vTe
  integer                                             :: i, l
  real(PR), dimension(:), allocatable                 :: F1,F2,F3
  !$omp  PARALLEL DO DEFAULT(NONE) &
  !$omp& SHARED(Nvx, Nx, dvx, vx0, fn, ne, je, ve, vTe) &
  !$omp& PRIVATE(i,l,F1,F2,F3) COLLAPSE(1)
  do i=-1,Nx+2,1
    allocate(F1(1:Nvx),F2(1:Nvx),F3(1:Nvx))
    F1(1:Nvx) = fn(i,1:Nvx)*dvx
    do l=1,Nvx,1
      F2(l) = - fn(i,l) * vx0(l) * dvx
      F3(l) = fn(i,l) * (vx0(l)**2._PR) * dvx
    end do
    ne(i)  = sum(F1(1:Nvx))
    je(i)  = sum(F2(1:Nvx))
    if (abs(ne(i)).lt.zero) then
      ve(i) = 0.
    else
      ve(i)  = - je(i) / ne(i)
    end if
    if (abs(ne(i)).ne.(ve(i)**2._PR)) then
      vTe(i) = ((sum(F3(1:Nvx))/ne(i))-(ve(i)**2._PR))**0.5_PR
    else
      vTe(i) = 0.
    end if
    deallocate(F1,F2,F3)
  end do
  !$omp END PARALLEL DO
end subroutine DENSITIES

subroutine ENERGIES(Nx, dx, &
                  & ne, ve, vTe, Exn, &
                  & dUK, dUT, dUE, &
                  & UK, UT, UE)
  implicit none
  integer                     , intent(in)    :: Nx
  real(PR)                    , intent(in)    :: dx
  real(PR), dimension(-1:Nx+2), intent(in)    :: ne, ve, vTe
  real(PR), dimension(-1:Nx+2), intent(in)    :: Exn
  real(PR), dimension(1:Nx)   , intent(inout) :: dUK, dUT, dUE
  real(PR)                    , intent(inout) :: UK, UT, UE
  integer                                     :: i
  !
  !$omp  PARALLEL DO DEFAULT(NONE) &
  !$omp& SHARED(Nx, ne, ve, vTe, Exn, dUK, dUT, dUE) &
  !$omp& PRIVATE(i) COLLAPSE(1)
  do i = 1,Nx,1
    dUK(i) = ne(i) * (ve(i)**2._PR)  / 2._PR
    dUT(i) = ne(i) * (vTe(i)**2._PR) / 2._PR
    dUE(i) = (Exn(i)**2._PR) / 2._PR
  end do
  !$omp END PARALLEL DO
  !
  UK = max(zero, sum(dUK(1:Nx)) * dx)
  UT = max(zero, sum(dUT(1:Nx)) * dx)
  UE = max(zero, sum(dUE(1:Nx)) * dx)
  !
end subroutine ENERGIES

subroutine MAXWELL_SOLVER(solver, bcond, &
                        & Nx, Nt, dt, dx, &
                        & je, ne, Exn, Exnp1, phin)
  implicit none
  integer,                      intent(in)    :: bcond, Nx
  integer,                      intent(in)    :: solver, Nt
  real(PR),                     intent(in)    :: dt, dx
  real(PR), dimension(-1:Nx+2), intent(in)    :: je, ne
  real(PR), dimension(-1:Nx+2), intent(inout) :: Exn
  real(PR), dimension(-1:Nx+2), intent(out)   :: Exnp1, phin
  if (solver == 1) then
    call AMPERE(bcond, Nx, Nt, dt, dx, &
              & je, ne, Exn, Exnp1, phin)
  else 
    call POISSON(bcond, Nx, dx, &
               & ne, Exn, phin)
  end if
end subroutine MAXWELL_SOLVER

subroutine FIELD_BOUNDARIES(bcond, Nx, Ex)
  implicit none
  integer, intent(in)                         :: bcond, Nx
  real(PR), dimension(-1:Nx+2), intent(inout) :: Ex
  !
  select case (bcond)
    ! absorbing
    case (1)
      Ex(0)    = Ex(1)
      Ex(-1)   = Ex(0)
      Ex(Nx+1) = Ex(Nx)
      Ex(Nx+2) = Ex(Nx+1)
    ! periodic
    case (2)
      Ex(0)    = Ex(Nx)
      Ex(-1)   = Ex(Nx-1)
      Ex(Nx+1) = Ex(1)
      Ex(Nx+2) = Ex(2)
  end select
end subroutine FIELD_BOUNDARIES

subroutine POISSON(bcond, Nx, dx, &
                 & ne, Exn, phin)
  implicit none
  integer,                      intent(in)  :: bcond, Nx
  real(PR),                     intent(in)  :: dx
  real(PR), dimension(-1:Nx+2), intent(in)  :: ne
  real(PR), dimension(-1:Nx+2), intent(out) :: Exn
  real(PR), dimension(-1:Nx+2), intent(out) :: phin
  integer                                   :: i
  real(PR), dimension(1:Nx+1)               :: a, b, c, d, e
  real(PR), dimension(1:Nx+1)               :: phi_temp 
  !
  !$omp  PARALLEL DO DEFAULT(NONE) &
  !$omp& SHARED(Nx, dx, ne, bcond, a, b, c, d, e, phi_temp) &
  !$omp& PRIVATE(i) COLLAPSE(1)
  do i=1,Nx,1
    a(i)        = -1._PR
    b(i)        =  2._PR
    c(i)        = -1._PR
    d(i)        = (1._PR - ne(i)) * (dx**2._PR)
    if (bcond.eq.2) then
      e(i)        = 0._PR
      phi_temp(i) = 0._PR
    end if
  end do
  !$omp END PARALLEL DO
  select case (bcond)
    ! absorbing boundary conditions
    case (1)
      call SOLVE_TRIDIAG(Nx, a(1:Nx),b(1:Nx),c(1:Nx),d(1:Nx), &
                       & phin(1:Nx))
      phin(Nx+1) = 0._PR
      phin(Nx+2) = phin(Nx) - phin(Nx-1)
      phin(0)    = 0._PR
      phin(-1)   = phin(1)  - phin(2)
    ! periodic boundary conditions
    case (2)
      call SOLVE_TRIDIAG(Nx-1, a(1:Nx-1),b(1:Nx-1),c(1:Nx-1),d(1:Nx-1), &
                       & phin(1:Nx-1))
      e(1)    = a(1)
      e(Nx-1) = c(Nx-1)
      call SOLVE_TRIDIAG(Nx-1, a(1:Nx-1),b(1:Nx-1),c(1:Nx-1),e(1:Nx-1), &
                       & phi_temp(1:Nx-1))
      phin(Nx) = (d(Nx)-(c(Nx)*phin(1))   -(a(Nx)*phin(Nx-1))) &
                &/ (b(Nx)+(c(Nx)*phi_temp(1))+(a(Nx)*phi_temp(Nx-1))) 
      phin(Nx+1) = phin(1)
      phin(Nx+2) = phin(2)
      phin(0)    = phin(Nx)
      phin(-1)   = phin(Nx-1)
  end select
  !$omp  PARALLEL DO DEFAULT(NONE) &
  !$omp& SHARED(Nx, dx, phin, Exn) &
  !$omp& PRIVATE(i) COLLAPSE(1)
  do i=1,Nx,1
    Exn(i) = (phin(i-1) - phin(i+1))/(2._PR*dx)
  end do
  !$omp END PARALLEL DO
  !
  call FIELD_BOUNDARIES(bcond, Nx, Exn)
end subroutine POISSON

subroutine AMPERE(bcond, Nx, Nt, dt, dx, &
                & je, ne, Exn, Exnp1, phin)
  implicit none
  integer,                      intent(in)    :: bcond, Nx
  integer,                      intent(in)    :: Nt
  real(PR),                     intent(in)    :: dt, dx
  real(PR), dimension(-1:Nx+2), intent(in)    :: je, ne
  real(PR), dimension(-1:Nx+2), intent(inout) :: Exn
  real(PR), dimension(-1:Nx+2), intent(out)   :: Exnp1, phin
  integer                                     :: i
  real(PR)                                    :: dx2
  !
  if (Nt.eq.1) then
    call POISSON(bcond, Nx, dx, ne, Exn, phin)
  else
    dx2 = 2._PR*dx
    !$omp  PARALLEL DO DEFAULT(NONE) &
    !$omp& SHARED(Nx, dt, dx2, je, Exn, phin, Exnp1) &
    !$omp& PRIVATE(i) COLLAPSE(1)
    do i=1,Nx,1
      Exnp1(i) = Exn(i)    - (dt *je(i)   )
      phin(i)  = phin(i-2) - (dx2*Exn(i-1))
    end do
    !$omp END PARALLEL DO
    !
    call FIELD_BOUNDARIES(bcond, Nx, Exnp1)
  end if
end subroutine AMPERE

subroutine DRIVE(Nx, dt, t0, x0, &
               & A1, Om0, K0, &
               & Exn, Exnp1, phin)
  implicit none
  integer,                      intent(in)  :: Nx
  real(PR),                     intent(in)  :: dt, t0
  real(PR),                     intent(in)  :: A1, Om0, K0
  real(PR), dimension(-1:Nx+2), intent(in)  :: x0
  real(PR), dimension(-1:Nx+2), intent(out) :: Exn, phin
  real(PR), dimension(-1:Nx+2), intent(out) :: Exnp1
  integer                                   :: i
  real(PR)                                  :: A2, kx0
  real(PR)                                  :: ot, ot2
  !
  A2   = A1 / K0
  ot  = Om0 * t0
  ot2 = Om0 * (t0+dt)
  !$omp  PARALLEL DO DEFAULT(NONE) &
  !$omp& SHARED(Nx, x0, A1, A2, K0, ot, ot2, Exn, Exnp1, phin) &
  !$omp& PRIVATE(i,kx0) COLLAPSE(1)
  do i=-1,Nx+2
    kx0      = K0 * x0(i)
    Exn(i)   = A1 * sin( ot  - kx0 )
    phin(i)  = A2 * cos( ot  - kx0 )
    Exnp1(i) = A1 * sin( ot2 - kx0 )
  end do
  !$omp END PARALLEL DO
end subroutine DRIVE

subroutine FLUXES(method, bVL, vx0, u_max, dt, d_mu, &
                & u_im2, u_im1, u_i, u_ip1, u_ip2, &
                & flux_l, flux_r) 
  implicit none
  integer,  intent(in)  :: method
  real(PR), intent(in)  :: bVL, vx0, u_max, dt, d_mu
  real(PR), intent(in)  :: u_im2, u_im1, u_i, u_ip1, u_ip2
  real(PR), intent(out) :: flux_l, flux_r
  real(PR)              :: eps_l, eps_r
  real(PR)              :: theta0, sigma_im1, sigma_i, sigma_ip1
  if (method.eq.NL_MUSCL1) then 
    if (vx0.ge.0._PR) then
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
      flux_r = vx0*(u_i  +(0.5_PR*eps_r*(u_ip1-u_i)))
      flux_l = vx0*(u_im1+(0.5_PR*eps_l*(u_i  -u_im1)))
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
      flux_r = vx0*(u_ip1-(0.5_PR*eps_r*(u_ip1-u_i)))
      flux_l = vx0*(u_i  -(0.5_PR*eps_l*(u_i-u_im1)))
    end if
  else if (method.eq.NL_MUSCL2) then
    if (vx0.ge.0._PR) then
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
      flux_r = vx0*(u_i  +(0.5_PR*eps_r*(u_ip1-u_i  )))
      flux_l = vx0*(u_im1+(0.5_PR*eps_l*(u_i  -u_im1)))
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
      flux_r = vx0*(u_ip1-(0.5_PR*eps_r*(u_ip1-u_i  )))
      flux_l = vx0*(u_i  -(0.5_PR*eps_l*(u_i  -u_im1)))
    end if
  else
    sigma_im1 = slope(method, bVL, u_im2, u_im1, u_i  )
    sigma_i   = slope(method, bVL, u_im1, u_i  , u_ip1)
    sigma_ip1 = slope(method, bVL, u_i  , u_ip1, u_ip2)
    theta0 = theta(vx0)
    flux_l = ( vx0 * ( ( (1._PR+theta0) * u_im1 ) + ( (1._PR-theta0) * u_i ) ) / 2._PR )    &
         & + ( theta0 * vx0 * ( 1._PR - ( theta0 * vx0 * dt / d_mu ) ) *                    &
         &   ( ( (1._PR+theta0) * sigma_im1 ) + ( (1._PR-theta0) * sigma_i   ) ) / 4._PR )
    flux_r = ( vx0 * ( ( (1._PR+theta0) * u_i ) + ( (1._PR-theta0) * u_ip1 ) ) / 2._PR )    &
         & + ( theta0 * vx0 * ( 1._PR - ( theta0 * vx0 * dt / d_mu ) ) *                    &
         &   ( ( (1._PR+theta0) * sigma_i )   + ( (1._PR-theta0) * sigma_ip1 ) ) / 4._PR )
  end if
end subroutine FLUXES

subroutine SOLVE_TRIDIAG(N0, a0, b0, c0, d0, x0)
  implicit none
  integer,intent(in)                 :: N0
  real(PR),dimension(N0),intent(in)  :: a0,b0,c0,d0
  real(PR),dimension(N0),intent(out) :: x0
  real(PR),dimension(N0)             :: cp,dp
  real(PR)                           :: yp
  integer                            :: i
  cp(1) = c0(1)/b0(1)
  dp(1) = d0(1)/b0(1)
  do i = 2,N0,1
     yp    = b0(i)-(cp(i-1)*a0(i))
     cp(i) = c0(i)/yp
     dp(i) = (d0(i)-(dp(i-1)*a0(i)))/yp
  end do
  x0(N0) = dp(N0)
  do i = N0-1,1,-1
     x0(i) = dp(i)-cp(i)*x0(i+1)
  end do 
end subroutine SOLVE_TRIDIAG

function slope(method, bVL, u_im1, u_i, u_ip1)
  implicit none
  integer,  intent(in) :: method
  real(PR), intent(in) :: bVL, u_im1, u_i, u_ip1
  real(PR)             :: slope
  select case (method)
    case default
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
      slope = minmod_3(bVL*(u_ip1 - u_i),0.5_PR*(u_ip1 - u_im1),bVL*(u_i - u_im1))
  end select
end function slope
  
function minmod(a0, b0)
  implicit none
  real(PR), intent(in) :: a0, b0
  real(PR)             :: minmod
  if ((a0*b0).le.0._PR) then
    minmod = 0._PR
  else
    if (abs(a0).gt.abs(b0)) then
      minmod = b0
    else
      minmod = a0
    end if
  end if
end function minmod
  
function minmod_3(a0, b0, c0)
  implicit none
  real(PR), intent(in) :: a0, b0, c0
  real(PR)             :: minmod_3
  minmod_3 = max(0._PR,min(a0,b0,c0)) + min(0._PR,max(a0,b0,c0))
end function minmod_3
  
function maxmod(a0, b0)
  implicit none
  real(PR), intent(in) :: a0, b0
  real(PR)             :: maxmod
  if ((a0*b0).le.0._PR) then
    maxmod = 0._PR
  else
    if (abs(a0).gt.abs(b0)) then
      maxmod = a0
    else
      maxmod = b0
    end if
  end if
end function maxmod
  
function theta(vx0)
  implicit none
  real(PR), intent(in) :: vx0
  real(PR)             :: theta
  if (vx0.ge.0) then
    theta =   1._PR
  else
    theta = - 1._PR
  end if
end function theta

end module library
