---
title: 'ESVM : An Electrostatic Vlasov-Maxwell Open-source Code'
tags:
  - Fortran
  - OpenMP
  - Python
  - Electrostatic Plasma
  - Poisson versus Maxwell-Ampere solver
  - Vlasov equation
  - donor cell
  - Lax-Wendroff
  - beam warming
  - Fromm 
  - minmod
  - superbee
  - Van Leer
  - MUSCL1
  - MUSCL2
  - finite volume numerical schemes for phase-space advection
  - Linear Landau damping
  - Non-linear Landau damping
  - Two-stream instability
  - Electrostatic wakefield
authors:
  - name: Michaël J-M R TOUATI
    orcid: 0000-0001-7590-0941
    affiliation: "1, 2, 3"
affiliations:
 - name: Department of Electrical Engineering, University of California Los Angeles, Los Angeles, CA 90095, USA
   index: 1
 - name: Group of Lasers and Plasmas, IPFN, IST, Universidade de Lisboa, Lisbon, Portugal
   index: 2
 - name: Centro de Láseres Pulsados de Salamanca (CLPU), Edificio M5, Parque Cientfico, C/ Adaja 8, 37185 Villamayor, Salamanca, Spain (current affiliation)
   index: 3
date: 5 August 2021
bibliography: paper.bib
---

# Summary

ESVM (ElectroStatic Vlasov-Maxwell) is a 1D-1V Vlasov-Maxwell Fortran code parallelized using OpenMP and developed 
to adapt simulations to specific Plasma Physics problems with explicit linear finite volume numerical advection schemes such as the donor cell, the Lax-Wendroff, the beam warming or the Fromm method and non-linear ones such as the minmod, the superbee, the Van Leer, the MUSCL1 or the MUSCL2 method coupled to the Maxwell-Gauss equation for the electrostatic field or the Maxwell-Ampere equation with Poisson equation computed at the first time step. Python scripts, using matplotlib and numpy packages, are provided to automatically extract the simulation results, to plot them and to save them. Compilation rules can be easily modified depending on the user compiler preferences using the provided makefile.

# Statement of need

`ESVM` allows for adapting simulations to a specific Plasma Physics problem : it is possible to chose a linear finite volume numerical advection schemes such as the donor cell, the Lax-Wendroff, the beam warming or the Fromm method or a non-linear one such as the minmod, the superbee, the Van Leer, the MUSCL1 or the MUSCL2 method in order to compute the Vlasov equation phase-space derivatives and to chose between computing the Poisson equation versus computing the Maxwell-Ampere equation (with Poisson equation computed at the first time step only). Well known academic Plasma Physics problems have been used to validate the code and are provided as examples (input deck and obtained simulation results). It is planned to :
- extend the code to relativistic 2D-2V phase-space simulations
- implement the Perfectly Matched Layer technique for Electromagnetic Fields absorption at the spatial simulation box boundaries
- implement its MPI parallelization
- implement the Belyaev-Budker relativistic collision operator

Four academic Plasma Physics cases are provided :
1) the linear Landau damping of an electrostatic wave; cf. \autoref{fig:linear-landau-damping}, 
2) the non-linear Landau damping of an electrostatic wave; cf. \autoref{fig:non-linear-landau-damping} and \autoref{fig:non-linear-landau-damping-2}, 
3) the two-stream instability; cf. \autoref{fig:two-stream-instability} and 
4) the emission of an electrostatic wakefield  by a Gaussian (in space and velocity-space) electron drifting at a mean velocity higher than the plasma electron thermal velocity; cf. \autoref{fig:electrostatic-wakefield}. For each Academic case, an example of input deck is provided together with the corresponding simulation result plots that the code typically generates.

# Mathematics

The equations computed explicitely by the code are the 1D-1V Vlasov equation for plasma electrons (ions are assumed to be fully ionized with an electrical charge $Z e$, $Z$ being the ion atomic number, and that they remain immobile with a density $n_i$): 
\begin{equation}
\label{eq:vlasov1d1v}
\displaystyle \frac{\partial f_e}{\partial t} (x,v_x,t) + \displaystyle \frac{\partial }{\partial x} \displaystyle \left ( v_x f_e(x,v_x,t) \right ) - \displaystyle \frac{\partial }{\partial v_x} \displaystyle \left ( \displaystyle \frac{e}{m_e} E_x (x,t) f_e (x,v_x,t)\right ) = 0
\end{equation}
and the coupled Poisson equation for the electrostatic field 
$$
\displaystyle \left \{ \begin{array}{l}
    \displaystyle \frac{\partial \Phi}{\partial x} (x,t) = - E_x (x,t)
\cr \displaystyle \frac{\partial E_x}{\partial x} (x,t) = 4 \pi \displaystyle \left ( Z e n_i - e \displaystyle \int_{-\infty}^\infty f_e (x,v_x,t) \, d v_x \right )
\end{array} \right .
$$
\begin{equation}
\label{eq:poisson}
\Rightarrow \displaystyle \frac{\partial^2 \Phi}{\partial x^2} (x,t) = - 4 \pi \displaystyle \left ( Z e n_i - e \displaystyle \int_{-\infty}^\infty f_e (x,v_x,t) \, d v_x\right )
\end{equation}
or equivalently, the coupled Maxwell-Ampere equation with Poisson equation computed at $t=0$ only
\begin{equation}
\label{eq:ampere}
\displaystyle \left \{ \begin{array}{l}
    \displaystyle \frac{\partial^2 \Phi}{\partial x^2} (x,t=0) = - 4 \pi \displaystyle \left ( Z e n_i - e \displaystyle \int_{-\infty}^\infty f_e (x,v_x,t=0) \, d v_x\right )
\cr  \displaystyle \frac{\partial E_x }{\partial t } (x,t) = 4 \pi e \displaystyle \int_{-\infty}^\infty f_e (x,v_x,t) v_x \, d v_x
\end{array} \right .
\end{equation}

The code units consist in the commonly used electrostatic units : the electron mass $m_e$ for masses, the elementary charge $e$ for electrical charges, the inverse of the Langmuir plasma electron angular frequency $\omega_{p} = \displaystyle \sqrt{ 4 \pi Z n_i e^2 / m_e}$ for times, the Debye electron screening length $\lambda_{\mathrm{Debye}} = \displaystyle \sqrt{k_B T_e / 4 \pi Z n_i e^2}$ where $k_B$ is the Boltzmann constant and $T_e$ the plasma electron temperature (and therefore the thermal plasma electron velocity $v_{T} = \lambda_{\mathrm{Debye}} \omega_{p}$ for velocities) and the constant ion density $n_i$ for densities ($\underline{f}_e = f_e v_{T_e} / n_i$). The resulting electrostatic field unit consequently reads $\underline{E}_x = e E_x / m_e \omega_{p} v_{T}$.

Obviously, the spatial grid cells $\Delta x$ must be chosen lower than the Debye length for the simulations to be Physical, the velocity bins size $\Delta v_x$ and extrema $[v_{\mathrm{min}},v_{\mathrm{max}}]$ in agreement with the considered Plasma Physics problem. The CFL stability criterium is taken into account inside the code so that the user just needs to specify in the input deck the scalar parameter $\mathrm{cfl}$ such that the simulation time step respects
\begin{equation}
\Delta t = \mathrm{cfl} \times F(\Delta x, \Delta v_x) < F(\Delta x, \Delta v_x)
\end{equation}
where $F(\Delta x, \Delta v_x)$ depends on the chosen numerical scheme and is implemented in a code subroutine.

# Figures

![Linear Landau damping test case : Electrostatic field energy and Plasma electron kinetic energy versus time.\label{fig:linear-landau-damping}](test-cases/Linear-Landau-Damping/figures-Poisson/energy.png)

![Non Linear Landau damping test case : Electrostatic field energy and Plasma electron kinetic energy versus time.\label{fig:non-linear-landau-damping}](test-cases/Non-Linear-Landau-Damping/figures-Poisson/energy.png)

![Non Linear Landau damping test case : Plasma electrons phase-space.\label{fig:non-linear-landau-damping-2}](test-cases/Non-Linear-Landau-Damping/figures-Poisson/f_log/f_log_69.png)

![Two stream instability test case : Plasma electrons phase-space.\label{fig:two-stream-instability}](test-cases/Two-Stream-Instability/figures-Poisson/f/f_81.png)

![Electrostatic wakefield test case : Electrostatic wakefield.\label{fig:electrostatic-wakefield}](test-cases/Wakefield-Emission/figures-Poisson/Ex/Ex_30.png)
