---
title: 'ESVM: an open-source finite volume Electrostatic Vlasov-Maxwell code'
tags:
  - Fortran
  - OpenMP
  - Python
  - Electrostatic 
  - Collisionless
  - Plasma
  - Poisson 
  - Maxwell-Gauss
  - Maxwell-Ampere
  - Vlasov
  - Advection
  - Finite volume
  - Donor-cell
  - Lax-Wendroff
  - Beam-Warming
  - Fromm 
  - minmod
  - superbee
  - Van Leer
  - MUSCL
  - Landau damping
  - Two-stream instability
  - Electrostatic wakefield
authors:
  - name: Michaël J TOUATI
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

A plasma is a set of charged particles consisting of electrons and ionized atoms whose quantity is sufficiently large to behave collectively through the long-distance electromagnetic fields they produce. It is thought that more than 99.9% of visible matter in the universe is in the plasma state. In a collisionless plasma  consisting in an ionized gas composed of electrons moving inbetween much heavier ions, any electrostatic field is rapidly screened by the plasma electrons over the Debye screening distance [@DebyeHuckel:1923]. When the number of electrons in these Debye spheres can be assumed to be infinite, the plasma electron population is correctly described by the Vlasov equation [@Vlasov:1938] that neglects all correlations between particles such as the binary Coulomb collisions between them. In addition to being simple, the resulting Vlasov-Maxwell set of equations is extremely rich in physics and has many applications ranging from astrophysics and theoretical plasma physics to intense laser-matter interaction experiments. [ESVM](https://github.com/michaeltouati/ESVM) (ElectroStatic Vlasov-Maxwell) is a Vlasov-Maxwell Fortran 95 standard-compliant code, parallelized with OpenMP and using Python 3 for post-processing, that allows for the study of these collisionless plasmas. Many finite volume advection schemes [@Godunov:1959] are implemented in order to discretize the Vlasov equation, namely:

- the donor-cell scheme, i.e., the downwind / upwind scheme [@Courant:1952] depending on the advection direction in each phase-space cell, 
- the Lax-Wendroff scheme [@LaxWendroff:1960], 
- the Fromm scheme [@Fromm:1968],
- the Beam-Warming scheme [@BeamWarming:1976],
- the Van Leer scheme [@VanLeerIII:1977],
- the minmod scheme [@Roe:1986], 
- the superbee scheme [@Roe:1986], and 
- two Monotonic Upwind-centered Scheme for Conservation Laws (MUSCL) [@VanLeerV:1977] schemes MUSCL2 [@Crouseilles:2004] and MUSCL1 [@Duclous:2009]. 

Unlike the linear second order Lax-Wendroff, Fromm, and Beam-Warming schemes, the non-linear second order minmod, superbee, Van Leer, and MUSCL schemes make use of a Total Variation Diminishing (TVD) non-linear flux limiter with the price of becoming a first order scheme in some phase-space cells to limit the numerical oscillations. The donor-cell scheme is a first order method and has the pros of limiting such eventual oscillations but the cons of being numerically less consistent and more diffusive. In ESVM, the discretized Vlasov equation is coupled with the self-consistent Maxwell-Gauss equation or equivalently with the Maxwell-Ampere equation with Maxwell-Gauss equation computed at the first time step, only. While the second order Maxwell-Gauss solver needs a computationally expensive inversion of a tridiagonal matrix for the computation of the Poisson equation, the Maxwell-Ampere equation solver makes use of a faster first order numerical scheme (in time). Both absorbing and periodic boundary conditions for both the particles and the fields are implemented. Python scripts, using the Matplotlib and Numpy packages, are provided to automatically extract and plot the stored simulation results. The simulation parameters are described in the [input deck](https://github.com/michaeltouati/ESVM/blob/master/input-deck) and they can be modified without having to recompile the code. Compilation rules can be modified in the [makefile](https://github.com/michaeltouati/ESVM/blob/master/makefile) depending on the user's compiler preferences. Classical plasma physics academic case simulations that need less than one CPU$\times$hour each, tools for testing the compilation of the code, and tools for checking the simulation parameters are provided. 

# Statement of need

[ESVM](https://github.com/michaeltouati/ESVM) has been developed in order to adapt simulations to specific plasma physics problems by chosing the more adequate finite volume numerical advection scheme in order to compute the Vlasov equation phase-space advection derivatives and to chose between computing the Maxwell-Gauss equation or the Maxwell-Ampere equation with Maxwell-Gauss equation computed at the first time step, only. The code aims at beeing used by the open-source highly parallel computing (HPC) plasma physics community, ranging from under or post-graduate students to teachers and researchers who usually use Particle-In-Cell (PIC) codes [@Dawson:1962] to study collisionless plasmas. Indeed, the PIC method may prohibit the study of plasma physical processes on large time scales and/or for very dense collisionless plasmas due to the statistical and numerical fluctuations of the computed quantities imposed by the use of a finite number of particles. Also, plasma instabilities naturally develop in PIC codes, seeded by the available fluctuations spatial spectrum k-vector for which the instability growth rate is maximum and some small amplitude plasma physical processes may be hidden under the fluctuactions level. Compared to the many open source PIC codes such as [Smilei](https://github.com/SmileiPIC/Smilei) [@Derouillat:2018] and semi-Lagrangian codes such as [vmf90](https://github.com/pdebuyl/vmf90) [@Debuyl:2014], there are no open source finite volume Vlasov-Maxwell codes in the literature that are not based on an expansion method such as [OSHUN](https://github.com/UCLA-Plasma-Simulation-Group/OSHUN) [@Tzoufras:2011],  [AMoRE](https://github.com/michaeltouati/AMoRE) [@Touati:2014], or [Vlapy](https://github.com/joglekara/VlaPy) [@Joglekar:2020]. Finally, since the Vlasov equation is a conservation equation of the probable number of particles in the phase-space, using a finite volume method in order to compute the Vlasov equation has the advantage of allowing for the use of numerical schemes that are numerically flux conserving and/or that ensure the distribution function positivity compared to other numerical methods. ESVM has already been used during courses for under and post-graduate students about the "numerical tools for laser-plasma interaction Physics" and it is currently used for theoretical Plasma Physics investigations.

# Future work

The author plans, in the near future, to:

1. provide another plasma physics academic simulation about one non-linear BGK electron plasma wave, from the name of its finder I. B. Bernstein, J. M. Greene, and M. D. Kruskal [@BernsteinGreenKruskal:1957]
2. provide another plasma physics academic simulation about the echo of two plasma electron waves [@Gould:1967]
3. implement non-equally spaced phase-space cells
4. implement high order Weighted Essentially Non-Oscillatory (WENO) advection schemes [@Liu:1994]
5. compute the plasma ion Vlasov equation to allow for the ions to be mobile 
6. store the simulation results in hdf5 files instead of text files
7. implement MPI parallelization
8. implement vectorization
9. extend the code to the relativistic regime: ESVM $\Rightarrow$ RESVM for open source Relativistic ElectroStatic Vlasov-Maxwell code
10. implement a BGK collision operator, from the name of its finder P. L. Bhatnagar, E. P. Gross,  and M. Krook [@BhatnagarGrossKrook:1954]
11. extend the code to 1D-2V and 1D-3V phase-space electrostatic plasma simulations
12. implement the Landau [@Landau:1936] and Belaiev-Budker [@BelaievBudker:1957] relativistic collision operators using the Rosenbluth potentials [@Rosenbluth:1957] and their relativistic Braams-Karney extension [@BraamsKarney:1987]: (R)ESVM $\Rightarrow$ (R)ESVFPM for open source (Relativistic) ElectroStatic Vlasov-Fokker-Planck-Maxwell code
13. extend the code to electrostatic 2D-1V, 2D-2V, and 2D-3V phase-space plasma simulations: (R)ESV(FP)M $\Rightarrow$ (R)ESV(FP)M2 for open source (Relativistic) ElectroStatic Vlasov-(Fokker-Planck-)Maxwell in 2D
14. extend the code with the second order finite difference Yee scheme [@Yee:1966] to electromagnetic 2D-1V, 2D-2V, and 2D-3V phase-space plasma simulations: (R)ESV(FP)M(2) $\Rightarrow$ (R)EMV(FP)M(2) for open source (Relativistic) ElectroMagnetic Vlasov-(Fokker-Planck)-Maxwell (in 2D)
15. implement the Perfectly Matched Layer (PML) technique [@Berenger:1994] to absorb the electromagnetic fields at the spatial simulation box boundaries
16. deploy the code to GPU architectures.

# References
