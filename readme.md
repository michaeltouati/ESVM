## Written by Dr Michaël J TOUATI

# Compiling the code

Modify the makefile as a function of the wished compilation options and the compilers installed on your computer and then type :

```sh
make
```
It is recommended to compile the code with the double-precision option (-fdefault-real-8 for the GNU compiler gfortran or -r8 for the INTEL compiler ifort). 

The compilation can be tested by typing
```sh
make test
```
The tests consist in performing a series of 'diff file1 file2' where :
* file1 is one test simulation terminal output performed with an input deck located in the directory test-cases/Test/ and
* file2 is the terminal output of the corresponding simulation also located in test-cases/Test/ 
concerning the simulation of a drifting plasma at Maxwell-Boltzmann equilibrium.

# Running a simulation

Fill the wished input-deck (all parameters are described inside), eventually check them by typing
```sh
make check
```
and type :
```sh
./esvm
```
or
```sh
make run
```
to execute the simulation.

# Plotting the simulation results

All simulation results are stored in files located in the folder 'results'. 
Python scripts allowing to extract and plot the simulation results are located in the folder 'sources/extract'.
They can be used by simply typing :
```sh
make extract
```
to plot all the results. If you want to plot the 1D1V plasma electron distribution function phase-space density maps in logarithmic scale, just type :
```sh
make extract_logfe  
```
The resulting plots will be stored in the folder 'figures'. It can also be plotted separately during the simulation run :
- the energies scalar plots by typing :
```sh
make extract_energy  
```
- the 1D plasma electron hydrodynamic moments space-time density maps by typing :
```sh
make extract_hydro2D  
```
- the 1D plasma electron hydrodynamic moments scalar plots by typing : 
```sh
make extract_hydro1D
```
- or the 1D1V plasma electron distribution function phase-space density maps by typing :
```sh
make extract_fe 
```

# Cleaning the directory

If you want to remove from the ESVM directory :
- the compilation files and ESVM executable, type :
```sh
make clean
```
- the directory figures/ containing all simulation result plots, type :
```sh
make clean_figures
```
- the directory results/ containing all simulation result data files, type :
```sh
make clean_results
```
- the three previous ones, type :
```sh
make clean_all
```
