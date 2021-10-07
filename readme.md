# Compiling the code

Modify the makefile as a function of the wished compilation options and the compilers installed on your computer and then type :

```sh
make
```

It is recommended to compile the code with the double-precision option (-fdefault-real-8 for the GNU compiler gfortran or -r8 for the INTEL compiler ifort). 

The compilation can be tested by typing

```sh
make check
```

The tests consist in performing a series of 'diff file1 file2' where :
* file1 is one test simulation terminal output performed with an input deck located in the directory test-cases/Test/ and
* file2 is the terminal output of the corresponding simulation also located in test-cases/Test/ 
concerning the simulation of a drifting plasma at Maxwell-Boltzmann equilibrium.

# Running a simulation

Fill the wished input-deck (all parameters are described inside) and type :

```sh
./esvm
```

The resulting simulation parameters can be checked by typing

```sh
make check
```

# Plotting the simulation results

All simulation results are stored in text files located in the folder 'results'. 
Python scripts allowing to extract and plot the simulation results are located in the folder 'sources'.
They can be used by simply typing :

```sh
make extract
```

If you want to plot 1D1V distribution function phase-space density maps at each dumped time iteration in logarithmic scale, just type :

```sh
make extract_logfe  
```

The resulting plots will be located in the folder 'figures'.

You can also plot separately during the simulation run :
- the energy plots by typing :
```sh
make extract_energy  
```
- the hydrodynamic moments space-time density maps by typing :
```sh
make extract_hydro2D  
```
- the hydrodynamic moments scalar plots by typing : 
```sh
make extract_hydro1D
```
- or the 1D1V distribution function phase-space density maps by typing :
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
