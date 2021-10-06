Written by Michael J TOUATI - https://github.com/michaeltouati

# Compiling the code

Modify the makefile as a function of the wished compilation options and the compilers installed on your computer and then type :

make

It is recommended to compile the code with the double-precision option (-fdefault-real-8 for the GNU compiler gfortran or -r8 for the INTEL compiler ifort). 

The compilation can be tested by typing

make check

The tests consist in performing a series of 'diff file1 file2' where :
* file1 is one test simulation terminal output performed with an input deck located in the directory test-cases/Test/ and
* file2 is the terminal output of the corresponding simulation also located in test-cases/Test/ 
concerning the simulation of a drifting plasma at Maxwell-Boltzmann equilibrium.

# Running a simulation

Fill the wished input-deck (all parameters are described inside) and type :

./esvm

The resulting simulation parameters can be checked by typing

make check

# Plotting the simulation results

All simulation results are stored in text files located in the folder 'results'. 
Python scripts allowing to extract and plot the simulation results are located in the folder 'sources'.
They can be used by simply typing :

make extract

If you want to plot the distribution function phase-space in logarithmic scale, just type :

make extract_logfe  

The resulting plots will be located in the folder 'figures'.

# Cleaning the directory

If you want to remove from the ESVM directory :

- the compilation files and ESVM executable, type :

make clean

- the directory figures/ containing all simulation result plots, type :

make clean_figures

- the directory results/ containing all simulation result data files, type :

make clean_results

- the three previous ones, type :

make clean_all
