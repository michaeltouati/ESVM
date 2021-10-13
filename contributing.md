## Written by Michaël J TOUATI

“Physics is like sex: sure, it may give some practical results, but that's not why we do it.”-Richard P. Feynman

“We need science education to produce scientists, but we need it equally to create literacy in the public. Man has a fundamental urge to comprehend the world about him, and science gives today the only world picture which we can consider as valid. It gives an understanding of the inside of the atom and of the whole universe, or the peculiar properties of the chemical substances and of the manner in which genes duplicate in biology. An educated layman can, of course, not contribute to science, but can enjoy and participate in many scientific discoveries which as constantly made. Such participation was quite common in the 19th century, but has unhappily declined. Literacy in science will enrich a person's life.”-Hans Bethe

# Contributing to ESVM

Thank you very much for your interest in contributing to [ESVM](https://github.com/michaeltouati/ESVM). Sharing Richard P. Feynmann enthusiasm to study Mother Nature Physics laws (and to play drum too 😄) while also deeply convinced by this Hans Bethe statement, it seemed to me natural to develop [ESVM](https://github.com/michaeltouati/ESVM) with Github in order to try to understand the complexity of plasmas, numerical analysis and High Parallel Computing (HPC) and to make accessible this small piece of Modern theoretical and computational Physics knowledge to the largest Public as possible at the same time.
In the following can be found how to :
- report a bug you faced when compiling or using [ESVM](https://github.com/michaeltouati/ESVM)
- propose new features that would make [ESVM](https://github.com/michaeltouati/ESVM) better
- fix a bug or submit a new feature

# Report a bug

For problems related with plotting, please make sure first LaTeX, dvipng and ghostscript softwares are each working and on your PATH for the Matplotlib Python package to be able to render tex fonts. If the bug persists or if it is related to another problem, follow these steps :
1) Go to ['Issues'](https://github.com/michaeltouati/ESVM/issues) on the [ESVM master branch](https://github.com/michaeltouati/ESVM) 
2) Click on 'New issue'
4) Describe the bug the more "clear and concise" as possible in the title starting with "Bug :"
5) Describe with details the bug providing the more pieces of information as possible such as input-deck, terminal output and/or screenshot, etc... but more importantly, your environment :
- OS: [e.g., Ubuntu 20.04]
- Fortran compiler [e.g., gfortran 11.2.0]
- Python version [e.g., Python 3.7.11]
- Matplotlib Python package version [e.g., Matplotlib 3.4.3]
- Numpy Python package version [e.g., Numpy 1.21.2]
6) Click on 'Submit new issue'

# Propose a new feature

If you want to propose a new feature in [ESVM](https://github.com/michaeltouati/ESVM) that is not already in the perspectives of the code mentioned in the last section of the documentation [esvm.pdf](https://github.com/michaeltouati/ESVM/blob/master/esvm.pdf), follow these steps :
1) Go to ['Issues'](https://github.com/michaeltouati/ESVM/issues) on the [ESVM master branch](https://github.com/michaeltouati/ESVM)
2) Click on 'New issue'
4) Describe the new feature request the more "clear and concise" as possible in the title starting with "Feature request :"
5) Describe with details the requested feature and
6) Click on 'Submit new issue'

It will be a pleasure to discuss about it.

# Fix a bug or submit a new feature

[ESVM](https://github.com/michaeltouati/ESVM) uses the [Fork and pull model](https://docs.github.com/en/github/collaborating-with-pull-requests/getting-started/about-collaborative-development-models). In order to fix a bug or to submit a new feature that you've added in [ESVM](https://github.com/michaeltouati/ESVM), follow these steps:

1) Fork the repo and create your own branch from the [ESVM 'master' branch](https://github.com/michaeltouati/ESVM).
2) Add your code. Please, keep the code Fortran 95 standard compliant with same coding style by : 
- respecting 2 spaces for indentation rather than tabs
- not using object oriented Fortran 2003 features
- not forgetting to deallocate allocated arrays
- not using allocatable arrays with assumed-shape dummy argument in subroutines or functions for bounds-preservation
- etc...
3) If you've added a new feature that needs new simulation parameters in the [input-deck](https://github.com/michaeltouati/ESVM/blob/master/input-deck) and updated correspondingly the source file [input.f90](https://github.com/michaeltouati/ESVM/blob/master/sources/input.f90), please add their descriptions in the [input-deck](https://github.com/michaeltouati/ESVM/blob/master/input-deck) following the same style
```sh
##                                                                   ##
## T  = electron temperature in eV                                   ##
##                                                                   ##
#######################################################################
#
#T  1000.
#
```
4) If you've added a new feature, let's call it 'my-new-feature', please, create a directory 'my-new-feature' in 'test-cases/Tests/New-features' and add a typical input deck inside that uses your new feature with the simulation name `#simu New-features` and parameters such as :
- `#x_min 0.`, `#x_max 5.`, `#d_x 0.25`, `#vx_min -5.`, `#vx_max 5.`, `#d_vx 0.1` and
- `#cfl 9.e-1`, `#L_t 5.`, `#dt_diag 0.25`
for the test to be fast and a file entitled 'output' containing the terminal output generated by the corresponding simulation. For this, you can just type from the directory that you've created 'test-cases/Tests/New-features/my-new-feature'
```sh
../../../.././esvm > output
```
and remove the directory generated 'results'. Finally, please, add in the makefile at the section ```## TESTING ##``` (after the line ```TEST_DIR=' Maxwell/Ampere'```) the following line
```sh
TEST_DIR += New-features/my-new-feature
```
5) Please, mention your contribution following the same style at the beginning of :
- the modified Fortran source file
```
!! Initial commit written by Michaël J TOUATI - Dec. 2015
!! my-new-feature or my-bug-fix commit by my-name - my-add-date
```
- the modified Python script
```
# Initial commit written by Michaël J Touati - Dec. 2015
# my-new-feature or my-bug-fix commit by my-name - my-add-date
```
- the input-deck if 3.
```sh
##   Initial commit written by Michaël J TOUATI - Dec. 2015          ##
##   my-new-feature commit by my-name - my-add-date                  ##
```
- the makefile if 4.
```sh
## Initial commit written by Michaël J TOUATI - Dec. 2015
## my-new-feature commit by my-name - my-add-date
```
6) Make sure your code passes all the tests by typing on your terminal:
```sh
make test
```
7) Ensure the code compiles with the following compiler debugging options:
- `-g -traceback -fopenmp -r8 -std95 -fpe0 -debug all -debug-parameters all -C` for the INTEL compiler ifort
- `-fdefault-real-8 -O -g -fopenmp -Wall -fcheck=all -fbacktrace -std=f95 -fall-intrinsics -ffpe-trap=invalid,zero,overflow` for the GNU compiler gfortran
8) SAVE YOUR IMPORTANT SIMULATIONS RESULTS AND/OR FIGURES somewhere else, clean the repo by typing
```sh
make clean
```
and issue that pull request! :relieved:

# License
When you submit code changes, your submissions are understood to be under the same [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) that covers [ESVM](https://github.com/michaeltouati/ESVM). 
