## Written by Dr Michaël J TOUATI

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

If you want to propose a new feature in [ESVM](https://github.com/michaeltouati/ESVM) that is not already in the perspectives of the code mentioned in the last section of the documentation [esvm.pdf](https://github.com/michaeltouati/ESVM/blob/main/esvm.pdf), follow these steps :
1) Go to ['Issues'](https://github.com/michaeltouati/ESVM/issues) on the [ESVM master branch](https://github.com/michaeltouati/ESVM)
2) Click on 'New issue'
4) Describe the new feature request the more "clear and concise" as possible in the title starting with "Feature request :"
5) Describe with details the requested feature and
6) Click on 'Submit new issue'

It will be a pleasure to discuss about it.

# Fix a bug or submit a new feature

[ESVM](https://github.com/michaeltouati/ESVM) uses the [Fork and pull model](https://docs.github.com/en/github/collaborating-with-pull-requests/getting-started/about-collaborative-development-models). In order to fix a bug or to submit a new feature that you've added in [ESVM](https://github.com/michaeltouati/ESVM), follow these steps:

1) Fork the repo and create your own branch from [ESVM 'main' branch](https://github.com/michaeltouati/ESVM).
2) Add your code. Please, keep the code Fortran 95 standard compliant by : 
- respecting 2 spaces for indentation rather than tabs
- not using object oriented Fortran 2003 features
- not forgetting to deallocate arrays
- etc ...
3) If you've added a new feature that needs new simulation parameters in the [input-deck](https://github.com/michaeltouati/ESVM/blob/main/input-deck) and updated correspondingly the source file [input.f90](https://github.com/michaeltouati/ESVM/blob/main/sources/input.f90), please add their descriptions in the [input-deck](https://github.com/michaeltouati/ESVM/blob/main/input-deck) following the same style
```sh
##                                                                   ##
## T  = electron temperature in eV                                   ##
##                                                                   ##
#######################################################################
#
#T  1000.
#
```
4) If you've added a new feature, let's call it 'my-new-feature', create a directory 'my-new-feature' in 'test-cases/Tests' and add a typical input deck inside that uses your new feature with parameters such as :
- `#x_min 0.`, `#x_max 5.`, `#d_x 0.25`, `#vx_min -5.`, `#vx_max 5.`, `#d_vx 0.1` for the discretized phase-space parameters and
- `#cfl 9.e-1`, `#L_t 5.`, `#dt_diag 0.25` for the simulation time parameters
(for the test to be fast) and a text file entitled output containing the terminal output generated by the corresponding simulation run using your new feature. You can type on the terminal :
```sh
./esvm > output
```
to generate such a file. Finally, add a section in the makefile such as
```sh
test_my-new-feature :
	@echo -n 'my-new-feature         : '
	@cp test-cases/Tests/my-new-feature/input-deck .
	@./esvm > test.output
	@tail -n +0 test.output | tail -r | tail -n +4 | tail -r > file1
	@tail -n +0 test-cases/Tests/my-new-feature/output | tail -r | tail -n +4 | tail -r > file2
	@diff file1 file2; \
	TST=$$?; \
	rm file1; rm file2; \
	if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \
```
that should be written if a new solver, a new boundary condition or ... is added. Finally, add your test at thye end of the makefile 
```sh
test :  test_start test_ampere test_poisson test_openMP test_periodic test_absorbing \
	test_donor_cell test_Lax_Wendroff test_Beam_Warming test_Fromm \
	test_minmod test_superbee test_Van_Leer test_MUSCL1 test_MUSCL2 \
	test_new-features_start test_my-new-feature test_end
```

5) Please, mention your contribution following the same style at the beginning of :
- the modified Fortran source file
```
!! Initial commit written by Dr Michaël J TOUATI - Dec. 2015         !!
!! my-new-feature or my-bug-fix commit by my-name - my-add-date.     !!
```
- the modified Python 3 script
```
# Initial commit written by Michaël J Touati
# my-new-feature or my-bug-fix commit by my-name - my-add-date
```
- input deck if 3.
```sh
##   Initial commit written by Dr Michaël J TOUATI - Dec. 2015       ##
##   my-new-feature commit by my-name - my-add-date                  ##
```
- makefile if 4.
```sh
## Initial commit written by Dr Michaël J TOUATI - Dec. 2015
## my-new-feature commit by my-name - my-add-date
```
- or at the beginning of a new source file
```
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!             ElectroStatic Vlasov-Maxwell (ESVM) code              !!
!!                                                                   !!
!! Initial commit written by commit by my-name - my-add-date         !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
```
6) Make sure your code passes all the tests by typing on your terminal:
```sh
make test
```
7) Ensure the code compiles with the following compiler debugging options:
- `-g -traceback -fopenmp -r8 -std95 -fpe0 -debug all -debug-parameters all -C` for the INTEL compiler ifort
- `-fdefault-real-8 -O -g -fopenmp -Wall -fcheck=all -fbacktrace -std=f95 -fall-intrinsics -ffpe-trap=invalid,zero,overflow` for the GNU compiler gfortran
8) If 3. and 4., update the tests generator bash script 'test-cases/Tests/test-gen.sh' by adding a command line that corresponds to your new feature to be tested
```sh
cp test-cases/Tests/my-new-feature/input-deck . && ./esvm > output && cp output test-cases/Tests/my-new-feature/ 
```
10) Issue that pull request!

# License
When you submit code changes, your submissions are understood to be under the same [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) that covers [ESVM](https://github.com/michaeltouati/ESVM). 
