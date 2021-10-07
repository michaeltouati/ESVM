# Contributing to ESVM

First of all, thank you very much for your interest in contributing to [ESVM](https://github.com/michaeltouati/ESVM)! 

Hans Bethe once said : “We need science education to produce scientists, but we need it equally to create literacy in the public. Man has a fundamental urge to comprehend the world about him, and science gives today the only world picture which we can consider as valid. It gives an understanding of the inside of the atom and of the whole universe, or the peculiar properties of the chemical substances and of the manner in which genes duplicate in biology. An educated layman can, of course, not contribute to science, but can enjoy and participate in many scientific discoveries which as constantly made. Such participation was quite common in the 19th century, but has unhappily declined. Literacy in science will enrich a person's life.”

But Richard P. Feynmann also added : “Physics is like sex: sure, it may give some practical results, but that's not why we do it.”

Being fully in agreement with Richard P. feynman while firmly believing in the recommendations of Hans Bethe, it seemed to me natural to develop ESVM with Github in order to try to understand the complexity of plasmas, numerical analysis and High Parallel Computing (HPC) and make accessible this small piece of Modern Physics knowledge to the largest public as possible at the same time.

# Report a bug you faced when compiling or using 

For problems related with plotting, please make sure :
1) LaTeX, dvipng and ghostscript are each working and on your PATH for the Matplotlib Python package to be able to render tex fonts
2) you're using the main repository branch [ESVM](https://github.com/michaeltouati/ESVM) by typing on your terminal from your local ESVM
```sh
git pull
```
If the bug persists or if it is related to another problem, follow these steps :
1) Go to 'Issues' on the ESVM main repository branch [issues](https://github.com/michaeltouati/ESVM/issues)
2) Click on 'New issue'
4) Describe the bug the more clear and concise as possible in the title starting with "Bug :"
5) Describe with the more details as possible the bug providing the more pieces of information as possible (input-deck, terminal output and/or screeshot) and more importantly, specify your environment :

OS: [e.g. Ubuntu 20.04]

Fortran compiler [e.g. gfortran 11.2.0]

Python version [e.g. Python 3.7.11]

Matplotlib Python package version [e.g., Matplotlib 3.4.3]

Numpy Python package version [e.g., Numpy 1.21.2]

# Propose a new feature

If you want to propose a new feature for ESVM, follow these steps :
1) Go to 'Issues' on the ESVM main repository branch [issues](https://github.com/michaeltouati/ESVM/issues)
2) Click on 'New issue'
4) Describe the new feature request the more clear and concise as possible in the title starting with "Feature request :"
5) Describe with the more details as possible the new feature request and
6) Click on 'Submit new issue'

We'll discuss about it.

# Fix a bug or submit a new feature

ESVM Uses the [Fork and pull model](https://docs.github.com/en/github/collaborating-with-pull-requests/getting-started/about-collaborative-development-models).
In order to fix a bug or submit a new feature in ESVM, follow these steps:

1. Fork the repo and create your own branch from [ESVM 'main' branch](https://github.com/michaeltouati/ESVM).
2. Add your code.

Please, try to keep the code Fortran 90 standard compliant by : 
- respecting 2 spaces for indentation rather than tabs
- not using object oriented Fortran 2003 features
- etc ...

3. If you've added parameters in the input-deck and updated the source file input.f90, please add their descriptions 
```javascript
##                                                                   ##
## T  = electron temperature in eV                                   ##
##                                                                   ##
```
in the input-deck following the same style:
```javascript
#
#T  1000.
#
```
4. Create a directory "my-new-feature" test-cases/Tests/"my-new-feature" by adding :
  1. a typical input deck that uses your new feature with:
- `#N_th 1` OpenMP threads, 
- discretized phase-space parameters `#x_min 0.`, `#x_max 5.`, `#d_x 0.25`, `#vx_min -5.`, `#vx_max 5.`, `#d_vx 0.1`
- simulation time parameters `#cfl 9.e-1`, `#L_t 5.`, `#dt_diag 0.25` and
- test case parameters `#perturb  0` and `#vd 1.`
for the test to be fast  
  2. a section in the makefile such as
```javascript
test_absorbing :
	@echo -n 'Absorbing bound. cond. : '
	@cp test-cases/Tests/Boundary-conditions/Absorbing/input-deck .
	@./esvm > test.output
	@diff test.output test-cases/Tests/Boundary-conditions/Absorbing/output; \
	TST=$$?;\
	if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \
```
that has been written to check 
if a new solver, a new boundary condition or... is added and make sure it passes the tests by typing on your terminal from your local branch cloned or downloaded from your forked branch:
```sh
make test
```
5. Ensure the code compiles with all compiler debugging options:
- `-g -traceback -fopenmp -r8 -std90 -fpe0 -debug all -debug-parameters all -C ` for the INTEL compiler ifort
- `-fdefault-real-8 -O -g -fopenmp -Wall -fcheck=all -fbacktrace -std=f95 -fall-intrinsics -ffpe-trap=invalid,zero,overflow` for the GNU compiler gfortran

6. Issue that pull request!

# License
When you submit code changes, your submissions are understood to be under the same [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html) that covers ESVM. 
