#######################################################################
##                                                                   ##
##             ElectroStatic Vlasov-Maxwell (ESVM) code              ##
##                                                                   ##
## Copyright © 2015 Michaël J TOUATI                                 ##
##                                                                   ##
## This file is part of ESVM.                                        ##
##                                                                   ##
## ESVM is free software: you can redistribute it and/or modify      ##
## it under the terms of the GNU General Public License as published ##
## by the Free Software Foundation, either version 3 of the License, ##
## or (at your option) any later version.                            ##
##                                                                   ##
## ESVM is distributed in the hope that it will be useful,           ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of    ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     ##
## GNU General Public License for more details.                      ##
##                                                                   ##
## You should have received a copy of the GNU General Public License ##
## along with ESVM. If not, see <https://www.gnu.org/licenses/>.     ##
##                                                                   ##
#######################################################################
## Initial commit written by Michaël J TOUATI - Dec. 2015

######################
######################
##  INTEL Compiler  ##
######################
######################

# F90 = ifort

##########
# openMP #
##########

# OPTS = -fopenmp -r8 -O3

#########
# debug #
#########

# OPTS = -g -traceback -fopenmp -r8 -std95 -fpe0 -debug all -debug-parameters all -C 

####################
####################
##  GNU Compiler  ##
####################
####################

F90 = gfortran

##########
# openMP #
##########

OPTS = -fopenmp -fdefault-real-8 -O3

#########
# debug #
#########

# OPTS = -fdefault-real-8 -O -g -fopenmp -Wall -fcheck=all -fbacktrace -std=f95 -fall-intrinsics -ffpe-trap=invalid,zero,overflow

#####################################
#####################################
##  Linkage and compilation rules  ##
#####################################
#####################################

SRC_PATH    = sources/

SRC_PATH_PY = ${SRC_PATH}/plot/

SRCS        = acuracy.f90 constants.f90 physics.f90 input.f90 diagnostics.f90 library.f90 esvm.f90

SRCS_CHK    = acuracy.f90 constants.f90 physics.f90 input.f90 input-chk.f90

OBJTS      := $(SRCS:%.f90=%.o)

OBJTS_CHK  := $(SRCS_CHK:%.f90=%.o)

%.o : $(SRC_PATH)%.f90
	$(F90) $(OPTS) -c $(SRC_PATH)$*.f90

all : check-input-deck esvm    \
	  remove-compilation-files

esvm : $(OBJTS)
	$(F90) $(OPTS) $(OBJTS) -o esvm

check-input-deck : $(OBJTS_CHK)
	$(F90) $(OPTS) $(OBJTS_CHK) -o check-input-deck

remove-compilation-files :
	@rm *.mod
	@rm *.o

########################
########################
##  Check input-deck  ##
########################
########################

check :
	@./check-input-deck

###############
###############
##  Running  ##
###############
###############

run :
	@./esvm

################
################
##  Cleaning  ##
################
################

distclean :
	@rm -f *.o *.mod esvm check-input-deck
	@rm -rf sources/plot/__pycache__

figclean :
	@rm -rf figures

resclean :
	@rm -rf results

clean : distclean figclean resclean

################
################
##  Plotting  ##
################
################

plot_energies : $(SRC_PATH_PY)plot_energies.py
	@python3 $(SRC_PATH_PY)plot_energies.py

plot_hydro2D : $(SRC_PATH_PY)plot_hydro2D.py
	@python3 $(SRC_PATH_PY)plot_hydro2D.py

plot_fe : $(SRC_PATH_PY)plot_fe.py
	@python3 $(SRC_PATH_PY)plot_fe.py

plot_hydro1D : $(SRC_PATH_PY)plot_hydro1D.py
	@python3 $(SRC_PATH_PY)plot_hydro1D.py

plot_logfe : $(SRC_PATH_PY)plot_logfe.py
	@python3 $(SRC_PATH_PY)plot_logfe.py

plot :  plot_energies \
		plot_hydro2D  \
		plot_fe       \
		plot_hydro1D

###############
###############
##  Testing  ##
###############
###############

TEST_DIR  = Maxwell/Ampere
TEST_DIR += Maxwell/Poisson
TEST_DIR += Parallelization/OpenMP
TEST_DIR += Boundary-conditions/Periodic
TEST_DIR += Boundary-conditions/Absorbing
TEST_DIR += Vlasov-linear/Donor-cell
TEST_DIR += Vlasov-linear/Lax-Wendroff
TEST_DIR += Vlasov-linear/Beam-Warming
TEST_DIR += Vlasov-linear/Fromm
TEST_DIR += Vlasov-nonlinear/Minmod
TEST_DIR += Vlasov-nonlinear/Superbee
TEST_DIR += Vlasov-nonlinear/Van-Leer
TEST_DIR += Vlasov-nonlinear/MUSCL1
TEST_DIR += Vlasov-nonlinear/MUSCL2
TEST_DIR += Academic-cases/Landau
TEST_DIR += Academic-cases/Wakefield
TEST_DIR += Academic-cases/Two-stream-insta
# TEST_DIR += New-features/Example

RED  =$(shell tput setaf 1)
GREEN=$(shell tput setaf 2)
RESET=$(shell tput sgr0)
TESTS:=$(sort ${TEST_DIR})
test :
	@mv input-deck input-deck-old
	@echo '---------------------------------'
	@echo '        TESTS DESCRIPTION        '
	@echo '---------------------------------'
	@echo '                                 '
	@echo 'The tests consist in performing  '
	@echo ' diff file1 file2 where :        '
	@echo ' * file1 is the test simulation  '
	@echo '   terminal output               '
	@echo ' * file2 the terminal output of  '
	@echo '   the corresponding simulation  '
	@echo '   performed by the developper   '
	@echo '   located in test-cases/Tests/  '
	@echo '                                 '
	@for tst in ${TESTS}; do \
		echo '---------------------------------' ; \
		echo $${tst}' :'  ; \
		cp test-cases/Tests/$${tst}/input-deck . ; \
		./esvm > test.output ; \
		if hash tac 2>/dev/null; then \
			tail -n +0 test.output | head -n -3 > file1 ; \
			tail -n +0 test-cases/Tests/$${tst}/output | head -n -3 > file2 ; \
		else \
			tail -n +0 test.output | tail -r | tail -n +4 | tail -r > file1 ; \
			tail -n +0 test-cases/Tests/$${tst}/output | tail -r | tail -n +4 | tail -r > file2 ; \
		fi ; \
		diff file1 file2; \
		TST=$$?; \
		rm file1; rm file2; \
		if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; \
	done
	@rm -f  test.output
	@rm -rf results/Academic-cases
	@rm -rf results/Boundary-conditions
	@rm -rf results/Maxwell
	@rm -rf results/Parallelization
	@rm -rf results/Vlasov-linear
	@rm -rf results/Vlasov-nonlinear
	@rm -rf results/New-features
	@mv input-deck-old input-deck
