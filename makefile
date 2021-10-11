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

SRC_PATH_PY = ${SRC_PATH}/extract/

SRCS        = acuracy.f90 constants.f90 physics.f90 input.f90 diagnostics.f90 library.f90 esvm.f90

SRCS_CHK    = acuracy.f90 constants.f90 physics.f90 input.f90 input-chk.f90

OBJTS := $(SRCS:%.f90=%.o)

OBJTS_CHK := acuracy.o constants.o physics.o input.o input-chk.o

all : ESVM

ESVM : $(OBJTS)
	$(F90) $(OPTS) $(OBJTS) -o esvm
	rm *.mod
	rm *.o

acuracy.o : $(SRC_PATH)acuracy.f90
	$(F90) $(OPTS) -c $(SRC_PATH)acuracy.f90
	
constants.o : $(SRC_PATH)constants.f90
	$(F90) $(OPTS) -c $(SRC_PATH)constants.f90

physics.o : $(SRC_PATH)physics.f90
	$(F90) $(OPTS) -c $(SRC_PATH)physics.f90

input.o : $(SRC_PATH)input.f90
	$(F90) $(OPTS) -c $(SRC_PATH)input.f90

input-chk.o : $(SRC_PATH)input-chk.f90
	$(F90) $(OPTS) -c $(SRC_PATH)input-chk.f90

diagnostics.o : $(SRC_PATH)diagnostics.f90
	$(F90) $(OPTS) -c $(SRC_PATH)diagnostics.f90

library.o : $(SRC_PATH)library.f90
	$(F90) $(OPTS) -c $(SRC_PATH)library.f90
	
esvm.o : $(SRC_PATH)esvm.f90
	$(F90) $(OPTS) -c $(SRC_PATH)esvm.f90

########################
########################
##  Check input-deck  ##
########################
########################

check-input-deck : $(OBJTS_CHK)
	$(F90) $(OPTS) $(OBJTS_CHK) -o check-input-deck
	@rm *.mod
	@rm *.o

check : check-input-deck
	@./check-input-deck
	@rm check-input-deck

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
	@rm -rf sources/extract/__pycache__ 
	@rm -f *.o *.mod esvm

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

extract_energy : $(SRC_PATH_PY)extract_energy.py
	@python3 $(SRC_PATH_PY)extract_energy.py

extract_hydro2D : $(SRC_PATH_PY)extract_hydro2D.py
	@python3 $(SRC_PATH_PY)extract_hydro2D.py

extract_fe : $(SRC_PATH_PY)extract_fe.py
	@python3 $(SRC_PATH_PY)extract_fe.py

extract_hydro1D : $(SRC_PATH_PY)extract_hydro1D.py
	@python3 $(SRC_PATH_PY)extract_hydro1D.py

extract_logfe : $(SRC_PATH_PY)extract_logfe.py
	@python3 $(SRC_PATH_PY)extract_logfe.py

extract :   extract_energy   \
			extract_hydro2D  \
			extract_fe       \
			extract_hydro1D

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
TEST_DIR += Academic-case/Landau
TEST_DIR += Academic-case/Wakefield
TEST_DIR += Academic-case/Two-stream-insta
# TEST_DIR += New-feature

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
		tail -n +0 test.output | tail -r | tail -n +4 | tail -r > file1 ; \
		tail -n +0 test-cases/Tests/$${tst}/output | tail -r | tail -n +4 | tail -r > file2 ; \
		diff file1 file2; \
		TST=$$?; \
		rm file1; rm file2; \
		if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; \
	done
	@rm -f test.output
	@rm -rf results/
	@mv input-deck-old input-deck
