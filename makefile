## Initial commit written by Dr Michaël J TOUATI - Dec. 2015
##
#############
#############
##  ifort  ##
#############
#############

# F90 = ifort

##########
# openMP #
##########

# OPTS = -fopenmp -r8 -O3

#########
# debug #
#########

# OPTS = -g -traceback -fopenmp -r8 -std95 -fpe0 -debug all -debug-parameters all -C 

################
################
##  gfortran  ##
################
################

F90 = gfortran

##########
# openMP #
##########

# OPTS = -fopenmp -fdefault-real-8 -O3

#########
# debug #
#########

OPTS = -fdefault-real-8 -O -g -fopenmp -Wall -fcheck=all -fbacktrace -std=f95 -fall-intrinsics -ffpe-trap=invalid,zero,overflow

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
##  CHECK INPUT DECK  ##
########################
########################

check-input-deck : $(OBJTS_CHK)
	$(F90) $(OPTS) $(OBJTS_CHK) -o check-input-deck
	@rm *.mod
	@rm *.o

check : check-input-deck
	./check-input-deck
	@rm check-input-deck

###############
###############
##  RUNNING  ##
###############
###############

run :
	./esvm

################
################
##  CLEANING  ##
################
################

clean :
	rm -rf *.o *.mod esvm

clean_figures :
	rm -rf figures

clean_results :
	rm -rf results

clean_all :
	rm -rf *.o *.mod esvm results figures

################
################
##  PLOTTING  ##
################
################

extract : $(SRC_PATH_PY)extract.py
	@python3 $(SRC_PATH_PY)extract.py

extract_hydro1D : $(SRC_PATH_PY)extract_hydro1D.py
	@python3 $(SRC_PATH_PY)extract_hydro1D.py

extract_energy : $(SRC_PATH_PY)extract_energy.py
	@python3 $(SRC_PATH_PY)extract_energy.py

extract_hydro2D : $(SRC_PATH_PY)extract_hydro2D.py
	@python3 $(SRC_PATH_PY)extract_hydro2D.py

extract_fe : $(SRC_PATH_PY)extract_fe.py
	@python3 $(SRC_PATH_PY)extract_fe.py

extract_logfe : $(SRC_PATH_PY)extract_logfe.py
	@python3 $(SRC_PATH_PY)extract_logfe.py

#############
#############
##  TESTS  ##
#############
#############

RED=$(shell tput setaf 1)
GREEN=$(shell tput setaf 2)
RESET=$(shell tput sgr0)

test_start : 
	@mv input-deck input-deck-old
	@echo '--------------------------------'
	@echo '        TESTS DESCRIPTION       '
	@echo '--------------------------------'
	@echo '                                '
	@echo 'The tests consist in performing '
	@echo ' diff file1 file2 where :       '
	@echo ' * file1 is the test simulation '
	@echo '   terminal output              '
	@echo ' * file2 the terminal output of '
	@echo '   the corresponding simulation '
	@echo '   performed by the developper  '
	@echo '   located in test-cases/Tests/ '
	@echo ' concerning a drifting plasma   '
	@echo ' simulation at Maxwell-Boltzmann'
	@echo ' equilibrium                    ' 
	@echo '                                '

test_ampere :
	@echo '--------------------------------'
	@echo '         Maxwell solver         '
	@echo '--------------------------------'
	@echo '                                '
	@echo -n 'Maxwell-Ampere solver  : '
	@cp test-cases/Tests/Ampere/input-deck .
	@./esvm > test.output
	@tail -n +0 test.output | tail -r | tail -n +4 | tail -r > file1
	@tail -n +0 test-cases/Tests/Ampere/output | tail -r | tail -n +4 | tail -r > file2
	@diff file1 file2; \
    TST=$$?; \
    rm file1; rm file2; \
    if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_poisson :
	@echo -n 'Maxwell-Poisson solver : '
	@cp test-cases/Tests/Poisson/input-deck .
	@./esvm > test.output
	@tail -n +0 test.output | tail -r | tail -n +4 | tail -r > file1
	@tail -n +0 test-cases/Tests/Poisson/output | tail -r | tail -n +4 | tail -r > file2
	@diff file1 file2; \
    TST=$$?; \
    rm file1; rm file2; \
    if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_openMP :
	@echo '--------------------------------'
	@echo '            OpenMP              '
	@echo '--------------------------------'
	@echo '                                '
	@echo -n 'OpenMP parallelization : '
	@cp test-cases/Tests/OpenMP/input-deck .
	@./esvm > test.output
	@tail -n +0 test.output | tail -r | tail -n +4 | tail -r > file1
	@tail -n +0 test-cases/Tests/OpenMP/output | tail -r | tail -n +4 | tail -r > file2
	@diff file1 file2; \
    TST=$$?; \
    rm file1; rm file2; \
    if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_periodic :
	@echo '--------------------------------'
	@echo '      Boundary conditions       '
	@echo '--------------------------------'
	@echo '                                '
	@echo -n 'Periodic bound. cond.  : '
	@cp test-cases/Tests/Boundary-conditions/Periodic/input-deck .
	@./esvm > test.output
	@tail -n +0 test.output | tail -r | tail -n +4 | tail -r > file1
	@tail -n +0 test-cases/Tests/Boundary-conditions/Periodic/output | tail -r | tail -n +4 | tail -r > file2
	@diff file1 file2; \
    TST=$$?; \
    rm file1; rm file2; \
    if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_absorbing :
	@echo -n 'Absorbing bound. cond. : '
	@cp test-cases/Tests/Boundary-conditions/Absorbing/input-deck .
	@./esvm > test.output
	@tail -n +0 test.output | tail -r | tail -n +4 | tail -r > file1
	@tail -n +0 test-cases/Tests/Boundary-conditions/Absorbing/output | tail -r | tail -n +4 | tail -r > file2
	@diff file1 file2; \
    TST=$$?; \
    rm file1; rm file2; \
    if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_donor_cell :
	@echo '--------------------------------'
	@echo '    Linear advection solvers    '
	@echo '--------------------------------'
	@echo '                                '
	@echo -n 'Donor cell solver      : '
	@cp test-cases/Tests/Linear-advection-schemes/Donor-cell/input-deck .
	@./esvm > test.output
	@tail -n +0 test.output | tail -r | tail -n +4 | tail -r > file1
	@tail -n +0 test-cases/Tests/Linear-advection-schemes/Donor-cell/output | tail -r | tail -n +4 | tail -r > file2
	@diff file1 file2; \
    TST=$$?; \
    rm file1; rm file2; \
    if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_Lax_Wendroff :
	@echo -n 'Lax-Wendroff solver    : '
	@cp test-cases/Tests/Linear-advection-schemes/Lax-Wendroff/input-deck .
	@./esvm > test.output
	@tail -n +0 test.output | tail -r | tail -n +4 | tail -r > file1
	@tail -n +0 test-cases/Tests/Linear-advection-schemes/Lax-Wendroff/output | tail -r | tail -n +4 | tail -r > file2
	@diff file1 file2; \
    TST=$$?; \
    rm file1; rm file2; \
    if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_Beam_Warming :
	@echo -n 'Beam-Warming solver    : '
	@cp test-cases/Tests/Linear-advection-schemes/Beam-Warming/input-deck .
	@./esvm > test.output
	@tail -n +0 test.output | tail -r | tail -n +4 | tail -r > file1
	@tail -n +0 test-cases/Tests/Linear-advection-schemes/Beam-Warming/output | tail -r | tail -n +4 | tail -r > file2
	@diff file1 file2; \
    TST=$$?; \
    rm file1; rm file2; \
	if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_Fromm :
	@echo -n 'Fromm solver           : '
	@cp test-cases/Tests/Linear-advection-schemes/Fromm/input-deck .
	@./esvm > test.output
	@tail -n +0 test.output | tail -r | tail -n +4 | tail -r > file1
	@tail -n +0 test-cases/Tests/Linear-advection-schemes/Fromm/output | tail -r | tail -n +4 | tail -r > file2
	@diff file1 file2; \
    TST=$$?; \
    rm file1; rm file2; \
	if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_minmod :
	@echo '--------------------------------'
	@echo '  Non-linear advection solvers  '
	@echo '--------------------------------'
	@echo '                                '
	@echo -n 'Minmod solver          : '
	@cp test-cases/Tests/Non-linear-advection-schemes/Minmod/input-deck .
	@./esvm > test.output
	@tail -n +0 test.output | tail -r | tail -n +4 | tail -r > file1
	@tail -n +0 test-cases/Tests/Non-linear-advection-schemes/Minmod/output | tail -r | tail -n +4 | tail -r > file2
	@diff file1 file2; \
    TST=$$?; \
    rm file1; rm file2; \
	if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_superbee :
	@echo -n 'Superbee solver        : '
	@cp test-cases/Tests/Non-linear-advection-schemes/Superbee/input-deck .
	@./esvm > test.output
	@tail -n +0 test.output | tail -r | tail -n +4 | tail -r > file1
	@tail -n +0 test-cases/Tests/Non-linear-advection-schemes/Superbee/output | tail -r | tail -n +4 | tail -r > file2
	@diff file1 file2; \
    TST=$$?; \
    rm file1; rm file2; \
	if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_Van_Leer :
	@echo -n 'Van Leer (b=1.5) solver: '
	@cp test-cases/Tests/Non-linear-advection-schemes/Van-Leer/input-deck .
	@./esvm > test.output
	@tail -n +0 test.output | tail -r | tail -n +4 | tail -r > file1
	@tail -n +0 test-cases/Tests/Non-linear-advection-schemes/Van-Leer/output | tail -r | tail -n +4 | tail -r > file2
	@diff file1 file2; \
    TST=$$?; \
    rm file1; rm file2; \
	if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_MUSCL1 :
	@echo -n 'MUSCL 1 solver         : '
	@cp test-cases/Tests/Non-linear-advection-schemes/MUSCL1/input-deck .
	@./esvm > test.output
	@tail -n +0 test.output | tail -r | tail -n +4 | tail -r > file1
	@tail -n +0 test-cases/Tests/Non-linear-advection-schemes/MUSCL1/output | tail -r | tail -n +4 | tail -r > file2
	@diff file1 file2; \
    TST=$$?; \
    rm file1; rm file2; \
	if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_MUSCL2 :
	@echo -n 'MUSCL 2 solver         : '
	@cp test-cases/Tests/Non-linear-advection-schemes/MUSCL2/input-deck .
	@./esvm > test.output
	@tail -n +0 test.output | tail -r | tail -n +4 | tail -r > file1
	@tail -n +0 test-cases/Tests/Non-linear-advection-schemes/MUSCL2/output | tail -r | tail -n +4 | tail -r > file2
	@diff file1 file2; \
    TST=$$?; \
    rm file1; rm file2; \
	if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_new-features_start :
	@echo '--------------------------------'
	@echo '          New Feature           '
	@echo '--------------------------------'

test_end :
	@rm -f test.output
	@rm -rf results/
	@mv input-deck-old input-deck

test :  test_start test_ampere test_poisson test_openMP test_periodic test_absorbing \
		test_donor_cell test_Lax_Wendroff test_Beam_Warming test_Fromm \
		test_minmod test_superbee test_Van_Leer test_MUSCL1 test_MUSCL2 \
		test_new-features_start test_end
