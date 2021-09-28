#############
#############
##  ifort  ##
#############
#############

#F90 = ifort

##########
# openMP #
##########

#OPTS = -openmp -r8 -i8 -O3

#########
# debug #
#########

#OPTS = -g -fpe0 -ftrapuv -debug all -debug-parameters all -C -traceback -openmp -openmp-report2 -r8

################
################
##  gfortran  ##
################
################

F90 = gfortran

##########
# openMP #
##########

OPTS = -fopenmp -fdefault-real-8 -fdefault-integer-8 -O3

#########
# debug #
#########

OPTS = -fdefault-integer-8 -fdefault-real-8 -O -g -fopenmp -Wall -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow

#####################################
#####################################
##  Linkage and compilation rules  ##
#####################################
#####################################

SRC_PATH    = sources/

SRC_PATH_PY = ${SRC_PATH}

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
	python3 $(SRC_PATH_PY)extract.py

extract_logfe : $(SRC_PATH_PY)extract_f_log.py
	python3 $(SRC_PATH_PY)extract_f_log.py

#############
#############
##  TESTS  ##
#############
#############

RED=$(shell tput setaf 1)
GREEN=$(shell tput setaf 2)
RESET=$(shell tput sgr0)

test_ampere :
	@echo '--------------------------------'
	@echo '            ${RED}WARNING${RESET}          '
	@echo '--------------------------------'
	@echo 'The code should be compiled with'
	@echo 'double floating point precisions'
	@echo '     for the tests to pass!     '
	@echo '--------------------------------'
	@echo '         Maxwell solver         '
	@echo '--------------------------------'
	@echo -n 'Maxwell-Ampere solver  : '
	@cp test-cases/Test/Ampere/input-deck .
	@./esvm > test.output
	@diff test.output test-cases/Test/Ampere/output; \
    TST=$$?;\
    if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_poisson :
	@echo -n 'Maxwell-Poisson solver : '
	@cp test-cases/Test/Poisson/input-deck .
	@./esvm > test.output
	@diff test.output test-cases/Test/Poisson/output; \
    TST=$$?;\
    if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_openMP :
	@echo '--------------------------------'
	@echo '            OpenMP              '
	@echo '--------------------------------'
	@echo -n 'OpenMP parallelization : '
	@cp test-cases/Test/OpenMP/input-deck .
	@./esvm > test.output
	@diff test.output test-cases/Test/OpenMP/output; \
    TST=$$?;\
    if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_donor_cell :
	@echo '--------------------------------'
	@echo '    Linear advection solvers    '
	@echo '--------------------------------'
	@echo -n 'Donor cell solver      : '
	@cp test-cases/Test/Linear-advection-schemes/Donor-cell/input-deck .
	@./esvm > test.output
	@diff test.output test-cases/Test/Linear-advection-schemes/Donor-cell/output; \
    TST=$$?;\
    if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_Lax_Wendroff :
	@echo -n 'Lax-Wendroff solver    : '
	@cp test-cases/Test/Linear-advection-schemes/Lax-Wendroff/input-deck .
	@./esvm > test.output
	@diff test.output test-cases/Test/Linear-advection-schemes/Lax-Wendroff/output; \
    TST=$$?;\
    if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_Beam_Warming :
	@echo -n 'Beam-Warming solver    : '
	@cp test-cases/Test/Linear-advection-schemes/Beam-Warming/input-deck .
	@./esvm > test.output
	@diff test.output test-cases/Test/Linear-advection-schemes/Beam-Warming/output; \
    TST=$$?;\
    if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_Fromm :
	@echo -n 'Fromm solver           : '
	@cp test-cases/Test/Linear-advection-schemes/Fromm/input-deck .
	@./esvm > test.output
	@diff test.output test-cases/Test/Linear-advection-schemes/Fromm/output; \
    TST=$$?;\
    if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_minmod :
	@echo '--------------------------------'
	@echo '  Non-linear advection solvers  '
	@echo '--------------------------------'
	@echo -n 'Minmod solver          : '
	@cp test-cases/Test/Non-linear-advection-schemes/Minmod/input-deck .
	@./esvm > test.output
	@diff test.output test-cases/Test/Non-linear-advection-schemes/Minmod/output; \
    TST=$$?;\
    if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_superbee :
	@echo -n 'Superbee solver        : '
	@cp test-cases/Test/Non-linear-advection-schemes/Superbee/input-deck .
	@./esvm > test.output
	@diff test.output test-cases/Test/Non-linear-advection-schemes/Superbee/output; \
    TST=$$?;\
    if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_Van_Leer :
	@echo -n 'Van Leer (b=1.5) solver: '
	@cp test-cases/Test/Non-linear-advection-schemes/Van-Leer/input-deck .
	@./esvm > test.output
	@diff test.output test-cases/Test/Non-linear-advection-schemes/Van-Leer/output; \
    TST=$$?;\
    if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_MUSCL1 :
	@echo -n 'MUSCL 1 solver         : '
	@cp test-cases/Test/Non-linear-advection-schemes/MUSCL1/input-deck .
	@./esvm > test.output
	@diff test.output test-cases/Test/Non-linear-advection-schemes/MUSCL1/output; \
    TST=$$?;\
    if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \

test_MUSCL2 :
	@echo -n 'MUSCL 2 solver         : '
	@cp test-cases/Test/Non-linear-advection-schemes/MUSCL2/input-deck .
	@./esvm > test.output
	@diff test.output test-cases/Test/Non-linear-advection-schemes/MUSCL2/output; \
    TST=$$?;\
    if [ $$TST -eq 0 ]; then echo "${GREEN}PASSED${RESET}"; else echo "${RED}NOT PASSED${RESET}"; fi; echo ' '; \
    rm -f test.output
	rm -rf results/

test : test_ampere test_poisson \
	   test_donor_cell test_Lax_Wendroff test_Beam_Warming test_Fromm \
	   test_minmod test_superbee test_Van_Leer test_MUSCL1 test_MUSCL2

