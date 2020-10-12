#############
#############
##  ifort  ##
#############
#############

#F90 = ifort

##########
# openMP #
##########

#OPTS = -prec-div -prec-sqrt -openmp -openmp-report2 -r8

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

OPTS = -fopenmp -O3 -ffixed-line-length-none  -fdefault-real-8

#########
# debug #
#########

#OPTS = -g -fopenmp -fbacktrace -ffpe-trap=zero,overflow,underflow

#####################################
#####################################
##  Linkage and compilation rules  ##
#####################################
#####################################

SRC_PATH    = sources/

SRC_PATH_PY = ${SRC_PATH}

SRCS        = acuracy.f90 constants.f90 physics.f90 input.f90 diagnostics.f90 library.f90 esvm.f90

OBJTS := $(SRCS:%.f90=%.o)

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

diagnostics.o : $(SRC_PATH)diagnostics.f90
	$(F90) $(OPTS) -c $(SRC_PATH)diagnostics.f90

library.o : $(SRC_PATH)library.f90
	$(F90) $(OPTS) -c $(SRC_PATH)library.f90
	
esvm.o : $(SRC_PATH)esvm.f90
	$(F90) $(OPTS) -c $(SRC_PATH)esvm.f90

clean_compilation :
	rm -rf *.o *.mod esvm

clean_figures :
	rm -rf figures

clean_results :
	rm -rf results

clean_all :
	rm -rf *.o *.mod esvm results figures

extract : $(SRC_PATH_PY)extract.py
	python3 $(SRC_PATH_PY)extract.py

extract_logfe : $(SRC_PATH_PY)extract_f_log.py
	python3 $(SRC_PATH_PY)extract_f_log.py
