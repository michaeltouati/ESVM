F90=nvfortran
OPTS=-mp -acc -fast -ta=tesla -Minfo=accel

SRCS = acuracy.f90 constants.f90 physics.f90 input.f90 diagnostics.f90 library.f90 esvm.f90

OBJTS := $(SRCS:%.f90=%.o)

all : ESVM

ESVM : $(OBJTS)
	$(F90) $(OPTS) $(OBJTS) -o esvm
	rm *.mod
	rm *.o

acuracy.o : acuracy.f90
	$(F90) $(OPTS) -c acuracy.f90
	
constants.o : constants.f90
	$(F90) $(OPTS) -c constants.f90

physics.o : physics.f90
	$(F90) $(OPTS) -c physics.f90

input.o : input.f90
	$(F90) $(OPTS) -c input.f90

diagnostics.o : diagnostics.f90
	$(F90) $(OPTS) -c diagnostics.f90

library.o : library.f90
	$(F90) $(OPTS) -c library.f90
	
esvm.o : esvm.f90
	$(F90) $(OPTS) -c esvm.f90

clean :
	rm -rf *.o *.mod esvm

