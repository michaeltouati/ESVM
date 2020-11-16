# Copyright (c) 2016, NVIDIA CORPORATION. All rights reserved.
COMPILER ?= NVHPC-tesla
FC = nvfortran
CFLAGS = -DUSE_DOUBLE
ifeq ($(COMPILER),GCC)
	CFLAGS += -march=native -O3 -lm
else ifeq ($(COMPILER),NVHPC)
	CFLAGS += -fast
else ifeq ($(COMPILER),NVHPC-tesla)
	CFLAGS += -Minfo=accel -O3 -acc=gpu -gpu=managed
else ifeq ($(COMPILER),NVHPC-manual)
	CFLAGS += -Minfo=accel -O3 -acc=gpu
else ifeq ($(COMPILER),NVHPC-multicore)
	CFLAGS += -Minfo=accel -O3 -acc=multicore
endif

JSC_SUBMIT_CMD ?= srun --partition=gpus --gres=gpu:1 

all: poisson2d

poisson2d: acuracy.f90 constants.f90 physics.f90 input.f90 diagnostics.f90 library.f90 esvm.f90 Makefile
	$(FC) $(CFLAGS) acuracy.f90 constants.f90 physics.f90 input.f90 diagnostics.f90 library.f90 esvm.f90 -o $@

.PHONY: clean run profile
clean:
	rm -f poisson2d poisson2d*.o

run: poisson2d
	$(JSC_SUBMIT_CMD) ./poisson2d

profile: poisson2d
	$(JSC_SUBMIT_CMD) nsys nvprof ./poisson2d 10
profile_tofile: poisson2d
	$(JSC_SUBMIT_CMD) nsys profile -o poisson2d.qdrep ./poisson2d 10

.PHONY: task1 task1_profile task4
task1: poisson2d.f90 poisson2d_reference.f90 Makefile
	gfortran -g -pg -DUSE_DOUBLE -lm acuracy.f90 constants.f90 physics.f90 input.f90 diagnostics.f90 library.f90 esvm.f90 -o poisson2d

task1_profile:
	-rm gmon.out
	$(JSC_SUBMIT_CMD) ./poisson2d 100
	gprof -p -l ./poisson2d gmon.out
