COMPILER ?= gfortran

ifeq ($(COMPILER), gfortran)
	FC 	:= gfortran
	FFLAGS	:= -Wall -O2 -std=f2008 -fopenmp -frealloc-lhs
	ifeq ($(DEBUG), 1)
		FFLAGS += -fcheck=all -g -pg -ffpe-trap=invalid,zero,overflow \
		-pedantic -finit-real=nan
	endif

else ifeq ($(COMPILER), ifort)
	FC 	:= ifort
	FFLAGS	:= -warn-all -O2 -stand f2008 -openmp -assume realloc-lhs
	ifeq ($(DEBUG), 1)
		FFLAGS += -check all -g -p -debug all
	endif
endif

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))