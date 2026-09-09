# Makefile for gvpt2.x, the Fortran translation of step2/gvpt2.py
#
#   make                 build with gfortran (LAPACK from Accelerate on macOS,
#                        -llapack -lblas elsewhere)
#   make FC=ifx LIBS="-qmkl=sequential"     Intel compiler + MKL
#   make clean
#
# Only LAPACK's dsyev is needed (polyad diagonalisation).

ifeq ($(origin FC),default)
  FC = gfortran
endif
FFLAGS ?= -O2

UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
  LIBS ?= -framework Accelerate
else
  LIBS ?= -llapack -lblas
endif

OBJ := gvpt2_constants.o gvpt2_io.o gvpt2_polyad.o gvpt2_core.o gvpt2.o

.PHONY: all clean

all: gvpt2.x

gvpt2.x: $(OBJ)
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS)

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

# module dependencies
gvpt2_io.o:     gvpt2_constants.o
gvpt2_polyad.o: gvpt2_constants.o gvpt2_io.o
gvpt2_core.o:   gvpt2_constants.o gvpt2_polyad.o
gvpt2.o:        gvpt2_constants.o gvpt2_io.o gvpt2_polyad.o gvpt2_core.o

clean:
	rm -f *.o *.mod gvpt2.x
