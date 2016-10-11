PROG =	NMR

SRCS =	PDB_coords.f90 main_dist.f90

OBJS =	PDB_coords.o main_dist.o

LIBS =	

CC = cc
CFLAGS = -O
CCC = c++
CCFLAGS = -O
FC = f77
FFLAGS = -O
F90 = gfortran
F90FLAGS = -O
LDFLAGS = -s
M2C = \#

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90 .F90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.F90.o:
	$(F90) $(F90FLAGS) -c $<

%.o: %.mod

coords.o:
main.o: 
