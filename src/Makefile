include ../make.inc

OBJS=\
atm_spec.o \
diag_direct.o \
gvec.o \
hamiltonian.o \
io_vloc.o \
k_point.o \
lobpcg.o \
pwdft.o \
solver.o \
stdin.o

.SUFFIXES :
.SUFFIXES : .o .F90

all:pwdft.x

.F90.o:
	$(FC) $(FFLAGS) -c $<

pwdft.x:$(OBJS)
	$(FC) $(OBJS) $(FFLAGS) $(LDFLAGS) -o $@

clean:
	rm -rf *.o *.mod *.x

include make.depend
