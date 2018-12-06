include ../make.inc

OBJS=\
atm_spec.o \
constant.o \
diag_direct.o \
fftw_wrapper.o \
gvec.o \
hamiltonian.o \
hartree.o \
io_vloc.o \
k_point.o \
lobpcg.o \
pp.o \
pwdft.o \
scf.o \
solver.o \
stdin.o

.SUFFIXES :
.SUFFIXES : .o .F90

all:pwdft.x

.F90.o:
	$(FC) $(FFLAGS) -I./libtetrabz/ -c $<

pwdft.x:$(OBJS) tetra
	$(FC) $(OBJS) $(FFLAGS) $(LDFLAGS) -L./libtetrabz/ -ltetrabz -o $@

tetra:
	cd libtetrabz; make

clean:
	rm -rf *.o *.mod *.x
	cd libtetrabz; make clean

include make.depend
