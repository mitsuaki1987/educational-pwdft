include ../make.inc

OBJS=\
atm_spec.o \
constant.o \
diag_direct.o \
energy.o \
fftw_wrapper.o \
griddata.o \
gvec.o \
hamiltonian.o \
k_point.o \
kohn_sham.o \
lobpcg.o \
plot.o \
pp.o \
pwdft.o \
rho_v.o \
scf.o \
stdin.o

.SUFFIXES :
.SUFFIXES : .o .F90

all:tetra pwdft.x

.F90.o:
	$(FC) $(FFLAGS) -I./libtetrabz/ -c $<

pwdft.x:$(OBJS) libtetrabz/libtetrabz.a
	$(FC) $(OBJS) $(FFLAGS) $(LDFLAGS) -L./libtetrabz/ -ltetrabz -o $@

tetra:
	cd libtetrabz; make

clean:
	rm -rf *.o *.mod *.x
	cd libtetrabz; make clean

include make.depend
