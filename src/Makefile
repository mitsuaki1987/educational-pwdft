include ../make.inc

include make.depend

OBJS=\
atm_spec.o \
k_point.o \
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
