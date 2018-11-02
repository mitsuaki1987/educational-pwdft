include ../make.inc

include make.depend

all:pwdft.x

pwdft.x:
	$(FC) $(OBJS) $(FFLAGS) $(LDFLAGS) -o $@

clean:
	rm -rf *.o *.mod *.x
