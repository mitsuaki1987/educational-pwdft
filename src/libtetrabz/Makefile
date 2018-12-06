include ../../make.inc

OBJS = \
libtetrabz.o \
libtetrabz_dbldelta_mod.o \
libtetrabz_dblstep_mod.o \
libtetrabz_dos_mod.o \
libtetrabz_fermigr_mod.o \
libtetrabz_occ_mod.o \
libtetrabz_polcmplx_mod.o \
libtetrabz_polstat_mod.o \
libtetrabz_common.o

.SUFFIXES :
.SUFFIXES : .o .F90

libtetrabz.a:$(OBJS)
	ar cr $@ $(OBJS)

.F90.o:
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.a *.o *.mod

libtetrabz.o:libtetrabz_polcmplx_mod.o
libtetrabz.o:libtetrabz_fermigr_mod.o
libtetrabz.o:libtetrabz_polstat_mod.o
libtetrabz.o:libtetrabz_dbldelta_mod.o
libtetrabz.o:libtetrabz_dblstep_mod.o
libtetrabz.o:libtetrabz_dos_mod.o
libtetrabz.o:libtetrabz_occ_mod.o
libtetrabz_dbldelta_mod.o:libtetrabz_common.o
libtetrabz_dblstep_mod.o:libtetrabz_common.o
libtetrabz_dos_mod.o:libtetrabz_common.o
libtetrabz_fermigr_mod.o:libtetrabz_common.o
libtetrabz_occ_mod.o:libtetrabz_common.o
libtetrabz_polcmplx_mod.o:libtetrabz_common.o
libtetrabz_polstat_mod.o:libtetrabz_common.o

