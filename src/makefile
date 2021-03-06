#
# Makefile for program '4hetddft-anisotropic'
# 3D dynamic (real-time evolution) 4-He Density Functional Theory for AN-isotropic PES's
#

COMP = ifort
CFLAGS = -c -O3 -xAVX -align array64byte -qopenmp\
		 -parallel -qopt-matmul -unroll0 -module ./modules
LD_FLAGS = -threads -I${MKLROOT}/include/fftw -mkl=parallel -qopt-matmul
PROGNAME = 4hetddft-anisotropic

all: $(PROGNAME)
	cd ..
	ln -sf src/$(PROGNAME) ../$(PROGNAME)

#   Fortran objects
objs = init_deriv_parallel.o	modules.o	V_impur.o	DFT4He3d.o	derden.o\
			dimen.o				energy.o	fforma.o	fft.o		initcg.o\
			morse.o				mates.o		poten.o		printoutc.o	r_cm.o\
			instates.o			diag.o		readenc.o	respar.o	term_alfa.o\
			timer.o				titols.o	vareps.o	varmu.o		tstgrid.o\
			s13adf.o			pdergc.o	steprk.o	steppc.o	potenimp.o\
			potenimpini.o		newder.o	redef.o

.SUFFIXES: .f90 .f .o
$(PROGNAME): $(objs)
	$(COMP) -o $(PROGNAME) $(objs) $(LD_FLAGS)
.f90.o:
	$(COMP) $(CFLAGS) -o $(@) $<;
.f.o:
	$(COMP) $(CFLAGS) -o $(@) $<;

clean:
	rm -f *.o *.bak *.lst modules/*.mod;
distclean:
	rm -f *.o *.bak *.lst modules/*.mod $(PROGNAME) ../$(PROGNAME);
