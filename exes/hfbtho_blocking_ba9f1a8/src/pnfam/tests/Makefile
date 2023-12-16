################################################################################
# Makefile
# 
# Build test cases for PNFAM. There is an `export` command in the main PNFAM
# Makefile which sets the variables needed here. You shouldn't make anything
# from this directory. Run `make tests` in the main PNFAM directory to build
# these tests. 
#
# M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
################################################################################


all: test_blockmatrix_type.x test_complex_quadrature.x test_fermi.x test_phasespace.x

clean:
	-rm -vf *.mod *.o

cleanall: clean
	-rm -vf test_blockmatrix_type.x test_complex_quadrature.x test_fermi.x test_phasespace.x

# Fortran source
%.o : %.f90
	$(FC) -c -o $@ $< $(FFLAGS)
%.o : ../%.f90
	$(FC) -c -o $@ $< $(FFLAGS)

# C++ source
polyfit.o : ../polyfit.cpp
	$(CC) -c -o $@ $< $(CFLAGS) $(GSLHEADER)

# Targets
test_blockmatrix_type.x : blockmatrix_type.o blockmatrix_type_test.o
	$(FC) -o $@ $^ $(FFLAGS) $(LAPACK)
test_complex_quadrature.x : complex_quadrature.o complex_quadrature_test.o
	$(FC) -o $@ $^ $(FFLAGS)
test_fermi.x : pnfam_tests.o constants.o fermi.o fermi_test.o
	$(FC) -o $@ $^ $(FFLAGS)
test_phasespace.x : pnfam_tests.o complex_quadrature.o fermi.o polyfit.o phasespace.o phasespace_test.o
	$(FC) -o $@ $^ $(FFLAGS) $(GSL)
