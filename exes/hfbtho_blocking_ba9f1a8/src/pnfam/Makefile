################################################################################
# Makefile
#
# This Makefile is structured fairly linearly --- at the top are some reasonable
# default settings.  Following, we override these defaults and set up specific
# configurations for the various machines we use.  At the bottom are the nuts
# and bolts of the compilation which likely needn't be touched frequently.
#
# TARGETS:
# nompi    => pnfam_nompi.x
# mpi      => pnfam_mpi.x
# beta     => betadecay.x
# hfbinfo  => hfbinfo.x
# tests    => tests/
#
# M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
################################################################################

# ADJUST THE VARIABLES IN THIS SECTION FOR YOUR PLATFORM AND COMPILER

# the Fortran and C++ compilers (Note: you can use non-MPI compilers if you're only
# building the non-MPI version, i.e. "make nompi")
FC = mpif90
CXX = mpicxx

# Default compiler flags
FFLAGS_OPT   =-O2
FFLAGS_DEBUG =-O0 -Wall -Wextra -Wconversion -pedantic -fbounds-check -fimplicit-none -g -fbacktrace
CXXFLAGS     =-O2
OPENMP       =-fopenmp

# Locations for dependencies. YOU WILL ALMOST CERTAINLY HAVE TO MODIFY THESE
# Basic Linear Algebra Subroutines (BLAS):
LAPACK       =-lblas
# Gnu Scientific Library (GSL):
GSL          =-L/usr/local/lib -lgsl -lcblas
GSLHEADER    =-I/usr/local/include
# Asynchronous Dynamic Load-Balancing library (ADLB)
# Note: in the default configuration, libfadlb.a, libadlb.a, and adlbf.h must be
# found in the current directory.
ADLB       =-L. -lfadlb -ladlb

### Compile-time options ---
# can be overridden in the make command, e.g. "make openmp=1 all"
debug     = 0
openmp    = 0

# suggested flags for gfortran on OS X
ifeq ($(shell uname), Darwin)
  LAPACK =-framework Accelerate
  OPENMP =-fopenmp
  FFLAGS_OPT =-O3 -std=f2008 -pedantic
endif

# suggested flags for Intel Fortran on Linux:
#FFLAGS_OPT =-O3 -diag-disable=remark -lstdc++ -stand f08
#CXXFLAGS   =-O3 -diag-disable=remark
#ADLB       =<ADLB -I... -L... ...>
#GSLHEADER  =-<GSL -I...>
#GSL        =-<GSL -L...> -lgsl -lgslcblas -cxxlib
#LAPACK     =<LAPACK/BLAS -L... ...>
#OPENMP     =-openmp

################################################################################

### Compile-time preprocessing ---
# Debugging on/off
ifeq ($(debug),1)
  FFLAGS = $(FFLAGS_DEBUG)
else
  FFLAGS = $(FFLAGS_OPT)
endif

# OpenMP disabled by default but can be switched on (e.g. make openmp=1 all)
ifneq ($(openmp),1)
  OPENMP =
endif

# Unless overridden in the custom setting above, the CFLAGS should mirror the FFLAGS
ifeq ($(CXXFLAGS),)
  CXXFLAGS = $(FFLAGS)
endif

# This allows tests/Makefile to use the same settings automatically
export FC CC FFLAGS CFLAGS LAPACK GSLHEADER GSL

### Phony Targets ---
# The targets are put together into various groups:
# serial             - builds pnfam_serial.x
# parallel           - builds pnfam_parallel.x
# beta               - builds betadecay.x
# hfbinfo            - builds hfbinfo.x
# all (default)      - runs serial, beta, and hfbinfo
# tests              - builds module tests
# clean              - removes *.o and *.mod
# cleanall           - runs clean and also removes *.x
.PHONY : all serial parallel beta hfbinfo tests clean cleanall


# Main targets
mpi     : pnfam_mpi.x
nompi   : pnfam_nompi.x

beta    : betadecay.x

hfbinfo : hfbinfo.x

all     : nompi mpi beta hfbinfo

tests :
	make -C tests

clean :
	-rm -vf ./*.o
	-rm -vf ./*.mod
	cd tests; make clean

cleanall : clean
	-rm -vf pnfam_nompi.x pnfam_mpi.x betadecay.x hfbinfo.x
	cd tests; make cleanall


### Objects and executables ---
# Require OpenMP
hamiltonian.o : hamiltonian.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(OPENMP)
pnfam_common.o : pnfam_common.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(OPENMP)

# Require ADLB and OpenMP
pnfam_parallel.o : pnfam_parallel.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(OPENMP) $(ADLB)

# Require GSL
polyfit.o : polyfit.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(GSLHEADER)

# Other Fortran source
%.o : %.f90
	$(FC) -c -o $@ $< $(FFLAGS)

# Other C++ source
%.o : %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)


# Executables
pnfam_nompi.x : constants.o logger.o blockmatrix_type.o hfbtho_basis.o extfield_type.o extfield.o interaction.o broyden_mixer.o density_set_type.o hamiltonian.o pnfam_solver.o fermi.o complex_quadrature.o polyfit.o phasespace.o pnfam_common.o pnfam_serial.o
	$(FC) -o $@ $^ $(FFLAGS) $(OPENMP) $(LAPACK) $(GSL)

pnfam_mpi.x : constants.o logger.o logger_parallel.o blockmatrix_type.o hfbtho_basis.o extfield_type.o extfield.o interaction.o broyden_mixer.o density_set_type.o hamiltonian.o pnfam_solver.o fermi.o complex_quadrature.o polyfit.o phasespace.o pnfam_common.o pnfam_parallel.o
	$(FC) -o $@ $^ $(FFLAGS) $(OPENMP) $(LAPACK) $(GSL) $(ADLB)

betadecay.x : constants.o blockmatrix_type.o logger.o hfbtho_basis.o complex_quadrature.o fermi.o polyfit.o phasespace.o betadecay.o betadecay_main.o
	$(FC) -o $@ $^ $(FFLAGS) $(OPENMP) $(LAPACK) $(GSL)

hfbinfo.x : constants.o blockmatrix_type.o logger.o hfbtho_basis.o hfbinfo.o
	$(FC) -o $@ $^ $(FFLAGS) $(LAPACK)
