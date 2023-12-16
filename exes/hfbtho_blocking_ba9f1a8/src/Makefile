
SHELL := /bin/bash
COMPILE_DIR = $(shell pwd)

# Basic options
#   - COMPUTER..............: current choice is DESKTOP, QUARTZ LASSEN; only used to automatically set LINEAR_ALGEBRA below
#   - COMPILER .............: type of compiler: GFORTRAN, IFORT, CRAY or IBM; used to automatically set compiler options
#   - FORTRAN_MPI ..........: actual name of the Fortran compiler. Standard choices are gfortran, mpif90, ifort, xlf2008_r, etc.
COMPUTER    = DESKTOP
COMPILER    = GFORTRAN
FORTRAN_MPI = gfortran
# Build options (TRUE or FALSE)
#   - VERBOSE ..............: 0: short display; 1: full command displayed
#   - DEBUG ................: activates debugging options for each compiler
#   - DOCUMENTATION ........: builds doxygen documentation
VERBOSE       = 0
DEBUG         = FALSE
DOCUMENTATION = FALSE
# Preprocessor options
#   - USE_ESSL .............: uses IBM ESSL library for LAPACK and BLAS (different routine calls)
#   - USE_OPENMP ...........: enables OpenMP multithreading
#   - USE_PETSC ............: uses Petsc library for ...
USE_ESSL   = 0
USE_OPENMP = 0
USE_PETSC  = 0
# Physics options
#   - USE_MPI .............: 1, ANL optimization
#                            2, parallel mass table, driplines or PES mode
#   - DRIP_LINES ..........: 1, calculates a mass table from dripline to dripline; 0, inactive
#   - DO_MASSTABLE ........: 1, calculates a section of the mass table; 0, inactive
#   - DO_PES ..............: 1, calculates a potential energy surface for a set of nuclei; 0, inactive
#   - READ_FUNCTIONAL .....: 1, reads the parameters of the functional from a file
#                            0, energy functional is defined in the code based on the value of the input keyword 'functional'
#   - GOGNY_SYMMETRIES.....: 1, for production runs; assumes several symmetries in the finite range matrix elements.
#                            0, for debugging the finite range matrix elements; makes no assumption of any symmetry in the finite range matrix elements.
#                               Requires compiler compatible with 2008 standard of up to 15 dimensional arrays.
#   - GOGNY_HYPER..........: 1, uses hypergeometric function to calculate finite range matrix elements; accurate for big and small basis size
#                            0, uses the direct Gogny transformation to calculate the finite range matrix elements; accurate for small basis size ONLY (N<20)
#   - USE_QRPA ............: 1, produces output for the QRPA-pnFAM code of UNC Chapel Hill; 0, inactive
USE_MPI          = 0
DRIP_LINES       = 0
DO_MASSTABLE     = 0
DO_PES           = 0
READ_FUNCTIONAL  = 0
GOGNY_SYMMETRIES = 1
GOGNY_HYPER      = 1
USE_QRPA         = 0

# Names and paths to BLAS and LAPACK Libraries
#   - LC CLUSTER (MKL)
#       * Quartz, OpenMP ..: -L/usr/tce/packages/mkl/$(VERSION_LIBRARY)/lib -I/usr/tce/packages/mkl/$(VERSION_LIBRARY)/include \
#                            -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_avx2  -liomp
#       * Quartz, serial ..: -L/usr/tce/packages/mkl/$(VERSION_LIBRARY)/lib -I/usr/tce/packages/mkl/$(VERSION_LIBRARY)/include \
#                            -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
#   - LINUX DESKTOP
#       * local ...........: -L$(HOME)/local/lib -llapack_LINUX -lblas_LINUX
ifeq ($(COMPUTER),DESKTOP)
      LINEAR_ALGEBRA = -L/usr/lib/x86_64-linux-gnu -llapack -L/usr/lib/x86_64-linux-gnu -lblas
endif
ifeq ($(COMPUTER),QUARTZ)
      VERSION_LIBRARY = mkl-2020.0
      LINEAR_ALGEBRA = -L/usr/tce/packages/mkl/$(VERSION_LIBRARY)/lib -I/usr/tce/packages/mkl/$(VERSION_LIBRARY)/include \
                       -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_avx2  -liomp5
endif
ifeq ($(COMPUTER),LASSEN)
      LINEAR_ALGEBRA = -I$(ESSLHEADERDIR) -L$(ESSLLIBDIR64) -R$(ESSLLIBDIR64) -lesslsmp -llapackforessl
endif

#======================================================================#
#  Nothing beyond this line should be changed, in principle            #
#======================================================================#

# Add in options from file
ifneq ("$(wildcard MakeIncl.${HFBTHO_ARCH})","")
       include MakeIncl.${HFBTHO_ARCH}
endif

# Check if compile directory exists, returns error code and exit otherwise
EXISTS = $(firstword $(wildcard $(COMPILE_DIR)) )
$(info $$COMPILE_DIR is [${COMPILE_DIR}] )
ifeq ($(EXISTS), )
      $(error ${COMPILE_DIR} directory does not exist)
endif

# Consistency checks between MPI option and dripline or mass table modes
ifeq ($(DRIP_LINES),1)
     ifeq ($(DO_MASSTABLE),1)
           $(error DRIP_LINES and DO_MASSTABLE modes are incompatible)
     endif
     ifneq ($(USE_MPI),2)
            $(error DRIP_LINES mode requires USE_MPI=2)
     endif
endif
ifeq ($(USE_MPI),2)
      ifeq ($(DRIP_LINES),0)
            ifeq ($(DO_MASSTABLE),0)
                  ifeq ($(DO_PES),0)
                        $(error USE_MPI=2 not allowed for regular HFBTHO calculation)
                  endif
            endif
      endif
endif

ifeq ($(DOCUMENTATION),TRUE)
      DOXYGEN = doc
endif

# Names
HFBTHO_EXE     = hfbtho_main
HFBTHO_SOURCE  = hfbtho_main.f90
HFBTHO_OBJ     = hfbtho_main.o

# Defining compiler options for: IFORT FORTRAN COMPILER (ifort)
ifeq ($(COMPILER),IFORT)
      FORMAT_F90   = -fPIC -free -extend_source
      PREPROCESSOR = -fpp -DUSE_OPENMP=$(USE_OPENMP) \
                          -DUSE_MPI=$(USE_MPI)  \
                          -DUSE_PETSC=$(USE_PETSC) \
                          -DDO_PES=$(DO_PES) \
                          -DDO_MASSTABLE=$(DO_MASSTABLE) \
                          -DDRIP_LINES=$(DRIP_LINES) \
                          -DREAD_FUNCTIONAL=$(READ_FUNCTIONAL) \
                          -DGOGNY_SYMMETRIES=$(GOGNY_SYMMETRIES) \
                          -DGOGNY_HYPER=$(GOGNY_HYPER) \
                          -DUSE_QRPA=$(USE_QRPA) \
                          -DUSE_ESSL=$(USE_ESSL) $(FFLAGS)
      ifeq ($(DEBUG),FALSE)
            OPTIONS = $(PREPROCESSOR) -O3
      else
            OPTIONS = $(PREPROCESSOR) -g -check all -warn nodeclarations -warn nounused
      endif
      ifeq ($(USE_OPENMP),1)
            OPTIONS += -qopenmp
      endif
endif

# Defining compiler options for: GNU FORTRAN COMPILER (gfortran)
ifeq ($(COMPILER),GFORTRAN)
      FORMAT_F90   = -fPIC -ffree-form -ffree-line-length-none
      PREPROCESSOR = -cpp -DUSE_OPENMP=$(USE_OPENMP) \
                          -DUSE_MPI=$(USE_MPI)  \
                          -DUSE_PETSC=$(USE_PETSC) \
                          -DDO_PES=$(DO_PES) \
                          -DDO_MASSTABLE=$(DO_MASSTABLE) \
                          -DDRIP_LINES=$(DRIP_LINES) \
                          -DREAD_FUNCTIONAL=$(READ_FUNCTIONAL) \
                          -DGOGNY_SYMMETRIES=$(GOGNY_SYMMETRIES) \
                          -DGOGNY_HYPER=$(GOGNY_HYPER) \
                          -DUSE_QRPA=$(USE_QRPA) \
                          -DUSE_ESSL=$(USE_ESSL) $(FFLAGS)
      ifeq ($(DEBUG),FALSE)
            OPTIONS = $(PREPROCESSOR) -O3
      else
            OPTIONS = $(PREPROCESSOR) -g -O0 -Wall \
                      -Warray-bounds -Wunderflow -Warray-temporaries \
                      -Wcharacter-truncation -Wtabs -Wintrinsic-shadow -Walign-commons -frange-check \
                      -fbounds-check -Wconversion -Wuninitialized -pedantic \
                      -finit-real=nan \
                      -ftrapv
      endif
      ifeq ($(USE_OPENMP),1)
            OPTIONS += -fopenmp
      endif
endif

# Defining compiler options for: CRAY COMPILER
ifeq ($(COMPILER),CRAY)
      FORMAT_F90   = -f free
      PREPROCESSOR = -e Z -DUSE_OPENMP=$(USE_OPENMP) \
                          -DUSE_MPI=$(USE_MPI)  \
                          -DUSE_PETSC=$(USE_PETSC) \
                          -DDO_PES=$(DO_PES) \
                          -DDO_MASSTABLE=$(DO_MASSTABLE) \
                          -DDRIP_LINES=$(DRIP_LINES) \
                          -DREAD_FUNCTIONAL=$(READ_FUNCTIONAL) \
                          -DGOGNY_SYMMETRIES=$(GOGNY_SYMMETRIES) \
                          -DGOGNY_HYPER=$(GOGNY_HYPER) \
                          -DUSE_QRPA=$(USE_QRPA) \
                          -DUSE_ESSL=$(USE_ESSL) $(FFLAGS)
      ifeq ($(DEBUG),FALSE)
            OPTIONS = $(PREPROCESSOR) -O3
      else
            OPTIONS = $(PREPROCESSOR) -e c -e D
      endif
      ifeq ($(USE_OPENMP),0)
            OPTIONS += -h noomp
      endif
endif

# Defining compiler options for: CRAY COMPILER
ifeq ($(COMPILER),IBM)
      FORMAT_F90   = -qstrict -qfree=f90 -qsuffix=cpp=f90
      PREPROCESSOR = '-WF,-DUSE_OPENMP=$(USE_OPENMP)' \
                     '-WF,-DUSE_MPI=$(USE_MPI)' \
                     '-WF,-DUSE_PETSC=$(USE_PETSC)' \
                     '-WF,-DDO_PES=$(DO_PES)' \
                     '-WF,-DDO_MASSTABLE=$(DO_MASSTABLE)' \
                     '-WF,-DDRIP_LINES=$(DRIP_LINES)' \
                     '-WF,-DREAD_FUNCTIONAL=$(READ_FUNCTIONAL)' \
                     '-WF,-DGOGNY_SYMMETRIES=$(GOGNY_SYMMETRIES)' \
                     '-WF,-DGOGNY_HYPER=$(GOGNY_HYPER)' \
                     '-WF,-DUSE_QRPA=$(USE_QRPA)' \
                     '-WF,-DUSE_ESSL=$(USE_ESSL)' $(FFLAGS)
      ifeq ($(DEBUG),FALSE)
            OPTIONS = $(PREPROCESSOR) -O2 -qhot
      else
            OPTIONS = $(PREPROCESSOR) -g -C -qflttrap
      endif
      ifeq ($(USE_OPENMP),1)
            OPTIONS += -qsmp=omp
      endif
endif

#=========================#
# Beginning of the action #
#=========================#

ANL_DIR = $(COMPILE_DIR)/anl
HFBTHO_DIR = $(COMPILE_DIR)/hfbtho
DOC_DIR = $(HFBTHO_DIR)/doc

# Export HFBTHO names
export HFBTHO_EXE
export HFBTHO_SOURCE
export HFBTHO_OBJ
export HFBTHO_DIR
export DOC_DIR

# Export main compiler options
export FORTRAN_MPI
export PRECISION
export FORMAT_F77
export FORMAT_F90
export OPTIONS
export VERBOSE
export LINEAR_ALGEBRA

# Options are: anl, $(HFBTHO_EXE), $(DOXYGEN)
all: $(HFBTHO_EXE) $(DOXYGEN)

anl: $(HFBTHO_EXE)
	( cd $(ANL_DIR); $(MAKE); cd .. )

# Doxygen documentation
$(DOXYGEN): $(HFBTHO_EXE)
	( cd $(HFBTHO_DIR); $(MAKE) -C $(DOC_DIR); cd .. )

# HFBTHO executable
$(HFBTHO_EXE):
	( cd $(HFBTHO_DIR); $(MAKE); cd .. )

# Cleaning
clean :
	( rm -f *.o *.oo *.ipo *.mod; cd $(ANL_DIR); $(MAKE) clean; cd $(HFBTHO_DIR); $(MAKE) clean; cd ..)
	$(MAKE) clean -C $(DOC_DIR)

mrproper : clean
	$(MAKE) mrproper -C $(HFBTHO_DIR)
