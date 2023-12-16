#------------------------------------------------------------------------------
# Makefile
#
# Description:
# This makefile contains recipes for creating 4 executable files:
# pnfam_nompi.x, pnfam_mpi.x, pnfam_betadecay.x, pnfam_hfbinfo.x
# It also contains targets for making mfam and pnfam tests.
#
# Contents:
# 1. Compile-time user inputs that need to be changed
#    - Compiler type: GFORTRAN, IFORT, CRAY or IBM
#      (This is used to automatically set compiler options)
#    - Compiler name: Actual name of the Fortran compiler
#      (including mpi wrappers, e.g. gfortran, ifort, mpif90, etc.)
#    - Verbose: 0=short diplay; 1=full command displayed
#    - Debug: Use debugging compiler options
#    - Openmp: Compile with openmp multi-threading
#    - Library and header file locations
# 2. Compiling details that should not need to be changed
#    - Override library locations for popular machines
#    - Definition of compiler options
#    - Definition of main targets
#    - Definition of object file names
#    - Definition of object file dependencies
#    - Compiling/assembling recipes
#
# M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
# E.M. Ney, UNC Chapel Hill, 2018-
# T. Li, Michigan State Univ, 2020-
#------------------------------------------------------------------------------
SHELL := /bin/bash
COMPILE_DIR = $(shell pwd)

#==============================================================================#
#                           Compile-Time Options                               #
#==============================================================================#
COMPILER_TYPE = GFORTRAN
COMPILER = mpif90

VERBOSE   = 0
DEBUG     = FALSE
USE_OPENMP= 0

USE_HBLAS = 0
USE_ADLB  = 0
READ_FUNCTIONAL = 0

# Library (L) and header/mod files (I)
L_LINEARALGEBRA = -L/usr/lib/x86_64-linux-gnu -lblas -llapack
L_ADLB = -L. -lfadlb -ladlb
I_ADLB = -I.
L_HFBTHO = -L../hfbtho_gitlab/src/hfbtho -lhfbtho
I_HFBTHO = -I../hfbtho_gitlab/src/hfbtho

# SCALAPACK only for mfam, GSL to be deprecated
L_SCALAPACK =
L_GSL =-L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas
I_GSL =-I/usr/include

#==============================================================================#
#              Override Library Locations for Specific Machines                #
#==============================================================================#
ifneq ($(shell hostname | grep dogwood),)
  L_LINEARALGEBRA = -mkl
  L_SCALAPACK =
  L_GSL =-L/nas/longleaf/apps/gsl/2.4/lib -lgsl -lgslcblas -cxxlib
  I_GSL =-I/nas/longleaf/apps/gsl/2.4/include
endif

ifeq ($(USE_ADLB),0)
  L_ADLB =
  I_ADLB =
endif

#==============================================================================#
#                           Typical compiler flags                             #
#==============================================================================#
TYPE_CHECK=0
# Defining compiler options for: IFORT FORTRAN COMPILER (ifort)
ifeq ($(COMPILER_TYPE),IFORT)
      CCOMPILER = icc
      PREPROCESSOR = -fpp -DUSE_ADLB=$(USE_ADLB) \
                          -DUSE_HBLAS=$(USE_HBLAS) \
                          -DREAD_FUNCTIONAL=$(READ_FUNCTIONAL)
      FORMAT_F90   = -fPIC -free -extend_source
      ifeq ($(DEBUG),FALSE)
            OPTIONS = $(PREPROCESSOR) -O3
      else
            OPTIONS = $(PREPROCESSOR) -g -check all -warn nodeclarations -warn nounused
      endif
      ifeq ($(USE_OPENMP),1)
            OPTIONS += -qopenmp
      endif
      TYPE_CHECK=1
endif

# Defining compiler options for: GNU FORTRAN COMPILER (gfortran)
ifeq ($(COMPILER_TYPE),GFORTRAN)
      CCOMPILER = gcc
      PREPROCESSOR = -cpp -DUSE_ADLB=$(USE_ADLB) \
                          -DUSE_HBLAS=$(USE_HBLAS) \
                          -DREAD_FUNCTIONAL=$(READ_FUNCTIONAL)
      FORMAT_F90   = -fPIC -ffree-form -ffree-line-length-none
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
      TYPE_CHECK=1
endif

# Defining compiler options for: CRAY COMPILER
ifeq ($(COMPILER_TYPE),CRAY)
      PREPROCESSOR = -e Z -DUSE_ADLB=$(USE_ADLB) \
                          -DUSE_HBLAS=$(USE_HBLAS) \
                          -DREAD_FUNCTIONAL=$(READ_FUNCTIONAL)
      FORMAT_F90   = -f free
      ifeq ($(DEBUG),FALSE)
            OPTIONS = $(PREPROCESSOR) -O3
      else
            OPTIONS = $(PREPROCESSOR) -e c -e D
      endif
      ifeq ($(USE_OPENMP),0)
            OPTIONS += -h noomp
      endif
      TYPE_CHECK=1
endif

# Defining compiler options for: CRAY COMPILER
ifeq ($(COMPILER_TYPE),IBM)
      PREPROCESSOR = '-WF,-DUSE_ADLB=$(USE_ADLB)' \
                     '-WF,-DUSE_HBLAS=$(USE_HBLAS)' \
                     '-WF,-DREAD_FUNCTIONAL=$(READ_FUNCTIONAL)'
      FORMAT_F90   = -qstrict -qfree=f90 -qsuffix=cpp=f90
      ifeq ($(DEBUG),FALSE)
            OPTIONS = $(PREPROCESSOR) -O2 -qhot
      else
            OPTIONS = $(PREPROCESSOR) -g -C -qflttrap
      endif
      ifeq ($(USE_OPENMP),1)
            OPTIONS += -qsmp=omp
      endif
      TYPE_CHECK=1
endif

ifneq ($(TYPE_CHECK),1)
      $(error Invalid compiler type)
endif
#==============================================================================#
#                      Setup main targets and variables                        #
#==============================================================================#
# Executables
PNFAM_SERIAL_EXE    = pnfam_nompi.x
PNFAM_SERIAL_OBJ    = pnfam_serial.o
PNFAM_PARALLEL_EXE  = pnfam_mpi.x
PNFAM_PARALLEL_OBJ  = pnfam_parallel.o
PNFAM_BETADECAY_EXE = pnfam_betadecay.x
PNFAM_BETADECAY_OBJ = betadecay_main.o
PNFAM_HFBINFO_EXE   = pnfam_hfbinfo.x
PNFAM_HFBINFO_OBJ   = hfbinfo.o
PNFAM_TESTS         = pnfam_tests

# Directories
DOXYGEN  = doc
MFAM     = mfam
TESTS    = tests
TEST_DIR = $(COMPILE_DIR)/$(TESTS)
MFAM_DIR = $(COMPILE_DIR)/$(MFAM)
DOC_DIR  = $(COMPILE_DIR)/$(DOXYGEN)

# Export names
export PNFAM_EXE
export PNFAM_SRC
export PNFAM_OBJ
export TEST_DIR
export MFAM_DIR
export DOC_DIR

# Export compiler info
export COMPILER
export FORMAT_F90
export OPTIONS
export VERBOSE

# Export library info
export L_LINEARALGEBRA
export L_ADLB
export I_ADLB
export L_SCALAPACK
export L_GSL
export I_GSL

.PHONY : all clean cleanall $(MFAM) $(TESTS) $(DOXYGEN) $(PNFAM_TESTS)

all : $(PNFAM_SERIAL_EXE) $(PNFAM_BETADECAY_EXE) $(PNFAM_HFBINFO_EXE)

$(MFAM) : $(PNFAM_SERIAL_EXE)
	@( make -C $(MFAM_DIR) )

$(DOXYGEN) : $(PNFAM_SERIAL_EXE)
	@( make -C $(DOC_DIR))

$(TESTS) : $(PNFAM_TESTS)
	@( make -C $(TEST_DIR)/modules )

version.inc : update-version.sh
	$(SHELL) update-version.sh

clean :
	@( rm -f *.o *.mod version.inc )
	@( cd $(TEST_DIR)/modules; make clean; cd $(MFAM_DIR); make clean; cd $(DOC_DIR); make clean)
	@( cd $(TEST_DIR); rm -f ./*/*.out ./*/*.ctr ./*/thoout.dat ./*/*.hel ./*/*.log; cd $(COMPILE_DIR))

cleanall : clean
	@( rm -f *.x )
	@( cd $(TEST_DIR)/modules; make cleanall; cd $(MFAM_DIR); make cleanall; cd $(DOC_DIR); make clean; cd $(COMPILE_DIR))

#==============================================================================#
#                        Enumerate Object Files                                #
#==============================================================================#
# MATH
# ------------------------------------------------------------------------------
BROYDEN_MIXER_OBJ      = broyden_mixer.o
RATIONAL_INTERP_OBJ    = rational_interp.o
COMPLEX_QUADRATURE_OBJ = complex_quadrature.o
POLYNOMIAL_FIT_OBJ     = polyfit.o

# UTILITIES
# ------------------------------------------------------------------------------
CONSTANTS_OBJ          = constants.o
LOGGER_OBJ             = logger.o
LOGGER_PARALLEL_OBJ    = logger_parallel.o
DENSITY_SET_TYPE_OBJ   = density_set_type.o
BLOCKMATRIX_TYPE_OBJ   = blockmatrix_type.o
BBLOCKMATRIX_TYPE_OBJ  = bblockmatrix_type.o
EXTFIELD_TYPE_OBJ      = extfield_type.o

# HFBTHO
# ------------------------------------------------------------------------------
HFBTHO_SETUP_OBJ       = hfbtho_setup.o
HFBTHO_INTERFACE_OBJ   = hfbtho_interface.o

# BETADECAY
# ------------------------------------------------------------------------------
FERMI_OBJ              = fermi.o
PHASESPACE_OBJ         = phasespace.o
BETADECAY_OBJ          = betadecay.o

# PNFAM
# ------------------------------------------------------------------------------
EXTFIELD_OBJ           = extfield.o
INTERACTION_OBJ        = interaction.o
ifeq ($(USE_HBLAS),1)
  HAMILTONIAN_OBJ      = hamiltonian_blas.o
else
  HAMILTONIAN_OBJ      = hamiltonian.o
endif
PNFAM_COMMON_OBJ       = pnfam_common.o
PNFAM_OUTPUT_OBJ       = pnfam_output.o
PNFAM_STORAGE_OBJ       = pnfam_storage.o
PNFAM_SOLVER_OBJ       = pnfam_solver.o

#==============================================================================#
#                            Define Dependencies                               #
#==============================================================================#
# MATH
# ------------------------------------------------------------------------------
$(BROYDEN_MIXER_OBJ) : $(BBLOCKMATRIX_TYPE_OBJ)
$(RATIONAL_INTERP_OBJ) :
$(COMPLEX_QUADRATURE_OBJ) :
$(POLYNOMIAL_FIT_OBJ) : polyfit.cpp
ifeq ($(VERBOSE),1)
	$(CCOMPILER) -c -o $@ $< $(OPTIONS) $(I_GSL)
else
	@echo "Compiling  \"$@\"..."
	@$(CCOMPILER) -c -o $@ $< $(OPTIONS) $(I_GSL)
endif

# UTILITIES
# ------------------------------------------------------------------------------
$(CONSTANTS_OBJ)        : version.inc
$(LOGGER_OBJ)           :
$(LOGGER_PARALLEL_OBJ)  : $(LOGGER_OBJ)
$(DENSITY_SET_TYPE_OBJ) :
$(BLOCKMATRIX_TYPE_OBJ) :
$(EXTFIELD_TYPE_OBJ)    : $(BLOCKMATRIX_TYPE_OBJ)

# HFBTHO
# ------------------------------------------------------------------------------
$(HFBTHO_SETUP_OBJ)     :
$(HFBTHO_INTERFACE_OBJ) : $(CONSTANTS_OBJ) $(HFBTHO_SETUP_OBJ)

# BETADECAY
# ------------------------------------------------------------------------------
$(FERMI_OBJ)      : $(CONSTANTS_OBJ)
$(PHASESPACE_OBJ) : $(RATIONAL_INTERP_OBJ) $(COMPLEX_QUADRATURE_OBJ) $(FERMI_OBJ) \
                    $(POLYNOMIAL_FIT_OBJ)
$(BETADECAY_OBJ)  : $(CONSTANTS_OBJ) $(LOGGER_OBJ) $(BLOCKMATRIX_TYPE_OBJ) \
                    $(HFBTHO_INTERFACE_OBJ) $(PHASESPACE_OBJ) $(FERMI_OBJ) \
                    $(PNFAM_COMMON_OBJ)

# PNFAM
# ------------------------------------------------------------------------------
$(EXTFIELD_OBJ)     : $(CONSTANTS_OBJ) $(LOGGER_OBJ) $(HFBTHO_INTERFACE_OBJ) \
                      $(EXTFIELD_TYPE_OBJ)
$(INTERACTION_OBJ)  : $(CONSTANTS_OBJ) $(LOGGER_OBJ) $(HFBTHO_INTERFACE_OBJ)
$(HAMILTONIAN_OBJ)  : $(CONSTANTS_OBJ) $(HFBTHO_INTERFACE_OBJ) $(INTERACTION_OBJ) \
                      $(BLOCKMATRIX_TYPE_OBJ) $(DENSITY_SET_TYPE_OBJ)
$(PNFAM_COMMON_OBJ) : $(CONSTANTS_OBJ) $(LOGGER_OBJ) $(HFBTHO_INTERFACE_OBJ) \
                      $(BLOCKMATRIX_TYPE_OBJ) $(EXTFIELD_TYPE_OBJ) \
                      $(INTERACTION_OBJ) $(HAMILTONIAN_OBJ) $(EXTFIELD_OBJ) \
                      $(BROYDEN_MIXER_OBJ) $(BBLOCKMATRIX_TYPE_OBJ)
$(PNFAM_SOLVER_OBJ) : $(CONSTANTS_OBJ) $(PNFAM_COMMON_OBJ) $(PNFAM_OUTPUT_OBJ) \
                      $(BROYDEN_MIXER_OBJ) $(BBLOCKMATRIX_TYPE_OBJ) $(PNFAM_STORAGE_OBJ)
$(PNFAM_OUTPUT_OBJ) : $(PNFAM_COMMON_OBJ)
$(PNFAM_STORAGE_OBJ): $(PNFAM_COMMON_OBJ) $(INTERACTION_OBJ) $(CONSTANTS_OBJ) $(HAMILTONIAN_OBJ)

# EXECUTABLES
# ------------------------------------------------------------------------------
$(PNFAM_SERIAL_OBJ)    : $(PNFAM_SOLVER_OBJ)
$(PNFAM_PARALLEL_OBJ)  : $(PNFAM_SOLVER_OBJ) $(LOGGER_PARALLEL_OBJ)
$(PNFAM_BETADECAY_OBJ) : $(BETADECAY_OBJ)
$(PNFAM_HFBINFO_OBJ)   : $(CONSTANTS_OBJ) $(LOGGER_OBJ) $(HFBTHO_INTERFACE_OBJ) \
                         $(BLOCKMATRIX_TYPE_OBJ) $(LOGGER_OBJ)
$(PNFAM_TESTS)         : $(BLOCKMATRIX_TYPE_OBJ) $(COMPLEX_QUADRATURE_OBJ) $(FERMI_OBJ) \
                         $(PHASESPACE_OBJ)

#==============================================================================#
#                              Program Builds                                  #
#==============================================================================#
SERIAL_OBJECTS = \
   $(CONSTANTS_OBJ) \
   $(LOGGER_OBJ) \
   $(HFBTHO_SETUP_OBJ) \
   $(HFBTHO_INTERFACE_OBJ) \
   $(BLOCKMATRIX_TYPE_OBJ) \
   $(BBLOCKMATRIX_TYPE_OBJ) \
   $(EXTFIELD_OBJ) \
   $(EXTFIELD_TYPE_OBJ) \
   $(INTERACTION_OBJ) \
   $(HAMILTONIAN_OBJ) \
   $(DENSITY_SET_TYPE_OBJ) \
   $(BROYDEN_MIXER_OBJ) \
   $(PNFAM_COMMON_OBJ) \
   $(PNFAM_OUTPUT_OBJ) \
   $(PNFAM_STORAGE_OBJ) \
   $(PNFAM_SOLVER_OBJ)

PARALLEL_OBJECTS = \
   $(SERIAL_OBJECTS) \
   $(LOGGER_PARALLEL_OBJ)

BETADECAY_OBJECTS = \
   $(CONSTANTS_OBJ) \
   $(LOGGER_OBJ) \
   $(BLOCKMATRIX_TYPE_OBJ) \
   $(RATIONAL_INTERP_OBJ) \
   $(COMPLEX_QUADRATURE_OBJ) \
   $(POLYNOMIAL_FIT_OBJ) \
   $(HFBTHO_SETUP_OBJ) \
   $(HFBTHO_INTERFACE_OBJ) \
   $(FERMI_OBJ) \
   $(PHASESPACE_OBJ) \
   $(BETADECAY_OBJ)

HFBINFO_OBJECTS = \
   $(CONSTANTS_OBJ) \
   $(LOGGER_OBJ) \
   $(BLOCKMATRIX_TYPE_OBJ) \
   $(HFBTHO_SETUP_OBJ) \
   $(HFBTHO_INTERFACE_OBJ)

LIBRARIES = $(L_HFBTHO) $(L_LINEARALGEBRA) $(L_ADLB) $(L_GSL)
INCLUDE   = $(I_HFBTHO)

$(PNFAM_SERIAL_EXE) : $(PNFAM_SERIAL_OBJ) $(SERIAL_OBJECTS)
ifeq ($(VERBOSE),1)
	$(COMPILER) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(LIBRARIES)
else
	@echo "Assembling \"$@\"..."
	@$(COMPILER) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(LIBRARIES)
endif

$(PNFAM_PARALLEL_EXE) : $(PNFAM_PARALLEL_OBJ) $(PARALLEL_OBJECTS)
ifeq ($(VERBOSE),1)
	$(COMPILER) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(LIBRARIES)
else
	@echo "Assembling \"$@\"..."
	@$(COMPILER) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(LIBRARIES)
endif

$(PNFAM_BETADECAY_EXE) : $(PNFAM_BETADECAY_OBJ) $(SERIAL_OBJECTS)
ifeq ($(VERBOSE),1)
	$(COMPILER) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(LIBRARIES)
else
	@echo "Assembling \"$@\"..."
	@$(COMPILER) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(LIBRARIES)
endif

$(PNFAM_HFBINFO_EXE) : $(PNFAM_HFBINFO_OBJ) $(HFBINFO_OBJECTS)
ifeq ($(VERBOSE),1)
	$(COMPILER) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(LIBRARIES)
else
	@echo "Assembling \"$@\"..."
	@$(COMPILER) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(LIBRARIES)
endif

# Compile object files (header/mod files need to be included here)
%.o : %.f90
ifeq ($(VERBOSE),1)
	$(COMPILER) -c $< $(FORMAT_F90) $(OPTIONS) $(INCLUDE)
else
	@echo "Compiling  \"$@\"..."
	@$(COMPILER) -c $< $(FORMAT_F90) $(OPTIONS) $(INCLUDE)
endif
