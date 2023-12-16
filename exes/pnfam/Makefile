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
COMPILER = gfortran

VERBOSE   = 0
DEBUG     = FALSE
USE_OPENMP= 1
# Turn above on with two-body currents

USE_HBLAS = 1
READ_FUNCTIONAL = 0

# Library (L) and header/mod files (I)
L_LINEARALGEBRA = -L/usr/lib/x86_64-linux-gnu -lblas -llapack
#L_HFBTHO = -L../hfbtho_gitlab/src/hfbtho -lhfbtho
#I_HFBTHO = -I../hfbtho_gitlab/src/hfbtho

L_HFBTHO = -L../hfbtho_blocking_ba9f1a8/src/hfbtho -lhfbtho
I_HFBTHO = -I../hfbtho_blocking_ba9f1a8/src/hfbtho


#==============================================================================#
#              Override Library Locations for Specific Machines                #
#==============================================================================#
ifneq ($(shell hostname | grep dogwood),)
  L_LINEARALGEBRA = -mkl
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
PNFAM_EXE = pnfam_main.x
PNFAM_LIB = libpnfam.a
CONTOUR_EXE = contour_main.x

# Directories
DOXYGEN  = doc
DOC_DIR  = $(COMPILE_DIR)/$(DOXYGEN)

# Export names
export PNFAM_EXE
export PNFAM_OBJ
export DOC_DIR

# Export compiler info
export COMPILER
export FORMAT_F90
export OPTIONS
export VERBOSE

# Export library info
export L_LINEARALGEBRA

.PHONY : all clean cleanall $(DOXYGEN)

all : $(PNFAM_EXE) $(PNFAM_LIB) $(DOXYGEN)

$(DOXYGEN) : $(PNFAM_EXE)
	@( make -C $(DOC_DIR))

clean :
	@( rm -f *.o *.mod )
	@( cd $(DOC_DIR); make clean )

cleanall : clean
	@( rm -f *.x *.a $(VERSION_INC) )

#==============================================================================#
#                        Enumerate Object Files                                #
#==============================================================================#
# Files to modify if not using hfbtho
HFB_INTERFACE_OBJ       = hfbtho_interface.o
HFB_SOLUTION_OBJ        = hfbtho_solution.o
PNFAM_EXTFIELD_OBJ      = pnfam_extfield.o
PNFAM_INTERACTION_OBJ   = pnfam_interaction.o
ifeq ($(USE_HBLAS),1)
  PNFAM_HAMILTONIAN_OBJ = pnfam_hamiltonian_blas.o
else
  PNFAM_HAMILTONIAN_OBJ = pnfam_hamiltonian.o
endif
PNFAM_SETUP_OBJ         = pnfam_setup.o

# pnFAM-specific files
VERSION_INC             = pnfam_version.inc
VERSION_SH              = pnfam_version.sh
TYPE_BLOCKMATRIX_OBJ    = pnfam_type_blockmatrix.o
TYPE_BBM_OBJ            = pnfam_type_bbm.o
TYPE_EXTFIELD_OBJ       = pnfam_type_extfield.o
TYPE_DENSITYSET_OBJ     = pnfam_type_densityset.o
PNFAM_LOGGER_OBJ        = pnfam_logger.o
PNFAM_CONSTANTS_OBJ     = pnfam_constants.o
PNFAM_BROYDEN_OBJ       = pnfam_broyden.o
PNFAM_TXTOUTPUT_OBJ     = pnfam_txtoutput.o
PNFAM_STORAGE_OBJ       = pnfam_storage.o
PNFAM_SOLVER_OBJ        = pnfam_solver.o
CONTOUR_SETUP_OBJ       = contour_setup.o

PNFAM_EXTFIELD_2BC_OBJ  = pnfam_extfield_2bc.o
TYPE_EXTFIELD_2BC_OBJ   = pnfam_type_extfield_2bc.o
TYPE_GAMDEL_2BC_OBJ     = pnfam_type_gamdel_2bc.o
PNFAM_SPATIAL_MTXELS_OBJ= pnfam_spatial_mtxels.o

# Programs
PNFAM_OBJ               = pnfam_prog.o
CONTOUR_OBJ             = contour_prog.o

#==============================================================================#
#                            Define Dependencies                               #
#==============================================================================#
$(VERSION_INC) : $(VERSION_SH)
ifeq ($(VERBOSE),1)
	$(SHELL) $(VERSION_SH) > $(VERSION_INC)
else
	@echo "Writing    \"$(VERSION_INC)\"..."
	@$(SHELL) $(VERSION_SH) > $(VERSION_INC)
endif

$(HFB_INTERFACE_OBJ)     :
$(HFB_SOLUTION_OBJ)      : $(PNFAM_CONSTANTS_OBJ) $(HFB_INTERFACE_OBJ)

$(PNFAM_EXTFIELD_OBJ)    : $(PNFAM_LOGGER_OBJ) $(PNFAM_CONSTANTS_OBJ) $(TYPE_EXTFIELD_OBJ) \
                               $(HFB_SOLUTION_OBJ) $(PNFAM_EXTFIELD_2BC_OBJ) $(TYPE_EXTFIELD_2BC_OBJ)
$(PNFAM_INTERACTION_OBJ) : $(PNFAM_LOGGER_OBJ) $(PNFAM_CONSTANTS_OBJ) $(HFB_SOLUTION_OBJ)
$(PNFAM_HAMILTONIAN_OBJ) : $(PNFAM_LOGGER_OBJ) $(PNFAM_CONSTANTS_OBJ) $(TYPE_BLOCKMATRIX_OBJ) \
                               $(TYPE_DENSITYSET_OBJ) $(HFB_SOLUTION_OBJ) $(PNFAM_INTERACTION_OBJ)

$(PNFAM_SETUP_OBJ)       : $(PNFAM_LOGGER_OBJ) $(PNFAM_CONSTANTS_OBJ) \
                               $(TYPE_BLOCKMATRIX_OBJ) $(TYPE_BBM_OBJ) $(TYPE_EXTFIELD_OBJ) \
                               $(HFB_SOLUTION_OBJ) $(PNFAM_INTERACTION_OBJ) $(PNFAM_HAMILTONIAN_OBJ) \
                               $(PNFAM_BROYDEN_OBJ)
$(TYPE_BLOCKMATRIX_OBJ)  :
$(TYPE_BBM_OBJ)          : $(PNFAM_LOGGER_OBJ) $(BLOCKMATRIX_TYPE_OBJ)
$(TYPE_EXTFIELD_OBJ)     : $(TYPE_BLOCKMATRIX_OBJ) $(TYPE_GAMDEL_2BC_OBJ)
$(TYPE_DENSITYSET_OBJ)   :
$(PNFAM_LOGGER_OBJ)      :
$(PNFAM_CONSTANTS_OBJ)   : $(PNFAM_LOGGER_OBJ) $(VERSION_INC)
$(PNFAM_BROYDEN_OBJ)     : $(PNFAM_LOGGER_OBJ) $(TYPE_BBM_OBJ)
$(PNFAM_TXTOUTPUT_OBJ)   : $(PNFAM_SETUP_OBJ) $(PNFAM_CONSTANTS_OBJ)
$(PNFAM_STORAGE_OBJ)     : $(PNFAM_CONSTANTS_OBJ) $(PNFAM_SETUP_OBJ) $(PNFAM_INTERACTION_OBJ) \
                           $(PNFAM_BROYDEN_OBJ) $(PNFAM_EXTFIELD_2BC_OBJ) $(TYPE_GAMDEL_2BC_OBJ)
$(PNFAM_SOLVER_OBJ)      : $(PNFAM_LOGGER_OBJ) $(PNFAM_CONSTANTS_OBJ) $(PNFAM_SETUP_OBJ) \
                               $(TYPE_BBM_OBJ) $(PNFAM_BROYDEN_OBJ) $(PNFAM_TXTOUTPUT_OBJ) \
                               $(TYPE_EXTFIELD_OBJ) $(PNFAM_EXTFIELD_OBJ) $(PNFAM_STORAGE_OBJ) \
                               $(TYPE_GAMDEL_2BC_OBJ)
$(CONTOUR_SETUP_OBJ)     : $(PNFAM_LOGGER_OBJ) $(PNFAM_SETUP_OBJ) $(PNFAM_TXTOUTPUT_OBJ)


$(PNFAM_EXTFIELD_2BC_OBJ)  : $(TYPE_BLOCKMATRIX_OBJ) $(TYPE_EXTFIELD_OBJ) $(TYPE_EXTFIELD_2BC_OBJ) \
                                 $(HFB_SOLUTION_OBJ) $(PNFAM_LOGGER_OBJ) $(TYPE_GAMDEL_2BC_OBJ)
$(TYPE_EXTFIELD_2BC_OBJ)   : $(PNFAM_CONSTANTS_OBJ) $(PNFAM_SPATIAL_MTXELS_OBJ) $(PNFAM_SETUP_OBJ) \
                                 $(PNFAM_LOGGER_OBJ)
$(TYPE_GAMDEL_2BC_OBJ)     :
$(PNFAM_SPATIAL_MTXELS_OBJ): $(HFB_SOLUTION_OBJ)


$(PNFAM_OBJ)             : $(PNFAM_SOLVER_OBJ)
$(CONTOUR_OBJ)           : $(CONTOUR_SETUP_OBJ) $(PNFAM_SOLVER_OBJ)

#==============================================================================#
#                              Program Builds                                  #
#==============================================================================#
GLOBAL_OBJECTS = \
    $(HFB_INTERFACE_OBJ) \
    $(HFB_SOLUTION_OBJ) \
    $(TYPE_DENSITYSET_OBJ) \
    $(TYPE_EXTFIELD_OBJ) \
    $(TYPE_BLOCKMATRIX_OBJ) \
    $(TYPE_BBM_OBJ) \
    $(PNFAM_LOGGER_OBJ) \
    $(PNFAM_CONSTANTS_OBJ) \
    $(PNFAM_BROYDEN_OBJ) \
    $(PNFAM_EXTFIELD_OBJ) \
    $(PNFAM_INTERACTION_OBJ) \
    $(PNFAM_HAMILTONIAN_OBJ) \
    $(PNFAM_SETUP_OBJ) \
    $(PNFAM_TXTOUTPUT_OBJ) \
    $(PNFAM_STORAGE_OBJ) \
    $(PNFAM_SOLVER_OBJ) \
    $(PNFAM_EXTFIELD_2BC_OBJ) \
    $(TYPE_EXTFIELD_2BC_OBJ) \
    $(TYPE_GAMDEL_2BC_OBJ) \
    $(PNFAM_SPATIAL_MTXELS_OBJ)


CTR_OBJECTS =\
    $(CONTOUR_SETUP_OBJ)


LIBRARIES = $(L_HFBTHO) $(L_LINEARALGEBRA)
INCLUDE   = $(I_HFBTHO)

$(PNFAM_EXE) : $(PNFAM_OBJ) $(GLOBAL_OBJECTS)
ifeq ($(VERBOSE),1)
	$(COMPILER) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(LIBRARIES)
else
	@echo "Assembling \"$@\"..."
	@$(COMPILER) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(LIBRARIES)
endif

$(PNFAM_LIB): $(GLOBAL_OBJECTS)
ifeq ($(VERBOSE),1)
	ar rcs $@ $(GLOBAL_OBJECTS)
else
	@echo "Building library \"$@\"..."
	@ar rcs $@ $(GLOBAL_OBJECTS)
endif

$(CONTOUR_EXE) : $(CONTOUR_OBJ) $(CTR_OBJECTS) $(GLOBAL_OBJECTS)
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


tbc_test.x : $(GLOBAL_OBJECTS) tbc_test.o
	$(COMPILER) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(LIBRARIES)
