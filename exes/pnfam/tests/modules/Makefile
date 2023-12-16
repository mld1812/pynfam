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
	@-rm -vf *.mod *.o

cleanall: clean
	@-rm -vf test_blockmatrix_type.x test_complex_quadrature.x test_fermi.x test_phasespace.x

# Dependency names
COMPLEX_QUADRATURE_OBJ = ../../complex_quadrature.o
BLOCKMATRIX_TYPE_OBJ   = ../../blockmatrix_type.o
FERMI_OBJ              = ../../fermi.o
CONSTANTS_OBJ          = ../../constants.o
PHASESPACE_OBJ         = ../../phasespace.o
RATIONAL_INTERP_OBJ    = ../../rational_interp.o
POLYNOMIAL_FIT_OBJ     = ../../polyfit.o

# Test names
blockmatrix_type_test_obj   = blockmatrix_type_test.o
complex_quadrature_test_obj = complex_quadrature_test.o
fermi_test_obj              = fermi_test.o
phasespace_test_obj         = phasespace_test.o
pnfam_tests_obj             = pnfam_tests.o

# "Global" object files for each test (should already be compiled via main Makefile)
blockmatrix_type_test_objs   = $(blockmatrix_type_test_obj) \
                               $(BLOCKMATRIX_TYPE_OBJ)
complex_quadrature_test_objs = $(complex_quadrature_test_obj) \
                               $(COMPLEX_QUADRATURE_OBJ)
fermi_test_objs              = $(fermi_test_obj) $(pnfam_tests_obj) \
                               $(FERMI_OBJ) $(CONSTANTS_OBJ)
phasespace_test_objs         = $(phasespace_test_obj) $(pnfam_tests_obj) \
                               $(PHASESPACE_OBJ) $(RATIONAL_INTERP_OBJ) \
                               $(COMPLEX_QUADRATURE_OBJ) $(FERMI_OBJ) \
                               $(POLYNOMIAL_FIT_OBJ)

# Local test dependencies (other dependencies already compiled)
#(blockmatrix_type_test_obj) :
#(complex_quadrature_test_obj) :
$(fermi_test_obj) : $(pnfam_tests_obj)
$(phasespace_test_obj) : $(pnfam_tests_obj)

# Targets
test_blockmatrix_type.x : $(blockmatrix_type_test_objs)
ifeq ($(VERBOSE),1)
	$(COMPILER) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(L_LINEARALGEBRA)
else
	@echo "Assembling \"$@\"..."
	@$(COMPILER) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(L_LINEARALGEBRA)
endif

test_complex_quadrature.x : $(complex_quadrature_test_objs)
ifeq ($(VERBOSE),1)
	$(COMPILER) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(L_LINEARALGEBRA)
else
	@echo "Assembling \"$@\"..."
	@$(COMPILER) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(L_LINEARALGEBRA)
endif

test_fermi.x : $(fermi_test_objs)
ifeq ($(VERBOSE),1)
	$(COMPILER) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(L_LINEARALGEBRA)
else
	@echo "Assembling \"$@\"..."
	@$(COMPILER) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(L_LINEARALGEBRA)
endif

test_phasespace.x : $(phasespace_test_objs)
ifeq ($(VERBOSE),1)
	$(COMPILER) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(L_LINEARALGEBRA) $(L_GSL)
else
	@echo "Assembling \"$@\"..."
	@$(COMPILER) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(L_LINEARALGEBRA) $(L_GSL)
endif

# Compile object files (header/mod files need to be included here)
%.o : %.f90
ifeq ($(VERBOSE),1)
	$(COMPILER) -c $< $(FORMAT_F90) $(OPTIONS) -I../../
else
	@echo "Compiling  \"$@\"..."
	@$(COMPILER) -c $< $(FORMAT_F90) $(OPTIONS) -I../../
endif
