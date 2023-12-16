################################################################################
# Makefile
# 
# Build the matrix FAM solver. There is an `export` command in the main PNFAM
# Makefile which sets the variables needed here. You shouldn't make anything
# from this directory. Run `make mfam` in the main PNFAM directory to build
# these tests. 
#
# M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
################################################################################

all: pnfam_matrix_ab.x pnfam_matrix_oneqp.x pnfam_matrix_twoqp.x pnfam_matrix_odd_me.x pnfam_matrix_truncate.x pnfam_matrix_basis.x pnfam_matrix_odd_info.x pnfam_matrix_beta.x pnfam_matrix_diag.x 

clean:
	@-rm -vf ./*.mod ./*.o

cleanall: clean
	@-rm -vf ./*.x

# Fortran source
%.o : %.f90
	$(FC) -I.. -c -o $@ $< $(FFLAGS)

# Targets
pnfam_matrix_ab.x: ../constants.o ../blockmatrix_type.o ../logger.o ../hfbtho_basis.o type_hfb_solution.o ../extfield_type.o ../extfield.o ../interaction.o ../density_set_type.o ../hamiltonian.o type_config.o mod_common.o main_abmatrix.o
	$(FC) -o $@ $^ $(FFLAGS) $(OPENMP) $(LAPACK)

pnfam_matrix_diag.x: ../constants.o ../blockmatrix_type.o ../logger.o ../hfbtho_basis.o ../extfield_type.o ../extfield.o type_hfb_solution.o ../density_set_type.o ../interaction.o ../hamiltonian.o type_config.o mod_common.o main_diagonalize.o
	$(FC) -o $@ $^ $(FFLAGS) $(OPENMP) $(LAPACK) $(SCALAPACK)

pnfam_matrix_twoqp.x: ../constants.o ../blockmatrix_type.o ../logger.o ../hfbtho_basis.o ../extfield_type.o ../extfield.o type_hfb_solution.o ../density_set_type.o ../interaction.o ../hamiltonian.o type_config.o mod_common.o main_transitions_twoqp.o
	$(FC) -o $@ $^ $(FFLAGS) $(OPENMP) $(LAPACK)

pnfam_matrix_oneqp.x: ../constants.o ../blockmatrix_type.o ../extfield_type.o ../logger.o ../hfbtho_basis.o ../extfield.o type_hfb_solution.o ../density_set_type.o ../interaction.o ../hamiltonian.o type_config.o mod_common.o main_transitions_oneqp.o
	$(FC) -o $@ $^ $(FFLAGS) $(OPENMP) $(LAPACK)

pnfam_matrix_odd_me.x: ../logger.o ../blockmatrix_type.o ../constants.o ../hfbtho_basis.o ../extfield_type.o ../extfield.o type_hfb_solution.o ../density_set_type.o ../interaction.o ../hamiltonian.o type_config.o mod_common.o main_odd_matrix_elements.o
	$(FC) -o $@ $^ $(FFLAGS) $(OPENMP) $(LAPACK)

pnfam_matrix_truncate.x: ../constants.o ../blockmatrix_type.o ../logger.o ../hfbtho_basis.o ../extfield_type.o ../extfield.o type_hfb_solution.o ../density_set_type.o ../interaction.o ../hamiltonian.o type_config.o mod_common.o aux_basis_cutdown.o
	$(FC) -o $@ $^ $(FFLAGS) $(OPENMP) $(LAPACK)

pnfam_matrix_basis.x: ../constants.o ../logger.o ../blockmatrix_type.o ../hfbtho_basis.o ../extfield_type.o ../extfield.o type_hfb_solution.o ../density_set_type.o ../interaction.o ../hamiltonian.o type_config.o mod_common.o aux_basis_sizes.o
	$(FC) -o $@ $^ $(FFLAGS) $(OPENMP) $(LAPACK)

pnfam_matrix_odd_info.x: ../constants.o ../logger.o ../blockmatrix_type.o ../hfbtho_basis.o ../interaction.o ../density_set_type.o ../hamiltonian.o ../extfield_type.o ../extfield.o type_hfb_solution.o type_config.o mod_common.o aux_odd_information.o
	$(FC) -o $@ $^ $(FFLAGS) $(OPENMP) $(LAPACK)

pnfam_matrix_beta.x: ../constants.o ../logger.o ../blockmatrix_type.o ../hfbtho_basis.o type_hfb_solution.o ../complex_quadrature.o ../fermi.o ../polyfit.o ../rational_interp.o ../phasespace.o ../extfield_type.o ../extfield.o ../density_set_type.o ../interaction.o ../hamiltonian.o type_config.o mod_common.o main_betadecay.o
	$(FC) -o $@ $^ $(FFLAGS) $(OPENMP) $(LAPACK) $(GSL)
