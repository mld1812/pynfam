all: anl_unedf0 anl_unedf1 anl_unedf2

VPATH = ..

HFBTHO_OBJS = $(HFBTHO_DIR)/hfbtho_bessel.o  $(HFBTHO_DIR)/hfbtho_unedf.o \
              $(HFBTHO_DIR)/hfbtho_elliptic_integrals.o  $(HFBTHO_DIR)/hfbtho_utilities.o \
              $(HFBTHO_DIR)/hfbtho_multipole_moments.o  $(HFBTHO_DIR)/hfbtho_collective.o \
              $(HFBTHO_DIR)/hfbtho_fission.o \
              $(HFBTHO_DIR)/hfbtho_gauss.o $(HFBTHO_DIR)/hfbtho_solver.o \
              $(HFBTHO_DIR)/hfbtho_variables.o $(HFBTHO_DIR)/hfbtho_linear_algebra.o  \
              $(HFBTHO_DIR)/hfbtho_tho.o $(HFBTHO_DIR)/hfbtho_version.o

anl_unedf0: anl_unedf0.o
	${FORTRAN_MPI} ${FORMAT_F90} ${OPTIONS} -o $@ $^ anl_utils.o ${HFBTHO_OBJS} ${LINEAR_ALGEBRA}
anl_unedf1: anl_unedf1.o
	${FORTRAN_MPI} ${FORMAT_F90} ${OPTIONS} -o $@ $^ anl_utils.o ${HFBTHO_OBJS} ${LINEAR_ALGEBRA}
anl_unedf2: anl_unedf2.o
	${FORTRAN_MPI} ${FORMAT_F90} ${OPTIONS} -o $@ $^ anl_utils.o ${HFBTHO_OBJS} ${LINEAR_ALGEBRA}
anl_unedf1_argonne: anl_unedf1_argonne.o
	${FORTRAN_MPI} ${FORMAT_F90} ${OPTIONS} -o $@ $^ anl_utils.o ${HFBTHO_OBJS} ${LINEAR_ALGEBRA}
anl_unedf1_argonne_e6: anl_unedf1_argonne.o 
	${FORTRAN_MPI} ${FORMAT_F90} ${OPTIONS} -o $(HFBTHO_DIR)/$(HFBTHO_ARCH)/bin/$@ $^ anl_utils.o ${HFBTHO_OBJS} ${LINEAR_ALGEBRA}
anl_unedf1_single: anl_unedf1_single.o 
	${FORTRAN_MPI} ${FORMAT_F90} ${OPTIONS} -o $(HFBTHO_DIR)/$(HFBTHO_ARCH)/bin/$@ $^ anl_utils.o ${HFBTHO_OBJS} ${LINEAR_ALGEBRA}

anl_utils.o: anl_utils.f90
	${FORTRAN_MPI} -I$(HFBTHO_DIR) ${FORMAT_F90} ${OPTIONS} -c $<

anl_unedf0.o: anl_unedf0.f90 anl_utils.o
	${FORTRAN_MPI} -I$(HFBTHO_DIR) ${FORMAT_F90} ${OPTIONS} -c $<

anl_unedf1.o: anl_unedf1.f90 anl_utils.o
	${FORTRAN_MPI} -I$(HFBTHO_DIR) ${FORMAT_F90} ${OPTIONS} -c $<
anl_unedf2.o: anl_unedf2.f90 anl_utils.o
	${FORTRAN_MPI} -I$(HFBTHO_DIR) ${FORMAT_F90} ${OPTIONS} -c $<

anl_unedf1_argonne.o: anl_unedf1_argonne.f90 anl_utils.o
	${FORTRAN_MPI} -I$(HFBTHO_DIR) ${FORMAT_F90} ${OPTIONS} -c $<
anl_unedf1_single.o: anl_unedf1_single.f90 anl_utils.o
	${FORTRAN_MPI} -I.. ${FORMAT_F90} ${OPTIONS} -c $<

# Cleaning
clean ::
	-rm -f *.o *.oo *.ipo *.mod
