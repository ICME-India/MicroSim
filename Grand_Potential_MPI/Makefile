FDIR=functions
SDIR=solverloop
TDIR=tdbs
CC=h5pcc #-O3
# CFLAGS=-I.
#LDIR=$(GSL_ROOT_DIR)/lib
#IDIR=$(GSL_ROOT_DIR)/include

LDIR=/share/apps/gsl/lib
IDIR=/share/apps/gsl/include
CFLAGS=-I $(IDIR)

_FUNCTION_DEPS = global_vars.h functions.h matrix.h utility_functions.h functionH.h functionF_01.h functionF_02.h functionF_03.h functionF_04.h functionF_05.h  functionF_06.h functionQ.h \
                 functionF_elast.h functionW_01.h functionW_02.h function_A_00.h function_A_01.h anisotropy_01.h functionTau.h \
                 functionD.h filling.h reading_input_parameters.h reading_input_parameters_new.h dynamic_input_parameters.h read_boundary_conditions.h print_boundary_conditions.h print_input_parameters.h \
                 initialize_variables.h free_variables.h fill_domain.h shift.h Temperature_gradient.h

DEPS = $(patsubst %,$(FDIR)/%,$(_FUNCTION_DEPS))

_SOLVERLOOP_DEPS = serialinfo_xy.h mpiinfo_xyz.h gradients.h simplex_projection.h calculate_gradients.h \
		   calculate_fluxes_concentration.h calculate_divergence_phasefield.h calculate_divergence_concentration.h calculate_divergence_stress.h calculate_LBM_distribution.h calculate_LBM_distribution3D.h\
		   file_writer.h file_writer_3D.h initialize_functions_solverloop.h solverloop.h boundary_mpi.h


DEPS += $(patsubst %,$(SDIR)/%,$(_SOLVERLOOP_DEPS))

_tdb_DEPS = Thermo.c Thermo.h

DEPS += $(patsubst %,$(TDIR)/%,$(_tdb_DEPS))

LIBS =-L $(LDIR) -lm -lgsl -lgslcblas

all : microsim_gp

# echo "Making microsim_gp solver (grand-potential based phase-field model)"

microsim_gp : microsim_gp.o
	$(CC) -o microsim_gp microsim_gp.o $(CFLAGS) $(LIBS) -Wall

microsim_gp.o : $(DEPS)

_HEADERS_XDMF_WRITER = global_vars.h functions.h matrix.h utility_functions.h reading_input_parameters.h

_HEADERS_RECONSTRUCT = global_vars.h functions.h matrix.h utility_functions.h reading_input_parameters.h

DEPS_XDMF_WRITER = $(patsubst %,$(FDIR)/%,$(_HEADERS_XDMF_WRITER))

DEPS_RECONSTRUCT = $(patsubst %,$(FDIR)/%,$(_HEADERS_RECONSTRUCT))


# echo "Making xdmf file write for .h5 files. Required for viewing in paraview"

write_xdmf : write_xdmf.o
	$(CC) -o write_xdmf write_xdmf.o $(CFLAGS) $(LIBS)

write_xdmf.o : $(DEPS_XDMF_WRITER)

# echo "Making reconstruct for collating different processor files into consolidated .vtk files. Valid for non .h5 files"

run: microsim_gp
	mpirun -np 2 ./microsim_gp Inp.md Fill.in outpuut 1 2 1

.PHONY : clean

clean :
	rm microsim_gp.o microsim_gp
