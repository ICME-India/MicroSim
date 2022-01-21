FDIR=functions
SDIR=solverloop
CC=h5pcc
CFLAGS=-I.

_FUNCTION_DEPS = global_vars.h functions.h matrix.h utility_functions.h functionH.h functionF_01.h functionQ.h \
                 functionW_01.h functionW_02.h function_A_00.h function_A_01.h anisotropy_01.h functionTau.h \
                 functionD.h filling.h reading_input_parameters.h read_boundary_conditions.h print_boundary_conditions.h print_input_parameters.h \
                 initialize_variables.h free_variables.h fill_domain.h shift.h Temperature_gradient.h
                 
DEPS = $(patsubst %,$(FDIR)/%,$(_FUNCTION_DEPS))

_SOLVERLOOP_DEPS = serialinfo_xy.h mpiinfo_xy.h gradients.h simplex_projection.h calculate_gradients.h \
		   calculate_fluxes_concentration.h calculate_divergence_phasefield.h calculate_divergence_concentration.h \
		   initialize_functions_solverloop.h solverloop.h boundary_mpi.h \
		   file_writer.h

DEPS += $(patsubst %,$(SDIR)/%,$(_SOLVERLOOP_DEPS))


LIBS =-lm

all : microsim_gp write_xdmf reconstruct 


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

reconstruct : reconstruct.o
	$(CC) -o reconstruct reconstruct.o $(CFLAGS) $(LIBS)
	
reconstruct.o : $(DEPS_RECONSTRUCT)


.PHONY : clean

clean :

	rm microsim_gp.o microsim_gp reconstruct.o reconstruct write_xdmf.o write_xdmf


