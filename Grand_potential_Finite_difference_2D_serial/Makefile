FDIR=functions
SDIR=solverloop
CC=gcc
CFLAGS=-I.

_FUNCTION_DEPS = global_vars.h functions.h matrix.h utility_functions.h functionH.h functionF_01.h functionQ.h \
                 functionW_01.h functionW_02.h function_A_00.h function_A_01.h anisotropy_01.h functionTau.h \
                 functionD.h filling.h reading_input_parameters.h read_boundary_conditions.h initialize_variables.h \
                 free_variables.h fill_domain.h shift.h Temperature_gradient.h
                 
DEPS = $(patsubst %,$(FDIR)/%,$(_FUNCTION_DEPS))

_SOLVERLOOP_DEPS = serialinfo_xy.h gradients.h simplex_projection.h calculate_gradients.h \
		   calculate_fluxes_concentration.h calculate_divergence_phasefield.h calculate_divergence_concentration.h \
		   initialize_functions_solverloop.h solverloop.h boundary_serial.h \
		   file_writer.h

DEPS += $(patsubst %,$(SDIR)/%,$(_SOLVERLOOP_DEPS))


LIBS =-lm

microsim_gp : microsim_gp.o 
	$(CC) -o microsim_gp microsim_gp.o $(CFLAGS) $(LIBS)

microsim_gp.o : $(DEPS)       

.PHONY : clean

clean :

	rm microsim_gp.o microsim_gp


