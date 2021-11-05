FDIR=functions
SDIR=solverloop
CC=gcc
CFLAGS=-I.

IDIR =-I/usr/local/cuda/include/

LDIR =-L/usr/local/cuda/lib64/

_FUNCTION_DEPS = global_vars.h functions.h matrix.h utility_functions.h  \
                 filling.h reading_input_parameters.h read_boundary_conditions.h  \
                 free_variables.h fill_domain.h
                 
DEPS = $(patsubst %,$(FDIR)/%,$(_FUNCTION_DEPS))

_SOLVERLOOP_DEPS = serialinfo_xy.h  file_writer.h

DEPS += $(patsubst %,$(SDIR)/%,$(_SOLVERLOOP_DEPS))

LIBS =-lm -lOpenCL

microsim_kks_opencl : microsim_kks_opencl.o 
	$(CC)  -o microsim_kks_opencl microsim_kks_opencl.o $(CFLAGS)  $(IDIR) $(LDIR) $(LIBS)

microsim_kks_opencl.o : $(DEPS)       

.PHONY : clean

clean :

	rm -f microsim_kks_opencl.o microsim_kks_opencl *.out *.o  *.dat *.gz *.vtk *.txt *.mv *.log *.err *.epsc *.eps *.mat  *.png *.op *.jpg *.jpeg \
	rm -rf DATA
