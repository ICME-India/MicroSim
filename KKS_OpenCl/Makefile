SHELL ?= /bin/sh

# Toolchains (Spack injects wrappers)
CC ?= mpicc
HC ?= h5pcc

# Directories
FDIR ?= functions
SDIR ?= solverloop
TDIR ?= tdbs

# Flags (do NOT overwrite – Spack extends these)
CFLAGS   ?=
CPPFLAGS ?=
LDFLAGS  ?=
LIBS     ?=

# Local includes
CPPFLAGS += -I.

# Libraries (Spack provides correct paths)
LIBS += -lm -lgsl -lgslcblas -lOpenCL

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------

_FUNCTION_DEPS = global_vars.h CL_global_vars.h CL_initialize_variables.h \
                 CL_device_kernel_build.h CL_buffer_allocation.h \
                 CL_create_kernel_args.h CL_Initialize_domain.h \
                 CL_kernel_init_temperature.h functions.h matrix.h \
                 utility_functions.h functionH.h functionF_01.h \
                 functionF_02.h functionF_03.h CL_functionF_03_HostUpdate.h \
                 functionF_04.h CL_functionF_04_HostUpdate.h \
                 initialize_variables.h filling.h \
                 reading_input_parameters.h read_boundary_conditions.h \
                 free_variables.h fill_domain.h CL_DeviceToHost.h \
                 Temperature_gradient.h shift.h CL_Shift.h

DEPS = $(patsubst %,$(FDIR)/%,$(_FUNCTION_DEPS))

_SOLVERLOOP_DEPS = serialinfo_xy.h file_writer.h \
                   initialize_functions_solverloop.h mpi_xy.h \
                   CL_Solve_phi_com_Function_F_2.h \
                   CL_Solve_phi_com_Function_F_3.h \
                   CL_Solve_phi_com_Function_F_4.h \
                   CL_Solve_phi_com.h CL_Update_Temperature.h \
                   defines.h defines1.h Function_Thermo_calls.h \
                   FunctionDefines_CL.h FunctionF_2_c_mu_NR_GEData2.h \
                   functionH_CL_5th.h functionW_02_CL.h \
                   GibbsEnergyData2.h matrix_nsolCL.h \
                   FunctionF_4_SplineCPU.h CL_Kim_Kernel.cl \
                   ThermoCL_2.c ThermoCL_2.h

DEPS += $(patsubst %,$(SDIR)/%,$(_SOLVERLOOP_DEPS))

_TDB_DEPS = Thermo.c Thermo.h
DEPS += $(patsubst %,$(TDIR)/%,$(_TDB_DEPS))

# -----------------------------------------------------------------------------
# Targets
# -----------------------------------------------------------------------------

TARGET = microsim_kks_opencl
RECONSTRUCT = reconstruct

all: $(TARGET) $(RECONSTRUCT)

$(TARGET): microsim_kks_opencl.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

microsim_kks_opencl.o: $(DEPS)

# -----------------------------------------------------------------------------
# Reconstruct utility
# -----------------------------------------------------------------------------

_HEADERS_RECONSTRUCT = global_vars.h CL_global_vars.h functions.h matrix.h \
                       utility_functions.h reading_input_parameters.h

DEPS_RECONSTRUCT = $(patsubst %,$(FDIR)/%,$(_HEADERS_RECONSTRUCT))

$(RECONSTRUCT): reconstruct.o
	$(HC) $(LDFLAGS) -o $@ $^ $(LIBS)

reconstruct.o: $(DEPS_RECONSTRUCT)

# -----------------------------------------------------------------------------
# Install (required by Spack)
# -----------------------------------------------------------------------------

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin

install: $(TARGET) $(RECONSTRUCT)
	mkdir -p $(DESTDIR)$(BINDIR)
	cp $(TARGET) $(RECONSTRUCT) $(DESTDIR)$(BINDIR)/

# -----------------------------------------------------------------------------
# Cleanup
# -----------------------------------------------------------------------------

clean:
	rm -f *.o $(TARGET) $(RECONSTRUCT)

distclean: clean
	rm -f *.out *.dat *.gz *.vtk *.txt *.log *.err *.png *.jpg *.jpeg

.PHONY: all install clean distclean
