# Directories
FDIR ?= functions
SDIR ?= solverloop
TDIR ?= tdbs

# Flags (Spack controls these)
CFLAGS   ?=
CPPFLAGS ?=
LDFLAGS  ?=
LIBS     ?=

# Local includes
CPPFLAGS += -I.

# Libraries (resolved by Spack)
LIBS += -lm -lgsl -lgslcblas

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------

_FUNCTION_DEPS = global_vars.h functions.h matrix.h utility_functions.h \
                 functionH.h functionF_01.h functionF_02.h functionF_03.h \
                 functionF_04.h functionF_05.h functionF_06.h functionQ.h \
                 functionF_elast.h functionW_01.h functionW_02.h \
                 function_A_00.h function_A_01.h anisotropy_01.h \
                 functionTau.h functionD.h filling.h \
                 reading_input_parameters.h reading_input_parameters_new.h \
                 dynamic_input_parameters.h read_boundary_conditions.h \
                 print_boundary_conditions.h print_input_parameters.h \
                 initialize_variables.h free_variables.h fill_domain.h \
                 shift.h Temperature_gradient.h

DEPS = $(patsubst %,$(FDIR)/%,$(_FUNCTION_DEPS))

_SOLVERLOOP_DEPS = serialinfo_xy.h mpiinfo_xyz.h gradients.h \
                   simplex_projection.h calculate_gradients.h \
                   calculate_fluxes_concentration.h \
                   calculate_divergence_phasefield.h \
                   calculate_divergence_concentration.h \
                   calculate_divergence_stress.h \
                   calculate_LBM_distribution.h \
                   calculate_LBM_distribution3D.h \
                   file_writer.h file_writer_3D.h \
                   initialize_functions_solverloop.h \
                   solverloop.h boundary_mpi.h

DEPS += $(patsubst %,$(SDIR)/%,$(_SOLVERLOOP_DEPS))

_TDB_DEPS = Thermo.c Thermo.h
DEPS += $(patsubst %,$(TDIR)/%,$(_TDB_DEPS))

# -----------------------------------------------------------------------------
# Main target
# -----------------------------------------------------------------------------

TARGET = microsim_gp

all: $(TARGET)

$(TARGET): microsim_gp.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

microsim_gp.o: $(DEPS)

# -----------------------------------------------------------------------------
# XDMF writer utility
# -----------------------------------------------------------------------------

XDMF = write_xdmf

_HEADERS_XDMF = global_vars.h functions.h matrix.h \
                utility_functions.h reading_input_parameters.h

DEPS_XDMF = $(patsubst %,$(FDIR)/%,$(_HEADERS_XDMF))

$(XDMF): write_xdmf.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

write_xdmf.o: $(DEPS_XDMF)

# -----------------------------------------------------------------------------
# Install (required for Spack)
# -----------------------------------------------------------------------------

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin

install: $(TARGET) $(XDMF)
	mkdir -p $(DESTDIR)$(BINDIR)
	cp $(TARGET) $(XDMF) $(DESTDIR)$(BINDIR)/

# -----------------------------------------------------------------------------
# Run helper (not used by Spack, but kept)
# -----------------------------------------------------------------------------

run: $(TARGET)
	mpirun -np 2 ./$(TARGET) Inp.md Fill.in output 1 2 1

# -----------------------------------------------------------------------------
# Cleanup
# -----------------------------------------------------------------------------

clean:
	rm -f *.o $(TARGET) $(XDMF)

distclean: clean
	rm -f *.out *.dat *.vtk *.log *.err

.PHONY: all install run clean distclean
