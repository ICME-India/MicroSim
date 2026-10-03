SHELL ?= /bin/sh

# Toolchain (Spack injects compiler wrapper)
CC ?= cc

# Directories
FDIR ?= functions
SDIR ?= solverloop

# Flags (do not overwrite Spack flags)
CFLAGS   ?=
CPPFLAGS ?=
LDFLAGS  ?=
LIBS     ?=

# Local includes
CPPFLAGS += -I.

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------

_FUNCTION_DEPS = global_vars.h functions.h matrix.h utility_functions.h \
                 filling.h reading_input_parameters.h free_variables.h \
                 fill_domain.h

DEPS = $(patsubst %,$(FDIR)/%,$(_FUNCTION_DEPS))

_SOLVERLOOP_DEPS = serialinfo_xy.h GibbsEnergyData.h \
                   functions_fftw.h file_writer.h

DEPS += $(patsubst %,$(SDIR)/%,$(_SOLVERLOOP_DEPS))

# ------------------------------------------------------------------------------
# Libraries (FFTW provided by Spack)
# ------------------------------------------------------------------------------
LIBS += -lfftw3 -lm

# ------------------------------------------------------------------------------
# Targets
# ------------------------------------------------------------------------------
TARGET = microsim_ch_fft

all: $(TARGET)

$(TARGET): microsim_ch_fft.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

microsim_ch_fft.o: $(DEPS)

# ------------------------------------------------------------------------------
# Install (required by Spack)
# ------------------------------------------------------------------------------
PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin

install: $(TARGET)
	mkdir -p $(DESTDIR)$(BINDIR)
	cp $(TARGET) $(DESTDIR)$(BINDIR)/

# ------------------------------------------------------------------------------
# Cleanup
# ------------------------------------------------------------------------------
clean:
	rm -f microsim_ch_fft.o $(TARGET)

.PHONY: all install clean
