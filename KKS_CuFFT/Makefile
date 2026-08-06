SHELL ?= /bin/sh

# Compilers (Spack injects these)
NVCC ?= nvcc

# Flags (Spack-friendly)
NVCCFLAGS ?=
CUDA_ARCH ?=
CPPFLAGS  ?=
LDFLAGS   ?=
LIBS      ?=

# Optimisation (can be overridden by Spack)
NVCCFLAGS += -Xptxas -O3

# CUDA architecture (Spack can pass this)
CUDA_ARCH ?= -arch=sm_70

# CUDA libraries (Spack provides correct paths via wrappers)
LIBS += -lcufft -lcudart -lcuda

# Directories
SRCDIR ?= .
BINDIR ?= DATA

# Sources
SRCFILES = $(SRCDIR)/KKS_cuFFT.cu
HEADERS  = functions/global_vars.h

TARGET = microsim_kks_cufft

all: build

build: $(TARGET)

$(TARGET): $(SRCFILES) $(HEADERS)
	$(NVCC) $(NVCCFLAGS) $(CUDA_ARCH) $(CPPFLAGS) $< -o $@ $(LDFLAGS) $(LIBS)

# Install support for Spack
PREFIX ?= /usr/local
BINDIR_INSTALL ?= $(PREFIX)/bin

install: $(TARGET)
	mkdir -p $(DESTDIR)$(BINDIR_INSTALL)
	cp $(TARGET) $(DESTDIR)$(BINDIR_INSTALL)/

# Optional runtime data directory
build_dir:
	mkdir -p $(BINDIR)

clean:
	rm -f $(TARGET)

distclean: clean
	rm -rf $(BINDIR)

.PHONY: all build install clean distclean
