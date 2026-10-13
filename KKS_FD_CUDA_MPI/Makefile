SHELL ?= /bin/sh
NPROCS ?= 1

# ------------------------------------------------------------------------------
# Toolchains (Spack injects wrappers)
# ------------------------------------------------------------------------------
MPICC ?= mpicc
MPILD ?= mpicc
H5CC  ?= h5pcc
NVCC  ?= nvcc

# ------------------------------------------------------------------------------
# Feature toggles (Spack variants)
# ------------------------------------------------------------------------------
ENABLE_HDF5     ?= 1
ENABLE_CUFFTMP  ?= 1
ENABLE_NVSHMEM  ?= 1

# ------------------------------------------------------------------------------
# Directories
# ------------------------------------------------------------------------------
SOLVERLOOP ?= solverloop
FUNCTIONS  ?= functions
TDBS       ?= tdbs

# ------------------------------------------------------------------------------
# Flags (Spack-controlled)
# ------------------------------------------------------------------------------
CPPFLAGS  ?=
CFLAGS    ?=
CXXFLAGS  ?=
NVCCFLAGS ?=
LDFLAGS   ?=
LIBS      ?=

CPPFLAGS += -I$(SOLVERLOOP) -I$(FUNCTIONS) -I$(TDBS)

# ------------------------------------------------------------------------------
# CUDA architecture (Spack sets cuda_arch)
# ------------------------------------------------------------------------------
CUDA_ARCH ?= sm_70
NVCCFLAGS += -Xcompiler -fPIC -gencode arch=compute_$(CUDA_ARCH),code=sm_$(CUDA_ARCH)

# ------------------------------------------------------------------------------
# Libraries (Spack resolves paths)
# ------------------------------------------------------------------------------
LIBS += -lm -lgsl -lgslcblas

ifeq ($(ENABLE_HDF5),1)
  CPPFLAGS += -DENABLE_HDF5=1
  LIBS     += -lhdf5
endif

ifeq ($(ENABLE_CUFFTMP),1)
  CPPFLAGS += -DENABLE_CUFFTMP=1
  LIBS     += -lcufftMp
else
  LIBS     += -lcufft
endif

ifeq ($(ENABLE_NVSHMEM),1)
  CPPFLAGS += -DENABLE_NVSHMEM=1
  LIBS     += -lnvshmem
endif

LIBS += -lcudart

# ------------------------------------------------------------------------------
# Output
# ------------------------------------------------------------------------------
TARGET = microsim_kks_fd_cuda_mpi

# ------------------------------------------------------------------------------
# Objects
# ------------------------------------------------------------------------------
HOST_OBJS = \
  $(SOLVERLOOP)/obj/fileWriter.o \
  $(FUNCTIONS)/obj/utilityFunctions.o \
  $(FUNCTIONS)/obj/initialize_variables.o \
  $(FUNCTIONS)/obj/filling.o \
  $(FUNCTIONS)/obj/reading_input_parameters.o \
  $(FUNCTIONS)/obj/fill_domain.o \
  $(FUNCTIONS)/obj/read_boundary_conditions.o \
  $(FUNCTIONS)/obj/calc_bn.o \
  MicroSim_KKS_FD_CUDA_MPI.o

DEVICE_OBJS = \
  $(TDBS)/obj/Thermo.o \
  $(FUNCTIONS)/obj/functionF.o \
  $(FUNCTIONS)/obj/functionH.o \
  $(FUNCTIONS)/obj/functionTau.o \
  $(FUNCTIONS)/obj/functionW_02.o \
  $(FUNCTIONS)/obj/anisotropy_01.o \
  $(FUNCTIONS)/obj/functionA_01.o \
  $(FUNCTIONS)/obj/functionA_02.o \
  $(FUNCTIONS)/obj/matrix.o \
  $(SOLVERLOOP)/obj/utilityKernels.o \
  $(SOLVERLOOP)/obj/boundary.o \
  $(SOLVERLOOP)/obj/calcPhaseComp.o \
  $(SOLVERLOOP)/obj/smooth.o \
  $(SOLVERLOOP)/obj/computeDrivingForce.o \
  $(SOLVERLOOP)/obj/computeElastic.o \
  $(SOLVERLOOP)/obj/updateComposition.o \
  $(SOLVERLOOP)/obj/updatePhi.o \
  $(FUNCTIONS)/obj/error_checks.o

# ------------------------------------------------------------------------------
# Build rules
# ------------------------------------------------------------------------------
all: $(TARGET)

%.o: %.cu
	$(NVCC) -dc $< -o $@ $(CPPFLAGS) $(NVCCFLAGS)

%.o: %.c
	$(MPICC) -c $< -o $@ $(CPPFLAGS) $(CFLAGS)

MicroSim_KKS_FD_CUDA_MPI.o: microsim_kks_fd_cuda_mpi.c
	$(MPICC) -c $< -o $@ $(CPPFLAGS) $(CFLAGS)

# Device link
device_link.o: $(DEVICE_OBJS)
	$(NVCC) -dlink $^ -o $@

# Final link
$(TARGET): $(HOST_OBJS) $(DEVICE_OBJS) device_link.o
	$(H5CC) $^ $(LDFLAGS) $(LIBS) -o $@

# ------------------------------------------------------------------------------
# Install (Spack requirement)
# ------------------------------------------------------------------------------
PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin

install: $(TARGET)
	mkdir -p $(DESTDIR)$(BINDIR)
	cp $(TARGET) $(DESTDIR)$(BINDIR)/

# ------------------------------------------------------------------------------
# Run helpers (not used by Spack)
# ------------------------------------------------------------------------------
run:
	mpirun -n $(NPROCS) ./$(TARGET) $(INPUT) $(FILLING) $(OUTPUT)

# ------------------------------------------------------------------------------
# Cleanup
# ------------------------------------------------------------------------------
clean:
	rm -rf $(SOLVERLOOP)/obj/*.o $(FUNCTIONS)/obj/*.o \
	       $(TDBS)/obj/*.o *.o $(TARGET)

.PHONY: all install run clean
