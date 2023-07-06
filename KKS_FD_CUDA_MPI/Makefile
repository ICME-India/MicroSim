# Build KKS_FD_CUDA_MPI
SHELL = /bin/sh
NPROCS = 2

HPCSDK_VERSION ?= 23.1

NVHPC_ROOT       ?= /opt/nvidia
NVHPC_HOME       ?= ${NVHPC_ROOT}/hpc_sdk/Linux_x86_64/${HPCSDK_VERSION}
NVHPC_CUDA_HOME  ?= ${NVHPC_HOME}/cuda
MPI_HOME         ?= ${NVHPC_HOME}/comm_libs/openmpi4/openmpi-4.0.5
HDF5_HOME        ?= /usr/local
CUFFT_LIB        ?= ${NVHPC_HOME}/math_libs/lib64
CUFFT_INC        ?= ${NVHPC_HOME}/math_libs/include/cufftmp
NVSHMEM_LIB      ?= ${NVHPC_HOME}/math_libs/12.0/lib64/compat/nvshmem_2.6.0-1
GSL_HOME         ?= /usr/local

# Arguments
INPUT = Input.in
FILLING = Filling.in
OUTPUT = Output
ENABLE_HDF5 = 1
ENABLE_CUFFTMP = 1

# Compilers
MPICC = $(PREP) ${MPI_HOME}/bin/mpicxx
MPILD = $(PREP) ${MPI_HOME}/bin/mpicxx
ifeq ($(ENABLE_HDF5), 1)
H5CC  = $(PREP) ${HDF5_HOME}/bin/h5c++
else
H5CC  = $(PREP) ${MPI_HOME}/bin/mpicxx
endif
NVCC = $(PREP) ${NVHPC_CUDA_HOME}/bin/nvcc

# Flags
CFLAGS = -lm 
GSLFLAGS = -I${GSL_HOME}/include -L${GSL_HOME}/lib -lgsl -lgslcblas
MPICFLAGS = -I${MPI_HOME}/include --diag_suppress declared_but_not_referenced --diag_suppress set_but_not_used
CUDACFLAGS = -I${NVHPC_CUDA_HOME}/include

GENCODE_SM70    := -gencode arch=compute_70,code=\"sm_70,compute_70\"
GENCODE_FLAGS   := $(GENCODE_SM70)

NVCCFLAGS = -Xcompiler -march=native -x cu $(GENCODE_FLAGS)

ifeq ($(ENABLE_CUFFTMP), 1)
CUDALDFLAGS = -L${NVHPC_CUDA_HOME}/lib64 -L${CUFFT_LIB} -lcudart -lm -lcufftMp
else
CUDALDFLAGS = -L${NVHPC_CUDA_HOME}/lib64 -lcudart -lm
endif
MPILDFLAGS = -L${MPI_HOME}/lib

# Description of binaries
KKS_OUT = ./microsim_kks_fd_cuda_mpi
SOLVERLOOP = ./solverloop
FUNCTIONS = ./functions
TDBS =./tdbs

ifeq ($(ENABLE_CUFFTMP), 1)
INCFLAGS = -I$(SOLVERLOOP) -I$(FUNCTIONS) -I$(TDBS) -I${CUFFT_INC} -I${MPI_HOME}/include
else
INCFLAGS = -I$(SOLVERLOOP) -I$(FUNCTIONS) -I$(TDBS) -I${MPI_HOME}/include
endif

host_objects = $(SOLVERLOOP)/obj/fileWriter.o $(FUNCTIONS)/defines.h $(FUNCTIONS)/obj/utilityFunctions.o $(FUNCTIONS)/obj/initialize_variables.o $(FUNCTIONS)/obj/filling.o $(FUNCTIONS)/obj/reading_input_parameters.o $(FUNCTIONS)/obj/fill_domain.o $(FUNCTIONS)/obj/read_boundary_conditions.o $(FUNCTIONS)/obj/calc_bn.o MicroSim_KKS_FD_CUDA_MPI.o Makefile

device_objects = $(TDBS)/obj/Thermo.o $(FUNCTIONS)/obj/functionF.o $(FUNCTIONS)/obj/functionH.o $(FUNCTIONS)/obj/functionTau.o $(FUNCTIONS)/obj/functionW_02.o $(FUNCTIONS)/obj/anisotropy_01.o $(FUNCTIONS)/obj/functionA_01.o $(FUNCTIONS)/obj/functionA_02.o $(FUNCTIONS)/obj/matrix.o $(SOLVERLOOP)/obj/utilityKernels.o $(SOLVERLOOP)/obj/boundary.o $(SOLVERLOOP)/obj/calcPhaseComp.o $(SOLVERLOOP)/obj/smooth.o $(SOLVERLOOP)/obj/computeDrivingForce.o $(FUNCTIONS)/obj/error_checks.o $(SOLVERLOOP)/obj/computeElastic.o $(SOLVERLOOP)/obj/updateComposition.o $(SOLVERLOOP)/obj/updatePhi.o

# Commands
build: $(KKS_OUT)

$(TDBS)/obj/Thermo.o: $(TDBS)/Thermo.cu $(TDBS)/Thermo.cuh Makefile
	$(NVCC) -dc $(TDBS)/Thermo.cu -o $(TDBS)/obj/Thermo.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(CUDACFLAGS)

$(SOLVERLOOP)/obj/fileWriter.o: $(SOLVERLOOP)/fileWriter.c $(SOLVERLOOP)/fileWriter.h $(FUNCTIONS)/structures.h $(FUNCTIONS)/defines.h Makefile
	$(H5CC) $(MPICFLAGS) $(CFLAGS) $(CUDACFLAGS) $(INCFLAGS) -DENABLE_HDF5=$(ENABLE_HDF5) -c $(SOLVERLOOP)/fileWriter.c -o $(SOLVERLOOP)/obj/fileWriter.o

$(FUNCTIONS)/obj/utilityFunctions.o: $(FUNCTIONS)/utilityFunctions.c $(FUNCTIONS)/utilityFunctions.h $(TDBS)/Thermo.cu $(FUNCTIONS)/structures.h $(FUNCTIONS)/defines.h Makefile
	$(MPICC) $(MPICFLAGS) $(CFLAGS) $(CUDACFLAGS) $(INCFLAGS) -c $(FUNCTIONS)/utilityFunctions.c -o $(FUNCTIONS)/obj/utilityFunctions.o

$(FUNCTIONS)/obj/functionF.o: $(FUNCTIONS)/functionF.cu $(FUNCTIONS)/functionF.cuh $(FUNCTIONS)/utilityFunctions.c $(FUNCTIONS)/structures.h $(FUNCTIONS)/defines.h Makefile
	$(NVCC) -dc $(FUNCTIONS)/functionF.cu -o $(FUNCTIONS)/obj/functionF.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(CUDACFLAGS) $(CFLAGS) $(GSLFLAGS)

$(FUNCTIONS)/obj/initialize_variables.o: $(FUNCTIONS)/initialize_variables.c $(FUNCTIONS)/initialize_variables.h $(FUNCTIONS)/utilityFunctions.c $(FUNCTIONS)/structures.h $(FUNCTIONS)/defines.h Makefile
	$(MPICC) $(MPICFLAGS) $(CFLAGS) $(CUDACFLAGS) $(INCFLAGS) -c $(FUNCTIONS)/initialize_variables.c -o $(FUNCTIONS)/obj/initialize_variables.o

$(FUNCTIONS)/obj/box_iterator.o: $(FUNCTIONS)/box_iterator.cpp $(FUNCTIONS)/box_iterator.hpp Makefile
	$(NVCC) $(CFLAGS) $(CUDACFLAGS) $(INCFLAGS) -DENABLE_CUFFTMP=$(ENABLE_CUFFTMP) -dc $(FUNCTIONS)/box_iterator.cpp -o $(FUNCTIONS)/obj/box_iterator.o

$(FUNCTIONS)/obj/calc_bn.o: $(FUNCTIONS)/calc_bn.c $(FUNCTIONS)/calc_bn.h $(FUNCTIONS)/structures.h $(FUNCTIONS)/defines.h Makefile
	$(MPICC) $(MPICFLAGS) $(CFLAGS) $(CUDACFLAGS) $(INCFLAGS) -DENABLE_CUFFTMP=$(ENABLE_CUFFTMP) -c $(FUNCTIONS)/calc_bn.c -o $(FUNCTIONS)/obj/calc_bn.o

$(FUNCTIONS)/obj/error_checks.o: $(FUNCTIONS)/error_checks.cpp $(FUNCTIONS)/error_checks.hpp Makefile
	$(NVCC) $(CFLAGS) $(CUDACFLAGS) $(INCFLAGS) -DENABLE_CUFFTMP=$(ENABLE_CUFFTMP) -dc $(FUNCTIONS)/error_checks.cpp -o $(FUNCTIONS)/obj/error_checks.o

$(FUNCTIONS)/obj/filling.o: $(FUNCTIONS)/filling.c $(FUNCTIONS)/filling.h $(FUNCTIONS)/utilityFunctions.c $(FUNCTIONS)/structures.h $(FUNCTIONS)/defines.h Makefile
	$(MPICC) $(MPICFLAGS) $(CFLAGS) $(CUDACFLAGS) $(INCFLAGS) -c $(FUNCTIONS)/filling.c -o $(FUNCTIONS)/obj/filling.o

$(FUNCTIONS)/obj/reading_input_parameters.o: $(FUNCTIONS)/reading_input_parameters.c $(FUNCTIONS)/inputReader.h $(FUNCTIONS)/structures.h $(FUNCTIONS)/defines.h $(FUNCTIONS)/utilityFunctions.c Makefile
	$(MPICC) $(MPICFLAGS) $(MPILDFLAGS) $(CFLAGS) $(CUDACFLAGS) $(INCFLAGS) -c $(FUNCTIONS)/reading_input_parameters.c -o $(FUNCTIONS)/obj/reading_input_parameters.o

$(FUNCTIONS)/obj/fill_domain.o: $(FUNCTIONS)/fill_domain.c $(FUNCTIONS)/inputReader.h $(FUNCTIONS)/structures.h $(FUNCTIONS)/defines.h $(FUNCTIONS)/utilityFunctions.c Makefile
	$(MPICC) $(MPICFLAGS) $(CFLAGS) $(CUDACFLAGS) $(INCFLAGS) -c $(FUNCTIONS)/fill_domain.c -o $(FUNCTIONS)/obj/fill_domain.o

$(FUNCTIONS)/obj/read_boundary_conditions.o: $(FUNCTIONS)/read_boundary_conditions.c $(FUNCTIONS)/inputReader.h $(FUNCTIONS)/structures.h $(FUNCTIONS)/defines.h $(FUNCTIONS)/utilityFunctions.c Makefile
	$(MPICC) $(MPICFLAGS) $(CFLAGS) $(CUDACFLAGS) $(INCFLAGS) -c $(FUNCTIONS)/read_boundary_conditions.c -o $(FUNCTIONS)/obj/read_boundary_conditions.o

$(FUNCTIONS)/obj/matrix.o: $(FUNCTIONS)/matrix.cu $(FUNCTIONS)/matrix.cuh $(FUNCTIONS)/structures.h $(FUNCTIONS)/defines.h Makefile
	$(NVCC) -dc $(FUNCTIONS)/matrix.cu -o $(FUNCTIONS)/obj/matrix.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(CUDACFLAGS)

$(SOLVERLOOP)/obj/utilityKernels.o: $(SOLVERLOOP)/utilityKernels.cu $(SOLVERLOOP)/utilityKernels.cuh $(FUNCTIONS)/structures.h $(FUNCTIONS)/defines.h Makefile
	$(NVCC) -dc $(SOLVERLOOP)/utilityKernels.cu -o $(SOLVERLOOP)/obj/utilityKernels.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(CUDACFLAGS)

$(FUNCTIONS)/obj/functionTau.o: $(FUNCTIONS)/functionTau.cu $(FUNCTIONS)/functionTau.cuh Makefile
	$(NVCC) -dc $(FUNCTIONS)/functionTau.cu -o $(FUNCTIONS)/obj/functionTau.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(CUDACFLAGS)

$(FUNCTIONS)/obj/functionH.o: $(FUNCTIONS)/functionH.cu $(FUNCTIONS)/functionH.cuh Makefile
	$(NVCC) -dc $(FUNCTIONS)/functionH.cu -o $(FUNCTIONS)/obj/functionH.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(CUDACFLAGS)

$(FUNCTIONS)/obj/functionW_02.o: $(FUNCTIONS)/functionW_02.cu $(FUNCTIONS)/functionW_02.cuh Makefile
	$(NVCC) -dc $(FUNCTIONS)/functionW_02.cu -o $(FUNCTIONS)/obj/functionW_02.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(CUDACFLAGS)

$(FUNCTIONS)/obj/anisotropy_01.o: $(FUNCTIONS)/anisotropy_01.cu $(FUNCTIONS)/anisotropy_01.cuh Makefile
	$(NVCC) -dc $(FUNCTIONS)/anisotropy_01.cu -o $(FUNCTIONS)/obj/anisotropy_01.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(CUDACFLAGS)

$(FUNCTIONS)/obj/functionA_01.o: $(FUNCTIONS)/functionA_01.cu $(FUNCTIONS)/functionA_01.cuh Makefile
	$(NVCC) -dc $(FUNCTIONS)/functionA_01.cu -o $(FUNCTIONS)/obj/functionA_01.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(CUDACFLAGS)

$(FUNCTIONS)/obj/functionA_02.o: $(FUNCTIONS)/functionA_02.cu $(FUNCTIONS)/functionA_02.cuh Makefile
	$(NVCC) -dc $(FUNCTIONS)/functionA_02.cu -o $(FUNCTIONS)/obj/functionA_02.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(CUDACFLAGS)

MicroSim_KKS_FD_CUDA_MPI.o: ./microsim_kks_fd_cuda_mpi.c Makefile
	$(MPICC) $(MPICFLAGS) $(CFLAGS) $(GSLFLAGS) $(CUDACFLAGS) $(INCFLAGS) -DENABLE_HDF5=$(ENABLE_HDF5) -DENABLE_CUFFTMP=$(ENABLE_CUFFTMP) -c ./microsim_kks_fd_cuda_mpi.c -o ./MicroSim_KKS_FD_CUDA_MPI.o $(CUDALDFLAGS)

$(SOLVERLOOP)/obj/boundary.o: $(SOLVERLOOP)/boundary.cu $(SOLVERLOOP)/boundary.cuh $(FUNCTIONS)/structures.h $(FUNCTIONS)/defines.h Makefile
	$(NVCC) -dc $(SOLVERLOOP)/boundary.cu -o $(SOLVERLOOP)/obj/boundary.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS)

$(SOLVERLOOP)/obj/calcPhaseComp.o: $(SOLVERLOOP)/calcPhaseComp.cu $(SOLVERLOOP)/calcPhaseComp.cuh $(TDBS)/Thermo.cu $(FUNCTIONS)/structures.h $(FUNCTIONS)/defines.h Makefile
	$(NVCC) -dc $(SOLVERLOOP)/calcPhaseComp.cu -o $(SOLVERLOOP)/obj/calcPhaseComp.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(CUDACFLAGS)

$(SOLVERLOOP)/obj/smooth.o: $(SOLVERLOOP)/smooth.cu $(SOLVERLOOP)/smooth.cuh $(TDBS)/Thermo.cu $(FUNCTIONS)/structures.h $(FUNCTIONS)/defines.h Makefile
	$(NVCC) -dc $(SOLVERLOOP)/smooth.cu -o $(SOLVERLOOP)/obj/smooth.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(CUDACFLAGS)

$(SOLVERLOOP)/obj/computeDrivingForce.o: $(SOLVERLOOP)/computeDrivingForce.cu $(SOLVERLOOP)/computeDrivingForce.cuh $(TDBS)/Thermo.cu $(FUNCTIONS)/structures.h $(FUNCTIONS)/defines.h Makefile
	$(NVCC) -dc $(SOLVERLOOP)/computeDrivingForce.cu -o $(SOLVERLOOP)/obj/computeDrivingForce.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(CUDACFLAGS)

$(SOLVERLOOP)/obj/computeElastic.o: $(SOLVERLOOP)/computeElastic.cu $(SOLVERLOOP)/computeElastic.cuh $(TDBS)/Thermo.cu $(FUNCTIONS)/structures.h $(FUNCTIONS)/defines.h Makefile
	$(NVCC) -DENABLE_CUFFTMP=$(ENABLE_CUFFTMP) -dc $(SOLVERLOOP)/computeElastic.cu -o $(SOLVERLOOP)/obj/computeElastic.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(CUDACFLAGS) -std=c++17

$(SOLVERLOOP)/obj/updateComposition.o: $(SOLVERLOOP)/updateComposition.cu $(SOLVERLOOP)/updateComposition.cuh $(TDBS)/Thermo.cu $(FUNCTIONS)/structures.h $(FUNCTIONS)/defines.h Makefile
	$(NVCC) -dc $(SOLVERLOOP)/updateComposition.cu -o $(SOLVERLOOP)/obj/updateComposition.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(CUDACFLAGS)

$(SOLVERLOOP)/obj/updatePhi.o: $(SOLVERLOOP)/updatePhi.cu $(SOLVERLOOP)/updatePhi.cuh $(TDBS)/Thermo.cu $(FUNCTIONS)/structures.h $(FUNCTIONS)/defines.h Makefile
	$(NVCC) -dc $(SOLVERLOOP)/updatePhi.cu -o $(SOLVERLOOP)/obj/updatePhi.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(CUDACFLAGS)

link: $(device_objects) $(FUNCTIONS)/defines.h
	$(NVCC) $(GENCODE_FLAGS) -dlink ./MicroSim_KKS_FD_CUDA_MPI.o $(device_objects) -o $(SOLVERLOOP)/obj/deviceCode.o

$(KKS_OUT): $(host_objects) link
	$(H5CC) ./MicroSim_KKS_FD_CUDA_MPI.o $(SOLVERLOOP)/obj/*.o $(FUNCTIONS)/obj/*.o $(TDBS)/obj/Thermo.o $(CUDALDFLAGS) $(CUDACFLAGS) $(MPILDFLAGS) $(CFLAGS) $(GSLFLAGS)  -o $(KKS_OUT)
	rm -rf MicroSim_KKS_FD_CUDA_MPI.o $(SOLVERLOOP)/obj/fileWriter.o $(SOLVERLOOP)/obj/computeElastic.o $(FUNCTIONS)/obj/error_checks.o
thermo: $(INPUT)
	python GEdata_writer.py $(INPUT)

write_xdmf: write_xdmf.c $(FUNCTIONS)/structures.h $(FUNCTIONS)/reading_input_parameters.c $(SOLVERLOOP)/fileWriter.c $(INPUT) Makefile
	$(H5CC) $(INCFLAGS) $(CFLAGS) -I${NVHPC_CUDA_HOME}/include -L${NVHPC_CUDA_HOME}/lib64 -lcudart -DTHERMO=0 write_xdmf.c -o write_xdmf --diag_suppress declared_but_not_referenced

reconstruct: reconstruct.c $(FUNCTIONS)/structures.h $(FUNCTIONS)/reading_input_parameters.c $(SOLVERLOOP)/fileWriter.c $(INPUT) Makefile
	$(MPICC) $(INCFLAGS) $(CFLAGS)  --diag_suppress declared_but_not_referenced -I${NVHPC_CUDA_HOME}/include -L${NVHPC_CUDA_HOME}/lib64 -lcudart -DTHERMO=0 reconstruct.c -o reconstruct

run:
	LD_LIBRARY_PATH="${NVSHMEM_LIB}:${CUFFT_LIB}:${NVHPC_CUDA_HOME}/lib64:${LD_LIBRARY_PATH}" ${MPI_HOME}/bin/mpiexec -n $(NPROCS) $(KKS_OUT) $(INPUT) $(FILLING) $(OUTPUT)

profile:
	LD_LIBRARY_PATH="${NVSHMEM_LIB}:${CUFFT_LIB}:${NVHPC_CUDA_HOME}/lib64:${LD_LIBRARY_PATH}" nvprof  --profile-child-processes mpiexec -n $(NPROCS) $(KKS_OUT) $(INPUT) $(FILLING) $(OUTPUT)

memcheck:
	LD_LIBRARY_PATH="${NVSHMEM_LIB}:${CUFFT_LIB}:${NVHPC_CUDA_HOME}/lib64:${LD_LIBRARY_PATH}" compute-sanitizer --target-processes all mpiexec -n $(NPROCS) $(KKS_OUT) $(INPUT) $(FILLING) $(OUTPUT)

clear:
	rm -rf DATA/ *.bd */*.bd

.PHONY: clean
clean:
	rm -rf $(SOLVERLOOP)/obj/*.o $(FUNCTIONS)/obj/*.o $(KKS_OUT) *.o $(TDBS)/obj/*.o *.out
