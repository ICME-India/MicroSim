# Build KKS_FD_CUDA_MPI
SHELL = /bin/sh
NPROCS = 2
CUDA_HOME ?= $(shell dirname $$(command -v nvcc))/..
CUB_HOME = ${CUDA_HOME}/targets/x86_64-linux/include/cub

# Change these if required
# Decrease to reduce memory usage
MAX_NUMPHASES = 5
MAX_NUMCOMPONENTS = 5
NUM = $$(($$((($(MAX_NUMCOMPONENTS) - 1))) * $(MAX_NUMPHASES)))
MAXFLAGS = -DMAX_NUM_PHASES=$(MAX_NUMPHASES) -DMAX_NUM_COMP=$(MAX_NUMCOMPONENTS) -DMAX_NUM_PHASE_COMP=$(NUM)

# Arguments
INPUT = ./Input.in
FILLING = ./Filling.in
OUTPUT = Output

# Compilers
MPICC = $(PREP) mpicxx
MPILD = $(PREP) mpicxx
NVCC = $(PREP) nvcc

# Flags
CFLAGS = -lm
GSLFLAGS = -lgsl -lcblas
MPICFLAGS = -I${MPI_HOME}/include
CUDACFLAGS = -I${CUDA_HOME}/include -I${CUB_HOME}

GENCODE_SM70    := -gencode arch=compute_70,code=\"sm_70,compute_70\"
GENCODE_FLAGS   := $(GENCODE_SM70)

NVCCFLAGS = -Xcompiler -march=native -x cu

CUDALDFLAGS = -L${CUDA_HOME}/lib64 -lcudart -lm
MPILDFLAGS = -L/usr/lib/openmpi

# Description of binaries
KKS_OUT = ./microsim_kks_fd_cuda_mpi
SOLVERLOOP = ./solverloop
FUNCTIONS = ./functions
TDBS =./tdbs
INCFLAGS = -I$(SOLVERLOOP) -I$(FUNCTIONS) -I$(TDBS)

objects = $(TDBS)/obj/Thermo.o $(SOLVERLOOP)/obj/fileWriter.o $(FUNCTIONS)/obj/utilityFunctions.o $(FUNCTIONS)/obj/functionF.o $(FUNCTIONS)/obj/initialize_variables.o $(FUNCTIONS)/obj/filling.o $(FUNCTIONS)/obj/inputReader.o MicroSim_KKS_FD_CUDA_MPI.o $(SOLVERLOOP)/obj/utilityKernels.o $(SOLVERLOOP)/obj/calcPhaseComp.o  $(SOLVERLOOP)/obj/computeDrivingForce.o $(SOLVERLOOP)/obj/updateComposition.o $(SOLVERLOOP)/obj/updatePhi.o Makefile

# Commands
build: $(KKS_OUT)

$(TDBS)/obj/Thermo.o: $(TDBS)/Thermo.cu $(TDBS)/Thermo.cuh Input.in Makefile
	$(NVCC) -dc $(TDBS)/Thermo.cu -o $(TDBS)/obj/Thermo.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(CUDACFLAGS) $(MAXFLAGS)

$(SOLVERLOOP)/obj/fileWriter.o: $(SOLVERLOOP)/fileWriter.c $(SOLVERLOOP)/fileWriter.h $(FUNCTIONS)/structures.h Makefile
	$(MPILD) $(MPICFLAGS) $(CFLAGS) $(CUDACFLAGS) $(INCFLAGS) -c $(SOLVERLOOP)/fileWriter.c -o $(SOLVERLOOP)/obj/fileWriter.o

$(FUNCTIONS)/obj/utilityFunctions.o: $(FUNCTIONS)/utilityFunctions.c $(FUNCTIONS)/utilityFunctions.h $(FUNCTIONS)/structures.h Makefile
	$(MPICC) $(MPICFLAGS) $(CFLAGS) $(CUDACFLAGS) $(INCFLAGS) -c $(FUNCTIONS)/utilityFunctions.c -o $(FUNCTIONS)/obj/utilityFunctions.o

$(FUNCTIONS)/obj/functionF.o: $(FUNCTIONS)/functionF.cu $(FUNCTIONS)/functionF.cuh $(FUNCTIONS)/utilityFunctions.h $(FUNCTIONS)/structures.h Makefile
	$(NVCC) -dc $(FUNCTIONS)/functionF.cu -o $(FUNCTIONS)/obj/functionF.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(CUDACFLAGS) $(CFLAGS) $(GSLFLAGS) $(MAXFLAGS)

$(FUNCTIONS)/obj/initialize_variables.o: $(FUNCTIONS)/initialize_variables.c $(FUNCTIONS)/initialize_variables.h $(FUNCTIONS)/utilityFunctions.h $(FUNCTIONS)/structures.h Makefile
	$(MPICC) $(MPICFLAGS) $(CFLAGS) $(CUDACFLAGS) $(INCFLAGS) -c $(FUNCTIONS)/initialize_variables.c -o $(FUNCTIONS)/obj/initialize_variables.o

$(FUNCTIONS)/obj/filling.o: $(FUNCTIONS)/filling.c $(FUNCTIONS)/filling.h $(FUNCTIONS)/utilityFunctions.h $(FUNCTIONS)/structures.h Makefile
	$(MPICC) $(MPICFLAGS) $(CFLAGS) $(CUDACFLAGS) $(INCFLAGS) -c $(FUNCTIONS)/filling.c -o $(FUNCTIONS)/obj/filling.o

$(FUNCTIONS)/obj/inputReader.o: $(FUNCTIONS)/inputReader.c $(FUNCTIONS)/inputReader.h $(FUNCTIONS)/structures.h $(FUNCTIONS)/utilityFunctions.h Makefile
	$(MPICC) $(MPICFLAGS) $(CFLAGS) $(CUDACFLAGS) $(INCFLAGS) -c $(FUNCTIONS)/inputReader.c -o $(FUNCTIONS)/obj/inputReader.o

$(SOLVERLOOP)/obj/utilityKernels.o: $(SOLVERLOOP)/utilityKernels.cu $(SOLVERLOOP)/utilityKernels.cuh $(FUNCTIONS)/structures.h Makefile
	$(NVCC) -dc $(SOLVERLOOP)/utilityKernels.cu -o $(SOLVERLOOP)/obj/utilityKernels.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(CUDACFLAGS)

MicroSim_KKS_FD_CUDA_MPI.o: ./MicroSim_KKS_FD_CUDA_MPI.c Makefile
	$(MPICC) $(MPICFLAGS) $(CFLAGS) $(GSLFLAGS) $(CUDACFLAGS) $(INCFLAGS) $(MAXFLAGS) -c ./MicroSim_KKS_FD_CUDA_MPI.c -o ./MicroSim_KKS_FD_CUDA_MPI.o $(CUDALDFLAGS)

$(SOLVERLOOP)/obj/calcPhaseComp.o: $(SOLVERLOOP)/calcPhaseComp.cu $(SOLVERLOOP)/calcPhaseComp.cuh $(FUNCTIONS)/structures.h Makefile
	$(NVCC) -dc $(SOLVERLOOP)/calcPhaseComp.cu -o $(SOLVERLOOP)/obj/calcPhaseComp.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(MAXFLAGS)

$(SOLVERLOOP)/obj/computeDrivingForce.o: $(SOLVERLOOP)/computeDrivingForce.cu $(SOLVERLOOP)/computeDrivingForce.cuh $(FUNCTIONS)/structures.h Makefile
	$(NVCC) -dc $(SOLVERLOOP)/computeDrivingForce.cu -o $(SOLVERLOOP)/obj/computeDrivingForce.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(MAXFLAGS)

$(SOLVERLOOP)/obj/updateComposition.o: $(SOLVERLOOP)/updateComposition.cu $(SOLVERLOOP)/updateComposition.cuh $(FUNCTIONS)/structures.h Makefile
	$(NVCC) -dc $(SOLVERLOOP)/updateComposition.cu -o $(SOLVERLOOP)/obj/updateComposition.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(MAXFLAGS)

$(SOLVERLOOP)/obj/updatePhi.o: $(SOLVERLOOP)/updatePhi.cu $(SOLVERLOOP)/updatePhi.cuh $(FUNCTIONS)/structures.h Makefile
	$(NVCC) -dc $(SOLVERLOOP)/updatePhi.cu -o $(SOLVERLOOP)/obj/updatePhi.o $(INCFLAGS) $(NVCCFLAGS) $(CUDALDFLAGS) $(MAXFLAGS)

link: $(objects)
	$(NVCC) -dlink ./MicroSim_KKS_FD_CUDA_MPI.o $(TDBS)/obj/Thermo.o $(FUNCTIONS)/obj/functionF.o $(SOLVERLOOP)/obj/utilityKernels.o $(SOLVERLOOP)/obj/calcPhaseComp.o $(SOLVERLOOP)/obj/computeDrivingForce.o $(SOLVERLOOP)/obj/updateComposition.o $(SOLVERLOOP)/obj/updatePhi.o -o $(SOLVERLOOP)/obj/deviceCode.o

$(KKS_OUT): link
	$(MPILD) ./MicroSim_KKS_FD_CUDA_MPI.o $(SOLVERLOOP)/obj/*.o $(FUNCTIONS)/obj/*.o $(TDBS)/obj/Thermo.o $(CUDALDFLAGS) $(CUDACFLAGS) $(MPILDFLAGS) $(CFLAGS) $(GSLFLAGS)  -o $(KKS_OUT)
	rm -rf MicroSim_KKS_FD_CUDA_MPI.o

reconstruct: reconstruct.c $(FUNCTIONS)/structures.h $(FUNCTIONS)/inputReader.c $(SOLVERLOOP)/fileWriter.c $(INPUT) Makefile
	$(MPILD) $(INCFLAGS) $(CFLAGS) $(CUDACFLAGS) $(CUDALDFLAGS) reconstruct.c -o reconstruct.out

run:
	mpirun -n $(NPROCS) $(KKS_OUT) $(INPUT) $(FILLING) $(OUTPUT)

profile:
	nvprof --profile-child-processes mpirun -n $(NPROCS) $(KKS_OUT) $(INPUT) $(FILLING) $(OUTPUT)

memcheck:
	cuda-memcheck mpirun -n $(NPROCS) $(KKS_OUT) $(INPUT) $(FILLING) $(OUTPUT)

clear:
	rm -rf DATA/

clean:
	rm -rf DATA/Processor* $(SOLVERLOOP)/obj/*.o $(FUNCTIONS)/obj/*.o $(KKS_OUT) *.o $(TDBS)/obj/*.o
