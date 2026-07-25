# Compiler and options
CXX := mpiicpx
CXXFLAGS := -fsycl -std=c++20 -g -Wall -O1
INCLUDES := -Iinclude -I/home/arunc/hdf5-parallel/include
OBJDIR := obj
TARGET := bin/PF_SYCL

LIBS := -L/home/arunc/hdf5-parallel/lib -lhdf5 -lhdf5_cpp
# Add CUDA backend flag if GPU_CUDA=1 is passed to make

# By Default it is divice aware mpi in CPU
# use only if GPU_AWARE_MPI is not supported by your MPI_BUILD

ifeq ($(CPU),1)
	CXXFLAGS += -DCPU_H
#	-DENABLE_CPU_SHARED_MPI
endif

# CUDA GPU + disable device-aware MPI
ifeq ($(GPU),1)
	CXXFLAGS += -DGPU 
	CXXFLAGS += -fsycl-targets=nvptx64-nvidia-cuda

# for device-aware MPI (if supported by MPI build)
#	CXXFLAGS += -DENABLE_GPU_AWARE_MPI
#	INCLUDES += -I$(HOME)/openmpi-cuda/include
#	LIBS += -L$(HOME)/openmpi-cuda/lib -lmpi -lopen-pal -lpmix -lprrte
endif


# Common Source Files
SRC := \
    src/main.cpp \
    src/comm/device_info.cpp \
    src/utilities/Logger.cpp \
    src/utilities/IO.cpp \
    src/utilities/microsim.cpp \
    src/comm/MPI_Info.cpp \
    src/physics/filling.cpp \
    src/thermo/Function_F4.cpp

# Common Object Files
OBJ := $(OBJDIR)/main.o \
       $(OBJDIR)/comm/device_info.o \
       $(OBJDIR)/utilities/Logger.o \
       $(OBJDIR)/utilities/IO.o \
       $(OBJDIR)/utilities/microsim.o \
       $(OBJDIR)/comm/MPI_Info.o \
       $(OBJDIR)/physics/filling.o \
       $(OBJDIR)/thermo/Function_F4.o

# Model selection
ifeq ($(KKS),1)
    # --- KKS Selected ---
    TARGET := bin/KKS_SYCL
    CXXFLAGS += -D__KKS

    SRC += src/physics/Kim_Kim_Suzuky.cpp
    OBJ += $(OBJDIR)/physics/Kim_Kim_Suzuky.o
    
else ifeq ($(FID),1)
    # --- FID Selected ---
    TARGET := bin/FID_SYCL
    CXXFLAGS += -D__FID
    
    SRC += src/physics/Finite_Interface_Dessipation.cpp
    OBJ += $(OBJDIR)/physics/Finite_Interface_Dessipation.o
    
else ifeq ($(RKB),1)
    # --- GP Selected Default ---
    #TARGET := bin/RKB_SYCL
    CXXFLAGS += -D__GP

    SRC += src/physics/kobayshi.cpp
    OBJ += $(OBJDIR)/physics/kobayshi.o

else
    # --- GP Selected Default ---
    CXXFLAGS += -D__GP
    #TARGET := bin/GP_SYCL 

    SRC += src/physics/Grand_Potential.cpp
    OBJ += $(OBJDIR)/physics/Grand_Potential.o
endif


ifdef DEBUG
    CXXFLAGS += -g -O0 -Wextra
endif


# Libraries
LDFLAGS := $(LIBS) -lgsl -lgslcblas -lm 

# Default rule
all: $(TARGET)

# Linking final executable
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)


# Compile
$(OBJDIR)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean_obj:
	rm -rf $(OBJDIR)

run: $(TARGET)
	mpirun -np 1 ./$(TARGET) ./sys/Input_aniso.in ./sys/Filling_aniso.in Dec 1 1

# Clean rule
clean:
	rm -rf $(OBJDIR) $(TARGET) rank*

.PHONY: all clean