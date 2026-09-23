cmake_minimum_required(VERSION 3.18)
project(microsim_kks_fd_cuda_mpi LANGUAGES C CXX CUDA)

# ── Dependencies ────────────────────────────────────────────────────────────
find_package(MPI         REQUIRED COMPONENTS C CXX)
find_package(HDF5        REQUIRED COMPONENTS C HL)
find_package(GSL         REQUIRED)
find_package(CUDAToolkit REQUIRED)

# ── Build options ─────────────────────────────────────────────────────────────
option(ENABLE_HDF5    "Enable HDF5 output"        ON)
option(ENABLE_CUFFTMP "Enable multi-GPU cuFFTMp"  OFF)

# GPU architecture – default sm_70; override via -DCUDA_ARCH=<arch>
if(NOT DEFINED CUDA_ARCH)
    set(CUDA_ARCH "70")
endif()

# ── Helper: convenience path variables ───────────────────────────────────────
set(FUNCTIONS  ${CMAKE_CURRENT_SOURCE_DIR}/functions)
set(SOLVERLOOP ${CMAKE_CURRENT_SOURCE_DIR}/solverloop)
set(TDBS       ${CMAKE_CURRENT_SOURCE_DIR}/tdbs)

# ── Sources ──────────────────────────────────────────────────────────────────
set(HOST_SOURCES
    microsim_kks_fd_cuda_mpi.c
    ${SOLVERLOOP}/fileWriter.c
    ${FUNCTIONS}/utilityFunctions.c
    ${FUNCTIONS}/initialize_variables.c
    ${FUNCTIONS}/filling.c
    ${FUNCTIONS}/reading_input_parameters.c
    ${FUNCTIONS}/fill_domain.c
    ${FUNCTIONS}/read_boundary_conditions.c
    ${FUNCTIONS}/calc_bn.c
    ${FUNCTIONS}/box_iterator.cpp
    ${FUNCTIONS}/error_checks.cpp
)

set(DEVICE_SOURCES
    ${TDBS}/Thermo.cu
    ${FUNCTIONS}/functionF.cu
    ${FUNCTIONS}/functionH.cu
    ${FUNCTIONS}/functionTau.cu
    ${FUNCTIONS}/functionW_02.cu
    ${FUNCTIONS}/anisotropy_01.cu
    ${FUNCTIONS}/functionA_01.cu
    ${FUNCTIONS}/functionA_02.cu
    ${FUNCTIONS}/matrix.cu
    ${SOLVERLOOP}/utilityKernels.cu
    ${SOLVERLOOP}/boundary.cu
    ${SOLVERLOOP}/calcPhaseComp.cu
    ${SOLVERLOOP}/smooth.cu
    ${SOLVERLOOP}/computeDrivingForce.cu
    ${SOLVERLOOP}/computeElastic.cu
    ${SOLVERLOOP}/updateComposition.cu
    ${SOLVERLOOP}/updatePhi.cu
)

# ── Target ───────────────────────────────────────────────────────────────────
add_executable(microsim_kks_fd_cuda_mpi
    ${HOST_SOURCES}
    ${DEVICE_SOURCES}
)

# CUDA separable compilation handles the -dc / device-link steps automatically.
set_target_properties(microsim_kks_fd_cuda_mpi PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES         ${CUDA_ARCH}
)

# computeElastic requires C++17
set_source_files_properties(${SOLVERLOOP}/computeElastic.cu
    PROPERTIES COMPILE_OPTIONS "$<$<COMPILE_LANGUAGE:CUDA>:-std=c++17>"
)

# ── Include directories ──────────────────────────────────────────────────────
target_include_directories(microsim_kks_fd_cuda_mpi PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${FUNCTIONS}
    ${SOLVERLOOP}
    ${TDBS}
    ${HDF5_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIRS}
    ${CUDAToolkit_INCLUDE_DIRS}
)

if(ENABLE_CUFFTMP)
    # cuFFTMp headers are shipped with the HPC SDK; users can override the path.
    if(NOT DEFINED CUFFTMP_INCLUDE_DIR)
        message(WARNING "ENABLE_CUFFTMP=ON but CUFFTMP_INCLUDE_DIR is not set. "
                        "Provide -DCUFFTMP_INCLUDE_DIR=<path> to the cuFFTMp headers.")
    else()
        target_include_directories(microsim_kks_fd_cuda_mpi PRIVATE ${CUFFTMP_INCLUDE_DIR})
    endif()
endif()

# ── Compile definitions ───────────────────────────────────────────────────────
target_compile_definitions(microsim_kks_fd_cuda_mpi PRIVATE
    ENABLE_HDF5=$<BOOL:${ENABLE_HDF5}>
    ENABLE_CUFFTMP=$<BOOL:${ENABLE_CUFFTMP}>
    ${HDF5_DEFINITIONS}
)

# Suppress common NVHPC/nvcc diagnostic noise
target_compile_options(microsim_kks_fd_cuda_mpi PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:
        --diag_suppress=declared_but_not_referenced
        --diag_suppress=set_but_not_used>
    $<$<COMPILE_LANGUAGE:C>:
        --diag_suppress=declared_but_not_referenced
        --diag_suppress=set_but_not_used>
)

# ── Link libraries ────────────────────────────────────────────────────────────
target_link_libraries(microsim_kks_fd_cuda_mpi PRIVATE
    MPI::MPI_CXX
    GSL::gsl
    GSL::gslcblas
    CUDA::cudart
    m
)

if(ENABLE_HDF5)
    target_link_libraries(microsim_kks_fd_cuda_mpi PRIVATE ${HDF5_C_LIBRARIES})
endif()

if(ENABLE_CUFFTMP)
    # cuFFTMp – library location can be overridden with -DCUFFTMP_LIB_DIR=<path>
    if(DEFINED CUFFTMP_LIB_DIR)
        target_link_directories(microsim_kks_fd_cuda_mpi PRIVATE ${CUFFTMP_LIB_DIR})
    endif()
    target_link_libraries(microsim_kks_fd_cuda_mpi PRIVATE cufftMp)
else()
    target_link_libraries(microsim_kks_fd_cuda_mpi PRIVATE CUDA::cufft)
endif()

# ── Utility: write_xdmf ──────────────────────────────────────────────────────
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/write_xdmf.c")
    add_executable(write_xdmf_kks write_xdmf.c)

    target_include_directories(write_xdmf_kks PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${FUNCTIONS}
        ${SOLVERLOOP}
        ${HDF5_INCLUDE_DIRS}
        ${CUDAToolkit_INCLUDE_DIRS}
    )

    target_compile_definitions(write_xdmf_kks PRIVATE
        THERMO=0
        ENABLE_HDF5=$<BOOL:${ENABLE_HDF5}>
        ${HDF5_DEFINITIONS}
    )

    target_link_libraries(write_xdmf_kks PRIVATE
        MPI::MPI_CXX
        ${HDF5_C_LIBRARIES}
        CUDA::cudart
        m
    )
endif()

# ── Utility: reconstruct ──────────────────────────────────────────────────────
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/reconstruct.c")
    add_executable(reconstruct_kks reconstruct.c)

    target_include_directories(reconstruct_kks PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${FUNCTIONS}
        ${SOLVERLOOP}
        ${CUDAToolkit_INCLUDE_DIRS}
    )

    target_compile_definitions(reconstruct_kks PRIVATE THERMO=0)

    target_link_libraries(reconstruct_kks PRIVATE
        MPI::MPI_CXX
        CUDA::cudart
        m
    )
endif()

install(TARGETS microsim_kks_fd_cuda_mpi RUNTIME DESTINATION bin)
