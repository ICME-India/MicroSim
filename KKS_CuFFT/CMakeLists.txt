cmake_minimum_required(VERSION 3.18)
project(microsim_kks_cufft LANGUAGES CUDA)

# ── Dependencies ────────────────────────────────────────────────────────────
find_package(CUDAToolkit REQUIRED)

# ── Target ───────────────────────────────────────────────────────────────────
add_executable(microsim_kks_cufft KKS_cuFFT.cu)

# Optimisation flags passed through to ptxas (matches -Xptxas -O3 in Makefile)
target_compile_options(microsim_kks_cufft PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:-Xptxas -O3>
)

# GPU architecture – default sm_70; override via -DCUDA_ARCH=<arch>
if(NOT DEFINED CUDA_ARCH)
    set(CUDA_ARCH "70")
endif()
set_target_properties(microsim_kks_cufft PROPERTIES
    CUDA_ARCHITECTURES ${CUDA_ARCH}
)

target_include_directories(microsim_kks_cufft PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/functions
)

target_link_libraries(microsim_kks_cufft PRIVATE
    CUDA::cufft
    CUDA::cuda_driver
    CUDA::cudart
)

# Create DATA output directory expected by the solver at run-time
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/DATA")

install(TARGETS microsim_kks_cufft RUNTIME DESTINATION bin)
