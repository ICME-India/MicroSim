cmake_minimum_required(VERSION 3.18)
project(microsim_kks_opencl LANGUAGES C)

# ── Dependencies ────────────────────────────────────────────────────────────
find_package(MPI    REQUIRED COMPONENTS C)
find_package(OpenCL REQUIRED)
find_package(HDF5   REQUIRED COMPONENTS C HL)
find_package(GSL    REQUIRED)

# ── Main solver target ───────────────────────────────────────────────────────
# Single-translation-unit C solver; functions/, solverloop/ and tdbs/ are
# header-only includes.
add_executable(microsim_kks_opencl microsim_kks_opencl.c)

target_include_directories(microsim_kks_opencl PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/functions
    ${CMAKE_CURRENT_SOURCE_DIR}/solverloop
    ${CMAKE_CURRENT_SOURCE_DIR}/tdbs
    ${HDF5_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIRS}
)

target_compile_definitions(microsim_kks_opencl PRIVATE ${HDF5_DEFINITIONS})

target_link_libraries(microsim_kks_opencl PRIVATE
    MPI::MPI_C
    OpenCL::OpenCL
    GSL::gsl
    GSL::gslcblas
    ${HDF5_C_LIBRARIES}
    m
)

# ── reconstruct utility ───────────────────────────────────────────────────────
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/reconstruct.c")
    add_executable(reconstruct_opencl reconstruct.c)

    target_include_directories(reconstruct_opencl PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/functions
        ${HDF5_INCLUDE_DIRS}
        ${GSL_INCLUDE_DIRS}
    )

    target_compile_definitions(reconstruct_opencl PRIVATE ${HDF5_DEFINITIONS})

    target_link_libraries(reconstruct_opencl PRIVATE
        MPI::MPI_C
        OpenCL::OpenCL
        GSL::gsl
        GSL::gslcblas
        ${HDF5_C_LIBRARIES}
        m
    )
endif()

install(TARGETS microsim_kks_opencl RUNTIME DESTINATION bin)
