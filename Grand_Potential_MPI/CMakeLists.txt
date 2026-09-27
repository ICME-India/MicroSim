cmake_minimum_required(VERSION 3.18)
project(microsim_gp_mpi LANGUAGES C)

# ── Dependencies ────────────────────────────────────────────────────────────
find_package(MPI  REQUIRED COMPONENTS C)
find_package(HDF5 REQUIRED COMPONENTS C HL)
find_package(GSL  REQUIRED)

# ── Main solver target ───────────────────────────────────────────────────────
# The build consists of a single .c translation unit; all other files
# in functions/, solverloop/ and tdbs/ are header-only includes.
add_executable(microsim_gp microsim_gp.c)

target_include_directories(microsim_gp PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/functions
    ${CMAKE_CURRENT_SOURCE_DIR}/solverloop
    ${CMAKE_CURRENT_SOURCE_DIR}/tdbs
    ${HDF5_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIRS}
)

target_compile_definitions(microsim_gp PRIVATE ${HDF5_DEFINITIONS})

target_link_libraries(microsim_gp PRIVATE
    MPI::MPI_C
    ${HDF5_C_LIBRARIES}
    GSL::gsl
    GSL::gslcblas
    m
)

# ── write_xdmf utility ───────────────────────────────────────────────────────
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/write_xdmf.c")
    add_executable(write_xdmf write_xdmf.c)

    target_include_directories(write_xdmf PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/functions
        ${CMAKE_CURRENT_SOURCE_DIR}/solverloop
        ${HDF5_INCLUDE_DIRS}
        ${GSL_INCLUDE_DIRS}
    )

    target_compile_definitions(write_xdmf PRIVATE ${HDF5_DEFINITIONS})

    target_link_libraries(write_xdmf PRIVATE
        MPI::MPI_C
        ${HDF5_C_LIBRARIES}
        GSL::gsl
        GSL::gslcblas
        m
    )
endif()

install(TARGETS microsim_gp RUNTIME DESTINATION bin)
