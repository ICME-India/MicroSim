cmake_minimum_required(VERSION 3.16)
project(MCCH_Solver VERSION 1.0.0 LANGUAGES CXX C)

# Prevent in-source builds to avoid polluting the source tree with non-portable CMakeCache.txt
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/CMakeCache.txt" "${CMAKE_BINARY_DIR}/CMakeFiles")
    message(FATAL_ERROR "In-source builds are disallowed. Please create a separate build directory (e.g. 'cmake -B build' or './install.sh') to preserve source tree portability.")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# -----------------------------------------------------------------------------
# Compiler & Optimization Flags
# -----------------------------------------------------------------------------
include(CheckCXXCompilerFlag)

option(ENABLE_NATIVE_TUNING "Optimize for the host CPU microarchitecture (-march=native / -mcpu=native)" OFF)

if(MSVC)
    add_compile_options(/O2 /W4 /EHsc /permissive-)
else()
    add_compile_options(-O3 -Wall -Wextra -Wno-unused-parameter)
    if(ENABLE_NATIVE_TUNING)
        check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
        if(COMPILER_SUPPORTS_MARCH_NATIVE)
            add_compile_options(-march=native)
            message(STATUS "Native CPU tuning ENABLED: -march=native")
        else()
            check_cxx_compiler_flag("-mcpu=native" COMPILER_SUPPORTS_MCPU_NATIVE)
            if(COMPILER_SUPPORTS_MCPU_NATIVE)
                add_compile_options(-mcpu=native)
                message(STATUS "Native CPU tuning ENABLED: -mcpu=native")
            else()
                message(WARNING "Native CPU tuning requested, but host flag not supported by compiler.")
            endif()
        endif()
    else()
        message(STATUS "Native CPU tuning: DISABLED (building portable generic binary)")
    endif()
endif()

# Math library (needed on Unix/Linux, built-in on Windows/MSVC)
if(UNIX AND NOT APPLE)
    set(MATH_LIBRARY m)
else()
    set(MATH_LIBRARY "")
endif()

# -----------------------------------------------------------------------------
# Dependencies: MPI & FFTW3
# -----------------------------------------------------------------------------
find_package(MPI REQUIRED)

find_package(FFTW3 REQUIRED COMPONENTS MPI)

# -----------------------------------------------------------------------------
# CUDA GPU Acceleration Support
# -----------------------------------------------------------------------------
option(ENABLE_CUDA "Enable CUDA GPU acceleration if available" ON)
option(REQUIRE_CUDA "Fail configuration if CUDA compiler is not found" OFF)

set(CUDA_SOURCES "")
set(EXTRA_CUDA_LIBS "")
set(CUDA_ACTIVE FALSE)

if(ENABLE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        set(CMAKE_CUDA_STANDARD 17)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
        if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
            if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.18")
                set(CMAKE_CUDA_ARCHITECTURES "all" CACHE STRING "CUDA architectures")
            endif()
        endif()
        add_compile_definitions(MCCH_ENABLE_CUDA)
        set(CUDA_SOURCES src/CudaKernels.cu)
        find_package(CUDAToolkit QUIET)
        if(TARGET CUDA::cudart)
            set(EXTRA_CUDA_LIBS CUDA::cudart)
        endif()
        set(CUDA_ACTIVE TRUE)
        message(STATUS "CUDA acceleration ENABLED with compiler: ${CMAKE_CUDA_COMPILER}")
    else()
        if(REQUIRE_CUDA)
            message(FATAL_ERROR "CUDA acceleration is required (REQUIRE_CUDA=ON), but nvcc was not found.")
        else()
            message(STATUS "CUDA compiler (nvcc) not found. Building CPU-only.")
        endif()
    endif()
else()
    message(STATUS "CUDA acceleration disabled by configuration. Building CPU-only.")
endif()

# -----------------------------------------------------------------------------
# Core Library: mcch_core
# -----------------------------------------------------------------------------
add_library(mcch_core
    src/Config.cpp
    tdbs/Thermo.c
    ${CUDA_SOURCES}
)

target_include_directories(mcch_core PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/tdbs>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/mcch>
)

target_link_libraries(mcch_core PUBLIC
    MPI::MPI_CXX
    FFTW3::FFTW3
    ${EXTRA_CUDA_LIBS}
    ${MATH_LIBRARY}
)

# -----------------------------------------------------------------------------
# Main Executable: mcch_solver
# -----------------------------------------------------------------------------
add_executable(mcch_solver
    src/main.cpp
)
target_link_libraries(mcch_solver PRIVATE mcch_core)

# -----------------------------------------------------------------------------
# Testing Suite
# -----------------------------------------------------------------------------
include(CTest)

if(BUILD_TESTING)
    # Check if MPI runner supports --oversubscribe (e.g. OpenMPI on low core counts)
    set(MPIEXEC_OVERSUBSCRIBE_FLAG "")
    if(MPIEXEC_EXECUTABLE)
        execute_process(
            COMMAND ${MPIEXEC_EXECUTABLE} --help
            OUTPUT_VARIABLE _MPIEXEC_HELP_TEXT
            ERROR_QUIET
        )
        if(_MPIEXEC_HELP_TEXT MATCHES "--oversubscribe")
            set(MPIEXEC_OVERSUBSCRIBE_FLAG "--oversubscribe")
        endif()
    endif()

    # Test 1: Ghost cell exchange & slab partitioning verification
    add_executable(test_slab_mpi tests/test_slab_mpi.cpp)
    target_link_libraries(test_slab_mpi PRIVATE mcch_core)
    add_test(NAME TestSlabMPI_1Rank COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_slab_mpi>)
    add_test(NAME TestSlabMPI_2Ranks COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_slab_mpi>)
    add_test(NAME TestSlabMPI_4Ranks COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_OVERSUBSCRIBE_FLAG} $<TARGET_FILE:test_slab_mpi>)

    # Test 2: Mass conservation & partition of unity
    add_executable(test_mass_conservation tests/test_mass_conservation.cpp)
    target_link_libraries(test_mass_conservation PRIVATE mcch_core)
    add_test(NAME TestMassConservation_2Ranks COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_mass_conservation>)

    # Test 3: Thermodynamic consistency (dF/dt <= 0)
    add_executable(test_free_energy_decrease tests/test_free_energy_decrease.cpp)
    target_link_libraries(test_free_energy_decrease PRIVATE mcch_core)
    add_test(NAME TestEnergyDecrease_2Ranks COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_free_energy_decrease>)

    # Test 4: Multi-component physics (ternary & quaternary)
    add_executable(test_multicomponent tests/test_multicomponent.cpp)
    target_link_libraries(test_multicomponent PRIVATE mcch_core)
    add_test(NAME TestMultiComponent_2Ranks COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_multicomponent>)

    # Test 5: Semi-implicit Fourier spectral solver with FFTW and MPI
    add_executable(test_semi_implicit tests/test_semi_implicit.cpp)
    target_link_libraries(test_semi_implicit PRIVATE mcch_core)
    add_test(NAME TestSemiImplicit_1Rank COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_semi_implicit>)
    add_test(NAME TestSemiImplicit_2Ranks COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_semi_implicit>)
    add_test(NAME TestSemiImplicit_4Ranks COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_OVERSUBSCRIBE_FLAG} $<TARGET_FILE:test_semi_implicit>)

    # Test 6: GPU execution verification (Explicit Euler, RK2, RK4)
    add_executable(test_gpu tests/test_gpu.cpp)
    target_link_libraries(test_gpu PRIVATE mcch_core)
    add_test(NAME TestGPU_1Rank COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_gpu>)

    # Test 7: Redlich-Kister thermodynamic function & solver verification
    add_executable(test_redlich_kister tests/test_redlich_kister.cpp)
    target_link_libraries(test_redlich_kister PRIVATE mcch_core)
    add_test(NAME TestRedlichKister_1Rank COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_redlich_kister>)
    add_test(NAME TestRedlichKister_2Ranks COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_redlich_kister>)
endif()

# -----------------------------------------------------------------------------
# Installation Rules (GNU Standard Directories)
# -----------------------------------------------------------------------------
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(MCCH_CONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/MCCH")

install(TARGETS mcch_solver mcch_core
    EXPORT MCCHTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mcch
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mcch
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
)

install(FILES tdbs/Thermo.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mcch
)

install(DIRECTORY examples/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/mcch/examples
    FILES_MATCHING PATTERN "*.json"
)

install(DIRECTORY tdbs/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/mcch/tdbs
    FILES_MATCHING PATTERN "*.TDB" PATTERN "*.h" PATTERN "*.c"
)

install(DIRECTORY python/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/mcch/python
    FILES_MATCHING PATTERN "*.py"
)

# Export CMake Package Configuration
install(EXPORT MCCHTargets
    FILE MCCHTargets.cmake
    NAMESPACE MCCH::
    DESTINATION ${MCCH_CONFIG_INSTALL_DIR}
)

configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/MCCHConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/MCCHConfig.cmake"
    INSTALL_DESTINATION ${MCCH_CONFIG_INSTALL_DIR}
    PATH_VARS CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR CMAKE_INSTALL_BINDIR
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/MCCHConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/MCCHConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/MCCHConfigVersion.cmake"
    "${PROJECT_SOURCE_DIR}/cmake/FindFFTW3.cmake"
    DESTINATION ${MCCH_CONFIG_INSTALL_DIR}
)

# -----------------------------------------------------------------------------
# CPack Packaging Support (.tar.gz, .deb, .rpm, .zip)
# -----------------------------------------------------------------------------
set(CPACK_PACKAGE_NAME "mcch-solver")
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Generalized Multi-Component Cahn-Hilliard Phase-Field Solver")
set(CPACK_PACKAGE_VENDOR "MCCH Team")
set(CPACK_PACKAGE_CONTACT "mcch-solver@users.noreply.github.com")
set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.md")
set(CPACK_GENERATOR "TGZ;ZIP")
if(UNIX AND NOT APPLE)
    find_program(DPKG_CMD dpkg)
    if(DPKG_CMD)
        list(APPEND CPACK_GENERATOR "DEB")
        set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
    endif()
    find_program(RPMBUILD_CMD rpmbuild)
    if(RPMBUILD_CMD)
        list(APPEND CPACK_GENERATOR "RPM")
    endif()
endif()

# Exclude build directories, caches, and VCS files from CPack source packages
set(CPACK_SOURCE_IGNORE_FILES
    "/build/"
    "/\\\\.git/"
    "/\\\\.github/"
    "/\\\\.venv/"
    "/venv/"
    "/__pycache__/"
    "/.*\\\\.egg-info/"
    "/results.*"
    "/Testing/"
    "/\\\\.cache/"
    "\\\\.swp$"
    "\\\\.DS_Store$"
    "CMakeCache\\\\.txt"
    "CMakeFiles"
)

include(CPack)

# -----------------------------------------------------------------------------
# Distclean Target (Purges CMake build cache and generated files)
# -----------------------------------------------------------------------------
add_custom_target(distclean
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}/CMakeFiles"
    COMMAND ${CMAKE_COMMAND} -E remove -f "${CMAKE_BINARY_DIR}/CMakeCache.txt"
    COMMAND ${CMAKE_COMMAND} -E remove -f "${CMAKE_BINARY_DIR}/cmake_install.cmake"
    COMMAND ${CMAKE_COMMAND} -E remove -f "${CMAKE_BINARY_DIR}/Makefile"
    COMMAND ${CMAKE_COMMAND} -E remove -f "${CMAKE_BINARY_DIR}/CTestTestfile.cmake"
    COMMAND ${CMAKE_COMMAND} -E remove -f "${CMAKE_BINARY_DIR}/CPackConfig.cmake"
    COMMAND ${CMAKE_COMMAND} -E remove -f "${CMAKE_BINARY_DIR}/CPackSourceConfig.cmake"
    COMMENT "Purging CMake build cache and generated files for clean portability..."
)
