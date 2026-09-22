#[=======================================================================[.rst:
FindFFTW3
---------

Find the FFTW3 (Fastest Fourier Transform in the West) library and its MPI
component.

Hints and Environment Variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This module accepts the following variables:
  ``FFTW3_ROOT`` or ``FFTW_ROOT``
  ``FFTW3_DIR`` or ``FFTW_DIR``
  ``ENV{CONDA_PREFIX}``
  ``ENV{SPACK_ROOT}``

Components
^^^^^^^^^^
The following components may be requested:
  ``MPI``     - Find FFTW3 MPI library (fftw3_mpi and fftw3-mpi.h)
  ``threads`` - Find FFTW3 pthreads library (fftw3_threads)
  ``omp``     - Find FFTW3 OpenMP library (fftw3_omp)

Result Variables
^^^^^^^^^^^^^^^^
This module defines the following variables:
  ``FFTW3_FOUND``        - True if FFTW3 and all requested components were found
  ``FFTW3_INCLUDE_DIRS`` - Include directories for FFTW3
  ``FFTW3_LIBRARIES``    - Libraries needed to link FFTW3
  ``FFTW3_VERSION``      - Detected FFTW3 version

Imported Targets
^^^^^^^^^^^^^^^^
This module provides the following imported targets:
  ``FFTW3::fftw3``      - Core FFTW3 library
  ``FFTW3::fftw3_mpi``  - FFTW3 MPI library (if MPI component requested and found)
  ``FFTW3::FFTW3``      - Interface target containing all found components
#]=======================================================================]

include(FindPackageHandleStandardArgs)

# Search hints / paths
set(_FFTW3_HINTS
    ${FFTW3_ROOT}
    ${FFTW_ROOT}
    ${FFTW3_DIR}
    ${FFTW_DIR}
    $ENV{FFTW3_ROOT}
    $ENV{FFTW_ROOT}
    $ENV{FFTW3_DIR}
    $ENV{FFTW_DIR}
    $ENV{CONDA_PREFIX}
    $ENV{SPACK_ROOT}
    $ENV{HOME}/.local
    /opt/homebrew
    /usr/local
    /opt/local
    /usr
)

set(_FFTW3_INCLUDE_HINTS)
set(_FFTW3_LIB_HINTS)
foreach(_hint ${_FFTW3_HINTS})
    if(_hint)
        list(APPEND _FFTW3_INCLUDE_HINTS "${_hint}/include")
        list(APPEND _FFTW3_LIB_HINTS "${_hint}/lib")
        if(CMAKE_LIBRARY_ARCHITECTURE)
            list(APPEND _FFTW3_LIB_HINTS "${_hint}/lib/${CMAKE_LIBRARY_ARCHITECTURE}")
        endif()
        list(APPEND _FFTW3_LIB_HINTS "${_hint}/lib64")
    endif()
endforeach()

# Try PkgConfig for additional hints
find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(_PC_FFTW3 QUIET fftw3)
    pkg_check_modules(_PC_FFTW3_MPI QUIET fftw3_mpi)
    if(_PC_FFTW3_INCLUDEDIR)
        list(APPEND _FFTW3_INCLUDE_HINTS ${_PC_FFTW3_INCLUDEDIR})
    endif()
    if(_PC_FFTW3_LIBDIR)
        list(APPEND _FFTW3_LIB_HINTS ${_PC_FFTW3_LIBDIR})
    endif()
    if(_PC_FFTW3_MPI_INCLUDEDIR)
        list(APPEND _FFTW3_INCLUDE_HINTS ${_PC_FFTW3_MPI_INCLUDEDIR})
    endif()
    if(_PC_FFTW3_MPI_LIBDIR)
        list(APPEND _FFTW3_LIB_HINTS ${_PC_FFTW3_MPI_LIBDIR})
    endif()
endif()

# Find core header fftw3.h
find_path(FFTW3_INCLUDE_DIR
    NAMES fftw3.h
    HINTS ${_FFTW3_INCLUDE_HINTS}
    PATH_SUFFIXES include
)

# Find core library fftw3
find_library(FFTW3_LIBRARY
    NAMES fftw3 libfftw3
    HINTS ${_FFTW3_LIB_HINTS}
    PATH_SUFFIXES lib lib64
)

# Extract version from fftw3.h
set(FFTW3_VERSION "")
if(FFTW3_INCLUDE_DIR AND EXISTS "${FFTW3_INCLUDE_DIR}/fftw3.h")
    file(STRINGS "${FFTW3_INCLUDE_DIR}/fftw3.h" _fftw_version_line REGEX "^#define[ \t]+FFTW_VERSION[ \t]+\".*\"")
    if(_fftw_version_line)
        string(REGEX REPLACE "^#define[ \t]+FFTW_VERSION[ \t]+\"([^\"]+)\".*" "\\1" FFTW3_VERSION "${_fftw_version_line}")
    endif()
endif()

set(_FFTW3_REQ_VARS FFTW3_LIBRARY FFTW3_INCLUDE_DIR)

# Component: MPI
set(FFTW3_MPI_FOUND FALSE)
if("MPI" IN_LIST FFTW3_FIND_COMPONENTS OR FFTW3_FIND_REQUIRED_MPI)
    find_path(FFTW3_MPI_INCLUDE_DIR
        NAMES fftw3-mpi.h
        HINTS ${FFTW3_INCLUDE_DIR} ${_FFTW3_INCLUDE_HINTS}
        PATH_SUFFIXES include
    )
    find_library(FFTW3_MPI_LIBRARY
        NAMES fftw3_mpi libfftw3_mpi
        HINTS ${_FFTW3_LIB_HINTS}
        PATH_SUFFIXES lib lib64
    )
    if(FFTW3_MPI_LIBRARY AND FFTW3_MPI_INCLUDE_DIR)
        set(FFTW3_MPI_FOUND TRUE)
    endif()
    list(APPEND _FFTW3_REQ_VARS FFTW3_MPI_LIBRARY FFTW3_MPI_INCLUDE_DIR)
endif()

# Component: threads
set(FFTW3_threads_FOUND FALSE)
if("threads" IN_LIST FFTW3_FIND_COMPONENTS)
    find_library(FFTW3_THREADS_LIBRARY
        NAMES fftw3_threads libfftw3_threads
        HINTS ${_FFTW3_LIB_HINTS}
        PATH_SUFFIXES lib lib64
    )
    if(FFTW3_THREADS_LIBRARY)
        set(FFTW3_threads_FOUND TRUE)
    endif()
    list(APPEND _FFTW3_REQ_VARS FFTW3_THREADS_LIBRARY)
endif()

# Component: omp
set(FFTW3_omp_FOUND FALSE)
if("omp" IN_LIST FFTW3_FIND_COMPONENTS)
    find_library(FFTW3_OMP_LIBRARY
        NAMES fftw3_omp libfftw3_omp
        HINTS ${_FFTW3_LIB_HINTS}
        PATH_SUFFIXES lib lib64
    )
    if(FFTW3_OMP_LIBRARY)
        set(FFTW3_omp_FOUND TRUE)
    endif()
    list(APPEND _FFTW3_REQ_VARS FFTW3_OMP_LIBRARY)
endif()

# Standard handling
find_package_handle_standard_args(FFTW3
    REQUIRED_VARS ${_FFTW3_REQ_VARS}
    VERSION_VAR FFTW3_VERSION
    HANDLE_COMPONENTS
)

if(FFTW3_FOUND)
    set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})
    if(FFTW3_MPI_INCLUDE_DIR AND NOT FFTW3_MPI_INCLUDE_DIR STREQUAL FFTW3_INCLUDE_DIR)
        list(APPEND FFTW3_INCLUDE_DIRS ${FFTW3_MPI_INCLUDE_DIR})
    endif()

    set(FFTW3_LIBRARIES ${FFTW3_LIBRARY})
    if(FFTW3_MPI_LIBRARY)
        list(APPEND FFTW3_LIBRARIES ${FFTW3_MPI_LIBRARY})
    endif()
    if(FFTW3_THREADS_LIBRARY)
        list(APPEND FFTW3_LIBRARIES ${FFTW3_THREADS_LIBRARY})
    endif()
    if(FFTW3_OMP_LIBRARY)
        list(APPEND FFTW3_LIBRARIES ${FFTW3_OMP_LIBRARY})
    endif()

    # Imported targets
    if(NOT TARGET FFTW3::fftw3)
        add_library(FFTW3::fftw3 UNKNOWN IMPORTED)
        set_target_properties(FFTW3::fftw3 PROPERTIES
            IMPORTED_LOCATION "${FFTW3_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}"
        )
    endif()

    if(FFTW3_MPI_FOUND AND NOT TARGET FFTW3::fftw3_mpi)
        add_library(FFTW3::fftw3_mpi UNKNOWN IMPORTED)
        set_target_properties(FFTW3::fftw3_mpi PROPERTIES
            IMPORTED_LOCATION "${FFTW3_MPI_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_MPI_INCLUDE_DIR}"
        )
    endif()

    if(NOT TARGET FFTW3::FFTW3)
        add_library(FFTW3::FFTW3 INTERFACE IMPORTED)
        set_target_properties(FFTW3::FFTW3 PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${FFTW3_LIBRARIES}"
        )
    endif()
endif()

mark_as_advanced(
    FFTW3_INCLUDE_DIR
    FFTW3_LIBRARY
    FFTW3_MPI_INCLUDE_DIR
    FFTW3_MPI_LIBRARY
    FFTW3_THREADS_LIBRARY
    FFTW3_OMP_LIBRARY
)
