# FindFFTW3.cmake
#
# Finds the FFTW3 double-precision library.
#
# Imported target:
#   FFTW3::fftw3
#
# Result variables:
#   FFTW3_FOUND         - TRUE if FFTW3 was found
#   FFTW3_INCLUDE_DIRS  - include directory
#   FFTW3_LIBRARIES     - libraries to link against
#
# Search hints (can be set by the user before calling find_package):
#   FFTW3_ROOT           - root installation prefix

# --- header ---
find_path(FFTW3_INCLUDE_DIR
    NAMES fftw3.h
    HINTS
        ${FFTW3_ROOT}
        ENV FFTW3_ROOT
        ENV FFTW3_DIR
    PATH_SUFFIXES include
)

# --- library ---
find_library(FFTW3_LIBRARY
    NAMES fftw3
    HINTS
        ${FFTW3_ROOT}
        ENV FFTW3_ROOT
        ENV FFTW3_DIR
    PATH_SUFFIXES lib lib64
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3
    REQUIRED_VARS FFTW3_LIBRARY FFTW3_INCLUDE_DIR
)

if(FFTW3_FOUND)
    set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})
    set(FFTW3_LIBRARIES    ${FFTW3_LIBRARY})

    if(NOT TARGET FFTW3::fftw3)
        add_library(FFTW3::fftw3 UNKNOWN IMPORTED)
        set_target_properties(FFTW3::fftw3 PROPERTIES
            IMPORTED_LOCATION             "${FFTW3_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}"
        )
    endif()
endif()

mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_LIBRARY)
