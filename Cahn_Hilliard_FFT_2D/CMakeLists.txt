cmake_minimum_required(VERSION 3.18)
project(microsim_ch_fft LANGUAGES C)

# Add the top-level cmake/ directory to the module path so that
# FindFFTW3.cmake is found whether this is built stand-alone or as part
# of the top-level project.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

# ── Dependencies ────────────────────────────────────────────────────────────
find_package(FFTW3 REQUIRED)   # provided by cmake/FindFFTW3.cmake

# ── Main solver target ───────────────────────────────────────────────────────
add_executable(microsim_ch_fft microsim_ch_fft.c)

target_include_directories(microsim_ch_fft PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/functions
    ${CMAKE_CURRENT_SOURCE_DIR}/solverloop
    ${FFTW3_INCLUDE_DIRS}
)

target_link_libraries(microsim_ch_fft PRIVATE
    FFTW3::fftw3
    m
)

install(TARGETS microsim_ch_fft RUNTIME DESTINATION bin)
