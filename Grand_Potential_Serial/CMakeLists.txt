cmake_minimum_required(VERSION 3.18)
project(microsim_gp_serial LANGUAGES C)

# Source files – the solver is a single translation unit; all other files
# are header-only and referenced via include paths.
add_executable(microsim_gp microsim_gp.c)

target_include_directories(microsim_gp PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/functions
    ${CMAKE_CURRENT_SOURCE_DIR}/solverloop
)

target_link_libraries(microsim_gp PRIVATE m)

install(TARGETS microsim_gp RUNTIME DESTINATION bin)
