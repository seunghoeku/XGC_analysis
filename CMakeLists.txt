cmake_minimum_required(VERSION 3.16)
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/CMake)

project (XGC_GB_Analysis)

add_subdirectory (diffusion)
add_subdirectory (heatload)
# add_subdirectory (middleman)
