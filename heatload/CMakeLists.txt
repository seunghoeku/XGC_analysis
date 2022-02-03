find_package(ADIOS2 2.7.1 REQUIRED COMPONENTS CXX)

add_executable (heatload main.cpp)
target_link_libraries(heatload adios2::adios2)