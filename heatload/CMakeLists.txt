set(CMAKE_CXX_STANDARD 14)

find_package (MPI COMPONENTS CXX)
find_package(ADIOS2 REQUIRED COMPONENTS CXX MPI)
find_package(OpenMP)

set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED OFF)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost COMPONENTS log log_setup program_options filesystem)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
    load.cpp
    heatload_calc.cpp
    initialize.cpp
    separate.cpp
    particles.cpp
    output.cpp
    heatload.cpp
)

add_library(heatloadlib ${SOURCES})
add_executable (heatload main.cpp)

if (OpenMP_FOUND)
    target_compile_definitions(heatloadlib PUBLIC USE_OMP)
    target_link_libraries(heatloadlib OpenMP::OpenMP_CXX)
    target_link_libraries(heatload OpenMP::OpenMP_CXX)
endif()

target_link_libraries(heatloadlib adios2::adios2 MPI::MPI_C)
target_link_libraries(heatloadlib ${Boost_LIBRARIES})

target_link_libraries(heatload heatloadlib adios2::adios2 MPI::MPI_C)
target_link_libraries(heatload ${Boost_LIBRARIES})
