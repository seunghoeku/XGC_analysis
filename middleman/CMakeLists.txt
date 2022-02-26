set(CMAKE_CXX_STANDARD 14)

find_package (MPI COMPONENTS CXX)
find_package (ADIOS2 REQUIRED COMPONENTS CXX MPI)
find_package (OpenMP)

set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED OFF)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost COMPONENTS log log_setup)
message(STATUS "Boost_INCLUDE_DIRS: " "${Boost_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES: " "${Boost_LIBRARIES}")
include_directories(${Boost_INCLUDE_DIRS})

add_executable (middleman middleman.cpp)
target_link_libraries (middleman adios2::adios2 MPI::MPI_C)
target_link_libraries (middleman ${Boost_LIBRARIES})

