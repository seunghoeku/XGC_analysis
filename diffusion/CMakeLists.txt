set(CMAKE_CXX_STANDARD 14)

find_package(MPI COMPONENTS CXX)
find_package(ADIOS2 REQUIRED COMPONENTS CXX MPI)
find_package(OpenMP)
find_package(CAMTIMERS)
message(STATUS ${CAMTIMERS_INCLUDE_DIRS})
message(STATUS ${CAMTIMERS_LIBRARIES})
if(CAMTIMERS_FOUND)
    include_directories(${CAMTIMERS_INCLUDE_DIRS})
    add_definitions(-DCAM_TIMERS)
    add_definitions(-DHAVE_MPI)
endif()

set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED OFF)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost COMPONENTS log log_setup program_options filesystem)
message(STATUS "Boost_INCLUDE_DIRS: " "${Boost_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES: " "${Boost_LIBRARIES}")
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
    diffusion_main.cpp
)

include_directories(${PROJECT_SOURCE_DIR}/middleman)
include_directories(${PROJECT_SOURCE_DIR}/util)

add_executable(diffusion diffusion_main.cpp ${PROJECT_SOURCE_DIR}/middleman/diffusion.cpp)
target_link_libraries(diffusion adios2::adios2 MPI::MPI_C)
target_link_libraries(diffusion ${Boost_LIBRARIES})

if(CAMTIMERS_FOUND)
    target_link_libraries(diffusion CAMTIMERS::CAMTIMERS)
endif()

if (OpenMP_FOUND)
    target_link_libraries(diffusion OpenMP::OpenMP_CXX)
endif()


