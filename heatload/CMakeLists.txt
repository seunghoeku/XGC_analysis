find_package(ADIOS2 2.7.1 REQUIRED COMPONENTS CXX)
find_package(OpenMP)

set(SOURCES
    main.cpp
    load.cpp
    heatload_calc.cpp
    initialize.cpp
    separate.cpp
    particles.cpp
)

add_executable (heatload ${SOURCES})
target_link_libraries(heatload adios2::adios2 OpenMP::OpenMP_CXX)
