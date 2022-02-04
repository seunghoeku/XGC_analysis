find_package(ADIOS2 2.7.1 REQUIRED COMPONENTS CXX)

set(SOURCES
    main.cpp
    load.cpp
    heatload_calc.cpp
    initialize.cpp
)
add_executable (heatload ${SOURCES})
target_link_libraries(heatload adios2::adios2)