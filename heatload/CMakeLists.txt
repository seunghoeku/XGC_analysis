find_package(ADIOS2 2.7.1 REQUIRED COMPONENTS CXX)
find_package(OpenMP)

set(SOURCES
    main.cpp
    load.cpp
    heatload_calc.cpp
    initialize.cpp
    separate.cpp
    particles.cpp
    output.cpp
)

add_executable (heatload ${SOURCES})
if (OpenMP_FOUND)
    target_compile_definitions(heatload PUBLIC USE_OMP)
    target_link_libraries(heatload adios2::cxx11 OpenMP::OpenMP_CXX)
else()
    target_link_libraries(heatload adios2::cxx11)
endif()
