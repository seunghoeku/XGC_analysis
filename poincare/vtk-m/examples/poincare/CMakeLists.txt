##============================================================================
##  Copyright (c) Kitware, Inc.
##  All rights reserved.
##  See LICENSE.txt for details.
##
##  This software is distributed WITHOUT ANY WARRANTY; without even
##  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##  PURPOSE.  See the above copyright notice for more information.
##============================================================================
cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)
project(Poincare CXX)

#Find the VTK-m package
find_package(VTKm REQUIRED QUIET)
find_package(ADIOS2 REQUIRED)

set(Poincare_srcs
 Poincare.cxx
 FindMaxR.cxx
 RunPoincare2.cxx
 EvalField.cxx
 perfstubs_api/timer.c
 )

add_executable(Poincare ${Poincare_srcs})

target_link_libraries(Poincare PRIVATE vtkm_cont vtkm_worklet vtkm_io adios2::adios2)
target_include_directories(Poincare PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

vtkm_add_target_information(Poincare
                            DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
                            DEVICE_SOURCES ${Poincare_srcs})

if(TARGET vtkm::tbb)
  target_compile_definitions(Poincare PRIVATE BUILDING_TBB_VERSION)
endif()
