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
#find_package(Fides REQUIRED)

#add_executable(Poincare Poincare.cxx io.cxx io.h)
#add_executable(Simple Simple.cxx)
#add_executable(Simple2 Simple2.cxx)
#add_executable(Simple3 Simple3.cxx)
#add_executable(Simple4 Simple4.cxx)
#add_executable(Simple2.2 Simple2.2.cxx)
add_executable(Simple2.3 Simple2.3.cxx)
#add_executable(SimpleFides SimpleFides.cxx)
#target_link_libraries(Poincare PRIVATE vtkm_filter vtkm_io fides)
#target_link_libraries(Simple PRIVATE vtkm_filter vtkm_io vtkm_rendering fides)
#target_link_libraries(Simple2 PRIVATE vtkm_filter vtkm_io vtkm_rendering fides)
#target_link_libraries(Simple3 PRIVATE vtkm_filter vtkm_io vtkm_rendering fides)
#target_link_libraries(Simple4 PRIVATE vtkm_filter vtkm_io vtkm_rendering fides)
#target_link_libraries(Simple2.2 PRIVATE vtkm_filter vtkm_io vtkm_rendering fides)
##target_link_libraries(Simple2.3 PRIVATE vtkm_filter vtkm_io vtkm_rendering fides)
target_link_libraries(Simple2.3 PRIVATE vtkm_filter vtkm_io adios2::adios2)
#target_link_libraries(SimpleFides PRIVATE vtkm_filter vtkm_io fides)
#vtkm_add_target_information(Poincare
#                            DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
#                            DEVICE_SOURCES Poincare.cxx)
#vtkm_add_target_information(Simple
#                            DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
#                            DEVICE_SOURCES Simple.cxx)
#vtkm_add_target_information(Simple2
#                            DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
#                            DEVICE_SOURCES Simple2.cxx)
#vtkm_add_target_information(Simple3
#                            DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
#                            DEVICE_SOURCES Simple3.cxx)
#vtkm_add_target_information(Simple4
#                            DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
#                            DEVICE_SOURCES Simple4.cxx)
#vtkm_add_target_information(Simple2.2
#                            DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
#                            DEVICE_SOURCES Simple2.2.cxx)
vtkm_add_target_information(Simple2.3
                            DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
                            DEVICE_SOURCES Simple2.3.cxx)
#vtkm_add_target_information(SimpleFides
#                            DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
#                            DEVICE_SOURCES SimpleFides.cxx)
if(TARGET vtkm::tbb)
  target_compile_definitions(Poincare PRIVATE BUILDING_TBB_VERSION)
  target_compile_definitions(Simple PRIVATE BUILDING_TBB_VERSION)
  target_compile_definitions(Simple2 PRIVATE BUILDING_TBB_VERSION)
  target_compile_definitions(Simple3 PRIVATE BUILDING_TBB_VERSION)
  target_compile_definitions(Simple4 PRIVATE BUILDING_TBB_VERSION)
endif()
