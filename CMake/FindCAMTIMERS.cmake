# Distributed under the OSI-approved BSD 3-Clause License.

#.rst:
# FindCAMTIMERS
# ---------
#
# Try to find CAMTIMERS
#
# Uses CAMTIMERS_ROOT in the cache variables or in the environment as a hint
# where to search
#
# IMPORTED Targets
# ^^^^^^^^^^^^^^^^
#
# This module defines :prop_tgt:`IMPORTED` target ``CAMTIMERS::CAMTIMERS``, if
# CAMTIMERS has been found.
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module defines the following variables:
#
# ``CAMTIMERS_FOUND``
#   system has CAMTIMERS
# ``CAMTIMERS_INCLUDE_DIRS``
#   the CAMTIMERS include directories
# ``CAMTIMERS_LIBRARIES``
#   Link these to use CAMTIMERS

find_path(CAMTIMERS_INCLUDE_DIRS NAMES gptl.h)
find_library(CAMTIMERS_LIBRARIES NAMES timers)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CAMTIMERS
                                  REQUIRED_VARS CAMTIMERS_LIBRARIES CAMTIMERS_INCLUDE_DIRS)

if (CAMTIMERS_FOUND)
    if(NOT TARGET CAMTIMERS::CAMTIMERS)
      add_library(CAMTIMERS::CAMTIMERS INTERFACE IMPORTED)
      set_property(TARGET CAMTIMERS::CAMTIMERS PROPERTY
        INTERFACE_INCLUDE_DIRECTORIES "${CAMTIMERS_INCLUDE_DIRS}")
      set_property(TARGET CAMTIMERS::CAMTIMERS PROPERTY
        INTERFACE_LINK_LIBRARIES "${CAMTIMERS_LIBRARIES};${CMAKE_DL_LIBS}")
    endif()
endif()
