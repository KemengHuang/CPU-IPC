# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2023 Google Inc. All rights reserved.
# http://ceres-solver.org/
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of Google Inc. nor the names of its contributors may be
#   used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# Author: alexs.mac@gmail.com (Alex Stewart)
#

#[=======================================================================[.rst:
FindSuiteSparse
===============

Module for locating SuiteSparse libraries and their dependencies.

This module defines the following variables:

``SuiteSparse_FOUND``
   ``TRUE`` if SuiteSparse and all dependencies have been found.

``SuiteSparse_VERSION``
   Extracted from ``SuiteSparse_config.h`` (>= v4).

``SuiteSparse_VERSION_MAJOR``
    Equal to 4 if ``SuiteSparse_VERSION`` = 4.2.1

``SuiteSparse_VERSION_MINOR``
    Equal to 2 if ``SuiteSparse_VERSION`` = 4.2.1

``SuiteSparse_VERSION_PATCH``
    Equal to 1 if ``SuiteSparse_VERSION`` = 4.2.1

The following variables control the behaviour of this module:

``SuiteSparse_NO_CMAKE``
  Do not attempt to use the native SuiteSparse CMake package configuration.


Targets
-------

The following targets define the SuiteSparse components searched for.

``SuiteSparse::AMD``
    Symmetric Approximate Minimum Degree (AMD)

``SuiteSparse::CAMD``
    Constrained Approximate Minimum Degree (CAMD)

``SuiteSparse::COLAMD``
    Column Approximate Minimum Degree (COLAMD)

``SuiteSparse::CCOLAMD``
    Constrained Column Approximate Minimum Degree (CCOLAMD)

``SuiteSparse::CHOLMOD``
    Sparse Supernodal Cholesky Factorization and Update/Downdate (CHOLMOD)

``SuiteSparse::Partition``
    CHOLMOD with METIS support

``SuiteSparse::SPQR``
    Multifrontal Sparse QR (SuiteSparseQR)

``SuiteSparse::Config``
    Common configuration for all but CSparse (SuiteSparse version >= 4).

Optional SuiteSparse dependencies:

``METIS::METIS``
    Serial Graph Partitioning and Fill-reducing Matrix Ordering (METIS)
]=======================================================================]

# Prefer an upstream SuiteSparse config package if available and allowed.
if(NOT SuiteSparse_NO_CMAKE)
  find_package(SuiteSparse NO_MODULE QUIET)
endif()

if(SuiteSparse_FOUND)
  return()
endif()

# Push CMP0057 to enable support for IN_LIST when cmake_minimum_required is < 3.3.
cmake_policy(PUSH)
cmake_policy(SET CMP0057 NEW)

if(NOT SuiteSparse_FIND_COMPONENTS)
  set(SuiteSparse_FIND_COMPONENTS
    AMD
    CAMD
    CCOLAMD
    CHOLMOD
    COLAMD
  )

  foreach(component IN LISTS SuiteSparse_FIND_COMPONENTS)
    set(SuiteSparse_FIND_REQUIRED_${component} TRUE)
  endforeach()
endif()

# Assume SuiteSparse was found and set it to false only if third-party
# dependencies could not be located. SuiteSparse components are handled by
# FindPackageHandleStandardArgs HANDLE_COMPONENTS option.
set(SuiteSparse_FOUND TRUE)

include(CheckLibraryExists)
include(CheckSymbolExists)
include(CMakePushCheckState)

# Config is a base component and thus always required.
set(SuiteSparse_IMPLICIT_COMPONENTS Config)

# CHOLMOD depends on AMD, CAMD, CCOLAMD, and COLAMD.
if(CHOLMOD IN_LIST SuiteSparse_FIND_COMPONENTS)
  list(APPEND SuiteSparse_IMPLICIT_COMPONENTS AMD CAMD CCOLAMD COLAMD)
endif()

# SPQR depends on CHOLMOD.
if(SPQR IN_LIST SuiteSparse_FIND_COMPONENTS)
  list(APPEND SuiteSparse_IMPLICIT_COMPONENTS CHOLMOD)
endif()

# Implicit components are always required.
foreach(component IN LISTS SuiteSparse_IMPLICIT_COMPONENTS)
  set(SuiteSparse_FIND_REQUIRED_${component} TRUE)
endforeach()

list(APPEND SuiteSparse_FIND_COMPONENTS ${SuiteSparse_IMPLICIT_COMPONENTS})

# Do not list components multiple times.
list(REMOVE_DUPLICATES SuiteSparse_FIND_COMPONENTS)

# Reset CALLERS_CMAKE_FIND_LIBRARY_PREFIXES to its value when
# FindSuiteSparse was invoked.
macro(SuiteSparse_RESET_FIND_LIBRARY_PREFIX)
  if(MSVC)
    set(CMAKE_FIND_LIBRARY_PREFIXES "${CALLERS_CMAKE_FIND_LIBRARY_PREFIXES}")
  endif()
endmacro()

# Called if we failed to find SuiteSparse or any of its required dependencies.
# Unsets all public variables and reports an error message depending on the
# REQUIRED/QUIET argument.
macro(SuiteSparse_REPORT_NOT_FOUND REASON_MSG)
  # Will be set to FALSE by find_package_handle_standard_args.
  unset(SuiteSparse_FOUND)

  # Do NOT unset SuiteSparse_REQUIRED_VARS here, as it is used by
  # FindPackageHandleStandardArgs() to generate the automatic error message on
  # failure which highlights which components are missing.

  suitesparse_reset_find_library_prefix()

  if(SuiteSparse_FIND_QUIETLY)
    message(STATUS "Failed to find SuiteSparse - " ${REASON_MSG} ${ARGN})
  elseif(SuiteSparse_FIND_REQUIRED)
    message(FATAL_ERROR "Failed to find SuiteSparse - " ${REASON_MSG} ${ARGN})
  else()
    message("-- Failed to find SuiteSparse - " ${REASON_MSG} ${ARGN})
  endif()

  # Do not call return() so that we keep processing if not called with REQUIRED
  # and report all missing components, rather than bailing after failing to find
  # the first.
endmacro()

# Handle possible presence of lib prefix for libraries on MSVC.
if(MSVC)
  set(CALLERS_CMAKE_FIND_LIBRARY_PREFIXES "${CMAKE_FIND_LIBRARY_PREFIXES}")
  # The empty string in this list represents shared libraries / DLLs with no prefix.
  set(CMAKE_FIND_LIBRARY_PREFIXES "lib" "" "${CMAKE_FIND_LIBRARY_PREFIXES}")
endif()

# Additional suffixes to try appending to each search path.
list(APPEND SuiteSparse_CHECK_PATH_SUFFIXES suitesparse)

# Wrappers to find_path/library that pass the SuiteSparse search hints/paths.
#
# suitesparse_find_component(<component> [FILES name1 [name2 ...]]
#                                        [LIBRARIES name1 [name2 ...]])
macro(suitesparse_find_component COMPONENT)
  include(CMakeParseArguments)
  set(MULTI_VALUE_ARGS FILES LIBRARIES)
  cmake_parse_arguments(SuiteSparse_FIND_COMPONENT_${COMPONENT}
    "" "" "${MULTI_VALUE_ARGS}" ${ARGN})

  set(SuiteSparse_${COMPONENT}_FOUND TRUE)
  if(SuiteSparse_FIND_COMPONENT_${COMPONENT}_FILES)
    find_path(SuiteSparse_${COMPONENT}_INCLUDE_DIR
      NAMES ${SuiteSparse_FIND_COMPONENT_${COMPONENT}_FILES}
      PATH_SUFFIXES ${SuiteSparse_CHECK_PATH_SUFFIXES})
    if(SuiteSparse_${COMPONENT}_INCLUDE_DIR)
      message(STATUS "Found ${COMPONENT} headers in: ${SuiteSparse_${COMPONENT}_INCLUDE_DIR}")
      mark_as_advanced(SuiteSparse_${COMPONENT}_INCLUDE_DIR)
    else()
      set(SuiteSparse_${COMPONENT}_FOUND FALSE)
      if(SuiteSparse_FIND_REQUIRED_${COMPONENT})
        suitesparse_report_not_found(
          "Did not find ${COMPONENT} header (required SuiteSparse component).")
      else()
        message(STATUS "Did not find ${COMPONENT} header (optional SuiteSparse component).")
        mark_as_advanced(SuiteSparse_${COMPONENT}_INCLUDE_DIR)
      endif()
    endif()
  endif()

  if(SuiteSparse_FIND_COMPONENT_${COMPONENT}_LIBRARIES)
    find_library(SuiteSparse_${COMPONENT}_LIBRARY
      NAMES ${SuiteSparse_FIND_COMPONENT_${COMPONENT}_LIBRARIES}
      PATH_SUFFIXES ${SuiteSparse_CHECK_PATH_SUFFIXES})
    if(SuiteSparse_${COMPONENT}_LIBRARY)
      message(STATUS "Found ${COMPONENT} library: ${SuiteSparse_${COMPONENT}_LIBRARY}")
      mark_as_advanced(SuiteSparse_${COMPONENT}_LIBRARY)
    else()
      set(SuiteSparse_${COMPONENT}_FOUND FALSE)
      if(SuiteSparse_FIND_REQUIRED_${COMPONENT})
        suitesparse_report_not_found(
          "Did not find ${COMPONENT} library (required SuiteSparse component).")
      else()
        message(STATUS "Did not find ${COMPONENT} library (optional SuiteSparse dependency)")
        mark_as_advanced(SuiteSparse_${COMPONENT}_LIBRARY)
      endif()
    endif()
  endif()

  if(SuiteSparse_FIND_REQUIRED_${COMPONENT})
    list(APPEND SuiteSparse_REQUIRED_VARS SuiteSparse_${COMPONENT}_INCLUDE_DIR)
    list(APPEND SuiteSparse_REQUIRED_VARS SuiteSparse_${COMPONENT}_LIBRARY)
  endif()

  # Define the target only if both the include directory and the library exist.
  if(SuiteSparse_${COMPONENT}_INCLUDE_DIR AND SuiteSparse_${COMPONENT}_LIBRARY)
    if(NOT TARGET SuiteSparse::${COMPONENT})
      add_library(SuiteSparse::${COMPONENT} IMPORTED UNKNOWN)
    endif()

    set_property(TARGET SuiteSparse::${COMPONENT} PROPERTY
      INTERFACE_INCLUDE_DIRECTORIES ${SuiteSparse_${COMPONENT}_INCLUDE_DIR})
    set_property(TARGET SuiteSparse::${COMPONENT} PROPERTY
      IMPORTED_LOCATION ${SuiteSparse_${COMPONENT}_LIBRARY})
  endif()
endmacro()

# Maintain a list of all variables that must be defined for SuiteSparse to be
# considered found. This helps the automatic failure message from
# FindPackageHandleStandardArgs().
unset(SuiteSparse_REQUIRED_VARS)

# BLAS is required by SuiteSparse.
find_package(BLAS QUIET)
if(NOT BLAS_FOUND)
  suitesparse_report_not_found(
    "Did not find BLAS library (required for SuiteSparse).")
endif()

# LAPACK is also required. Try OpenBLAS first because vcpkg's OpenBLAS package
# provides LAPACK symbols and may not expose a separate LAPACK target.
if(NOT LAPACK_FOUND)
  find_package(OpenBLAS QUIET)
  if(OpenBLAS_FOUND)
    set(LAPACK_FOUND TRUE)
    if(TARGET OpenBLAS::OpenBLAS)
      set(LAPACK_LIBRARIES OpenBLAS::OpenBLAS)
    else()
      set(LAPACK_LIBRARIES ${OpenBLAS_LIBRARIES})
    endif()
    message(STATUS "Using OpenBLAS for LAPACK.")
  else()
    find_package(LAPACK QUIET)
    if(NOT LAPACK_FOUND)
      suitesparse_report_not_found(
        "Did not find LAPACK library (required for SuiteSparse). "
        "Please install OpenBLAS or a LAPACK implementation.")
    endif()
  endif()
endif()

# Find each requested component.
foreach(component IN LISTS SuiteSparse_FIND_COMPONENTS)
  if(component STREQUAL Partition)
    # Partition is a meta component that neither provides additional headers nor
    # a separate library. It is strictly part of CHOLMOD.
    continue()
  endif()
  string(TOLOWER ${component} component_library)

  if(component STREQUAL "Config")
    set(component_header SuiteSparse_config.h)
    set(component_library suitesparseconfig)
  elseif(component STREQUAL "SPQR")
    set(component_header SuiteSparseQR.hpp)
  else()
    set(component_header ${component_library}.h)
  endif()

  suitesparse_find_component(${component}
    FILES ${component_header}
    LIBRARIES ${component_library})
endforeach()

# SuiteSparseQR may be compiled with Intel TBB. Assume that if TBB is installed,
# SuiteSparseQR was compiled with support for it.
if(TARGET SuiteSparse::SPQR)
  find_package(TBB QUIET NO_MODULE)
  if(TBB_FOUND)
    message(STATUS "Found Intel Thread Building Blocks (TBB) library "
      "(${TBB_VERSION_MAJOR}.${TBB_VERSION_MINOR} / ${TBB_INTERFACE_VERSION}). "
      "Assuming SuiteSparseQR was compiled with TBB.")
    set_property(TARGET SuiteSparse::SPQR APPEND PROPERTY
      INTERFACE_LINK_LIBRARIES TBB::tbb)
  else()
    message(STATUS "Did not find Intel TBB library, assuming SuiteSparseQR was "
      "not compiled with TBB.")
  endif()
endif()

check_library_exists(rt shm_open "" HAVE_LIBRT)

if(TARGET SuiteSparse::Config)
  # SuiteSparse_config requires librt for timing on Linux/Unix but not on macOS.
  if(HAVE_LIBRT)
    message(STATUS "Adding librt to SuiteSparse_config libraries "
      "(required on Linux & Unix [not OSX] if SuiteSparse is compiled with timing).")
    set_property(TARGET SuiteSparse::Config APPEND PROPERTY
      INTERFACE_LINK_LIBRARIES $<LINK_ONLY:rt>)
  else()
    message(STATUS "Could not find librt, but found SuiteSparse_config, "
      "assuming that SuiteSparse was compiled without timing.")
  endif()

  # Add BLAS and LAPACK as dependencies of SuiteSparse::Config for convenience
  # because all components depend on it.
  if(BLAS_FOUND)
    if(TARGET BLAS::BLAS)
      set_property(TARGET SuiteSparse::Config APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES $<LINK_ONLY:BLAS::BLAS>)
    else()
      set_property(TARGET SuiteSparse::Config APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES ${BLAS_LIBRARIES})
    endif()
  endif()

  if(LAPACK_FOUND)
    if(TARGET LAPACK::LAPACK)
      set_property(TARGET SuiteSparse::Config APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES $<LINK_ONLY:LAPACK::LAPACK>)
    else()
      set_property(TARGET SuiteSparse::Config APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES ${LAPACK_LIBRARIES})
    endif()
  endif()

  # Extract SuiteSparse version from SuiteSparse_config.h (>= v4).
  set(SuiteSparse_VERSION_FILE
    ${SuiteSparse_Config_INCLUDE_DIR}/SuiteSparse_config.h)
  if(NOT EXISTS ${SuiteSparse_VERSION_FILE})
    suitesparse_report_not_found(
      "Could not find file: ${SuiteSparse_VERSION_FILE} containing version "
      "information for >= v4 SuiteSparse installs, but SuiteSparse_config was "
      "found (only present in >= v4 installs).")
  else()
    file(READ ${SuiteSparse_VERSION_FILE} Config_CONTENTS)

    string(REGEX MATCH "#define SUITESPARSE_MAIN_VERSION[ \t]+([0-9]+)"
      SuiteSparse_VERSION_LINE "${Config_CONTENTS}")
    set(SuiteSparse_VERSION_MAJOR ${CMAKE_MATCH_1})

    string(REGEX MATCH "#define SUITESPARSE_SUB_VERSION[ \t]+([0-9]+)"
      SuiteSparse_VERSION_LINE "${Config_CONTENTS}")
    set(SuiteSparse_VERSION_MINOR ${CMAKE_MATCH_1})

    string(REGEX MATCH "#define SUITESPARSE_SUBSUB_VERSION[ \t]+([0-9]+)"
      SuiteSparse_VERSION_LINE "${Config_CONTENTS}")
    set(SuiteSparse_VERSION_PATCH ${CMAKE_MATCH_1})

    unset(SuiteSparse_VERSION_LINE)

    # Keep this on a single line so CMake does not interpret it as a list and
    # insert ';' separators.
    set(SuiteSparse_VERSION
      "${SuiteSparse_VERSION_MAJOR}.${SuiteSparse_VERSION_MINOR}.${SuiteSparse_VERSION_PATCH}")

    if(SuiteSparse_VERSION MATCHES "[0-9]+\\.[0-9]+\\.[0-9]+")
      set(SuiteSparse_VERSION_COMPONENTS 3)
    else()
      message(WARNING "Could not parse SuiteSparse_config.h: SuiteSparse "
        "version will not be available")

      unset(SuiteSparse_VERSION)
      unset(SuiteSparse_VERSION_MAJOR)
      unset(SuiteSparse_VERSION_MINOR)
      unset(SuiteSparse_VERSION_PATCH)
    endif()
  endif()
endif()

# CHOLMOD requires AMD, CAMD, CCOLAMD, and COLAMD.
if(TARGET SuiteSparse::CHOLMOD)
  foreach(component IN ITEMS AMD CAMD CCOLAMD COLAMD)
    if(TARGET SuiteSparse::${component})
      set_property(TARGET SuiteSparse::CHOLMOD APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES SuiteSparse::${component})
    else()
      # Consider CHOLMOD not found if a dependency cannot be found.
      set(SuiteSparse_CHOLMOD_FOUND FALSE)
    endif()
  endforeach()
endif()

# SPQR requires CHOLMOD.
if(TARGET SuiteSparse::SPQR)
  if(TARGET SuiteSparse::CHOLMOD)
    set_property(TARGET SuiteSparse::SPQR APPEND PROPERTY
      INTERFACE_LINK_LIBRARIES SuiteSparse::CHOLMOD)
  else()
    # Consider SPQR not found if CHOLMOD cannot be found.
    set(SuiteSparse_SPQR_FOUND FALSE)
  endif()
endif()

# Add SuiteSparse::Config as a dependency to all components.
if(TARGET SuiteSparse::Config)
  foreach(component IN LISTS SuiteSparse_FIND_COMPONENTS)
    if(component STREQUAL Config)
      continue()
    endif()

    if(TARGET SuiteSparse::${component})
      set_property(TARGET SuiteSparse::${component} APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES SuiteSparse::Config)
    endif()
  endforeach()
endif()

# Check whether CHOLMOD was compiled with METIS support.
if(TARGET SuiteSparse::CHOLMOD)
  # If SuiteSparse was compiled as a static library we need to link against
  # METIS already during the check. Otherwise the check can fail due to
  # undefined references even though SuiteSparse was compiled with METIS.
  find_package(METIS)

  if(TARGET METIS::METIS)
    cmake_push_check_state(RESET)
    set(CMAKE_REQUIRED_LIBRARIES SuiteSparse::CHOLMOD METIS::METIS)
    check_symbol_exists(cholmod_metis cholmod.h SuiteSparse_CHOLMOD_USES_METIS)
    cmake_pop_check_state()

    if(SuiteSparse_CHOLMOD_USES_METIS)
      set_property(TARGET SuiteSparse::CHOLMOD APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES $<LINK_ONLY:METIS::METIS>)

      # Provide the SuiteSparse::Partition component whose availability indicates
      # that CHOLMOD was compiled with the Partition module.
      if(NOT TARGET SuiteSparse::Partition)
        add_library(SuiteSparse::Partition IMPORTED INTERFACE)
      endif()

      set_property(TARGET SuiteSparse::Partition APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES SuiteSparse::CHOLMOD)
    endif()
  endif()
endif()

if(TARGET SuiteSparse::Partition)
  set(SuiteSparse_Partition_FOUND TRUE)
else()
  set(SuiteSparse_Partition_FOUND FALSE)
endif()

suitesparse_reset_find_library_prefix()

# Handle REQUIRED and QUIET arguments to find_package().
include(FindPackageHandleStandardArgs)
if(SuiteSparse_FOUND)
  find_package_handle_standard_args(SuiteSparse
    REQUIRED_VARS ${SuiteSparse_REQUIRED_VARS}
    VERSION_VAR SuiteSparse_VERSION
    FAIL_MESSAGE "Failed to find some/all required components of SuiteSparse."
    HANDLE_COMPONENTS)
else()
  # Do not pass VERSION_VAR if we failed to find SuiteSparse to avoid a
  # confusing autogenerated failure message stating 'not found (missing: FOO)
  # (found version: x.y.z)'.
  find_package_handle_standard_args(SuiteSparse
    REQUIRED_VARS ${SuiteSparse_REQUIRED_VARS}
    FAIL_MESSAGE "Failed to find some/all required components of SuiteSparse."
    HANDLE_COMPONENTS)
endif()

# Pop CMP0057.
cmake_policy(POP)
