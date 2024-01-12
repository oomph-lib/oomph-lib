# cmake-format: off
# ==============================================================================
#                                          888             888 d8b 888
#                                          888             888 Y8P 888
#                                          888             888     888
#  .d88b.   .d88b.  88888b.d88b.  88888b.  88888b.         888 888 88888b.
# d88""88b d88""88b 888 "888 "88b 888 "88b 888 "88b        888 888 888 "88b
# 888  888 888  888 888  888  888 888  888 888  888 888888 888 888 888  888
# Y88..88P Y88..88P 888  888  888 888 d88P 888  888        888 888 888 d88P
#  "Y88P"   "Y88P"  888  888  888 88888P"  888  888        888 888 88888P"
#                                 888
#                                 888
#                                 888
# ------------------------------------------------------------------------------
# oomph-lib is an object-oriented, open-source finite-element library for the
# simulation of multi-physics problems. It is developed and maintained by
# Matthias Heil and Andrew Hazel of the School of Mathematics at The University
# of Manchester, along with many other contributors.
#
# This library is distributed with a GNU LGPL License; see LICENSE.txt.
# ------------------------------------------------------------------------------
# The version number follows the symantic versioning strategy here:
#                           https://semver.org/
# ------------------------------------------------------------------------------
# See https://cmake-basis.github.io/standard/fhs.html for a description of
# recommended source, build, and installation layouts.
#
# Additional notes for post-renovation:
# See: https://nvie.com/posts/a-successful-git-branching-model/ for a Git
# branching model.
# ------------------------------------------------------------------------------
# The "cmake-format" commands littered around the project are used to
# disable/enable the use of the cmake-format extension used to automatically
# format CMake files. Do not delete them.
# ==============================================================================
# cmake-format: on
cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(
  oomphlib
  DESCRIPTION "A finite-element library for simulating multi-physics problems."
  VERSION 2.9.0
  HOMEPAGE_URL "https://oomph-lib.github.io/oomph-lib/doc/html/"
  LANGUAGES C CXX Fortran)

# ---------------------------[ INITIAL PROJECT SETUP ]--------------------------

# Specify the location of non-standard CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/")

# Specify the C++ standard
set(CMAKE_C_STANDARD 17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# Build in Release (fully optimised) mode by default
if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# Sort out config files
include(OomphConfigureCMakeInfoForProject)

# Select level of verbosity (e.g. STATUS, VERBOSE, etc.). Only valid from CMake
# 3.17 onwards. For earlier CMake versions, setting this variable has no effect
set(CMAKE_MESSAGE_LOG_LEVEL VERBOSE)

# ----------------------[ FIND REQUIRED/DESIRED PACKAGES ]----------------------

# Built-in modules to find out whether certain packages are installed
find_package(Git)
find_package(Python3 REQUIRED)

# ------------------------------[ BUILD OPTIONS ]-------------------------------
# cmake-format: off
include(OomphCustomOptions)

# Regular CMake option; doesn't influence the build so we won't print its value at the end
option(OOMPH_PRINT_SETTINGS_AFTER_INSTALL "Display the oomph-lib settings at the end of the install step" ON)

# FIXME: Delete
oomph_option(OOMPH_USE_OLD_SUPERLU_DIST "Revert to SuperLU Dist 3.0" OFF)

# General options
oomph_option(OOMPH_DONT_SILENCE_USELESS_WARNINGS "Display (harmless) warnings from external_src/ and src/ that are silenced" OFF)
oomph_option(OOMPH_ENABLE_MPI "Enable the use of MPI for parallel processing" OFF)
oomph_option(OOMPH_ENABLE_PARANOID "Enable the PARANOID flag in Debug" OFF)
oomph_option(OOMPH_ENABLE_RANGE_CHECKING "Enable RANGE_CHECKING flag in Debug" OFF)

# Fine-grained control of build of external_srcs.
# FIXME: Make sure the SUPPRESS variables below actually work
# TODO: These should be set on the fly after checking whether the appropriate variable WITH_* has been set...
oomph_option(OOMPH_SUPPRESS_TRIANGLE_LIB "Suppress build of oomph-lib's copy of the triangle library" OFF)
oomph_option(OOMPH_SUPPRESS_TETGEN_LIB "Suppress build of oomph-lib's copy of the tetgen library" OFF)

# Paths to packages that the user has installed themselves
oomph_path_option(
  FLAG OOMPH_USE_OPENBLAS_FROM
  DOCSTRING "Path to OpenBLAS installation."
  REQUIRED
)
oomph_path_option(
  FLAG OOMPH_USE_GMSH_FROM
  DOCSTRING "Path to 'gmsh' executable"
)
oomph_path_option(
  FLAG OOMPH_USE_MPI_FROM
  DOCSTRING "Path to MPI installation"
)
oomph_path_option(
  FLAG OOMPH_USE_GMP_FROM
  DOCSTRING "Path to GMP installation"
)
oomph_path_option(
  FLAG OOMPH_USE_MPFR_FROM
  DOCSTRING "Path to MPFR installation"
)
oomph_path_option(
  FLAG OOMPH_USE_BOOST_FROM
  DOCSTRING "Path to Boost installation"
)
oomph_path_option(
  FLAG OOMPH_USE_CGAL_FROM
  DOCSTRING "Path to CGAL installation"
)
oomph_path_option(
  FLAG OOMPH_USE_MUMPS_FROM
  DOCSTRING "Path to MUMPS installation"
)
oomph_path_option(
  FLAG OOMPH_USE_HYPRE_FROM
  DOCSTRING "Path to Hypre installation"
)
oomph_path_option(
  FLAG OOMPH_USE_TRILINOS_FROM
  DOCSTRING "Path to Trilinos installation"
)

# cmake-format: on
# -------------------------[ PROCESS SELECTED OPTIONS ]-------------------------
# Introduce oomph_add_cxx_compile_definitions(...) to add compile definitions to
# C++ files only
include(OomphAddCXXCompileDefinitions)

# FIXME: XCode 15.0 has a new linker implementation which breaks things; setting
# this values tells CMake to revert to the old linker
if(APPLE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 15.0)
      add_link_options("-Wl,-ld_classic")
    endif()
  endif()
endif()

# Let CMake take care of the famous -fPIC flag
if(BUILD_SHARED_LIBS)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# The user needs to provide at least v0.3.25 of OpenBLAS
find_package(
  OpenBLAS 0.3.25 REQUIRED
  PATHS ${OOMPH_USE_OPENBLAS_FROM}
  NO_DEFAULT_PATH)

# FIXME: This oddity (of looking for BLAS even though we've already found
# OpenBLAS) is to handle an issue in HYPRE. Although we built HYPRE with
# OpenBLAS, HYPRE searches (using CMake) for the BLAS library (not OpenBLAS). In
# some cases it won't be able to find it, resulting in an unexpected error
# during the configuration of the main oomph-lib library. To handle this, we'll
# tell CMake to locate the OpenBLAS library when looking for BLAS and find it
# ourselves. When HYPRE later calls find_dependency(BLAS), CMake will tell it
# that it has already found BLAS (the OpenBLAS version) and will make everything
# (hopefully) work!
set(BLA_VENDOR OpenBLAS)
set(BACKUP_CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}")
set(CMAKE_PREFIX_PATH ${OOMPH_USE_OPENBLAS_FROM})
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
set(CMAKE_PREFIX_PATH ${BACKUP_CMAKE_PREFIX_PATH})

# Do we want to use Gmsh?
if(OOMPH_USE_GMSH_FROM)
  include(OomphFindGmsh)
endif()

# We silence several useless warnings but if you want to see them, enable the
# OOMPH_DONT_SILENCE_USELESS_WARNINGS flag
if(NOT OOMPH_DONT_SILENCE_USELESS_WARNINGS)
  include(OomphSilenceWarnings)
endif()

# Set up MPI functionality (see cmake/OomphMPI.cmake)
if(OOMPH_ENABLE_MPI)
  include(OomphMPI)
endif()

# Add PARANOID and RANGE_CHECKING for Debug configurations (ONLY!) if desired
if((OOMPH_ENABLE_PARANOID OR OOMPH_ENABLE_RANGE_CHECKING)
   AND (NOT CMAKE_BUILD_TYPE STREQUAL "Debug"))
  message(
    WARNING
      "You can only enable\n\tOOMPH_ENABLE_PARANOID\nand\n\tOOMPH_ENABLE_RANGE_CHECKING\nwhen building with\n\tCMAKE_BUILD_TYPE=\"Debug\"\nbut you're building with\n\tCMAKE_BUILD_TYPE=\"${CMAKE_BUILD_TYPE}\"\nI'm therefore going to ignore this flag..."
  )
endif()

# Only set the PARANOID or RANGE_CHECKING compiler definitions if building the
# Debug configuration. It doesn't make sense to build with full optimisation and
# to enable these flags
if(OOMPH_ENABLE_PARANOID)
  oomph_add_cxx_compile_definitions($<$<CONFIG:DEBUG,>:PARANOID>)
endif()
if(OOMPH_ENABLE_RANGE_CHECKING)
  oomph_add_cxx_compile_definitions($<$<CONFIG:DEBUG,>:RANGE_CHECKING>)
endif()

# --------------------------[ FIND REQUIRED HEADERS ]---------------------------

include(OomphCheckForRequiredHeaders)

# Enable @rpath in the install name for any shared library being built note:
# this is overwritten by INSTALL_NAME_DIR. See: https://tinyurl.com/feksedtz
set(CMAKE_MACOSX_RPATH TRUE)

# -------------------------[ ADD CONFIG. HEADER FLAG ]--------------------------

set(OOMPH_ADD_CONFIG_H TRUE)
if(OOMPH_ADD_CONFIG_H)
  oomph_add_cxx_compile_definitions(HAVE_CONFIG_H)
endif()

# ------------------------------[ BUILD LIBRARY ]-------------------------------

message(VERBOSE "<====== Starting oomph-lib library configuration ======>")

# Configure the external sources (i.e. third-party libraries that oomph-lib
# always ships with)
add_subdirectory(external_src)

# FIXME: Having to duplicate CXX definitions for C code to make SuperLUDist code
# work...

# Add C++ definitions for all of the defined third-party sources that we build
# with oomph-lib, e.g. -DOOMPH_HAS_SUPERLU, -DOOMPH_HAS_TRIANGLE, etc.
oomph_add_c_compile_definitions(${EXTERNAL_SRC_CXX_DEFINITIONS})
oomph_add_cxx_compile_definitions(${EXTERNAL_SRC_CXX_DEFINITIONS})

# Locate the third-party libraries that the user has provided via the
# OOMPH_USE_<LIBRARY>_FROM flags
include(OomphLocateOptionalThirdPartyLibs)

# Add C++ definitions for all of the defined (optional) third-party libraries,
# e.g. -DOOMPH_HAS_CGAL, -DOOMPH_HAS_MUMPS, etc.
oomph_add_c_compile_definitions(${EXTERNAL_DIST_CXX_DEFINITIONS})
oomph_add_cxx_compile_definitions(${EXTERNAL_DIST_CXX_DEFINITIONS})

# Configure all of the main oomph-lib library sources
add_subdirectory(src)

message(VERBOSE "<====== Finished oomph-lib library configuration ======>")

# --------------------------------[ EXPORTING ]---------------------------------

# Now export the package so it can be found using find_package(...)
include(OomphInstallLibrary)

# ---------[ CONFIGURE oomph-lib-config.h AND oomphlibUninstall.cmake ]---------

# Construct a oomph-lib-specific config file. Currently
if(OOMPH_ADD_CONFIG_H)
  configure_file(${CMAKE_SOURCE_DIR}/cmake/oomph-lib-config.h.in
                 ${CMAKE_SOURCE_DIR}/oomph-lib-config.h @ONLY)
  configure_file(${CMAKE_SOURCE_DIR}/cmake/oomph-lib-config.h.in
                 ${CMAKE_BINARY_DIR}/src/oomph-lib-config.h @ONLY)
  install(FILES ${CMAKE_BINARY_DIR}/src/oomph-lib-config.h
          DESTINATION "${OOMPH_INSTALL_INCLUDE_DIR}")
endif()

# Configure the oomphlibUninstall.cmake file to populate the file with the list
# of auto-generated combined headers
configure_file(
  "${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Uninstall.cmake.in"
  "${CMAKE_BINARY_DIR}/${PROJECT_NAME}Uninstall.cmake" IMMEDIATE @ONLY)

# Define the 'uninstall' command to handle uninstalling the installed targets,
# to be used, e.g. with the commands 'make uninstall' or 'ninja install'
# (depending on the user-selected generator)
add_custom_target(
  oomph_uninstall COMMAND ${CMAKE_COMMAND} -P
                          ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Uninstall.cmake)

# --------------------[ PRINT PROJECT-SPECIFIC SETTINGS ]-----------------------

include(OomphDocProjectSettings)
if(OOMPH_PRINT_SETTINGS_AFTER_INSTALL)
  oomph_doc_project_settings(ENABLE_ALSO_PRINT_SETTINGS_AFTER_INSTALL)
else()
  oomph_doc_project_settings()
endif()

# ------------------------------------------------------------------------------
