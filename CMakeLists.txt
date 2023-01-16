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
cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
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
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# Build in Release (fully optimised) mode by default
if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# Sort out config files
include(OomphConfigureProject)

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

# General options
oomph_option(OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY "Build tests with library build" OFF)
oomph_option(OOMPH_DONT_SILENCE_USELESS_WARNINGS "Display (harmless) warnings from external_src/ and src/ that are silenced" OFF)
oomph_option(OOMPH_ENABLE_MPI "Enable the use of MPI for parallel processing" OFF)
oomph_option(OOMPH_ENABLE_PARANOID "Enable the PARANOID flag in Debug" OFF)
oomph_option(OOMPH_ENABLE_RANGE_CHECKING "Enable RANGE_CHECKING flag in Debug" OFF)
oomph_option(OOMPH_TRANSITION_TO_VERSION_3 "Try to build with up-to-date external sources" OFF)
oomph_option(OOMPH_USE_DEPRECATED_SUPERLU "Use oomph-lib's deprecated version of SuperLU (4.3)" OFF)

# Fine-grained control of build of external_srcs.
# FIXME: Make sure the SUPPRESS variables below actually work
# TODO: These should be set on the fly after checking whether the appropriate variable WITH_* has been set...
oomph_option(OOMPH_SUPPRESS_TRIANGLE_LIB "Suppress build of oomph-lib's copy of the triangle library" OFF)
oomph_option(OOMPH_SUPPRESS_TETGEN_LIB "Suppress build of oomph-lib's copy of the tetgen library" OFF)

# FIXME: Override OOMPH_WANT_<LIB> variables if OOMPH_USE_<LIB>_FROM is specified

# Third-party libraries. [WORK IN PROGRESS.]
oomph_option(OOMPH_WANT_CGAL "Do we want to build the CGAL library?" OFF)
oomph_option(OOMPH_WANT_HYPRE "Use Hypre library" OFF)
oomph_option(OOMPH_WANT_MUMPS "Use MUMPS library" OFF)
oomph_option(OOMPH_WANT_TRILINOS "Use Trilinos library" OFF)
oomph_option(OOMPH_WANT_NLOHMANN_JSON "Use the nlohmann JSON library" OFF)

# Enable the use of code coverage flags so that coverage data can be collected
oomph_option(OOMPH_ENABLE_CODE_COVERAGE "Enable collection of code coverage results" OFF)

# Paths to packages that the user has installed themselves
oomph_path_option(OOMPH_USE_MPI_FROM "Path to MPI installation" "${OOMPH_USE_MPI_FROM}")
oomph_path_option(OOMPH_USE_GMSH_FROM "Path to 'gmsh' executable" "${OOMPH_USE_GMSH_FROM}")
oomph_path_option(OOMPH_USE_BLAS_FROM "Path to BLAS installation" "${OOMPH_USE_BLAS_FROM}")
oomph_path_option(OOMPH_USE_LAPACK_FROM "Path to LAPACK installation" "${OOMPH_USE_LAPACK_FROM}")
oomph_path_option(OOMPH_USE_GMP_FROM "Path to GMP installation" "${OOMPH_USE_GMP_FROM}")
oomph_path_option(OOMPH_USE_MPFR_FROM "Path to MPFR installation" "${OOMPH_USE_MPFR_FROM}")
oomph_path_option(OOMPH_USE_BOOST_FROM "Path to Boost installation" "${OOMPH_USE_BOOST_FROM}")
oomph_path_option(OOMPH_USE_CGAL_FROM "Path to CGAL installation" "${OOMPH_USE_CGAL_FROM}")
oomph_path_option(OOMPH_USE_MUMPS_FROM "Path to MUMPS installation" "${OOMPH_USE_MUMPS_FROM}")
oomph_path_option(OOMPH_USE_HYPRE_FROM "Path to Hypre installation" "${OOMPH_USE_HYPRE_FROM}")
oomph_path_option(OOMPH_USE_TRILINOS_FROM "Path to Trilinos installation" "${OOMPH_USE_TRILINOS_FROM}")

# cmake-format: on
# -------------------------[ PROCESS SELECTED OPTIONS ]-------------------------

# If OOMPH_USE_GMSH_FROM is the path to an executable, make sure it is a
# program. If it's a directory, check that we can find 'gmsh' in the directory
if(OOMPH_USE_GMSH_FROM)
  if(NOT EXISTS "${OOMPH_USE_GMSH_FROM}")
    message(
      FATAL_ERROR
        "You said that 'gmsh' can be found here:\n\t${OOMPH_USE_GMSH_FROM}\nbut this file/directory does not exist!"
    )
  endif()

  # Extract path components
  set(OOMPH_GMSH_DIR ${OOMPH_USE_GMSH_FROM})
  set(OOMPH_GMSH_PROGRAM)
  if(NOT IS_DIRECTORY ${OOMPH_USE_GMSH_FROM})
    cmake_path(GET OOMPH_USE_GMSH_FROM PARENT_PATH OOMPH_GMSH_DIR)
    cmake_path(GET OOMPH_USE_GMSH_FROM STEM OOMPH_GMSH_PROGRAM)
  endif()

  # Search for a program;
  #
  # TODO: Add a 'VALIDATOR' argument to check that gmsh works
  find_program(
    OOMPH_GMSH_COMMAND
    NAMES gmsh ${OOMPH_GMSH_PROGRAM}
    PATHS ${OOMPH_GMSH_DIR} REQUIRED
    NO_DEFAULT_PATH)

  # Notify user
  message(STATUS "Found 'gmsh' at: ${OOMPH_USE_GMSH_FROM}")
endif()

# The demo_driver self-tests use shell scripts and symbolic links, both of which
# require a Unix-based system
if(NOT UNIX)
  set(OOMPH_ENABLE_CODE_COVERAGE OFF CACHE BOOL
      "Enable collection of code coverage results" FORCE)
  set(OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY OFF)
endif()

# Let CMake take care of the famous -fPIC flag
if(BUILD_SHARED_LIBS)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

if(NOT OOMPH_DONT_SILENCE_USELESS_WARNINGS)
  include(OomphSilenceWarnings)
endif()

# Set up MPI functionality (see cmake/OomphMPI.cmake)
if(OOMPH_ENABLE_MPI)
  include(OomphMPI)

  # Sanity check
  if(NOT OOMPH_MPI_NUM_PROC MATCHES "^[0-9]+$")
    message(
      FATAL_ERROR
        "The flag 'OOMPH_MPI_NUM_PROC' must be set an integer, not ${OOMPH_MPI_NUM_PROC}!"
    )
  endif()
endif()

# Add PARANOID and RANGE_CHECKING for Debug configurations (ONLY!) if desired.
# We have to use generator expressions so that it works correctly for multi-
# configuration generators like Xcode, Visual Studio or Ninja Multi-Config
if(OOMPH_ENABLE_PARANOID)
  add_compile_definitions("$<$<CONFIG:DEBUG>:PARANOID>")
endif()
if(OOMPH_ENABLE_RANGE_CHECKING)
  add_compile_definitions("$<$<CONFIG:DEBUG>:RANGE_CHECKING>")
endif()

# If we're using up-to-date external sources/dists then make sure we have a
# preprocessor definition to handle changes that need to be made from within
# source files (e.g. API changes)
if(OOMPH_TRANSITION_TO_VERSION_3)
  add_compile_definitions(OOMPH_TRANSITION_TO_VERSION_3)
endif()
if(OOMPH_USE_DEPRECATED_SUPERLU)
  add_compile_definitions(OOMPH_USE_DEPRECATED_SUPERLU)
endif()

if(BLA_VENDOR
   OR OOMPH_USE_BLAS_FROM
   OR OOMPH_USE_LAPACK_FROM)
  # Temporarily set CMAKE_PREFIX_PATH to the values of OOMPH_USE_BLAS_FROM and
  # OOMPH_USE_LAPACK_FROM if either is provided
  if(OOMPH_USE_BLAS_FROM OR OOMPH_USE_LAPACK_FROM)
    set(BACKUP_CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}" CACHE INTERNAL "")
    set(CMAKE_PREFIX_PATH "${OOMPH_USE_BLAS_FROM};${OOMPH_USE_LAPACK_FROM}"
        CACHE INTERNAL "" FORCE)
  endif()

  # Search for the BLAS/LAPACK libraries
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)

  # Reset the value of CMAKE_PREFIX_PATH
  if(OOMPH_USE_BLAS_FROM OR OOMPH_USE_LAPACK_FROM)
    set(CMAKE_PREFIX_PATH "${BACKUP_CMAKE_PREFIX_PATH}" CACHE INTERNAL "" FORCE)
  endif()

  # Store a list of directories containing the BLAS libraries
  set(BLAS_DIRS)
  foreach(BLAS_LIBRARY IN LISTS BLAS_LIBRARIES)
    cmake_path(GET BLAS_LIBRARY PARENT_PATH BLAS_DIR)
    list(APPEND BLAS_DIRS ${BLAS_DIR})
  endforeach()

  # Store a list of directories containing the LAPACK libraries
  set(LAPACK_DIRS)
  foreach(LAPACK_LIBRARY IN LISTS LAPACK_LIBRARIES)
    cmake_path(GET LAPACK_LIBRARY PARENT_PATH LAPACK_DIR)
    list(APPEND LAPACK_DIRS ${LAPACK_DIR})
  endforeach()
endif()

if(OOMPH_ENABLE_CODE_COVERAGE)
  # We can only get coverage for the library sources if the demo drivers are
  # built with the library itself, which requires
  # OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY to be enabled
  if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
    message(
      WARNING
        "Code coverage was requested which requires building the demo drivers with the library\n"
        "but you have not enabled it. I'm going to enable it now.")
    set(OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY ON CACHE BOOL
        "Build tests with library build" FORCE)
  endif()

  include(OomphEnableCodeCoverage)
  oomph_enable_code_coverage()
endif()

# --------------------------[ FIND REQUIRED HEADERS ]---------------------------

include(OomphCheckForRequiredHeaders)

# Enable @rpath in the install name for any shared library being built note:
# this is overwritten by INSTALL_NAME_DIR. See: https://tinyurl.com/feksedtz
set(CMAKE_MACOSX_RPATH TRUE)

# -------------------------[ ADD CONFIG. HEADER FLAG ]--------------------------

set(OOMPH_ADD_CONFIG_H TRUE)
if(OOMPH_ADD_CONFIG_H)
  add_compile_definitions(HAVE_CONFIG_H)
endif()

# ------------------------------[ BUILD LIBRARY ]-------------------------------

message(VERBOSE "<====== Starting oomph-lib library configuration ======>")

# Process the oomph-lib directories containing the library files
add_subdirectory(external_src)
add_compile_definitions(${EXTERNAL_SRC_CXX_DEFINITIONS})
add_subdirectory(external_distributions)
add_compile_definitions(${EXTERNAL_DIST_CXX_DEFINITIONS})
add_subdirectory(src)

message(VERBOSE "<====== Finished oomph-lib library configuration ======>")

# --------------------------------[ EXPORTING ]---------------------------------

# Now export the package so it can be found using find_package(...)
include(OomphInstallLibrary)

# ---------------------------------[ TESTING ]----------------------------------

# Enable testing for the current directory and all child directories but ONLY
# when this is the root project, i.e. it isn't being built as part of some
# larger project hierarchy
if(OOMPH_IS_MAIN_PROJECT AND OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  include(CTest)
  include(OomphAddExecutable)
  include(OomphAddTest)
  include(OomphAddPureCppTest)
  include(OomphDefineTestData)
  include(OomphGetHashedTargetName)
  add_subdirectory(demo_drivers)
  # FIXME: Need to integrate self_test files into the self-test framework
  # add_subdirectory(self_test)
endif()

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

include(OomphPrintProjectSettings)
oomph_print_project_settings()

# ------------------------------------------------------------------------------
