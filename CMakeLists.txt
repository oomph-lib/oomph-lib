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
cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(
  oomphlib
  DESCRIPTION "A finite-element library for simulating multi-physics problems."
  VERSION 1.6.0
  HOMEPAGE_URL "http://oomph-lib.maths.man.ac.uk/doc/html/index.html"
  LANGUAGES C CXX Fortran)

# Specify the location of non-standard CMake modules and the C++ standard
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# ---------------------------[ INITIAL PROJECT SETUP ]--------------------------

if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "" FORCE)
endif()

# Sort out config files
include(OomphConfigureProjectEnvironment)

# Select level of verbosity (e.g. STATUS, VERBOSE, etc.). Only valid from CMake
# 3.17 onwards. For earlier CMake versions, setting this variable has no effect
set(CMAKE_MESSAGE_LOG_LEVEL VERBOSE)

# ----------------------[ FIND REQUIRED/DESIRED PACKAGES ]----------------------

# Built-in modules to find out whether certain packages are installed
find_package(Doxygen)
find_package(Git)
find_package(Python3 REQUIRED)

# ------------------------------[ BUILD OPTIONS ]-------------------------------
# cmake-format: off

# We silence a few harmless warnings in external_src/ and in src/. If you want,
# you can choose to have them displayed during the build
option(OOMPH_DONT_SILENCE_USELESS_WARNINGS "Display all warnings" OFF)

# Attempt to build and use up-to-date external sources. [WORK IN PROGRESS]
# CHANGES WITH 'OOMPH_TRANSITION_TO_VERSION_3'
option(OOMPH_TRANSITION_TO_VERSION_3 "Try to build with up-to-date external sources" OFF)
option(OOMPH_ENABLE_64BIT_INTS "Use 64-bit ints for matrix indices" OFF)

# Enable the Clang static code analyser/formatter. [WORK IN PROGRESS]
# WARNING: If you enable clang-tidy, you MUST NOT(!) run the Makefile/Ninja
# build in parallel as it will likely lead to data races and thus the output
# will be mangled. ALSO, you should enable as many options as possible for the
# build so that most, or all, of the code is visible to the formatter (which
# formats the code visible during the build).
#
# This functionality should really be reserved for a CI server to run...
option(OOMPH_ENABLE_CLANG_TIDY "Enable clang-tidy C++ static code analysis" OFF)
option(OOMPH_ENABLE_CLANG_FORMAT "Enable Clang-based C++ code formatting" OFF)
option(OOMPH_ENABLE_INCLUDE_WHAT_YOU_USE "Enable 'Include What You Use' sanitisation" OFF)

# Extra debugging compile options; only used by Debug configurations
option(OOMPH_ENABLE_PARANOID "Enable the PARANOID flag in Debug" ON)
option(OOMPH_ENABLE_RANGE_CHECKING "Enable RANGE_CHECKING flag in Debug" ON)

# General options
option(BUILD_SHARED_LIBS "Build using shared libraries; static otherwise" OFF)
option(SUPPRESS_DOC "Suppress Doxygen creation of API documentation" ON)
option(SYMBOLIC_LINKS_FOR_HEADERS "Replace headers by symbolic links" ON)

# Testing [PM: not sure if this is needed; should build self-tests separately]
option(OOMPH_ENABLE_TESTING "Build unit tests" OFF)

# Multiprocessing/GPU usage. Note: CUDA functionality is current very(!) limited
# and only accessible through certain third-party libraries
option(OOMPH_ENABLE_MPI "Enable the use of MPI for parallel processing" OFF)
option(OOMPH_ENABLE_CUDA "Enable the use of CUDA where possible" OFF)

# Optionally build OpenBLAS/LAPACK
option(OOMPH_USE_OPENBLAS "Use the OpenBLAS implementation of BLAS/LAPACK" ON)
# option(OOMPH_USE_MKL_BLAS_AND_LAPACK "Use the MKL implementation of BLAS and LAPACK" OFF)

# TODO: Discuss CGAL with MH: "CGAL 5.0 is the first release of CGAL that
# requires a C++ compiler with the support of C++14"

# FIXME: spdlog and Backward need to be "exported"(?!) with 'generic'. Not sure
# what this means....

# Optionally build third-party libraries. [WORK IN PROGRESS.]
option(OOMPH_WANT_BACKWARD "Use the Backward library for stack trace pretty-printing" OFF)
option(OOMPH_WANT_CGAL "Do we want to build the CGAL library?" OFF)
option(OOMPH_WANT_HYPRE "Use Hypre library" OFF)
option(OOMPH_WANT_MUMPS "Use MUMPS library" OFF)
option(OOMPH_WANT_NLOHMANN_JSON "Use the nlohmann JSON library" OFF)
option(OOMPH_WANT_TRILINOS "Use Trilinos library" OFF)
option(OOMPH_WANT_SPDLOG "Use the spdlog logging library" OFF)

# Fine-grained control of build of external_srcs. NOTE: These should be set on
# the fly after checking whether the appropriate variable WITH_* has been set!
# TODO: ADD OOMPH_ PREFIXES TO THESE VARIABLES!
# option(ALREADY_HAVE_BLAS "Use BLAS library specified via absolute path" OFF)
# option(ALREADY_HAVE_LAPACK "Use LAPACK library specified via absolute path" OFF)
# option(ALREADY_HAVE_METIS "Suppress build of oomph-lib's own copy of metis" OFF)
# option(ALREADY_HAVE_PARMETIS "Suppress build of oomph-lib's copy of parmetis" OFF)
# option(ALREADY_HAVE_SUPERLU "Suppress build of oomph-lib's copy of serial SuperLU" OFF)
# option(ALREADY_HAVE_SUPERLU_DIST "Suppress build of oomph-lib's copy of serial SuperLU_dist" OFF)
option(SUPPRESS_TRIANGLE_LIB "Suppress build of oomph-lib's copy of the triangle library" OFF)
option(SUPPRESS_TETGEN_LIB "Suppress build of oomph-lib's copy of the tetgen library" OFF)

# Support for deprecated packages:
option(OOMPH_USE_DEPRECATED_SUPERLU "Use oomph-lib's deprecated version of SuperLU" OFF)

# cmake-format: on
# -------------------------[ PROCESS SELECTED OPTIONS ]-------------------------

# Run clang-tidy during the build if the user requires
if(OOMPH_ENABLE_CLANG_TIDY)
  include(OomphStaticAnalysers)
endif()

# Run clang-format during the build if the user requires
if(OOMPH_ENABLE_CLANG_FORMAT)
  include(OomphClangFormat)
endif()

# Set the matrix integer index size
if(OOMPH_ENABLE_64BIT_INTS)
  set(OOMPH_INTSIZE 64)
else()
  set(OOMPH_INTSIZE 32)
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

# Let CMake take care of the famous -fPIC flag
if(BUILD_SHARED_LIBS)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# The demo_driver self-tests use shell scripts and symbolic links, both of which
# require a Unix-based system
if(NOT UNIX)
  set(OOMPH_ENABLE_TESTING OFF)
  set(SYMBOLIC_LINKS_FOR_HEADERS OFF)
endif()

# Set up MPI functionality (see cmake/OomphMPI.cmake)
if(OOMPH_ENABLE_MPI)
  include(OomphMPI)
endif()

# Set up CUDA functionality. Nothing in the library uses CUDA yet. However, the
# plan is to support SuperLU Dist 7.0.0 in the future, which is capable of
# harnessing CUDA functionality for (hopefully damn fast) computations.
if(OOMPH_ENABLE_CUDA)
  find_package(CUDA REQUIRED)
endif()

if(OOMPH_USE_OPENBLAS)
  # Search for OpenBLAS/LAPACK libraries. Ideally, installed via 'sudo apt-get'
  # or 'brew install' and the user simply sets 'CMAKE_PREFIX_PATH' in their
  # CMakeUserPresets.json.
  list(APPEND CMAKE_PREFIX_PATH "/usr/local/Cellar/openblas/0.3.17")
  set(BLA_VENDOR OpenBLAS)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)

  # NOTE: For BLAS and LAPACK, there should be just a single library artifact
  # but we'll generically loop over the *_LIBRARIES variables set by the Find*
  # modules to be on the safe side...

  # Store a list of directories containing the BLAS libraries
  set(BLAS_DIRS)
  foreach(BLAS_LIBRARY IN LISTS BLAS_LIBRARIES)
    get_filename_component(BLAS_DIR ${BLAS_LIBRARY} DIRECTORY)
    list(APPEND BLAS_DIRS ${BLAS_DIR})
  endforeach()

  # Store a list of directories containing the LAPACK libraries
  set(LAPACK_DIRS)
  foreach(LAPACK_LIBRARY IN LISTS LAPACK_LIBRARIES)
    get_filename_component(LAPACK_DIR ${LAPACK_LIBRARY} DIRECTORY)
    list(APPEND LAPACK_DIRS ${LAPACK_DIR})
  endforeach()

  # PM: Disable Trilinos for now; not sure how to link it against OpenBLAS
  set(OOMPH_WANT_TRILINOS OFF CACHE BOOL
      "Disable Trilinos when OpenBLAS is used" FORCE)
endif()

# TODO: FINISH OFF INTEGRATION OF THIS VARIABLES INTO BUILD
if(DEFINED WITH_BLACS)
  set(OOMPH_ALREADY_HAVE_BLACS TRUE)
endif()
if(DEFINED WITH_SCALAPACK)
  set(OOMPH_ALREADY_HAVE_SCALAPACK TRUE)
endif()
if(DEFINED WITH_PORD)
  set(OOMPH_ALREADY_HAVE_PORD TRUE)
endif()
if(DEFINED WITH_MUMPS)
  set(OOMPH_ALREADY_HAVE_MUMPS TRUE)
endif()
if(DEFINED WITH_BLAS)
  set(OOMPH_ALREADY_HAVE_BLAS TRUE)
endif()
if(DEFINED WITH_LAPACK)
  set(OOMPH_ALREADY_HAVE_LAPACK TRUE)
endif()
if(DEFINED WITH_HYPRE)
  set(OOMPH_ALREADY_HAVE_HYPRE TRUE)
endif()
if(DEFINED WITH_TRILINOS)
  set(OOMPH_ALREADY_HAVE_TRILINOS TRUE)
endif()
if(DEFINED WITH_BOOST)
  set(OOMPH_ALREADY_HAVE_BOOST TRUE)
endif()
if(DEFINED WITH_GMP)
  set(OOMPH_ALREADY_HAVE_GMP TRUE)
endif()
if(DEFINED WITH_MPFR)
  set(OOMPH_ALREADY_HAVE_MPFR TRUE)
endif()
if(DEFINED WITH_CGAL)
  set(OOMPH_ALREADY_HAVE_CGAL TRUE)
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

message(VERBOSE "<====== Starting oomph-lib library build ======>")

# Process the oomph-lib directories containing the library files
add_subdirectory(scripts)
add_subdirectory(external_src)
add_compile_definitions(${EXTERNAL_SRC_CXX_DEFINITIONS})
add_subdirectory(external_distributions)
add_compile_definitions(${EXTERNAL_DIST_CXX_DEFINITIONS})
add_subdirectory(src)

message(VERBOSE "<====== Finished oomph-lib library build ======>")

# --------------------------------[ EXPORTING ]---------------------------------

# Now export the package so it can be found using find_package(...)
include(OomphExportPackage)

# ---------------------------------[ TESTING ]----------------------------------

# The last part to the build is compiling tests. The demo_drivers currently work
# when run *outside* the library build. I haven't tried integrating their build
# into the build of the library itself yet. I imagine it should be okay for the
# most part. There might be an issue with calling find_package(oomphlib) from
# the demo_drivers CMakeLists.txt files though. Also, note that private/ and
# user_src/ and user_drivers/ are not built with the library.

# Enable testing for the current directory and all child directories but only
# when this is the root project, i.e. it isn't being built as part of some
# larger project hierarchy.
if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME) AND OOMPH_ENABLE_TESTING)
  add_subdirectory(demo_drivers)
  add_subdirectory(self_test)
  include(OomphGetDoctest)
  enable_testing()
  add_subdirectory(tests)
endif()

# -----------------------[ CONFIGURE oomph-lib-config.h ]-----------------------

# Construct a oomph-lib-specific config file. Currently
if(OOMPH_ADD_CONFIG_H)
  configure_file(${CMAKE_SOURCE_DIR}/cmake/oomph-lib-config.h.in
                 ${CMAKE_SOURCE_DIR}/oomph-lib-config.h @ONLY)
  configure_file(${CMAKE_SOURCE_DIR}/cmake/oomph-lib-config.h.in
                 ${CMAKE_BINARY_DIR}/src/oomph-lib-config.h @ONLY)
  install(FILES ${CMAKE_BINARY_DIR}/src/oomph-lib-config.h
          DESTINATION "${OOMPH_INSTALL_INCLUDE_DIR}")
endif()

# -------------------[ DEBUG SHIT: DELETE AFTER OVERHAUL ]----------------------

include(OomphPrintCompilerDefinitions)

# --------------------------------[ PACKAGING ]---------------------------------

# Work in progress...

# ------------------------------------------------------------------------------
