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
  VERSION 1.5.0
  HOMEPAGE_URL "http://oomph-lib.maths.man.ac.uk/doc/html/index.html"
  LANGUAGES C CXX Fortran)

# Specify the location of non-standard CMake modules and the C++ standard
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# ---------------------------[ Initial project setup ]--------------------------

if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "" FORCE)
endif()

# Sort out config files
include(OomphConfigureProjectEnvironment)

# TODO: Change first and delete second once finished overhauling
set(CMAKE_MESSAGE_LOG_LEVEL VERBOSE)
include(OomphSilenceWarnings)

# ----------------------[ Find required/desired packages ]----------------------

# Built-in modules to find out whether certain packages are installed
find_package(BLAS)
find_package(Doxygen)
find_package(Git)
find_package(LAPACK)
find_package(Python3 REQUIRED)

# ------------------------------[ Build options ]-------------------------------
# Enable the clang-tidy static code analyser. NOTE: You should NOT run the
# Makefile/Ninja build in parallel as it may lead to data races. This
# functionality should probably be reserved for a CI server to run
option(OOMPH_ENABLE_CLANG_TIDY "Enable clang-tidy C++ static code analysis" OFF)
option(OOMPH_ENABLE_CLANG_FORMAT "Enable Clang-based C++ code formatting" OFF)

# Attempt to build and use up-to-date external sources. [WORK IN PROGRESS]
option(OOMPH_BLEEDING_EDGE "Try to build with up-to-date external sources" OFF)
option(OOMPH_ENABLE_64BIT_INTS "Use 64-bit ints for matrix indices" OFF)

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
option(OOMPH_ENABLE_MPI "Enable the use of MPI for parallel processing" ON)
option(OOMPH_ENABLE_CUDA "Enable the use of CUDA where possible" OFF)

# Optionally build third-party libraries
option(OOMPH_WANT_CGAL "Do we want to build the CGAL library?" OFF)
option(OOMPH_WANT_HYPRE "Use Hypre library" OFF)
option(OOMPH_WANT_MUMPS "Use MUMPS library" OFF)
option(OOMPH_WANT_TRILINOS "Use Trilinos library" OFF)
option(OOMPH_WANT_NLOHMANN_JSON "Use the nlohmann JSON library" OFF)

# Fine-grained control of build of external_srcs
option(ALREADY_HAVE_BLAS "Use BLAS library specified via absolute path" OFF)
option(ALREADY_HAVE_LAPACK "Use LAPACK library specified via absolute path" OFF)
option(ALREADY_HAVE_METIS "Suppress build of oomph-lib's own copy of metis" OFF)
option(ALREADY_HAVE_PARMETIS "Suppress build of oomph-lib's copy of parmetis"
       OFF)
option(ALREADY_HAVE_SUPERLU
       "Suppress build of oomph-lib's copy of serial SuperLU" OFF)
option(ALREADY_HAVE_SUPERLU_DIST
       "Suppress build of oomph-lib's copy of serial SuperLU_dist" OFF)
option(SUPPRESS_TRIANGLE_LIB
       "Suppress build of oomph-lib's copy of the triangle library" OFF)
option(SUPPRESS_TETGEN_LIB
       "Suppress build of oomph-lib's copy of the tetgen library" OFF)

# Support for deprecated packages:

option(OOMPH_USE_DEPRECATED_SUPERLU
       "Use oomph-lib's deprecated version of SuperLU" OFF)

# -------------------------[ Process selected options ]-------------------------

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

# Add PARANOID and RANGE_CHECKING for Debug configurations if desired. Have to
# use generator expressions to make sure it works correctly for multi-config.
# generators like Xcode, Visual Studio or Ninja Multi-Config
if(OOMPH_ENABLE_PARANOID)
  add_compile_definitions("$<$<CONFIG:DEBUG>:PARANOID>")
endif()
if(OOMPH_ENABLE_RANGE_CHECKING)
  add_compile_definitions("$<$<CONFIG:DEBUG>:RANGE_CHECKING>")
endif()

# If we're using up-to-date external sources/dists then make sure we have a
# preprocessor definition to handle changes that need to be made from within
# source files (e.g. API changes)
if(OOMPH_BLEEDING_EDGE)
  add_compile_definitions(OOMPH_BLEEDING_EDGE)
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

# Set up CUDA functionality (nothing in the library uses CUDA at the moment, so
# I'm just doing the bare minimum to incorporate it for now)
if(OOMPH_ENABLE_CUDA)
  find_package(CUDA REQUIRED)
endif()

# --------------------------[ Find required headers ]---------------------------

include(OomphCheckForRequiredHeaders)

# Enable @rpath in the install name for any shared library being built note:
# this is overwritten by INSTALL_NAME_DIR. See:
# https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling#:~:text=Since%20version%202.8.,overridden%20for%20the%20install%20tree.
set(CMAKE_MACOSX_RPATH TRUE)

# --------------------------[ Add oomph-lib-config.h ]--------------------------

set(OOMPH_ADD_CONFIG_H TRUE)
if(OOMPH_ADD_CONFIG_H)
  add_compile_definitions(HAVE_CONFIG_H)
endif()

# -------------------------[ Build & install library ]--------------------------

message(VERBOSE "<====== Starting oomph-lib library build ======>")

# Process the oomph-lib directories containing the library files
add_subdirectory(scripts)
add_subdirectory(external_src)
add_compile_definitions(${EXTERNAL_SRC_CXX_DEFINITIONS})
add_subdirectory(external_distributions)
add_compile_definitions(${EXTERNAL_DIST_CXX_DEFINITIONS})
add_subdirectory(src)

message(VERBOSE "<====== Finished oomph-lib library build ======>")

# Now export the package so it can be found using find_package(...)
include(OomphExportPackage)

# ---------------------------------[ Testing ]----------------------------------

# The last part to the build is compiling tests. The demo_drivers are currently
# working when run *outside* the library build. The integration of their build
# with the build of the library itself has not yet been tested. I think they
# should be fine, but there might be an issue with find_package(oomphlib). Note
# that private/ and user_src/ and user_drivers/ are not built with the library.

# Enable testing for the current directory and all child directories
if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME) AND OOMPH_ENABLE_TESTING)
  add_subdirectory(demo_drivers)
  add_subdirectory(self_tests)
  find_package(Doctest)
  enable_testing()
  add_subdirectory(tests)
endif()

# -----------------------[ Configure oomph-lib-config.h ]-----------------------

# Construct a oomph-lib-specific config file. Currently
if(OOMPH_ADD_CONFIG_H)
  configure_file(${CMAKE_SOURCE_DIR}/cmake/oomph-lib-config.h.in
                 ${CMAKE_SOURCE_DIR}/oomph-lib-config.h @ONLY)
  configure_file(${CMAKE_SOURCE_DIR}/cmake/oomph-lib-config.h.in
                 ${CMAKE_BINARY_DIR}/src/oomph-lib-config.h @ONLY)
  install(FILES ${CMAKE_BINARY_DIR}/src/oomph-lib-config.h
          DESTINATION "${OOMPH_INCLUDE_INSTALL_DIR}/${PROJECT_NAME}")
endif()

# -------------------[ Debug shit: delete after overhaul ]----------------------

include(OomphPrintCompilerDefinitions)

# --------------------------------[ Packaging ]---------------------------------

# Work in progress...

# ------------------------------------------------------------------------------
