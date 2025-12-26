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
# Matthias Heil and Andrew Hazel of the Department of Mathematics at The
# University of Manchester, along with many other contributors.
# ==============================================================================
# cmake-format: on

# ------------------------------------------------------------------------------
# OOMPH auto-inserted: warn once if user configures inside an existing build.
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt"
   AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(
    WARNING "=====================================================\n"
            "You're re-configuring an existing build tree. \n"
            "Delete the whole build folder and start with a fresh \n"
            "directory if your attempt to re-configure leads to \n"
            "strange error messages below here. \n"
            "=====================================================\n")
endif()
# ------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(
  oomphlib
  DESCRIPTION "A finite-element library for simulating multi-physics problems."
  VERSION 3.0.1
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
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build mode for performance."
      FORCE)
endif()

# Flag that we'll use to enable Debug-specific options
set(HAVE_DEBUG_BUILD FALSE)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(HAVE_DEBUG_BUILD TRUE)
endif()

# Limits the parallelism in ninja? Note this needs to go before any
# targets are created, otherwise ninja gets confused. Here this means
# it needs to go before the inclusion of OomphConfigureCMakeInfoForProject
if (CMAKE_GENERATOR MATCHES "Ninja")
  if (OOMPH_DO_NOT_LIMIT_NINJA_PARALLEL_JOBS)
    message(VERBOSE "Not limiting the number of parallel ninja jobs when building oomph-lib")
  else()
    message(VERBOSE "Limiting the number of parallel ninja jobs when building oomph-lib")
    include(OomphLimitNinjaJobs)
    include(ProcessorCount)
    ProcessorCount(NPROC)
    if (NPROC GREATER 2)
        # By default, ninja uses nproc+2! Let's go
        # the other way for safety
        math(EXPR oomph_NJOBS "${NPROC} - 2")
    else()
        set(oomph_NJOBS 1)
    endif()

    # Do it (for oomph-lib main build only!)
    oomph_limit_ninja_jobs(${oomph_NJOBS})
  endif()
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

# Flag to bypass our limits on ninja parallel jobs: default: do the
# limiting because we don't want to melt innocent users' underpowered laptops.
# Note that this limit is not automatically exported to
# stand-alone projects that use oomph-lib. Users have to
# explicitly call the relevant function to limit the number of ninja jobs
# in their own CMakeLists.txt files. (We do this in our demo drivers
# so have a look there to see how it's done.)
option(OOMPH_DO_NOT_LIMIT_NINJA_PARALLEL_JOBS "Allow ninja to decide on the number of parallel jobs" OFF)

# Regular CMake option; doesn't influence the build so we won't print its value at the end
option(OOMPH_PRINT_SETTINGS_AFTER_INSTALL "Display the oomph-lib settings at the end of the install step" ON)



# Limit the number of cores used by ninja to (number of processors - 2)
# (rather than (number of processors + 2) which is the default and
# causes under-powered computers to die).
include(ProcessorCount)
ProcessorCount(NPROC)
if (NPROC GREATER 2)
    math(EXPR NJOBS "${NPROC} - 2")
else()
    set(NJOBS 1)
endif()
set_property(GLOBAL PROPERTY JOB_POOLS limit_pool=${NJOBS})
set(CMAKE_JOB_POOL_COMPILE limit_pool)
set(CMAKE_JOB_POOL_LINK    limit_pool)
set(CMAKE_JOB_POOL_CUSTOM  limit_pool)
message(STATUS "Detected ${NPROC} cores → limiting Ninja to ${NJOBS} parallel jobs")

# Pass it on to sub-directories/sub-projects
set(NJOBS "${NJOBS}" CACHE STRING "Number of Ninja parallel jobs")


# Control the installation
oomph_option(OOMPH_ALLOW_INSTALL_AS_SUPERUSER "Allow the user to install to the default system install path (if CMAKE_INSTALL_PREFIX is not set)" OFF)
oomph_option(OOMPH_INSTALL_HEADERS_AS_SYMLINKS "Install symlinks to the oomph-lib headers instead of copying them" OFF)

# General build options
oomph_option(OOMPH_DONT_SILENCE_USELESS_WARNINGS "Display (harmless) warnings from external_src/ and src/ that are silenced" OFF)
oomph_option(OOMPH_ENABLE_MPI "Enable the use of MPI for parallel processing" OFF)
oomph_option(OOMPH_ENABLE_MPI_OVERSUBSCRIPTION "Enable oversubscription with MPI" OFF)
oomph_option(OOMPH_ENABLE_PARANOID "Enable the PARANOID flag in Debug" ${HAVE_DEBUG_BUILD})
oomph_option(OOMPH_ENABLE_RANGE_CHECKING "Enable RANGE_CHECKING flag in Debug" ${HAVE_DEBUG_BUILD})
oomph_option(OOMPH_ENABLE_SANITIZER_SUPPORT "Enable sanitizer support in Debug builds" OFF)

# Choose the default linear solver. By default, we'll disable
# Mumps because it not only requires mpi (which would be easy to
# check) but also requires MPI to have been initialised in the driver
# code before it can be used; not good for all the serial self-tests.
# The user can set -DOOMPH_ENABLE_MUMPS_AS_DEFAULT_LINEAR_SOLVER=ON
# to override this setting.
oomph_option(OOMPH_ENABLE_MUMPS_AS_DEFAULT_LINEAR_SOLVER "Use MumpsSolver as the default linear solver (when available)" OFF)

# Fine-grained control of build of external_srcs.
# TODO: Discuss with MH; do we actually want these options anymore?
oomph_option(OOMPH_SUPPRESS_TRIANGLE_LIB "Suppress build of oomph-lib's copy of the triangle library" OFF)
oomph_option(OOMPH_SUPPRESS_TETGEN_LIB "Suppress build of oomph-lib's copy of the tetgen library" OFF)

# Extra flags for the user to provide
oomph_option(OOMPH_EXTRA_COMPILE_DEFINES "Additional C/C++ compile defines to pass to the oomph-lib build." "")

# Paths to packages that the user has installed themselves
oomph_path_option(
  FLAG OOMPH_USE_OPENBLAS_FROM
  DOCSTRING "Path to OpenBLAS installation."
  REQUIRED
)
oomph_path_option(
  FLAG OOMPH_USE_GKLIB_FROM
  DOCSTRING "Path to GKlib installation"
  REQUIRED
)
oomph_path_option(
  FLAG OOMPH_USE_METIS_FROM
  DOCSTRING "Path to METIS installation"
  REQUIRED
)
oomph_path_option(
  FLAG OOMPH_USE_SUPERLU_FROM
  DOCSTRING "Path to SuperLU installation"
  REQUIRED
)
oomph_path_option(
  FLAG OOMPH_USE_PARMETIS_FROM
  DOCSTRING "Path to ParMETIS installation"
)
oomph_path_option(
  FLAG OOMPH_USE_SUPERLU_DIST_FROM
  DOCSTRING "Path to SuperLU_DIST installation"
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
else()
  if(OOMPH_USE_SUPERLU_DIST_FROM)
    message(FATAL_ERROR "You can only enable SuperLU_DIST when MPI is enabled!")
  endif()
endif()

# Add PARANOID and RANGE_CHECKING for Debug configurations (ONLY!) if desired
if(OOMPH_ENABLE_PARANOID OR OOMPH_ENABLE_RANGE_CHECKING)
  if(NOT (CMAKE_BUILD_TYPE STREQUAL "Debug"))
    message(
      FATAL_ERROR
        "You can only enable OOMPH_ENABLE_PARANOID and OOMPH_ENABLE_RANGE_CHECKING when CMAKE_BUILD_TYPE is set to \"Debug\"; you're building with CMAKE_BUILD_TYPE set to \"${CMAKE_BUILD_TYPE}\"."
    )
  endif()
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

# Enable Address, UndefinedBehaviour, Leak sanitizer support, if desired
if(OOMPH_ENABLE_SANITIZER_SUPPORT)
  include(OomphAddSanitizerSupport)
endif()

# Do we want Mumps to override SuperLU as the default solver?
if(OOMPH_ENABLE_MUMPS_AS_DEFAULT_LINEAR_SOLVER)
  oomph_add_cxx_compile_definitions(OOMPH_ENABLE_MUMPS_AS_DEFAULT_LINEAR_SOLVER)
endif()

# If the user provided extra compile defines to use in the build
if(OOMPH_EXTRA_COMPILE_DEFINES)
  oomph_add_cxx_compile_definitions(${OOMPH_EXTRA_COMPILE_DEFINES})
endif()


# Install ninja limiter script even if we're not building with ninja
install(FILES cmake/OomphLimitNinjaJobs.cmake
        DESTINATION lib/cmake/oomphlib)



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

# FIXME: Having to duplicate CXX definitions for C code to make SuperLUDist code
# work...

message(VERBOSE "<====== Starting oomph-lib library configuration ======>")

# Locate the third-party libraries that the user has provided via the
# OOMPH_USE_<LIBRARY>_FROM flags
include(OomphLocateThirdPartyLibs)

# Add C++ definitions for all of the defined (optional) third-party libraries,
# e.g. -DOOMPH_HAS_CGAL, -DOOMPH_HAS_MUMPS, etc.
oomph_add_c_compile_definitions(${EXTERNAL_DIST_CXX_DEFINITIONS})
oomph_add_cxx_compile_definitions(${EXTERNAL_DIST_CXX_DEFINITIONS})

# Configure the external sources (i.e. third-party libraries that oomph-lib
# always ships with)
add_subdirectory(external_src)

# Add C++ definitions for all of the defined third-party sources that we build
# with oomph-lib, e.g. -DOOMPH_HAS_SUPERLU, -DOOMPH_HAS_TRIANGLE, etc.
oomph_add_c_compile_definitions(${EXTERNAL_SRC_CXX_DEFINITIONS})
oomph_add_cxx_compile_definitions(${EXTERNAL_SRC_CXX_DEFINITIONS})

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
