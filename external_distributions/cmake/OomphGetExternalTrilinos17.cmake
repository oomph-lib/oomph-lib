# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
#
# Build and install Trilinos for use by oomph-lib.
#
# This version targets Trilinos 17.0.0 and deliberately enables only the
# packages required for the initial migration of oomph-lib's iterative linear
# solver wrappers from Epetra/AztecOO to Tpetra/Belos.
#
# =============================================================================
include_guard()

set(
  TRILINOS_URL
  https://github.com/trilinos/Trilinos/archive/refs/tags/trilinos-release-17-0-0.tar.gz)

set(TRILINOS_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/trilinos17")

# Match CMake's default extension mode instead of passing an empty cache value
# to Trilinos, which prevents Kokkos from exporting its extension support.
if(NOT DEFINED CMAKE_CXX_EXTENSIONS OR "${CMAKE_CXX_EXTENSIONS}" STREQUAL "")
  set(TRILINOS_CXX_EXTENSIONS ON)
else()
  set(TRILINOS_CXX_EXTENSIONS "${CMAKE_CXX_EXTENSIONS}")
endif()

# Trilinos 17 requires CMake 3.27 or newer and a C++20-compatible compiler.
if(CMAKE_VERSION VERSION_LESS 3.27)
  message(FATAL_ERROR "Trilinos 17 requires CMake 3.27 or newer.")
endif()

# On Ubuntu, Trilinos doesn't appear to link to gfortran when using OpenBLAS,
# resulting in the error described here:
#           https://github.com/trilinos/Trilinos/issues/8632
# so we won't run the tests on Linux for now. Ideally we should provide
# Trilinos with the libgfortran library when we pass OpenBLAS. It only breaks
# the Trilinos tests for now though.
set(ENABLE_TRILINOS_TESTS OFF)
if(OOMPH_ENABLE_THIRD_PARTY_LIBRARY_TESTS)
  if(APPLE)
    set(ENABLE_TRILINOS_TESTS ON)
  endif()
endif()

set(
  TRILINOS_CMAKE_CONFIGURE_ARGS
  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  -DCMAKE_CXX_STANDARD=20
  -DCMAKE_CXX_EXTENSIONS=${TRILINOS_CXX_EXTENSIONS}
  -DCMAKE_INSTALL_PREFIX=${TRILINOS_INSTALL_DIR}
  -DTrilinos_ENABLE_TESTS=${ENABLE_TRILINOS_TESTS}
  -DTrilinos_ENABLE_Fortran=ON
  -DTrilinos_ENABLE_EXAMPLES=OFF
  -DTrilinos_ENABLE_ALL_PACKAGES=OFF
  -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF
  -DTrilinos_INSTALL_LIBRARIES_AND_HEADERS=ON
  -DTrilinos_ENABLE_INSTALL_CMAKE_CONFIG_FILES=ON
  -DKokkos_ENABLE_DEPRECATED_CODE_4=OFF)

# Minimal Trilinos 17 stack for oomph-lib's new iterative solver wrappers:
#
#   Teuchos  : Core utilities and parameter lists.
#
#   Tpetra   : Distributed sparse matrices, vectors and operators. This
#              replaces the role played by Epetra in the old wrappers.
#
#   Belos    : Iterative linear solvers. This provides the new implementations
#              of GMRES, CG and BiCGStab that will replace AztecOO.
#
#   Ifpack2  : Successor to Ifpack. This will ultimately be required for
#              Trilinos-native preconditioners. The first Belos wrapper can
#              still use an oomph-lib Preconditioner through a Tpetra::Operator
#              adapter.
#
# Deliberately omitted for now:
#
#   Amesos2  : Add when the direct-solver wrappers are migrated.
#
#   MueLu    : Add when the old ML-based multigrid wrappers are migrated.
#
#   Anasazi  : Eigensolvers; unrelated to the initial iterative linear-solver
#              migration.
#
# Kokkos and KokkosKernels are dependencies of Tpetra/Ifpack2 and are enabled
# through Trilinos' dependency handling.
set(
  DESIRED_TRILINOS_PACKAGES
  Teuchos
  Tpetra
  Belos
  Ifpack2)

foreach(TRILINOS_PACKAGE IN LISTS DESIRED_TRILINOS_PACKAGES)
  list(
    APPEND
    TRILINOS_CMAKE_CONFIGURE_ARGS
    -DTrilinos_ENABLE_${TRILINOS_PACKAGE}=ON)
endforeach()

# Point Trilinos at our OpenBLAS installation if oomph-lib is building one.
if(OOMPH_BUILD_OPENBLAS)
  list(
    APPEND
    TRILINOS_CMAKE_CONFIGURE_ARGS
    -DTPL_ENABLE_BLAS=ON
    -DTPL_ENABLE_LAPACK=ON
    -DBLAS_LIBRARY_DIRS=${CMAKE_INSTALL_PREFIX}/openblas/lib
    -DLAPACK_LIBRARY_DIRS=${CMAKE_INSTALL_PREFIX}/openblas/lib
    -DBLAS_LIBRARY_NAMES=openblas
    -DLAPACK_LIBRARY_NAMES=openblas)
endif()

# Configure MPI support consistently with the remainder of oomph-lib's
# third-party-library build.
if(OOMPH_ENABLE_MPI)
  list(
    APPEND
    TRILINOS_CMAKE_CONFIGURE_ARGS
    -DTPL_ENABLE_MPI=ON
    -DMPI_USE_COMPILER_WRAPPERS=ON
    -DMPI_C_COMPILER=${MPI_C_COMPILER}
    -DMPI_CXX_COMPILER=${MPI_CXX_COMPILER}
    -DMPI_Fortran_COMPILER=${MPI_Fortran_COMPILER})

  # Find the base MPI installation directory from the MPI C++ include
  # directories. Trilinos can use this to locate the remainder of the MPI
  # installation.
  set(MPI_BASE_DIR "")

  foreach(MPI_INCLUDE_DIR IN LISTS MPI_CXX_INCLUDE_DIRS)
    get_filename_component(MPI_DIR "${MPI_INCLUDE_DIR}" DIRECTORY)

    if(EXISTS "${MPI_DIR}/bin")
      set(MPI_BASE_DIR ${MPI_DIR})
      message(STATUS "Found base MPI directory: ${MPI_BASE_DIR}")
      if(EXISTS "${MPI_DIR}/bin")
        message(STATUS "Yay! It also contains a bin/ directory!")
      endif()
      break()
    else()
      message(
        STATUS
        "Couldn't find root MPI directory from MPI_CXX_INCLUDE_DIRS.")
      message(
        STATUS
        "For Trilinos I will default to: -DMPI_BASE_DIR=${MPI_BASE_DIR}")
    endif()
  endforeach()

  # Now append to the full list of arguments.
  list(
    APPEND
    TRILINOS_CMAKE_CONFIGURE_ARGS
    -DMPI_BASE_DIR=${MPI_BASE_DIR})
endif()

# TODO: The TeuchosCore_TypeConversions_UnitTest and *_MPI_4 tests die on
# macOS. Not sure how to fix this yet, so filter them out for now.
set(CTEST_EXCLUDE_REGEX_STRING "")
if(APPLE)
  set(
    CTEST_EXCLUDE_REGEX_STRING
    "(TeuchosCore_TypeConversions_UnitTest|.*_MPI_4)")
endif()

# Define how to configure/build/install the project.
#
# IMPORTANT:
#
# Use the GitHub archive for the fixed Trilinos release rather than cloning
# the git repository. The Trilinos 17.0.0 git tree contains historical
# submodule metadata for packages/framework/GenConfig.bck which causes
# CMake's ExternalProject git-submodule initialisation to fail. None of those
# framework/CI repositories is required for our Trilinos library build.
oomph_get_external_project_helper(
  PROJECT_NAME trilinos17
  URL ${TRILINOS_URL}
  INSTALL_DIR ${TRILINOS_INSTALL_DIR}
  CONFIGURE_COMMAND
    ${CMAKE_COMMAND}
    -G=${CMAKE_GENERATOR}
    ${TRILINOS_CMAKE_CONFIGURE_ARGS}
    -B=build
  BUILD_COMMAND
    ${CMAKE_COMMAND}
    --build build
  INSTALL_COMMAND
    ${CMAKE_COMMAND}
    --install build
  TEST_COMMAND
    ${CMAKE_CTEST_COMMAND}
    --test-dir build
    -j ${OOMPH_NUM_JOBS}
    --output-on-failure
    --verbose
    -E "${CTEST_EXCLUDE_REGEX_STRING}")

# Trilinos depends on OpenBLAS. If we're building OpenBLAS ourselves then we
# need to make sure that it gets built before Trilinos.
if(TARGET openblas)
  add_dependencies(trilinos17 openblas)
endif()

# ---------------------------------------------------------------------------------
# cmake-format: on
