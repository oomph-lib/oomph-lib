# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
#
# NOTE: The OpenBLAS installation automatically runs self-tests but it's hard to
# extract stats from them (partly because I don't know how a failed test would
# be reported; there's no executive summary.
#
# USAGE:
# ------
#
# ...to be filled in...
#
# EXAMPLE:
# --------
#
# ...to be filled in...
#
# =============================================================================
include_guard()

set(TRILINOS_TARBALL_URL
    https://github.com/trilinos/Trilinos/archive/refs/tags/trilinos-release-16-0-0.tar.gz
)
set(TRILINOS_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/trilinos")

# On Ubuntu, Trilinos doesn't appear to link to gfortran when using OpenBLAS,
# resulting in error described here:
#           https://github.com/trilinos/Trilinos/issues/8632
# so we won't run the tests on Linux for now. Ideally we should provide Trilinos with the
# libgfortran library when we pass OpenBLAS. It only breaks the Trilinos tests for now though.
set(ENABLE_TRILINOS_TESTS OFF)
if(OOMPH_ENABLE_THIRD_PARTY_LIBRARY_TESTS)
  if(APPLE)
    set(ENABLE_TRILINOS_TESTS ON)
  endif()
endif()

set(TRILINOS_CMAKE_CONFIGURE_ARGS
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
    -DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}
    -DTrilinos_ENABLE_TESTS=${ENABLE_TRILINOS_TESTS}
    -DTrilinos_ENABLE_Fortran=ON
    -DTrilinos_ENABLE_EXAMPLES=OFF
    -DTrilinos_ENABLE_ALL_PACKAGES=OFF
    -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF
    -DTrilinos_INSTALL_LIBRARIES_AND_HEADERS=ON
    -DTrilinos_ENABLE_INSTALL_CMAKE_CONFIG_FILES=ON
    -DKOKKOS_USE_CXX_EXTENSIONS=ON
    -DTPL_ENABLE_BLAS=ON
    -DTPL_ENABLE_LAPACK=ON
    -DTPL_BLAS_LIBRARIES=${OpenBLAS_LIBRARIES}
    -DTPL_LAPACK_LIBRARIES=${OpenBLAS_LIBRARIES}
    -DTPL_ENABLE_MPI=${OOMPH_ENABLE_MPI})

set(DESIRED_TRILINOS_PACKAGES
    Amesos
    Anasazi
    AztecOO
    Epetra
    EpetraExt
    Ifpack
    ML
    Teuchos
    Triutils)

# TODO: Handle deprecated packages more gracefully
message(
  WARNING
    "Disabling deprecated Trilinos package warnings. Do not move past v16.0.0 without making the necessary oomph-lib changes to handle this."
)

# Enable the package but disable the deprecated warning
foreach(TPL_PACKAGE IN LISTS DESIRED_TRILINOS_PACKAGES)
  list(APPEND TRILINOS_CMAKE_CONFIGURE_ARGS -DTrilinos_ENABLE_${TPL_PACKAGE}=ON)
  list(APPEND TRILINOS_CMAKE_CONFIGURE_ARGS
       -D${TPL_PACKAGE}_SHOW_DEPRECATED_WARNINGS=OFF)
endforeach()

if(OOMPH_ENABLE_MPI)
  if(NOT MPI_CXX_INCLUDE_DIRS)
    message(FATAL_ERROR "Requested MPI but MPI_CXX_INCLUDE_DIRS is not set!")
  endif()

  # Need to do some filtering on MPI_CXX_INCLUDE_DIRS if there are multiple
  # include paths as Trilinos doesn't accept a list of arguments to
  # MPI_BASE_DIR. The approach here is to just take the first path from that
  # list that has lib/ and include/ subfolders. If we can't find a match then
  # we'll just take the first path in the list

  # Default to first entry of MPI_CXX_INCLUDE_DIRS as base directory
  list(GET MPI_CXX_INCLUDE_DIRS 0 MPI_BASE_DIR)

  # Loop over the include directories, look at its parent directory to see
  # whether it has the required bin/, include/, and lib/ dirs. Would prefer to
  # just pass MPI_CXX_INCLUDE_DIRS but there appear to be issues with making
  # sure multiple path values are interpreted correctly.
  foreach(MPI_INCLUDE_DIR IN LISTS MPI_CXX_INCLUDE_DIRS)
    cmake_path(GET MPI_INCLUDE_DIR PARENT_PATH MPI_DIR)
    message(STATUS "Checking if ${MPI_DIR} is the root MPI directory")

    # See if it has an include/ and lib/ directory (and optionally a bin/
    # directory)
    if((EXISTS "${MPI_DIR}/include") AND (EXISTS "${MPI_DIR}/lib"))
      set(MPI_BASE_DIR ${MPI_DIR})
      message(STATUS "Found base MPI directory: ${MPI_BASE_DIR}")
      if(EXISTS "${MPI_DIR}/bin")
        message(STATUS "Yay! It also contains a bin/ directory!")
      endif()
      break()
    else()
      message(
        STATUS "Couldn't find root MPI directory from MPI_CXX_INCLUDE_DIRS.")
      message(
        STATUS "For Trilinos I will default to: -DMPI_BASE_DIR=${MPI_BASE_DIR}")
    endif()
  endforeach()

  # Now append to the full list of arguments
  list(APPEND TRILINOS_CMAKE_CONFIGURE_ARGS -DMPI_BASE_DIR=${MPI_BASE_DIR})
endif()


set(CTEST_EXCLUDE_REGEX_STRING "")
if(APPLE)
  # TODO: The TeuchosCore_TypeConversions_UnitTest and *_MPI_4 tests die on macOS. Not sure
  # how to fix this yet why so we'll just filter it out for now. Need to come back to this.
  set(CTEST_EXCLUDE_REGEX_STRING "(TeuchosCore_TypeConversions_UnitTest|.*_MPI_4)")
endif()

# There's an issue with one of the test files; we need to patch it
set(PATCH_FILE ${CMAKE_CURRENT_LIST_DIR}/patches/0001-Patch-packages-amesos-test-Test_Basic-Amesos_TestDri.patch)

# Look for 'patch'
find_program(PATCH_EXECUTABLE NAMES patch REQUIRED)

# Throw an error if it's not found
if(NOT PATCH_EXECUTABLE)
  message(FATAL_ERROR "Error: 'patch' command not found! Please install patch (e.g., 'apt install patch' or 'brew install patch').")
endif()

# Define how to configure/build/install the project
oomph_get_external_project_helper(
  PROJECT_NAME trilinos
  URL "${TRILINOS_TARBALL_URL}"
  INSTALL_DIR "${TRILINOS_INSTALL_DIR}"
  PATCH_COMMAND ${PATCH_EXECUTABLE} -p1 -i ${PATCH_FILE}
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -G=${CMAKE_GENERATOR} ${TRILINOS_CMAKE_CONFIGURE_ARGS} -B=build
  BUILD_COMMAND ${CMAKE_COMMAND} --build build -j ${NUM_JOBS}
  INSTALL_COMMAND ${CMAKE_COMMAND} --install build
  TEST_COMMAND ${CMAKE_CTEST_COMMAND} --test-dir build -j ${NUM_JOBS} --output-on-failure --verbose -E "${CTEST_EXCLUDE_REGEX_STRING}")

# Trilinos depends on OpenBLAS. If we're building OpenBLAS ourselves then we
# need to make sure that it gets built before Trilinos
if(TARGET openblas)
  add_dependencies(trilinos openblas)
endif()

# ---------------------------------------------------------------------------------
# cmake-format: on
