# ==============================================================================
# Build the Hypre library.
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered hypre subdirectory")

# Look for an already-installed version of Hypre. If OOMPH_USE_HYPRE_FROM is
# defined then we'll look under there
if(OOMPH_USE_HYPRE_FROM)
  find_package(
    HYPRE 2.27.0 REQUIRED
    PATHS "${OOMPH_USE_HYPRE_FROM}"
    NO_DEFAULT_PATH)
endif()
find_package(HYPRE 2.27.0 QUIET)

# If we couldn't find a pre-installed version, lets build it ourselves
if(HYPRE_FOUND)
  message(STATUS "Found an existing installation of HYPRE at: ${HYPRE_DIR}")
else()
  # Hypre build options
  set(HYPRE_ENABLE_SHARED OFF)
  set(HYPRE_ENABLE_BIGINT OFF)
  set(HYPRE_ENABLE_MIXEDINT OFF)
  set(HYPRE_ENABLE_SINGLE OFF)
  set(HYPRE_ENABLE_LONG_DOUBLE OFF)
  set(HYPRE_ENABLE_COMPLEX OFF)
  set(HYPRE_WITH_MPI ${OOMPH_ENABLE_MPI})
  set(HYPRE_WITH_OPENMP OFF)
  set(HYPRE_ENABLE_HOPSCOTCH OFF)
  set(HYPRE_WITH_DSUPERLU OFF)
  set(HYPRE_PRINT_ERRORS OFF)
  set(HYPRE_TIMING OFF)
  set(HYPRE_BUILD_EXAMPLES OFF)
  set(HYPRE_BUILD_TESTS OFF)
  # set(HYPRE_WITH_EXTRA_CFLAGS) set(HYPRE_WITH_EXTRA_CXXFLAGS)

  # Should Hypre just use its own internal BLAS/LAPACK libraries?
  set(HYPRE_ENABLE_HYPRE_BLAS OFF)
  set(HYPRE_ENABLE_HYPRE_LAPACK OFF)

  # Specify the absolute paths to the BLAS/LAPACK libraries. If we're using
  # OpenBLAS, specify the path to the libraries that we found. If not, and we're
  # just going to build our own, then use a generator expression to define the
  # path to the library (the value will be populated properly at build-time,
  # when the actual library artifact is generated)

  # BLAS/LAPACK
  set(TPL_ENABLE_BLAS ON)
  set(TPL_ENABLE_LAPACK ON)
  if(OOMPH_ENABLE_USE_OPENBLAS AND BLAS_FOUND)
    set(TPL_BLAS_LIBRARIES "${BLAS_LIBRARIES}")
    set(TPL_LAPACK_LIBRARIES "${LAPACK_LIBRARIES}")
  else()
    message(
      FATAL_ERROR "Currently unable to link to oomph-lib BLAS/LAPACK libraries."
    )
    set(TPL_BLAS_LIBRARIES oomph_blas)
    set(TPL_LAPACK_LIBRARIES oomph_lapack)
  endif()

  # Make it look like we've changed the installation path so that Hypre won't
  # try to set its own custom installation path...
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT FALSE)
  endif()

  # TODO: Fork Hypre repo and update install paths so they are in the lib/hypre,
  # include/hypre directories...

  # Download HYPRE repo
  include(FetchContent)
  FetchContent_Declare(
    hypre
    GIT_REPOSITORY https://github.com/hypre-space/hypre.git
    GIT_TAG v2.27.0
    GIT_PROGRESS TRUE
    GIT_SHALLOW TRUE
    SOURCE_SUBDIR "src")
  FetchContent_MakeAvailable(hypre)
endif()

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and make the updated
# variable values visible to the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_HYPRE)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_HYPRE TRUE CACHE INTERNAL "")

message(STATUS "Leaving hypre subdirectory")
# ------------------------------------------------------------------------------
