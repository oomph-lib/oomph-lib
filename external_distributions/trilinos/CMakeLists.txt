# ==============================================================================
# Build the Trilinos library
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered Trilinos subdirectory")

# Look for an already-installed version of Trilinos. If OOMPH_USE_TRILINOS_FROM
# is defined then we'll look under there
if(OOMPH_USE_TRILINOS_FROM)
  find_package(
    Trilinos 13.5 REQUIRED
    PATHS ${OOMPH_USE_TRILINOS_FROM}
    NO_DEFAULT_PATH)
endif()
find_package(Trilinos 13.5 QUIET)

# If we couldn't find a pre-installed version, lets build it ourselves
if(Trilinos_FOUND)
  message(
    STATUS "Found an existing installation of Trilinos at: ${Trilinos_DIR}")
else()
  # Configuration flags
  set(Trilinos_ENABLE_TESTS OFF)
  set(Trilinos_ENABLE_EXAMPLES OFF)
  set(Trilinos_ENABLE_ALL_PACKAGES OFF)
  set(Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES OFF)
  set(Trilinos_ENABLE_Amesos ON)
  set(Trilinos_ENABLE_Anasazi ON)
  set(Trilinos_ENABLE_AztecOO ON)
  set(Trilinos_ENABLE_Epetra ON)
  set(Trilinos_ENABLE_EpetraExt ON)
  set(Trilinos_ENABLE_Ifpack ON)
  set(Trilinos_ENABLE_ML ON)
  set(Trilinos_ENABLE_Teuchos ON)
  set(Trilinos_ENABLE_Triutils ON)
  set(Trilinos_INSTALL_LIBRARIES_AND_HEADERS ON)
  set(Trilinos_ENABLE_INSTALL_CMAKE_CONFIG_FILES ON)

  # BLAS/LAPACK
  set(TPL_ENABLE_BLAS ON)
  set(TPL_ENABLE_LAPACK ON)
  if(BLAS_FOUND AND LAPACK_FOUND)
    set(TPL_BLAS_LIBRARIES "${BLAS_LIBRARIES}")
    set(TPL_LAPACK_LIBRARIES "${LAPACK_LIBRARIES}")
  else()
    # Leads to a circular dependency. Think it might work if we make sure to
    # build oomph_blas/oomph_lapack before Trilinos and specify the explicit
    # path to the build path and install path (using $<BUILD_INTERFACE:...> and
    # $<INSTALL_INTERFACE:...>)
    message(
      FATAL_ERROR "Currently unable to link to oomph-lib BLAS/LAPACK libraries."
    )
    set(TPL_BLAS_LIBRARIES "$<TARGET_LINKER_FILE:oomph_blas>")
    set(TPL_LAPACK_LIBRARIES "$<TARGET_LINKER_FILE:oomph_lapack>")
  endif()

  # MPI
  set(TPL_ENABLE_MPI ${OOMPH_ENABLE_MPI})
  if(OOMPH_ENABLE_MPI AND MPI_FOUND)
    set(MPI_BASE_DIR "${MPI_CXX_INCLUDE_DIRS}")
  endif()

  # NOTE (1): Grabbing a specific commit where Trilinos ships libraries using
  # modern CMake (see https://github.com/trilinos/Trilinos/pull/10614). It's not
  # a fully-fledged release yet so we have to use a specific commit hash
  #
  # NOTE (2): "GIT_SHALLOW TRUE" won't work while GIT_TAG corresponds to a
  # specific commit (that doesn't correspond to a tag or branch HEAD). Once
  # they've made an official 14.0 release, we can use that with a shallow clone
  include(FetchContent)
  FetchContent_Declare(
    Trilinos
    GIT_REPOSITORY https://github.com/trilinos/Trilinos.git
    GIT_TAG 0270d8bd4cd3c62fb37a30535f78ee16b9c77574
    GIT_PROGRESS TRUE)
  # GIT_SHALLOW TRUE)
  FetchContent_MakeAvailable(Trilinos)
endif()

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and make the updated
# variable values visible to the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_TRILINOS)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_TRILINOS TRUE CACHE INTERNAL "")

message(VERBOSE "Leaving trilinos subdirectory")
# ------------------------------------------------------------------------------
