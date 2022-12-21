# ==============================================================================
# Build the Trilinos library
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered Trilinos subdirectory")

# Build flags
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

# BLAS/LaPACK
set(TPL_ENABLE_BLAS ON)
set(TPL_ENABLE_LAPACK ON)
if(OOMPH_ENABLE_USE_OPENBLAS AND BLAS_FOUND)
  set(TPL_BLAS_LIBRARIES "${BLAS_LIBRARIES}")
  set(TPL_LAPACK_LIBRARIES "${LAPACK_LIBRARIES}")
else()
  # TODO: Try replacing this with the targets themselves
  set(TPL_BLAS_LIBRARIES "$<TARGET_LINKER_FILE:oomph::oomph_blas>")
  set(TPL_LAPACK_LIBRARIES "$<TARGET_LINKER_FILE:oomph::oomph_lapack>")
endif()

# MPI
set(TPL_ENABLE_MPI ${OOMPH_ENABLE_MPI})
if(OOMPH_ENABLE_MPI AND MPI_FOUND)
  set(MPI_BASE_DIR "${MPI_CXX_INCLUDE_DIRS}")
endif()

# NOTE (1): Grabbing a specific commit that adds modern CMake-ification to the
# Trilinos library (see https://github.com/trilinos/Trilinos/pull/10614). It's
# not a fully-fledged release yet so we have to use a specific commit hash
#
# NOTE (2): "GIT_SHALLOW TRUE" won't work while GIT_TAG corresponds to a commit.
# Once they've made an official 14.0 release, we can use that with a shallow
# clone
include(FetchContent)
FetchContent_Declare(
  Trilinos
  GIT_REPOSITORY https://github.com/trilinos/Trilinos.git
  GIT_TAG 963e937799ebb77f1f59f92bcccc0236d4e3dbfb
  GIT_PROGRESS TRUE)
FetchContent_MakeAvailable(Trilinos)

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
