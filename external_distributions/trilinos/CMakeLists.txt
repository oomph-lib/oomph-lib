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

# BLAS/LaPACK
set(TPL_ENABLE_BLAS ON)
set(TPL_ENABLE_LAPACK ON)
if(OOMPH_ENABLE_USE_OPENBLAS AND BLAS_FOUND)
  set(TPL_BLAS_LIBRARIES "${BLAS_LIBRARIES}")
  set(TPL_LAPACK_LIBRARIES "${LAPACK_LIBRARIES}")
else()
  set(TPL_BLAS_LIBRARIES "$<TARGET_LINKER_FILE:oomph::oomph_blas>")
  set(TPL_LAPACK_LIBRARIES "$<TARGET_LINKER_FILE:oomph::oomph_lapack>")
endif()

# MPI
set(TPL_ENABLE_MPI ${OOMPH_ENABLE_MPI})
if(OOMPH_ENABLE_MPI AND MPI_FOUND)
  set(MPI_BASE_DIR "${MPI_CXX_INCLUDE_DIRS}")
endif()

include(FetchContent)
FetchContent_Declare(
  Trilinos
  GIT_REPOSITORY https://github.com/trilinos/Trilinos.git
  GIT_TAG trilinos-release-13-4-0
  GIT_PROGRESS TRUE
  GIT_SHALLOW TRUE)
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
