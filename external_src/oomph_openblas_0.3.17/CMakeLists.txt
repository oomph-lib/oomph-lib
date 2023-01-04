# ==============================================================================
# Blas sources
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oomph_blas subdirectory")

set(BUILD_TESTING OFF)
set(C_LAPACK OFF)
set(BUILD_WITHOUT_CBLAS ON)
set(DYNAMIC_ARCH OFF)
set(DYNAMIC_OLDER OFF)
set(BUILD_RELAPACK OFF)
set(USE_LOCKING OFF)
set(USE_PERL OFF)
set(NO_AFFINITY ON)
set(CPP_THREAD_SAFETY_TEST OFF)
set(CPP_THREAD_SAFETY_GEMV OFF)
set(BUILD_STATIC_LIBS ON)
if(NOT BUILD_SHARED_LIBS)
  set(BUILD_STATIC_LIBS OFF)
endif()

cmake_policy(SET CMP0077 NEW)
include(FetchContent)
FetchContent_Declare(
  openblas
  GIT_REPOSITORY https://github.com/xianyi/OpenBLAS.git
  GIT_TAG v0.3.21
  GIT_PROGRESS TRUE
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(openblas)

# PM: Decide to how/where to define our BLAS/LAPACK libraries
set(BLAS_LIBRARIES BLAS::BLAS)
set(LAPACK_LIBRARIES LAPACK::LAPACK)

# ------------------------------------------------------------------------------
# BOILERPLATE: Update variables in parent scope

# Add the library to the list of external source libraries and make the updated
# variable value visible to the parent scope
list(APPEND EXTERNAL_SRC_LIBS ${LIBNAME})
set(EXTERNAL_SRC_LIBS ${EXTERNAL_SRC_LIBS} PARENT_SCOPE)

message(VERBOSE "Leaving oomph_blas subdirectory")
# ------------------------------------------------------------------------------
