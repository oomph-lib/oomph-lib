# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
#
# NOTE: The OpenBLAS installation automatically runs self-tests but it's hard
# to extract stats from them (partly because I don't know how a failed test
# would be reported; there's no executive summary.
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

# Where to get the code from and where to install it to
set(OPENBLAS_GIT_URL https://github.com/OpenMathLib/OpenBLAS.git)
set(OPENBLAS_GIT_TAG v0.3.29)
set(OPENBLAS_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/openblas")

# MUMPS build options
set(OPENBLAS_CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
  -DCMAKE_BUILD_TYPE=Release
  -DBUILD_SHARED_LIBS=OFF
)

# Define how to configure/build/install the project
oomph_get_external_project_helper(
  PROJECT_NAME openblas
  GIT_REPOSITORY ${OPENBLAS_GIT_URL}
  GIT_TAG ${OPENBLAS_GIT_TAG}
  INSTALL_DIR "${OPENBLAS_INSTALL_DIR}"
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -G=${CMAKE_GENERATOR} -B=build ${OPENBLAS_CMAKE_ARGS}
  BUILD_COMMAND ${CMAKE_COMMAND} --build build
  INSTALL_COMMAND ${CMAKE_COMMAND} --install build)

# Define the global variables OpenBLAS_ROOT and OpenBLAS_LIBRARIES for MUMPS,
# HYPRE and Trilinos to use
set(OpenBLAS_LIBNAME "${CMAKE_STATIC_LIBRARY_PREFIX}openblas${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(OpenBLAS_LIBRARIES "${OPENBLAS_INSTALL_DIR}/lib/${OpenBLAS_LIBNAME}" CACHE PATH "" FORCE)
set(OpenBLAS_ROOT ${OPENBLAS_INSTALL_DIR} CACHE PATH "" FORCE)

# ---------------------------------------------------------------------------------
# cmake-format: on
