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

# TODO: Upload to oomph-lib repo
set(OPENBLAS_TARBALL_URL ${OOMPH_THIRD_PARTY_TAR_FILE_URL}/OpenBLAS-0.3.25.tar.gz)
set(OPENBLAS_INSTALL_DIR "${OOMPH_THIRD_PARTY_INSTALL_DIR}/openblas")

# Define how to configure/build/install the project
oomph_get_external_project_helper(
  PROJECT_NAME openblas
  URL "${OPENBLAS_TARBALL_URL}"
  INSTALL_DIR "${OPENBLAS_INSTALL_DIR}"
  CONFIGURE_HANDLED_BY_BUILD
  BUILD_COMMAND ${MAKE_EXECUTABLE} --jobs=${NUM_JOBS} CXX=${CMAKE_CXX_COMPILER}
                CC=${CMAKE_C_COMPILER} FC=${CMAKE_Fortran_COMPILER}
  INSTALL_COMMAND ${MAKE_EXECUTABLE} --jobs=${NUM_JOBS}
                  PREFIX=${OPENBLAS_INSTALL_DIR} install)

# Define the global variables OpenBLAS_ROOT and OpenBLAS_LIBRARIES for MUMPS,
# HYPRE and Trilinos to use
set(OpenBLAS_LIBNAME "${CMAKE_STATIC_LIBRARY_PREFIX}openblas${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(OpenBLAS_LIBRARIES "${OPENBLAS_INSTALL_DIR}/lib/${OpenBLAS_LIBNAME}" CACHE PATH "" FORCE)

# ---------------------------------------------------------------------------------
# cmake-format: on
