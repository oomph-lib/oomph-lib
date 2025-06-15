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
set(HYPRE_GIT_URL https://github.com/hypre-space/hypre.git)
set(HYPRE_GIT_TAG v2.32.0)
set(HYPRE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/hypre")

# Hypre build options
set(HYPRE_CMAKE_CONFIGURE_ARGS
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${HYPRE_INSTALL_DIR}
    -DHYPRE_ENABLE_SHARED=OFF
    -DHYPRE_ENABLE_BIGINT=OFF
    -DHYPRE_ENABLE_MIXEDINT=OFF
    -DHYPRE_ENABLE_SINGLE=OFF
    -DHYPRE_ENABLE_LONG_DOUBLE=OFF
    -DHYPRE_ENABLE_COMPLEX=OFF
    -DHYPRE_WITH_MPI=${OOMPH_ENABLE_MPI}
    -DHYPRE_WITH_OPENMP=OFF
    -DHYPRE_ENABLE_HOPSCOTCH=OFF
    -DHYPRE_WITH_DSUPERLU=OFF
    -DHYPRE_PRINT_ERRORS=OFF
    -DHYPRE_TIMING=OFF
    -DHYPRE_BUILD_EXAMPLES=OFF
    -DHYPRE_BUILD_TESTS=OFF
    -DHYPRE_ENABLE_HYPRE_BLAS=OFF
    -DHYPRE_ENABLE_HYPRE_LAPACK=OFF
    -DTPL_ENABLE_BLAS=ON
    -DTPL_ENABLE_LAPACK=ON
    -DTPL_BLAS_LIBRARIES=${OpenBLAS_LIBRARIES}
    -DTPL_LAPACK_LIBRARIES=${OpenBLAS_LIBRARIES})

# Define how to configure/build/install the project
oomph_get_external_project_helper(
  PROJECT_NAME hypre
  GIT_REPOSITORY ${HYPRE_GIT_URL}
  GIT_TAG ${HYPRE_GIT_TAG}
  INSTALL_DIR ${HYPRE_INSTALL_DIR}
  CONFIGURE_COMMAND ${CMAKE_COMMAND} ${HYPRE_CMAKE_CONFIGURE_ARGS} -G=${CMAKE_GENERATOR} -S=src -B=src/cmbuild
  BUILD_COMMAND ${CMAKE_COMMAND} --build src/cmbuild -j ${OOMPH_NUM_JOBS}
  INSTALL_COMMAND ${CMAKE_COMMAND} --install src/cmbuild
  # TEST_COMMAND    ./src/cmbuild/test/ij
)

# Hypre depends on OpenBLAS being built. If we're building OpenBLAS ourselves
# then we need to make sure that it gets built before Hypre
if(TARGET openblas)
  add_dependencies(hypre openblas)
endif()

# ---------------------------------------------------------------------------------
# cmake-format: on
