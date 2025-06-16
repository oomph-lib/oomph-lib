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
set(OPENBLAS_GIT_URL https://github.com/OpenMathLib/OpenBLAS.git)
set(OPENBLAS_GIT_TAG v0.3.29)
set(OPENBLAS_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/openblas")

if(APPLE)
  set(OPENBLAS_CMAKE_ARGS
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DBUILD_SHARED_LIBS=OFF
    -DBUILD_WITHOUT_LAPACK=OFF
    -DBUILD_WITHOUT_LAPACKE=OFF
    -DBUILD_WITHOUT_CBLAS=OFF
    -DBUILD_LAPACK_DEPRECATED=ON
    -DBUILD_TESTING=${OOMPH_ENABLE_THIRD_PARTY_LIBRARY_TESTS}
    -DBUILD_BENCHMARKS=OFF
    -DC_LAPACK=OFF
    -DDYNAMIC_ARCH=ON
  )

  # Define how to configure/build/install the project
  oomph_get_external_project_helper(
    PROJECT_NAME openblas
    GIT_REPOSITORY ${OPENBLAS_GIT_URL}
    GIT_TAG ${OPENBLAS_GIT_TAG}
    INSTALL_DIR "${OPENBLAS_INSTALL_DIR}"
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> ${OPENBLAS_CMAKE_ARGS} -G ${CMAKE_GENERATOR} -B build
    BUILD_COMMAND ${CMAKE_COMMAND} --build build -j ${OOMPH_NUM_JOBS}
    INSTALL_COMMAND ${CMAKE_COMMAND} --install build)

else()

  # Define how to configure/build/install the project
  oomph_get_external_project_helper(
    PROJECT_NAME openblas
    GIT_REPOSITORY ${OPENBLAS_GIT_URL}
    GIT_TAG ${OPENBLAS_GIT_TAG}
    INSTALL_DIR "${OPENBLAS_INSTALL_DIR}"
    CONFIGURE_HANDLED_BY_BUILD
    BUILD_COMMAND ${MAKE_EXECUTABLE} --jobs=${OOMPH_NUM_JOBS} CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER} FC=${CMAKE_Fortran_COMPILER}
    INSTALL_COMMAND ${MAKE_EXECUTABLE} --jobs=${OOMPH_NUM_JOBS} PREFIX=${OPENBLAS_INSTALL_DIR} install)

endif()

# Define the global variables OpenBLAS_ROOT and OpenBLAS_LIBRARIES for MUMPS,
# HYPRE and Trilinos to use
set(OpenBLAS_LIBNAME "${CMAKE_STATIC_LIBRARY_PREFIX}openblas${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(OpenBLAS_LIBRARIES "${OPENBLAS_INSTALL_DIR}/lib/${OpenBLAS_LIBNAME}" CACHE PATH "" FORCE)

# ---------------------------------------------------------------------------------
# cmake-format: on
