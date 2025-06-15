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
set(MUMPS_GIT_URL https://github.com/puneetmatharu/mumps.git)
set(MUMPS_GIT_TAG v5.6.2.5)
set(MUMPS_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/mumps")

# Should we run the tests?
set(ENABLE_MUMPS_TESTS OFF)
if(OOMPH_ENABLE_THIRD_PARTY_LIBRARY_TESTS)
  set(ENABLE_MUMPS_TESTS ON)
endif()

# MUMPS build options
set(MUMPS_CMAKE_CONFIGURE_ARGS
    -DMUMPS_UPSTREAM_VERSION=5.6.2
    -Dparallel=${OOMPH_ENABLE_MPI}
    # -Dfind_static=ON # FIXME: Doesn't work on macOS!
    -Dgemmt=ON
    -Dintsize64=OFF
    -Dscotch=OFF
    -Dparmetis=OFF
    -Dmetis=OFF
    -Dopenmp=OFF
    -Dmatlab=OFF
    -Doctave=OFF
    -Dfind=OFF
    -DBUILD_TESTING=${ENABLE_MUMPS_TESTS}
    -DBUILD_SHARED_LIBS=OFF
    -DBUILD_SINGLE=ON
    -DBUILD_DOUBLE=ON
    -DBUILD_COMPLEX=OFF
    -DBUILD_COMPLEX16=OFF
    -DLAPACK_VENDOR=OpenBLAS
    -DLAPACK_ROOT=${OpenBLAS_ROOT}
    -DLAPACK_INCLUDE_DIR=${OpenBLAS_ROOT}/include
)

# Define how to configure/build/install the project
oomph_get_external_project_helper(
  PROJECT_NAME mumps
  GIT_REPOSITORY ${MUMPS_GIT_URL}
  GIT_TAG ${MUMPS_GIT_TAG}
  INSTALL_DIR "${MUMPS_INSTALL_DIR}"
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> ${MUMPS_CMAKE_CONFIGURE_ARGS} -G=${CMAKE_GENERATOR} -B=build
  BUILD_COMMAND ${CMAKE_COMMAND} --build build
  INSTALL_COMMAND ${CMAKE_COMMAND} --install build
  TEST_COMMAND ${CMAKE_CTEST_COMMAND} --test-dir build -j ${NUM_JOBS} --output-on-failure)

# If we're building OpenBLAS, make sure we build it before we get around to
# building MUMPS
if(TARGET openblas)
  add_dependencies(mumps openblas)
endif()

# ---------------------------------------------------------------------------------
# cmake-format: on
