# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
#
# ...to be filled in...
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

# Where to locate Boost/CGAL and the versions we need
set(BOOST_GIT_URL https://github.com/boostorg/boost.git)
set(BOOST_GIT_TAG boost-1.83.0)
set(CGAL_GIT_URL https://github.com/CGAL/cgal.git)
set(CGAL_GIT_TAG v6.0.1)

# Set the default installation paths
set(BOOST_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/boost")
set(CGAL_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/cgal")

# If we've already been given Boost, we'll use the files from there
if(OOMPH_USE_BOOST_FROM)
  set(BOOST_INSTALL_DIR "${OOMPH_USE_BOOST_FROM}")
endif()

# ----------------------------------------
# BOOST
# ----------------------------------------
if(NOT OOMPH_USE_BOOST_FROM)
  oomph_get_external_project_helper(
    PROJECT_NAME   boost
    GIT_REPOSITORY ${BOOST_GIT_URL}
    GIT_TAG ${BOOST_GIT_TAG}
    GIT_SUBMODULES_RECURSE TRUE
    INSTALL_DIR "${BOOST_INSTALL_DIR}"
    CONFIGURE_COMMAND ./bootstrap.sh --prefix=<INSTALL_DIR> --with-libraries=thread,system,program_options CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER}
    BUILD_COMMAND ./b2 install --jobs=${NUM_JOBS}
    INSTALL_COMMAND ""
  )
endif()

# ----------------------------------------
# CGAL
# ----------------------------------------
set(
  CGAL_CMAKE_ARGS
  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  -DCGAL_CMAKE_EXACT_NT_BACKEND=BOOST_BACKEND  # Use Boost.Multiprecision instead of GMP/MPFR
  -DCGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE=ON
  -DCGAL_DISABLE_GMP=ON
  -DCMAKE_DISABLE_FIND_PACKAGE_GMP=ON
)

set(
  CGAL_SELFTEST_ARGS
  ${CGAL_CMAKE_ARGS}
  -DBoost_NO_SYSTEM_PATHS=ON
  -DBoost_ROOT=${BOOST_INSTALL_DIR}
  -DWITH_CGAL_Qt5=OFF
)

oomph_get_external_project_helper(
  PROJECT_NAME cgal
  GIT_REPOSITORY ${CGAL_GIT_URL}
  GIT_TAG ${CGAL_GIT_TAG}
  INSTALL_DIR ${CGAL_INSTALL_DIR}
  PATCH_COMMAND ${CMAKE_CURRENT_LIST_DIR}/patches/patch_cgal.sh <SOURCE_DIR>
  CONFIGURE_COMMAND ${CMAKE_COMMAND} --install-prefix=<INSTALL_DIR> -G=${CMAKE_GENERATOR} ${CGAL_CMAKE_ARGS} -B=build
  BUILD_COMMAND ${CMAKE_COMMAND} --build build --parallel ${NUM_JOBS} --verbose
  INSTALL_COMMAND ${CMAKE_COMMAND} --install build
  TEST_COMMAND ${CMAKE_CURRENT_LIST_DIR}/scripts/run_cgal_self_test.sh <SOURCE_DIR> <LOG_DIR> ${CGAL_SELFTEST_ARGS}
)

# ----------------------------------------
# Set order of dependencies
# ----------------------------------------
# Make sure Boost gets built before CGAL if we're building it
if(TARGET boost)
  add_dependencies(cgal boost)
endif()

# -----------------------------------------------------------------------------
# cmake-format: on
