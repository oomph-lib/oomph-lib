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

# The tarballs for each library
set(BOOST_TARBALL_URL ${OOMPH_THIRD_PARTY_TAR_FILE_URL}/boost_1_83_0.tar.gz)
set(CGAL_TARBALL_URL https://github.com/CGAL/cgal/releases/download/v6.0.1/CGAL-6.0.1.tar.xz)

# Set the default installation paths
set(BOOST_INSTALL_DIR "${OOMPH_THIRD_PARTY_INSTALL_DIR}/boost")
set(CGAL_INSTALL_DIR "${OOMPH_THIRD_PARTY_INSTALL_DIR}/cgal")

# If we've already been given Boost, we'll use the files from there
if(OOMPH_USE_BOOST_FROM)
  set(BOOST_INSTALL_DIR "${OOMPH_USE_BOOST_FROM}")
endif()

# ----------------------------------------
# BOOST
# ----------------------------------------
if(NOT OOMPH_USE_BOOST_FROM)
  oomph_get_external_project_helper(
    PROJECT_NAME boost
    URL "${BOOST_TARBALL_URL}"
    INSTALL_DIR "${BOOST_INSTALL_DIR}"
    CONFIGURE_COMMAND ./bootstrap.sh --prefix=<INSTALL_DIR> --with-libraries=thread,system,program_options CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER}
    BUILD_COMMAND ./b2 install --jobs=${NUM_JOBS}
    INSTALL_COMMAND "")
endif()

# ----------------------------------------
# CGAL
# ----------------------------------------
set(
  CGAL_CMAKE_CONFIGURE_ARGS
  -DCMAKE_BUILD_TYPE=Release
  -DCGAL_CMAKE_EXACT_NT_BACKEND=BOOST_BACKEND # New: Use Boost.Multiprecision instead of GMP/MPFR
  -DCMAKE_VERBOSE_MAKEFILE=ON)

set(
  CGAL_CMAKE_CONFIGURE_ARGS_FOR_SELF_TEST
  ${CGAL_CMAKE_CONFIGURE_ARGS}
  -DCMAKE_BUILD_TYPE=Release
  -DBoost_NO_SYSTEM_PATHS=ON
  -DBoost_ROOT=${BOOST_INSTALL_DIR}
  -DWITH_CGAL_Qt5=OFF
)
oomph_get_external_project_helper(
  PROJECT_NAME cgal
  URL "${CGAL_TARBALL_URL}"
  INSTALL_DIR ${CGAL_INSTALL_DIR}
  PATCH_COMMAND ${CMAKE_CURRENT_LIST_DIR}/patches/patch_cgal.sh <SOURCE_DIR>
  CONFIGURE_COMMAND ${CMAKE_COMMAND} --install-prefix=<INSTALL_DIR> -G=${CMAKE_GENERATOR} ${CGAL_CMAKE_CONFIGURE_ARGS} -B=build
  BUILD_COMMAND ${CMAKE_COMMAND} --build build --parallel ${NUM_JOBS} --verbose
  INSTALL_COMMAND ${CMAKE_COMMAND} --install build
  TEST_COMMAND ${CMAKE_CURRENT_LIST_DIR}/scripts/run_cgal_self_test.sh <SOURCE_DIR> <LOG_DIR> ${CGAL_CMAKE_CONFIGURE_ARGS_FOR_SELF_TEST})

# -----------------------------------------------------------------------------

# Make sure Boost gets built before CGAL if we're building it
if(TARGET boost)
  add_dependencies(cgal boost)
endif()

# -----------------------------------------------------------------------------
# cmake-format: on
