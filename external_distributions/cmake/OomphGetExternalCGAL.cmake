# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
#
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
set(GMP_TARBALL_URL ${OOMPH_THIRD_PARTY_TAR_FILE_URL}/gmp-6.3.0.tar.xz)
set(MPFR_TARBALL_URL ${OOMPH_THIRD_PARTY_TAR_FILE_URL}/mpfr-4.2.1.tar.xz)
set(BOOST_TARBALL_URL ${OOMPH_THIRD_PARTY_TAR_FILE_URL}/boost_1_83_0.tar.gz)
set(CGAL_TARBALL_URL ${OOMPH_THIRD_PARTY_TAR_FILE_URL}/CGAL-5.6.tar.xz)

# Set the default installation paths
set(GMP_INSTALL_DIR "${OOMPH_THIRD_PARTY_INSTALL_DIR}/gmp")
set(MPFR_INSTALL_DIR "${OOMPH_THIRD_PARTY_INSTALL_DIR}/mpfr")
set(BOOST_INSTALL_DIR "${OOMPH_THIRD_PARTY_INSTALL_DIR}/boost")
set(CGAL_INSTALL_DIR "${OOMPH_THIRD_PARTY_INSTALL_DIR}/cgal")

# If we've already been given GMP, MPFR and/or Boost, we'll use the files from
# there
if(OOMPH_USE_GMP_FROM)
  set(GMP_INSTALL_DIR "${OOMPH_USE_GMP_FROM}")
endif()
if(OOMPH_USE_MPFR_FROM)
  set(MPFR_INSTALL_DIR "${OOMPH_USE_MPFR_FROM}")
endif()
if(OOMPH_USE_BOOST_FROM)
  set(BOOST_INSTALL_DIR "${OOMPH_USE_BOOST_FROM}")
endif()

# ----------------------------------------
# GMP
# ----------------------------------------
# Expected library path and include directory
set(GMP_C_LIBNAME ${CMAKE_STATIC_LIBRARY_PREFIX}gmp${CMAKE_STATIC_LIBRARY_SUFFIX})
set(GMP_C_LIBRARIES ${GMP_INSTALL_DIR}/lib/${GMP_C_LIBNAME} CACHE PATH "Path to GMP C libraries")
set(GMP_C_INCLUDE_DIR ${GMP_INSTALL_DIR}/include CACHE PATH "Path to GMP C include directory")

# If we need to build GMP
if(NOT OOMPH_USE_GMP_FROM)
  oomph_get_external_project_helper(
    PROJECT_NAME gmp
    URL "${GMP_TARBALL_URL}"
    INSTALL_DIR "${GMP_INSTALL_DIR}"
    CONFIGURE_COMMAND ./configure --prefix=<INSTALL_DIR> CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER}
    BUILD_COMMAND ${MAKE_EXECUTABLE} --jobs=${NUM_JOBS}
    INSTALL_COMMAND ${MAKE_EXECUTABLE} --jobs=${NUM_JOBS} install
    TEST_COMMAND ${MAKE_EXECUTABLE} check)
    # FIXME: Can't use INSTALL_BYPRODUCTS until CMake v3.26
    # INSTALL_BYPRODUCTS ${GMP_C_LIBRARIES} ${GMP_C_INCLUDE_DIR}/gmp.h)
endif()

# ----------------------------------------
# MPFR
# ----------------------------------------
# Expected library path and include directory
set(MPFR_LIBNAME ${CMAKE_STATIC_LIBRARY_PREFIX}mpfr${CMAKE_STATIC_LIBRARY_SUFFIX})
set(MPFR_LIBRARIES ${MPFR_INSTALL_DIR}/lib/${MPFR_LIBNAME} CACHE PATH "Path to GMP libraries")
set(MPFR_INCLUDE_DIR ${MPFR_INSTALL_DIR}/include CACHE PATH "Path to GMP include directory")

# If we need to build MPFR
if(NOT OOMPH_USE_MPFR_FROM)
  # Define how to configure/build/install the project
  oomph_get_external_project_helper(
    PROJECT_NAME mpfr
    URL "${MPFR_TARBALL_URL}"
    INSTALL_DIR "${MPFR_INSTALL_DIR}"
    CONFIGURE_COMMAND ./configure --prefix=<INSTALL_DIR> --with-gmp=${GMP_INSTALL_DIR} CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER}
    BUILD_COMMAND ${MAKE_EXECUTABLE} --jobs=${NUM_JOBS}
    INSTALL_COMMAND ${MAKE_EXECUTABLE} --jobs=${NUM_JOBS} install
    TEST_COMMAND ${MAKE_EXECUTABLE} check)
    # FIXME: Can't use INSTALL_BYPRODUCTS until CMake v3.26
    # INSTALL_BYPRODUCTS ${MPFR_LIBRARIES} "${MPFR_INCLUDE_DIR}/mpfr.h")
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
  CGAL_CMAKE_BUILD_ARGS_FOR_SELF_TEST
  -DCMAKE_BUILD_TYPE=Release
  -DGMP_INCLUDE_DIR=${GMP_C_INCLUDE_DIR}
  -DGMP_LIBRARIES=${GMP_C_LIBRARIES}
  -DMPFR_INCLUDE_DIR=${MPFR_INCLUDE_DIR}
  -DMPFR_LIBRARIES=${MPFR_LIBRARIES}
  -DBoost_NO_SYSTEM_PATHS=ON
  -DBoost_ROOT=${BOOST_INSTALL_DIR}
  -DWITH_CGAL_Qt5=OFF
)
oomph_get_external_project_helper(
  PROJECT_NAME cgal
  URL "${CGAL_TARBALL_URL}"
  INSTALL_DIR ${CGAL_INSTALL_DIR}
  PATCH_COMMAND ${CMAKE_CURRENT_LIST_DIR}/patches/patch_cgal.sh <SOURCE_DIR>
  CONFIGURE_COMMAND ${CMAKE_COMMAND} --install-prefix=<INSTALL_DIR> -G=${CMAKE_GENERATOR} -DCMAKE_BUILD_TYPE=Release -B=build
  BUILD_COMMAND ${CMAKE_COMMAND} --build build --parallel ${NUM_JOBS}
  INSTALL_COMMAND ${CMAKE_COMMAND} --install build
  TEST_COMMAND ${CMAKE_CURRENT_LIST_DIR}/scripts/run_cgal_self_test.sh <SOURCE_DIR> <LOG_DIR> ${CGAL_CMAKE_BUILD_ARGS_FOR_SELF_TEST})

# -----------------------------------------------------------------------------

# MPFR relies on GMP, and CGAL relies on GMP, MPFR and Boost. We have to be
# careful to only define target dependencies if we're building the libraries
# ourselves
if((TARGET gmp) AND (TARGET mpfr))
  add_dependencies(mpfr gmp)
  add_dependencies(cgal gmp mpfr)
endif()
if(TARGET boost)
  add_dependencies(cgal boost)
endif()

# -----------------------------------------------------------------------------
# cmake-format: on
