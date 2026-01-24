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
# cmake-format: on
include_guard()

# Where to get the code from and where to install it to
set(MUMPS_GIT_URL https://github.com/puneetmatharu/mumps.git)
set(MUMPS_GIT_TAG v5.8.1.3-pm)
set(MUMPS_UPSTREAM_VERSION 5.8.1)
set(MUMPS_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/mumps")

# MUMPS build options
set(MUMPS_CMAKE_CONFIGURE_ARGS
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${MUMPS_INSTALL_DIR}
    -DMUMPS_UPSTREAM_VERSION=${MUMPS_UPSTREAM_VERSION}
    -DMUMPS_parallel=${OOMPH_ENABLE_MPI}
    -DMUMPS_scalapack=ON
    -DMUMPS_find_SCALAPACK=OFF
    # -Dfind_static=OFF # NOTE: using ON doesn't work on macOS!
    -DMUMPS_find_static=OFF # NOTE: using ON doesn't work on macOS!
    -Dgemmt=ON
    -DMUMPS_intsize64=OFF
    -DMUMPS_scotch=OFF
    -DMUMPS_parmetis=OFF
    -DMUMPS_metis=OFF
    -DMUMPS_openmp=OFF
    -DMUMPS_matlab=OFF
    -DSCALAPACK_BUILD_TESTING=${OOMPH_ENABLE_THIRD_PARTY_LIBRARY_TESTS}
    -DMUMPS_BUILD_TESTING=${OOMPH_ENABLE_THIRD_PARTY_LIBRARY_TESTS}
    -DBUILD_SINGLE=ON
    -DBUILD_DOUBLE=ON
    -DBUILD_COMPLEX=OFF
    -DBUILD_COMPLEX16=OFF
    -DLAPACK_VENDOR=OpenBLAS
    -DLAPACK_ROOT=${OpenBLAS_ROOT})

if(DEFINED BUILD_SHARED_LIBS)
  list(APPEND MUMPS_CMAKE_CONFIGURE_ARGS
       -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS})
else()
  list(APPEND MUMPS_CMAKE_CONFIGURE_ARGS -DBUILD_SHARED_LIBS=OFF)
endif()

# If we're build as a shared lib, we need to make sure transitive dependencies
# get propagated to the runtime linker. Here we ensure the rpath argument is
# used (with --disable-new-dtags) instead of runpath to ensure this
if(DEFINED BUILD_SHARED_LIBS AND BUILD_SHARED_LIBS)
  set(_OOMPH_MUMPS_RPATHS "${MUMPS_INSTALL_DIR}/lib" "${OpenBLAS_ROOT}/lib")
  if(_OOMPH_MUMPS_RPATHS)
    list(
      APPEND
      MUMPS_CMAKE_CONFIGURE_ARGS
      -DCMAKE_BUILD_RPATH=${_OOMPH_MUMPS_RPATHS}
      -DCMAKE_INSTALL_RPATH=${_OOMPH_MUMPS_RPATHS}
      -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON)
    if(UNIX AND NOT APPLE)
      list(APPEND MUMPS_CMAKE_CONFIGURE_ARGS
           -DCMAKE_SHARED_LINKER_FLAGS=-Wl,--disable-new-dtags
           -DCMAKE_EXE_LINKER_FLAGS=-Wl,--disable-new-dtags)
    endif()
  endif()
endif()

if(NOT OOMPH_MUMPS_TARBALL_PATH STREQUAL "")
  if(EXISTS "${OOMPH_MUMPS_TARBALL_PATH}")
    list(APPEND MUMPS_CMAKE_CONFIGURE_ARGS
         "-DMUMPS_url=${OOMPH_MUMPS_TARBALL_PATH}")
  else()
    message(
      FATAL_ERROR
        "OOMPH_MUMPS_TARBALL_PATH was set but the file does not exist:\n"
        "  ${OOMPH_MUMPS_TARBALL_PATH}")
  endif()
endif()

# Define how to configure/build/install the project
oomph_get_external_project_helper(
  PROJECT_NAME mumps
  GIT_REPOSITORY ${MUMPS_GIT_URL}
  GIT_TAG ${MUMPS_GIT_TAG}
  INSTALL_DIR ${MUMPS_INSTALL_DIR}
  CONFIGURE_COMMAND ${CMAKE_COMMAND} ${MUMPS_CMAKE_CONFIGURE_ARGS}
                    -G=${CMAKE_GENERATOR} -B=build
  BUILD_COMMAND ${CMAKE_COMMAND} --build build -j ${OOMPH_NUM_JOBS}
  INSTALL_COMMAND ${CMAKE_COMMAND} --install build
  TEST_COMMAND ${CMAKE_CTEST_COMMAND} --test-dir build -j ${OOMPH_NUM_JOBS}
               --output-on-failure)

# If we're building OpenBLAS, make sure we build it before we get around to
# building MUMPS
if(TARGET openblas)
  add_dependencies(mumps openblas)
endif()

# ---------------------------------------------------------------------------------
