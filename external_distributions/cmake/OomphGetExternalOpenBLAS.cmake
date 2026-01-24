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
  message(
    FATAL_ERROR
    "Building OpenBLAS is not currently supported on macOS. Instead, you need to "
    "install it with a package manager, e.g. 'brew install openblas' then set "
    "-DOOMPH_USE_OPENBLAS_FROM=$(brew --prefix openblas) during the project configuration "
    "step. Adding OpenBLAS installation support on macOS is a work in progress."
  )

  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # We can make OpenBLAS install with the CMake build system on macOS (and I
  # think(!)) on Linux, but the issue is with MUMPS/ScaLAPACK correctly locating
  # the lib/include dirs afterwards. It's a bit of mess.
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  # Relies on OpenBLAS 0.3.29:

  # set(OPENBLAS_CMAKE_ARGS
  #   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  #   -DCMAKE_INSTALL_PREFIX=${OPENBLAS_INSTALL_DIR}
  #   -DBUILD_SHARED_LIBS=OFF
  #   -DBUILD_WITHOUT_LAPACK=OFF
  #   -DBUILD_WITHOUT_LAPACKE=OFF
  #   -DBUILD_WITHOUT_CBLAS=OFF
  #   -DBUILD_LAPACK_DEPRECATED=ON
  #   -DBUILD_TESTING=${OOMPH_ENABLE_THIRD_PARTY_LIBRARY_TESTS}
  #   -DBUILD_BENCHMARKS=OFF
  #   -DC_LAPACK=OFF
  #   -DDYNAMIC_ARCH=ON
  # )

  # # Define how to configure/build/install the project
  # oomph_get_external_project_helper(
  #   PROJECT_NAME openblas
  #   GIT_REPOSITORY ${OPENBLAS_GIT_URL}
  #   GIT_TAG ${OPENBLAS_GIT_TAG}
  #   INSTALL_DIR ${OPENBLAS_INSTALL_DIR}
  #   CONFIGURE_COMMAND ${CMAKE_COMMAND} ${OPENBLAS_CMAKE_ARGS} -G "${CMAKE_GENERATOR}" -B build
  #   BUILD_COMMAND ${CMAKE_COMMAND} --build build -j ${OOMPH_NUM_JOBS}
  #   INSTALL_COMMAND ${CMAKE_COMMAND} --install build)

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

if(BUILD_SHARED_LIBS)
  set(OpenBLAS_LIBNAME "${CMAKE_SHARED_LIBRARY_PREFIX}openblas${CMAKE_SHARED_LIBRARY_SUFFIX}")
else()
  set(OpenBLAS_LIBNAME "${CMAKE_STATIC_LIBRARY_PREFIX}openblas${CMAKE_STATIC_LIBRARY_SUFFIX}")
endif()
set(OpenBLAS_LIBRARIES "${OPENBLAS_INSTALL_DIR}/lib/${OpenBLAS_LIBNAME}" CACHE PATH "" FORCE)

# ---------------------------------------------------------------------------------
# cmake-format: on
