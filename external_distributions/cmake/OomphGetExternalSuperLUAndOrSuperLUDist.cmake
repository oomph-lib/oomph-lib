# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
#
# ...to be filled in...
#
# =============================================================================
include_guard()

# Locate 'make'; required for building METIS
find_program(MAKE_EXECUTABLE NAMES make REQUIRED)

# Where to clone the repos from
set(GKLIB_GIT_URL https://github.com/KarypisLab/GKlib.git)
set(METIS_GIT_URL https://github.com/KarypisLab/METIS.git)
set(PARMETIS_GIT_URL https://github.com/KarypisLab/ParMETIS.git)
set(SUPERLU_GIT_URL https://github.com/xiaoyeli/superlu.git)
set(SUPERLU_DIST_GIT_URL https://github.com/xiaoyeli/superlu_dist.git)

# Set the default installation paths
set(GKLIB_INSTALL_DIR "${OOMPH_THIRD_PARTY_INSTALL_DIR}/gklib")
set(METIS_INSTALL_DIR "${OOMPH_THIRD_PARTY_INSTALL_DIR}/metis")
set(PARMETIS_INSTALL_DIR "${OOMPH_THIRD_PARTY_INSTALL_DIR}/parmetis")
set(SUPERLU_INSTALL_DIR "${OOMPH_THIRD_PARTY_INSTALL_DIR}/superlu")
set(SUPERLU_DIST_INSTALL_DIR "${OOMPH_THIRD_PARTY_INSTALL_DIR}/superlu_dist")

# If we've already been given GKLIB, METIS and/or Boost, we'll use the files from
# there
if(OOMPH_USE_GKLIB_FROM)
  set(GKLIB_INSTALL_DIR "${OOMPH_USE_GKLIB_FROM}")
endif()
if(OOMPH_USE_METIS_FROM)
  set(METIS_INSTALL_DIR "${OOMPH_USE_METIS_FROM}")
endif()
if(OOMPH_USE_PARMETIS_FROM)
  set(PARMETIS_INSTALL_DIR "${OOMPH_USE_PARMETIS_FROM}")
endif()

# -------------
# GKLIB:
# -------------
# Expected library path and include directory
set(GKLIB_LIBNAME ${CMAKE_STATIC_LIBRARY_PREFIX}GKlib${CMAKE_STATIC_LIBRARY_SUFFIX})
set(GKLIB_LIBRARIES ${GKLIB_INSTALL_DIR}/lib/${GKLIB_LIBNAME} CACHE PATH "Path to GKlib libraries")
set(GKLIB_INCLUDE_DIR ${GKLIB_INSTALL_DIR}/include CACHE PATH "Path to GKlib include directory")

# If we need to build GKlib
if(NOT OOMPH_USE_GKLIB_FROM)
  set(CONFIGURE_GKLIB_WITH_NO_X86 1)
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86|X86|amd64|AMD64")
    set(CONFIGURE_GKLIB_WITH_NO_X86 0)
  endif()

  set(TEST_COMMAND)

  ExternalProject_Add(
    gklib
    GIT_REPOSITORY ${GKLIB_GIT_URL}
    GIT_TAG master
    GIT_SHALLOW TRUE
    INSTALL_DIR "${GKLIB_INSTALL_DIR}"
    LOG_DIR "${CMAKE_BINARY_DIR}/logs"
    BUILD_IN_SOURCE TRUE
    LOG_PATCH TRUE
    LOG_UPDATE TRUE
    LOG_DOWNLOAD TRUE
    LOG_CONFIGURE TRUE
    LOG_BUILD TRUE
    LOG_INSTALL TRUE
    LOG_TEST TRUE
    LOG_MERGED_STDOUTERR TRUE
    LOG_OUTPUT_ON_FAILURE TRUE
    UPDATE_DISCONNECTED TRUE
    BUILD_ALWAYS FALSE
    CONFIGURE_HANDLED_BY_BUILD FALSE
    CONFIGURE_COMMAND ${CMAKE_COMMAND} --install-prefix=<INSTALL_DIR> -G=${CMAKE_GENERATOR} -DNO_X86=${CONFIGURE_GKLIB_WITH_NO_X86} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -B=build
    BUILD_COMMAND ${CMAKE_COMMAND} --build build --parallel ${NUM_JOBS}
    INSTALL_COMMAND ${CMAKE_COMMAND} --install build
    TEST_COMMAND ${TEST_COMMAND})
endif()

# ------
# METIS:
# ------
# Expected library path and include directory
# NOTE: This does assume that metis will be built as a static library...
set(METIS_LIBNAME ${CMAKE_STATIC_LIBRARY_PREFIX}metis${CMAKE_STATIC_LIBRARY_SUFFIX})
set(METIS_LIBRARIES ${METIS_INSTALL_DIR}/lib/${METIS_LIBNAME} CACHE PATH "Path to METIS libraries")
set(METIS_INCLUDE_DIR ${METIS_INSTALL_DIR}/include CACHE PATH "Path to METIS include directory")

# If we need to build METIS
if(NOT OOMPH_USE_METIS_FROM)
  ExternalProject_Add(
    metis
    GIT_REPOSITORY ${METIS_GIT_URL}
    GIT_TAG master
    GIT_SHALLOW TRUE
    INSTALL_DIR "${METIS_INSTALL_DIR}"
    LOG_DIR "${CMAKE_BINARY_DIR}/logs"
    BUILD_IN_SOURCE TRUE
    LOG_PATCH TRUE
    LOG_UPDATE TRUE
    LOG_DOWNLOAD TRUE
    LOG_CONFIGURE TRUE
    LOG_BUILD TRUE
    LOG_INSTALL TRUE
    LOG_TEST TRUE
    LOG_MERGED_STDOUTERR TRUE
    LOG_OUTPUT_ON_FAILURE TRUE
    UPDATE_DISCONNECTED TRUE
    BUILD_ALWAYS FALSE
    CONFIGURE_HANDLED_BY_BUILD TRUE
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${MAKE_EXECUTABLE} config cc=${CMAKE_C_COMPILER} prefix=<INSTALL_DIR> gklib_path=${GKLIB_INSTALL_DIR}
    INSTALL_COMMAND ${MAKE_EXECUTABLE} install)
endif()

# --------
# SUPERLU:
# --------
if (OOMPH_BUILD_SUPERLU)
  # Build options
  set(TPL_METIS_INCLUDE_DIRS_FOR_SUPERLU ${METIS_INCLUDE_DIR} ${GKLIB_INCLUDE_DIR})
  set(TPL_METIS_LIBRARIES_FOR_SUPERLU ${METIS_LIBRARIES} ${GKLIB_LIBRARIES})

  # Create a list with an alternate separator e.g. pipe symbol due to the weird way
  # that SuperLU parses args
  string(REPLACE ";" "|" TPL_METIS_INCLUDE_DIRS_FOR_SUPERLU "${TPL_METIS_INCLUDE_DIRS_FOR_SUPERLU}")
  string(REPLACE ";" "|" TPL_METIS_LIBRARIES_FOR_SUPERLU "${TPL_METIS_LIBRARIES_FOR_SUPERLU}")
  string(REPLACE ";" "|" TPL_BLAS_LIBRARIES_FOR_SUPERLU "${OpenBLAS_LIBRARIES}")

  set(SUPERLU_CMAKE_CONFIGURE_ARGS
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DTPL_ENABLE_METISLIB=ON
    -DTPL_METIS_INCLUDE_DIRS=${TPL_METIS_INCLUDE_DIRS_FOR_SUPERLU}
    -DTPL_METIS_LIBRARIES=${TPL_METIS_LIBRARIES_FOR_SUPERLU}
    -DTPL_BLAS_LIBRARIES=${TPL_BLAS_LIBRARIES_FOR_SUPERLU}
  )

  set(TEST_COMMAND)
  if(NOT OOMPH_DISABLE_THIRD_PARTY_LIBRARY_TESTS)
    set(TEST_COMMAND ${CMAKE_CTEST_COMMAND} --test-dir build --output-on-failure)
  endif()

  ExternalProject_Add(
    superlu
    GIT_REPOSITORY ${SUPERLU_GIT_URL}
    GIT_TAG v6.0.1
    GIT_SHALLOW TRUE
    INSTALL_DIR "${SUPERLU_INSTALL_DIR}"
    LOG_DIR "${CMAKE_BINARY_DIR}/logs"
    BUILD_IN_SOURCE TRUE
    LOG_PATCH TRUE
    LOG_UPDATE TRUE
    LOG_DOWNLOAD TRUE
    LOG_CONFIGURE TRUE
    LOG_BUILD TRUE
    LOG_INSTALL TRUE
    LOG_TEST TRUE
    LOG_MERGED_STDOUTERR TRUE
    LOG_OUTPUT_ON_FAILURE TRUE
    UPDATE_DISCONNECTED TRUE
    BUILD_ALWAYS FALSE
    CONFIGURE_HANDLED_BY_BUILD FALSE
    LIST_SEPARATOR | # Use the alternate list separator
    CONFIGURE_COMMAND ${CMAKE_COMMAND} --install-prefix=<INSTALL_DIR> ${SUPERLU_CMAKE_CONFIGURE_ARGS} -G=${CMAKE_GENERATOR} -B=build
    BUILD_COMMAND ${CMAKE_COMMAND} --build build -j ${NUM_JOBS}
    INSTALL_COMMAND ${CMAKE_COMMAND} --install build
    TEST_COMMAND ${TEST_COMMAND}
  )
endif()

# We can only build ParMETIS or SuperLUDist if MPI is enabled
if (OOMPH_ENABLE_MPI)
  if(NOT MPI_C_COMPILER)
    message(FATAL_ERROR "Something went wrong; MPI_C_COMPILER was not populated!")
  endif()

  # ---------
  # PARMETIS:
  # ---------
  # Expected library path and include directory
  set(PARMETIS_LIBNAME ${CMAKE_STATIC_LIBRARY_PREFIX}parmetis${CMAKE_STATIC_LIBRARY_SUFFIX})
  set(PARMETIS_LIBRARIES ${PARMETIS_INSTALL_DIR}/lib/${PARMETIS_LIBNAME} CACHE PATH "Path to ParMETIS libraries")
  set(PARMETIS_INCLUDE_DIR ${PARMETIS_INSTALL_DIR}/include CACHE PATH "Path to ParMETIS include directory")

  # If we need to build ParMETIS
  if(NOT OOMPH_USE_PARMETIS_FROM)
    ExternalProject_Add(
      parmetis
      GIT_REPOSITORY ${PARMETIS_GIT_URL}
      GIT_TAG main
      GIT_SHALLOW TRUE
      INSTALL_DIR "${PARMETIS_INSTALL_DIR}"
      LOG_DIR "${CMAKE_BINARY_DIR}/logs"
      BUILD_IN_SOURCE TRUE
      LOG_PATCH TRUE
      LOG_UPDATE TRUE
      LOG_DOWNLOAD TRUE
      LOG_CONFIGURE TRUE
      LOG_BUILD TRUE
      LOG_INSTALL TRUE
      LOG_TEST TRUE
      LOG_MERGED_STDOUTERR TRUE
      LOG_OUTPUT_ON_FAILURE TRUE
      UPDATE_DISCONNECTED TRUE
      BUILD_ALWAYS FALSE
      CONFIGURE_HANDLED_BY_BUILD TRUE
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ${MAKE_EXECUTABLE} config cc=${MPI_C_COMPILER} prefix=<INSTALL_DIR> gklib_path=${GKLIB_INSTALL_DIR} metis_path=${METIS_INSTALL_DIR}
      INSTALL_COMMAND ${MAKE_EXECUTABLE} install)
  endif()

  # -------------
  # SUPERLU_DIST:
  # -------------
  # Build options
  set(TPL_PARMETIS_INCLUDE_DIRS_FOR_SUPERLU_DIST ${PARMETIS_INCLUDE_DIR} ${METIS_INCLUDE_DIR} ${GKLIB_INCLUDE_DIR})
  set(TPL_PARMETIS_LIBRARIES_FOR_SUPERLU_DIST ${PARMETIS_LIBRARIES} ${METIS_LIBRARIES} ${GKLIB_LIBRARIES})

  # Create a list with an alternate separator e.g. pipe symbol
  string(REPLACE ";" "|" TPL_PARMETIS_INCLUDE_DIRS_FOR_SUPERLU_DIST "${TPL_PARMETIS_INCLUDE_DIRS_FOR_SUPERLU_DIST}")
  string(REPLACE ";" "|" TPL_PARMETIS_LIBRARIES_FOR_SUPERLU_DIST "${TPL_PARMETIS_LIBRARIES_FOR_SUPERLU_DIST}")
  string(REPLACE ";" "|" TPL_BLAS_LIBRARIES_FOR_SUPERLU_DIST "${OpenBLAS_LIBRARIES}")
  string(REPLACE ";" "|" MPIEXEC_PREFLAGS_FOR_SUPERLU_DIST "--oversubscribe")

  # Build options
  set(SUPERLU_DIST_CMAKE_CONFIGURE_ARGS
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -Denable_openmp=OFF
      -DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS_FOR_SUPERLU_DIST}
      -DTPL_BLAS_LIBRARIES=${TPL_BLAS_LIBRARIES_FOR_SUPERLU_DIST}
      -DTPL_LAPACK_LIBRARIES=${TPL_BLAS_LIBRARIES_FOR_SUPERLU_DIST}
      -DTPL_PARMETIS_INCLUDE_DIRS=${TPL_PARMETIS_INCLUDE_DIRS_FOR_SUPERLU_DIST}
      -DTPL_PARMETIS_LIBRARIES=${TPL_PARMETIS_LIBRARIES_FOR_SUPERLU_DIST}
      -DXSDK_ENABLE_Fortran=OFF
      -Denable_examples=OFF
      -Denable_tests=${OOMPH_DISABLE_THIRD_PARTY_LIBRARY_TESTS}
      -Denable_python=OFF)

  if(OOMPH_DISABLE_THIRD_PARTY_LIBRARY_TESTS)
    list(APPEND SUPERLU_DIST_CMAKE_CONFIGURE_ARGS -Denable_tests=OFF)
  else()
    list(APPEND SUPERLU_DIST_CMAKE_CONFIGURE_ARGS -Denable_tests=ON)
  endif()

  set(TEST_COMMAND)
  if(NOT OOMPH_DISABLE_THIRD_PARTY_LIBRARY_TESTS)
    set(TEST_COMMAND ${CMAKE_CTEST_COMMAND} --test-dir build --output-on-failure)
  endif()

  ExternalProject_Add(
    superlu_dist
    GIT_REPOSITORY ${SUPERLU_DIST_GIT_URL}
    GIT_TAG v9.1.0
    GIT_SHALLOW TRUE
    INSTALL_DIR "${SUPERLU_DIST_INSTALL_DIR}"
    LOG_DIR "${CMAKE_BINARY_DIR}/logs"
    BUILD_IN_SOURCE TRUE
    LOG_PATCH TRUE
    LOG_UPDATE TRUE
    LOG_DOWNLOAD TRUE
    LOG_CONFIGURE TRUE
    LOG_BUILD TRUE
    LOG_INSTALL TRUE
    LOG_TEST TRUE
    LOG_MERGED_STDOUTERR TRUE
    LOG_OUTPUT_ON_FAILURE TRUE
    UPDATE_DISCONNECTED TRUE
    BUILD_ALWAYS FALSE
    CONFIGURE_HANDLED_BY_BUILD FALSE
    LIST_SEPARATOR | # Use the alternate list separator
    CONFIGURE_COMMAND ${CMAKE_COMMAND} --install-prefix=<INSTALL_DIR> ${SUPERLU_DIST_CMAKE_CONFIGURE_ARGS} -G=${CMAKE_GENERATOR} -B=build
    BUILD_COMMAND ${CMAKE_COMMAND} --build build -j ${NUM_JOBS}
    INSTALL_COMMAND ${CMAKE_COMMAND} --install build
    TEST_COMMAND ${TEST_COMMAND})
endif()

# -----------------------------------------------------------------------------

# If we're building OpenBLAS, make sure it gets built before SuperLU or SuperLU_DIST
if(TARGET openblas)
  if(TARGET superlu)
    add_dependencies(superlu openblas)
  endif()
  if(TARGET superlu_dist)
    add_dependencies(superlu_dist openblas)
  endif()
endif()

# METIS depends on GKlib. ParMETIS depends on both METIS and GKlib
if((TARGET gklib) AND (TARGET metis))
  add_dependencies(metis gklib)
endif()
if((TARGET gklib) AND (TARGET parmetis))
  add_dependencies(parmetis gklib)
endif()
if((TARGET metis) AND (TARGET parmetis))
  add_dependencies(parmetis metis)
endif()

# SuperLU depends on GKlib and METIS
if(TARGET superlu)
  if (TARGET gklib)
    add_dependencies(superlu gklib)
  endif()
  if (TARGET metis)
    add_dependencies(superlu metis)
  endif()
endif()

# SuperLUDist depends on GKlib, METIS and ParMETIS
if(TARGET superlu_dist)
  if (TARGET gklib)
    add_dependencies(superlu_dist gklib)
  endif()
  if (TARGET metis)
    add_dependencies(superlu_dist metis)
  endif()
  if (TARGET parmetis)
    add_dependencies(superlu_dist parmetis)
  endif()
endif()

# ---------------------------------------------------------------------------------
# cmake-format: on
