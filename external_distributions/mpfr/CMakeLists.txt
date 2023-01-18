# ==============================================================================
# Install the MPFR library required by CGAL for multiprecision numbers.
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered mpfr subdirectory")

# Look for an already-installed version of MPFR using FindMPFR.cmake. If
# OOMPH_USE_MPFR_FROM is defined then FindMPFR.cmake will look under there for
# the required files.
find_package(MPFR 2.2.1 QUIET)
if(MPFR_FOUND)
  message(STATUS "Found an existing installation of MPFR: ${MPFR_LIBRARIES}")
else()
  message(STATUS "Third-party: creating target 'MPFR::MPFR'")

  # cmake-format: off
  set(PREFIX ${FETCHCONTENT_BASE_DIR}/mpfr)
  set(MPFR_LIBNAME ${CMAKE_STATIC_LIBRARY_PREFIX}mpfr${CMAKE_STATIC_LIBRARY_SUFFIX})
  set(MPFR_LIBRARIES ${CMAKE_INSTALL_PREFIX}/lib/${MPFR_LIBNAME} CACHE PATH "Path to GMP libraries")
  set(MPFR_INCLUDES ${CMAKE_INSTALL_PREFIX}/include CACHE PATH "Path to GMP include directory")
  # cmake-format: on

  # Extract library directory from full path to library
  cmake_path(GET GMP_C_LIBRARIES PARENT_PATH GMP_LIBRARY_DIR)

  include(ProcessorCount)
  ProcessorCount(NUM_CPU)
  ExternalProject_Add(
    mpfr
    URL https://ftp.gnu.org/gnu/mpfr/mpfr-4.1.0.tar.xz
    URL_MD5 bdd3d5efba9c17da8d83a35ec552baef
    PREFIX ${PREFIX}
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    UPDATE_DISCONNECTED true # need this to avoid constant rebuild
    CONFIGURE_HANDLED_BY_BUILD ON # avoid constant reconfigure
    CONFIGURE_COMMAND
      ${PREFIX}/src/mpfr/configure --disable-shared --with-pic
      --prefix=<INSTALL_DIR> --with-gmp-include=${GMP_C_INCLUDES}
      --with-gmp-lib=${GMP_LIBRARY_DIR}
    BUILD_COMMAND make -s -j${NUM_CPU} || make -j${NUM_CPU}
    INSTALL_COMMAND make -s -j${NUM_CPU} install || make -j${NUM_CPU} install
    TEST_COMMAND ""
    BUILD_BYPRODUCTS ${MPFR_LIBRARIES})
  ExternalProject_Get_Property(mpfr SOURCE_DIR)

  # Define CMake target to link to
  set(MPFR_LIBRARIES ${MPFR_LIBRARIES})
  add_library(MPFR::MPFR INTERFACE IMPORTED GLOBAL)
  file(MAKE_DIRECTORY ${MPFR_INCLUDES}) # avoid race condition
  target_include_directories(MPFR::MPFR INTERFACE ${MPFR_INCLUDES})
  target_link_libraries(
    MPFR::MPFR INTERFACE "${MPFR_LIBRARIES}" "${GMP_C_LIBRARIES}"
                         "${GMP_CXX_LIBRARIES}")
  add_dependencies(MPFR::MPFR mpfr)

  # If we're building gmp ourselves, make sure to build it before mpfr
  if(TARGET gmp)
    add_dependencies(mpfr gmp)
  endif()

  # Specify the path to find MPFR after installation
  message(STATUS "MPFR will be installed to: ${CMAKE_INSTALL_PREFIX}")
  set(OOMPH_USE_MPFR_FROM "${CMAKE_INSTALL_PREFIX}" CACHE PATH
      "Path to MPFR installation" FORCE)
endif()

if(NOT TARGET MPFR::MPFR)
  message(FATAL_ERROR "Creation of target 'MPFR::MPFR' failed")
endif()

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and make the updated
# variable values visible to the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_MPFR)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_MPFR TRUE CACHE INTERNAL "")

message(STATUS "Leaving mpfr subdirectory")
# ------------------------------------------------------------------------------
