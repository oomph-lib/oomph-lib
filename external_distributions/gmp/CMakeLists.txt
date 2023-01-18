# ==============================================================================
# Install the GMP library required by CGAL for multiprecision numbers.
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered gmp subdirectory")

# Look for an already-installed version of GMP using FindGMP.cmake. If
# OOMPH_USE_GMP_FROM is defined then FindGMP.cmake will look under there for the
# required files.
set(GMP_FOUND FALSE)
if(GMP_FOUND)
  find_package(GMP 4.2 QUIET)
  message(STATUS "Found an existing installation of GMP: ${GMP_C_LIBRARIES}")
else()
  message(STATUS "Third-party: creating target 'GMP::GMP'")

  # cmake-format: off
  set(PREFIX ${FETCHCONTENT_BASE_DIR}/gmp)
  set(GMP_C_LIBNAME ${CMAKE_STATIC_LIBRARY_PREFIX}gmp${CMAKE_STATIC_LIBRARY_SUFFIX})
  set(GMP_CXX_LIBNAME ${CMAKE_STATIC_LIBRARY_PREFIX}gmpxx${CMAKE_STATIC_LIBRARY_SUFFIX})
  set(GMP_C_LIBRARIES ${CMAKE_INSTALL_PREFIX}/lib/${GMP_C_LIBNAME} CACHE PATH "Path to GMP C libraries")
  set(GMP_CXX_LIBRARIES ${CMAKE_INSTALL_PREFIX}/lib/${GMP_CXX_LIBNAME} CACHE PATH "Path to GMP C++ libraries")
  set(GMP_C_INCLUDES ${CMAKE_INSTALL_PREFIX}/include CACHE PATH "Path to GMP C include directory")
  set(GMP_CXX_INCLUDES ${GMP_C_INCLUDES} CACHE PATH "Path to GMP C++ include directory")
  # cmake-format: on

  include(ProcessorCount)
  ProcessorCount(NUM_CPU)
  ExternalProject_Add(
    gmp
    URL https://gmplib.org/download/gmp/gmp-6.2.1.tar.xz
    URL_MD5 0b82665c4a92fd2ade7440c13fcaa42b
    PREFIX ${PREFIX}
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    UPDATE_DISCONNECTED true # need this to avoid constant rebuild
    CONFIGURE_HANDLED_BY_BUILD ON # avoid constant reconfigure
    CONFIGURE_COMMAND ${PREFIX}/src/gmp/configure --disable-shared --enable-cxx
                      --with-pic --prefix=<INSTALL_DIR>
    BUILD_COMMAND make -s -j${NUM_CPU} || make -j${NUM_CPU}
    INSTALL_COMMAND make -s -j${NUM_CPU} install || make -j${NUM_CPU} install
    TEST_COMMAND ""
    BUILD_BYPRODUCTS ${GMP_C_LIBRARIES} ${GMP_CXX_LIBRARIES})
  ExternalProject_Get_Property(gmp SOURCE_DIR)

  # Define CMake target to link to
  add_library(GMP::GMP INTERFACE IMPORTED GLOBAL)
  file(MAKE_DIRECTORY ${GMP_C_INCLUDES}) # avoid race condition
  target_include_directories(GMP::GMP INTERFACE ${GMP_C_INCLUDES}
                                                ${GMP_CXX_INCLUDES})
  target_link_libraries(GMP::GMP INTERFACE ${GMP_C_LIBRARIES}
                                           ${GMP_CXX_LIBRARIES})
  add_dependencies(GMP::GMP gmp)

  # Specify the path to find GMP after installation
  message(STATUS "GMP will be installed to: ${CMAKE_INSTALL_PREFIX}")
  set(OOMPH_USE_GMP_FROM "${CMAKE_INSTALL_PREFIX}" CACHE PATH
      "Path to GMP installation" FORCE)
endif()

if(NOT TARGET GMP::GMP)
  message(FATAL_ERROR "Creation of target 'GMP::GMP' failed")
endif()

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and make the updated
# variable values visible to the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_GMP)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_GMP TRUE CACHE INTERNAL "")

message(STATUS "Leaving gmp subdirectory")
# ------------------------------------------------------------------------------
