# ==============================================================================
# Build/install the CGAL library
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered cgal subdirectory")

# HANDLE WHAT HAPPENS IF WE WANT FIND BOOST...
find_package(Boost)

# Specify build options
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_LIST_DIR}/cgal_default_installation")
set(BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_C_COMPILER ${CMAKE_C_COMPILER})
set(CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER})
set(CMAKE_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
# set(OOMPH_WITH_GMP TRUE) set(GMP_INCLUDE_DIR "${GMP_INCLUDES}")
# set(GMP_LIBRARIES
# "${GMP_LINK}"))

# Download the CGAL library
FetchContent_Declare(
  CGAL_project
  GIT_REPOSITORY https://github.com/CGAL/cgal.git
  GIT_TAG v5.2.2
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(CGAL_project)

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and list of libraries
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_CGAL)
list(APPEND EXTERNAL_DIST_LIBS CGAL::CGAL)

# Make the updated variable values visible to the parent scope
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)
set(EXTERNAL_DIST_LIBS ${EXTERNAL_DIST_LIBS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_CGAL TRUE CACHE INTERNAL "")

message(STATUS "Leaving cgal subdirectory")
# ------------------------------------------------------------------------------
