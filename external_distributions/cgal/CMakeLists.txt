# ==============================================================================
# Build/install the CGAL library
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered cgal subdirectory")

# HANDLE WHAT HAPPENS IF WE WANT FIND BOOST...
find_package(Boost)

# Specify build options
set(BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
set(CMAKE_VERBOSE_MAKEFILE OFF)
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

# Update external distributions preprocessor definitions and make the updated
# variable values visible to the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_CGAL)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_CGAL TRUE CACHE INTERNAL "")

message(STATUS "Leaving cgal subdirectory")
# ------------------------------------------------------------------------------
