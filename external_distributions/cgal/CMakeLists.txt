# ==============================================================================
# Build/install the CGAL library
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered cgal subdirectory")

# Find dependencies
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
find_package(Boost REQUIRED COMPONENTS Thread System Program_options)

# Configuration flags
set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE CACHE INTERNAL "")

# Download the CGAL library
FetchContent_Declare(
  cgal_project
  GIT_REPOSITORY https://github.com/CGAL/cgal.git
  GIT_TAG v5.5.1 # releases/CGAL-4.11.3 # v5.5.1
  GIT_PROGRESS TRUE
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(cgal_project)

find_package(CGAL REQUIRED PATHS "${CGAL_BINARY_DIR}")

# FIXME: BUILD_TESTING is disable through CGAL and so whenever we import CGAL,
# testing gets disable. This is a particular pain when we want to run our demo
# drivers. Find a better way to do this...
#
# unset(BUILD_TESTING)

# if(NOT TARGET CGAL::CGAL) add_library(CGAL::CGAL ALIAS CGAL) endif()

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
