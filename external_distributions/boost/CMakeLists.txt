# ==============================================================================
# Build the Boost library. Be warned, it's a (really) huge library...
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered boost subdirectory")

# Download the Boost library
FetchContent_Declare(
  boost_project
  GIT_REPOSITORY https://github.com/boostorg/boost.git
  GIT_TAG boost-1.65.1
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(boost_project)

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and list of libraries
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_BOOST)
list(APPEND EXTERNAL_DIST_LIBS Boost::Boost)

# Make the updated variable values visible to the parent scope
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)
set(EXTERNAL_DIST_LIBS ${EXTERNAL_DIST_LIBS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_BOOST TRUE CACHE INTERNAL "")

message(STATUS "Leaving boost subdirectory")
# ------------------------------------------------------------------------------
