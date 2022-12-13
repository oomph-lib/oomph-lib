# ==============================================================================
# Build the Boost library. Be warned, it's a (really) huge library...
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered boost subdirectory")

# Look for an already-installed version of Boost
find_package(Boost 1.80.0 EXACT COMPONENTS Thread System Program_options)

# If we couldn't find the required version, lets build it ourselves
if(Boost_FOUND)
  message(
    STATUS
      "Found an existing installation of Boost so we're not going to build it ourselves."
  )
else()
  set(BOOST_INCLUDE_LIBRARIES thread system program_options)
  set(BOOST_ENABLE_CMAKE ON)

  # Download the Boost library (v1.80 only needs C++14)
  include(FetchContent)
  FetchContent_Declare(
    boost_project
    GIT_REPOSITORY https://github.com/boostorg/boost.git
    GIT_TAG boost-1.80.0
    GIT_PROGRESS TRUE
    GIT_SHALLOW TRUE)
  FetchContent_MakeAvailable(boost_project)
endif()

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and make the updated
# variable values visible to the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_BOOST)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_BOOST TRUE CACHE INTERNAL "")

message(STATUS "Leaving boost subdirectory")
# ------------------------------------------------------------------------------
