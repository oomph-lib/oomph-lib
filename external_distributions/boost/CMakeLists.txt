# ==============================================================================
# Build the Boost library. Be warned, it's a (really) huge library...
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered boost subdirectory")

# The list of Boost libraries we require
set(REQUIRED_COMPONENTS Thread System Program_options)

# Look for an already-installed version of Boost. If OOMPH_USE_BOOST_FROM is
# defined then we'll look under there
if(OOMPH_USE_BOOST_FROM)
  find_package(
    Boost 1.66.0
    COMPONENTS ${REQUIRED_COMPONENTS}
    REQUIRED
    PATHS "${OOMPH_USE_BOOST_FROM}"
    NO_DEFAULT_PATH)
endif()
find_package(Boost 1.66.0 QUIET COMPONENTS ${REQUIRED_COMPONENTS})

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
