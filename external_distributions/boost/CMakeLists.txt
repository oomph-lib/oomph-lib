# ==============================================================================
# Build the Boost library. Be warned, it's a (really) huge library...
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered boost subdirectory")

# The list of Boost libraries we require
set(REQUIRED_COMPONENTS Thread System Program_options)

if(OOMPH_USE_BOOST_FROM)
  # If the user specified OOMPH_USE_BOOST_FROM, look there
  find_package(
    Boost 1.83.0
    COMPONENTS ${REQUIRED_COMPONENTS}
    REQUIRED
    PATHS "${OOMPH_USE_BOOST_FROM}"
    NO_DEFAULT_PATH)
else()
  # Now search system wide
  find_package(Boost 1.83.0 REQUIRED COMPONENTS ${REQUIRED_COMPONENTS})
endif()

# If we couldn't find the required version, lets build it ourselves
if(Boost_FOUND)
  list(JOIN Boost_LIBRARIES ", " BOOST_COMPONENTS)
  message(
    STATUS
      "Found an existing installation of Boost! Components found: (${BOOST_COMPONENTS})"
  )
else()
  message(FATAL_ERROR "Could not find an existing installation of Boost!")
endif()

# TODO: Change OOMPH_HAS_<PACKAGE> to OOMPH_FOUND_<PACKAGE>...

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
