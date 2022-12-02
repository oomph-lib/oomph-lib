# ==============================================================================
# Collection of Bessel functions and the like, implemented in C++ by C. Bond
# (http://www.crbond.com/), using the algorithms in  Zhang and Jin's book
# "Computation of Special Functions". John Wiley and Sons, 1996. This is based
# on the version dated 06/04 which appears to have corrected errors in complex
# Bessel functions.
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oomph_crbond_bessel subdirectory")

# Define the sources
set(SOURCES crbond_bessel.cc)

# Define the headers
set(HEADERS crbond_bessel.h)

# Set the name of the library
set(LIBNAME oomph_crbond_bessel)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

# ------------------------------------------------------------------------------
# BOILERPLATE: Update variables in parent scope

# Add the library to the list of external source libraries and make the updated
# variable value visible to the parent scope
list(APPEND EXTERNAL_SRC_LIBS ${LIBNAME})
set(EXTERNAL_SRC_LIBS ${EXTERNAL_SRC_LIBS} PARENT_SCOPE)

message(VERBOSE "Leaving oomph_crbond_bessel subdirectory")
# ------------------------------------------------------------------------------
