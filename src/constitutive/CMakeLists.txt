# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered constitutive subdirectory")

# Define the headers
set(HEADERS constitutive_laws.h)

# Define the sources
set(SOURCES constitutive_laws.cc)

# Set the name of the library
set(LIBNAME constitutive)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving constitutive subdirectory")
# ------------------------------------------------------------------------------
