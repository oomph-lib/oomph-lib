# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered darcy subdirectory")

# Define the headers
set(HEADERS darcy_elements.h Tdarcy_elements.h darcy_face_elements.h)

# Define the sources
set(SOURCES darcy_elements.cc Tdarcy_elements.cc)

# Set the name of the library
set(LIBNAME darcy)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving darcy subdirectory")
# ------------------------------------------------------------------------------
