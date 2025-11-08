# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered axisym_displ_based_foeppl_von_karman subdirectory")

# Define the headers
set(HEADERS axisym_displ_based_fvk_elements.h)

# Define the sources
set(SOURCES axisym_displ_based_fvk_elements.cc)

# Set the name of the library
set(LIBNAME axisym_displ_based_foeppl_von_karman)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving axisym_displ_based_foeppl_von_karman subdirectory")
# ------------------------------------------------------------------------------
