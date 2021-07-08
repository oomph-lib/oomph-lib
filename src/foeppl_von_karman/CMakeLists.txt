# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered foeppl_von_karman subdirectory")

# Define the headers
set(HEADERS
    foeppl_von_karman_elements.h
    Tfoeppl_von_karman_elements.h
    foeppl_von_karman_volume_constraint_element.h
    displacement_based_foeppl_von_karman_elements.h
    Tdisplacement_based_foeppl_von_karman_elements.h)

# Define the sources
set(SOURCES foeppl_von_karman_elements.cc Tfoeppl_von_karman_elements.cc
    displacement_based_foeppl_von_karman_elements.cc
    Tdisplacement_based_foeppl_von_karman_elements.cc)

# Set the name of the library
set(LIBNAME foeppl_von_karman)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving foeppl_von_karman subdirectory")
# ------------------------------------------------------------------------------
