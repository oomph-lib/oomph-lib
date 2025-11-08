# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered space_time_block_preconditioner subdirectory")

# Define the headers
set(HEADERS general_purpose_space_time_block_preconditioner.h
    general_purpose_space_time_subsidiary_block_preconditioner.h
    general_purpose_space_time_block_preconditionable_elements.h)

# Define the sources
set(SOURCES general_purpose_space_time_block_preconditioner.cc
    general_purpose_space_time_subsidiary_block_preconditioner.cc)

# Set the name of the library
set(LIBNAME space_time_block_preconditioner)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving space_time_block_preconditioner subdirectory")
# ------------------------------------------------------------------------------
