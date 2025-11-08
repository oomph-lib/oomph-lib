# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE
        "Entered space_time_unsteady_heat_equal_order_galerkin subdirectory")

# Define the headers
set(HEADERS space_time_unsteady_heat_elements.h
    refineable_space_time_unsteady_heat_elements.h
    space_time_unsteady_heat_block_preconditionable_elements.h)

# Define the sources
set(SOURCES space_time_unsteady_heat_elements.cc
    refineable_space_time_unsteady_heat_elements.cc
    space_time_unsteady_heat_block_preconditionable_elements.cc)

# Set the name of the library
set(LIBNAME space_time_unsteady_heat_equal_order_galerkin)

# We depend on the generic library
set(LINKLIBS generic space_time_block_preconditioner)

# FIXME: Sort this out
#
# Unlike most libraries, the name of this library is not the same as the name of
# this subdirectory. Therefore we need to provide the directory name to the
# user. It's a bit clunky for now...
set(INCLUDE_SUBDIRECTORY galerkin_equal_order)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS}
  INCLUDE_SUBDIRECTORY ${INCLUDE_SUBDIRECTORY})

message(VERBOSE
        "Leaving space_time_unsteady_heat_equal_order_galerkin subdirectory")
# ------------------------------------------------------------------------------
