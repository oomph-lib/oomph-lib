# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered biharmonic subdirectory")

# Define the headers
set(HEADERS biharmonic_elements.h biharmonic_flux_elements.h
    biharmonic_problem.h biharmonic_preconditioner.h)

# Define the sources
set(SOURCES biharmonic_elements.cc biharmonic_flux_elements.cc
    biharmonic_problem.cc biharmonic_preconditioner.cc)

# Set the name of the library
set(LIBNAME biharmonic)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving biharmonic subdirectory")
# ------------------------------------------------------------------------------
