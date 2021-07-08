# ==============================================================================
# Triangle installed as library
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oomph_triangle subdirectory")

# Define the sources
set(SOURCES triangle.c)

# Define the headers
set(HEADERS triangle.h)

# If we can't find fpu_control.h, symlink the dummy version
if(NOT OOMPH_HAS_FPUCONTROLH)
  message(VERBOSE "Unable to find fpu_control.h")
  message(VERBOSE "Symlinking the dummy version instead.")

  # Now add the header to the list of headers required to build triangle.c
  list(APPEND HEADERS dummy_fpu_control.h)
  file(CREATE_LINK dummy_fpu_control.h
       "${CMAKE_CURRENT_LIST_DIR}/fpu_control.h" SYMBOLIC)
endif()

# Define the name of the library
set(LIBNAME oomph_triangle)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

# Update list of external sources preprocessor definitions and list of libraries
list(APPEND EXTERNAL_SRC_CXX_DEFINITIONS OOMPH_HAS_TRIANGLE_LIB)
list(APPEND EXTERNAL_SRC_LIBS ${LIBNAME})

# Make the updated variable values visible to the parent scope
set(EXTERNAL_SRC_CXX_DEFINITIONS ${EXTERNAL_SRC_CXX_DEFINITIONS} PARENT_SCOPE)
set(EXTERNAL_SRC_LIBS ${EXTERNAL_SRC_LIBS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_TRIANGLE TRUE CACHE INTERNAL "")

message(VERBOSE "Leaving oomph_triangle subdirectory")
# ------------------------------------------------------------------------------
