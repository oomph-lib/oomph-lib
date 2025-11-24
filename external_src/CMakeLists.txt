# ==============================================================================
# Handle the build of external sources, i.e. sources required by the project
# that are not developed as part of the library itself.
# ==============================================================================
# Helpful info. for the user
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered external_src subdirectory")

# Initialise the accumulated compiler flags and the external source libraries;
# these variables will/should be updated from each subdirectory added via an
# add_subdirectory(...) call
set(EXTERNAL_SRC_CXX_DEFINITIONS "")
set(EXTERNAL_SRC_LIBS)

# ZLib for gzip decompression
add_subdirectory(oomph_zlib)

# Add HSL and C.R. Bond Bessel functions
add_subdirectory(oomph_hsl)
add_subdirectory(oomph_crbond_bessel)

# Add Triangle
if(NOT OOMPH_SUPPRESS_TRIANGLE_LIB)
  add_subdirectory(oomph_triangle)
endif()

# Add Tetgen
if(NOT OOMPH_SUPPRESS_TETGEN_LIB)
  # TODO: Move to the latest tetgen
  if(OOMPH_USE_NEW_TETGEN)
    # NOTE: Builds fine but breaks the meshing/adaptive_tet_meshes/ demo
    # drivers... Needs further investigation
    add_subdirectory(oomph_tetgen_1.6)
  else()
    add_subdirectory(oomph_tetgen)
  endif()
endif()

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Make the compiler flags and external libraries visible in the parent scope
set(EXTERNAL_SRC_CXX_DEFINITIONS ${EXTERNAL_SRC_CXX_DEFINITIONS} PARENT_SCOPE)
set(EXTERNAL_SRC_LIBS ${EXTERNAL_SRC_LIBS} PARENT_SCOPE)

# Add suppressed Triangle C++ definition
if(OOMPH_SUPPRESS_TRIANGLE_LIB)
  list(APPEND EXTERNAL_SRC_CXX_DEFINITIONS SUPPRESS_TRIANGLE_LIB)
endif()

# Add suppressed Tetgen C++ definition
if(OOMPH_SUPPRESS_TETGEN_LIB)
  list(APPEND EXTERNAL_SRC_CXX_DEFINITIONS SUPPRESS_TETGEN_LIB)
endif()

message(VERBOSE "Leaving external_src subdirectory")
# ------------------------------------------------------------------------------
