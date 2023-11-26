# ==============================================================================
# Handle the build of external sources, i.e. sources required by the project
# that are not developed as part of the library itself.
# ==============================================================================
# Helpful info. for the user
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered external_src subdirectory")

if(ALREADY_HAVE_BLAS)
  message(FATAL_ERROR "The variable 'ALREADY_HAVE_BLAS' has been deprecated!")
endif()
if(ALREADY_HAVE_LAPACK)
  message(FATAL_ERROR "The variable 'ALREADY_HAVE_LAPACK' has been deprecated!")
endif()
if(ALREADY_HAVE_SUPERLU)
  message(
    FATAL_ERROR "The variable 'ALREADY_HAVE_SUPERLU' has been deprecated!")
endif()
if(ALREADY_HAVE_METIS)
  message(FATAL_ERROR "The variable 'ALREADY_HAVE_METIS' has been deprecated!")
endif()
if(ALREADY_HAVE_PARMETIS)
  message(
    FATAL_ERROR "The variable 'ALREADY_HAVE_PARMETIS' has been deprecated!")
endif()
if(ALREADY_HAVE_SUPERLU_DIST)
  message(
    FATAL_ERROR "The variable 'ALREADY_HAVE_SUPERLU_DIST' has been deprecated!")
endif()

# Initialise the accumulated compiler flags and the external source libraries;
# these variables will/should be updated from each subdirectory added via an
# add_subdirectory(...) call
set(EXTERNAL_SRC_CXX_DEFINITIONS "")
set(EXTERNAL_SRC_LIBS)

# ZLib for gzip decompression
add_subdirectory(oomph_zlib)

# Add HSL, C.R. Bond Bessel functions, and ARPACK
add_subdirectory(oomph_hsl)
add_subdirectory(oomph_crbond_bessel)
add_subdirectory(oomph_arpack)

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

# Add GKlib for METIS and METIS
add_subdirectory(oomph_gklib_from_metis_from_parmetis_4.0.3)
add_subdirectory(oomph_metis_from_parmetis_4.0.3)

# Add SuperLU
add_subdirectory(oomph_superlu_6.0.1)

# If compiled with MPI
if(OOMPH_HAS_MPI)
  # Provide up-to-date sources if the user wants a "bleeding edge" oomph-lib
  add_subdirectory(oomph_parmetis_4.0.3)
  add_subdirectory(oomph_superlu_dist_8.2.1)
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
