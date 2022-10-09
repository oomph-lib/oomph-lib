# ==============================================================================
# Handle the build of external sources, i.e. sources required by the project
# that are not developed as part of the library itself.
# ==============================================================================
# Helpful info. for the user
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered external_src subdirectory")

# cmake-format: off
# TODO:
# 1. Consider using cmake_add_fortran_subdirectory() for the fortran-only libs
# 2. Ask about whether the following would create any issues:
# # Build malloc.h in parmetis directory if needed (do it from here and before
# # descending into subdir because it may be build after library itself (in which
# # case the compilation will obviously fail.
# if(NOT ALREADY_HAVE_METIS)
#   # all-am: cd oomph_parmetis_3.1.1; make dummy; pwd; ls -l malloc*; cd ..
# endif()
# cmake-format: on

# Initialise the accumulated compiler flags and the external source libraries;
# these variables will/should be updated from each subdirectory added via an
# add_subdirectory(...) call
set(EXTERNAL_SRC_CXX_DEFINITIONS "")
set(EXTERNAL_SRC_LIBS)

# Add HSL, C.R. Bond Bessel functions, and ARPACK
add_subdirectory(oomph_hsl)
add_subdirectory(oomph_crbond_bessel)
add_subdirectory(oomph_arpack)

# Add Triangle
if(NOT SUPPRESS_TRIANGLE_LIB)
  add_subdirectory(oomph_triangle)
endif()

# Add Tetgen
if(NOT SUPPRESS_TETGEN_LIB)
  add_subdirectory(oomph_tetgen)
endif()

# Add BLAS
if(NOT OOMPH_ENABLE_USE_OPENBLAS)
  if(NOT ALREADY_HAVE_BLAS)
    add_subdirectory(oomph_blas)
  endif()

  # Add LAPACK/FORTRAN LAPACK
  if(NOT ALREADY_HAVE_LAPACK)
    add_subdirectory(oomph_lapack)
    add_subdirectory(oomph_flapack)
  endif()
else()

  # Need to use the 'develop' branch of: https://github.com/xianyi/OpenBLAS OR
  # USE: 'brew install openblas' & 'sudo apt-get install libopenblas-dev'
  #
  # add_subdirectory(oomph_openblas_0.3.17)

endif()

# Add SuperLU
if(NOT ALREADY_HAVE_SUPERLU)
  if(NOT OOMPH_USE_DEPRECATED_SUPERLU)
    add_subdirectory(oomph_superlu_5.2.2)
  else()
    add_subdirectory(oomph_superlu_4.3)
  endif()
endif()

# Provide up-to-date sources if the user wants a "bleeding edge" oomph-lib
if(OOMPH_TRANSITION_TO_VERSION_3)
  # Add GKlib for METIS and METIS
  if(NOT ALREADY_HAVE_METIS)
    add_subdirectory(oomph_gklib_from_metis_from_parmetis_4.0.3)
    add_subdirectory(oomph_metis_from_parmetis_4.0.3)
  endif()

  # If compiled with MPI
  if(OOMPH_HAS_MPI)
    # Add ParMETIS
    if(NOT ALREADY_HAVE_PARMETIS)
      add_subdirectory(oomph_parmetis_4.0.3)
    endif()

    # Add SuperLU Dist
    if(NOT ALREADY_HAVE_SUPERLU_DIST)
      add_subdirectory(oomph_superlu_dist_7.0.0)
    endif()
  endif()
else()
  # Add METIS
  if(NOT ALREADY_HAVE_METIS)
    add_subdirectory(oomph_metis_from_parmetis_3.1.1)
  endif()

  # If compiled with MPI
  if(OOMPH_HAS_MPI)
    # Add ParMETIS
    if(NOT ALREADY_HAVE_PARMETIS)
      add_subdirectory(oomph_parmetis_3.1.1)
    endif()

    # Add SuperLU Dist
    if(NOT ALREADY_HAVE_SUPERLU_DIST)
      add_subdirectory(oomph_superlu_dist_3.0)
    endif()
  endif()
endif()

# Make the compiler flags and external libraries visible in the parent scope
set(EXTERNAL_SRC_CXX_DEFINITIONS ${EXTERNAL_SRC_CXX_DEFINITIONS} PARENT_SCOPE)
set(EXTERNAL_SRC_LIBS ${EXTERNAL_SRC_LIBS} PARENT_SCOPE)

message(VERBOSE "Leaving external_src subdirectory")
# ------------------------------------------------------------------------------
