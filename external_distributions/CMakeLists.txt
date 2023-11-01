# ==============================================================================
# Build the third-party libraries
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered external_distributions subdirectory")

# CMake helpers to download/build/install external libraries
include(ExternalProject)
include(FetchContent)

# Targets to store the accumulated compiler flags and external libs to link
set(EXTERNAL_DIST_CXX_DEFINITIONS "")

if(OOMPH_USE_GMP_FROM
   OR OOMPH_USE_MPFR_FROM
   OR OOMPH_USE_BOOST_FROM
   OR OOMPH_USE_CGAL_FROM)
  # Sanity check
  if(NOT
     (OOMPH_USE_GMP_FROM
      AND OOMPH_USE_MPFR_FROM
      AND OOMPH_USE_BOOST_FROM
      AND OOMPH_USE_CGAL_FROM))
    message(
      FATAL_ERROR
        "If you want CGAL, you need to specify 'OOMPH_USE_GMP_FROM', 'OOMPH_USE_MPFR_FROM', 'OOMPH_USE_BOOST_FROM' AND 'OOMPH_USE_CGAL_FROM'!"
    )
  endif()

  add_subdirectory(gmp)
  add_subdirectory(mpfr)
  add_subdirectory(boost)
  add_subdirectory(cgal)
endif()

if(OOMPH_USE_HYPRE_FROM)
  message(FATAL_ERROR "Broken")
  add_subdirectory(hypre)
endif()

if(OOMPH_USE_MUMPS_FROM)
  message(FATAL_ERROR "Broken")
  add_subdirectory(mumps_and_scalapack)
endif()

if(OOMPH_USE_NLOHMANN_FROM_JSON)
  message(FATAL_ERROR "Broken")
  add_subdirectory(nlohmann_json)
endif()

if(OOMPH_USE_TRILINOS_FROM)
  message(FATAL_ERROR "Broken")
  add_subdirectory(trilinos)
endif()

# Make the accumulated compiler flags visible in the parent scope
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

message(VERBOSE "Leaving external_distributions subdirectory")
# ------------------------------------------------------------------------------
