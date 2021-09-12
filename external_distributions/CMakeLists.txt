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
set(EXTERNAL_DIST_LIBS)

if(OOMPH_WANT_BACKWARD)
  add_subdirectory(backward)
endif()

if(OOMPH_WANT_CGAL AND (NOT DONT_BUILD_CGAL))
  # FIXME: if(NOT DONT_BUILD_GMP) add_subdirectory(gmp) endif() if(NOT
  # DONT_BUILD_MPFR) add_subdirectory(mpfr) endif()
  if(NOT DONT_BUILD_BOOST)
    add_subdirectory(boost)
  endif()
  add_subdirectory(cgal)
endif()

if(OOMPH_WANT_HYPRE)
  add_subdirectory(hypre)
endif()

if(OOMPH_WANT_MUMPS)
  add_subdirectory(mumps_and_scalapack)
endif()

if(OOMPH_WANT_NLOHMANN_JSON)
  add_subdirectory(nlohmann_json)
endif()

if(OOMPH_WANT_SPDLOG)
  add_subdirectory(spdlog)
endif()

if(OOMPH_WANT_TRILINOS)
  add_subdirectory(trilinos)
endif()

# Make the accumulated compiler flags visible in the parent scope
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Make the accumulated external libraries target visible in the parent scope
set(EXTERNAL_DIST_LIBS ${EXTERNAL_DIST_LIBS} PARENT_SCOPE)

message(VERBOSE "Leaving external_distributions subdirectory")
# ------------------------------------------------------------------------------
