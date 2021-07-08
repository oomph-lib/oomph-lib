# ==============================================================================
# Build the Hypre library.
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered hypre subdirectory")

# Download HYPRE repo
include(FetchContent)
FetchContent_Declare(
  hypre_project
  GIT_REPOSITORY https://github.com/hypre-space/hypre.git
  GIT_TAG v2.21.0
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(hypre_project)

# Add OOMPH_HAS_HYPRE to the compiler flags to define the symbol and make the
# accumulated compiler flags visible in the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_HYPRE)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Add the HYPRE libraries to the list of external distribution libraries and
# make the updated variable value visible to the parent scope
list(APPEND EXTERNAL_DIST_LIBS HYPRE::HYPRE)
set(EXTERNAL_DIST_LIBS ${EXTERNAL_DIST_LIBS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_HYPRE TRUE CACHE INTERNAL "")

message(STATUS "Leaving hypre subdirectory")
# ------------------------------------------------------------------------------
