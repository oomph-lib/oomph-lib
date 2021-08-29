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

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and list of libraries
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_HYPRE)
list(APPEND EXTERNAL_DIST_LIBS HYPRE::HYPRE)

# Make the updated variable values visible to the parent scope
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)
set(EXTERNAL_DIST_LIBS ${EXTERNAL_DIST_LIBS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_HYPRE TRUE CACHE INTERNAL "")

message(STATUS "Leaving hypre subdirectory")
# ------------------------------------------------------------------------------
