# ==============================================================================
# Build the nlohmann JSON library
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered Backward subdirectory")

# Backward build option; enable testing (initially!)
set(BACKWARD_TESTS ON)

# Download and populate the Backward target
FetchContent_Declare(
  backward_project
  GIT_REPOSITORY https://github.com/bombela/backward-cpp.git
  GIT_SHALLOW TRUE
  GIT_TAG v1.6)
FetchContent_MakeAvailable(backward_project)

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and make the updated
# variable values visible to the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_BACKWARD)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_BACKWARD TRUE CACHE INTERNAL "")

message(STATUS "Leaving Backward subdirectory")
# ------------------------------------------------------------------------------
