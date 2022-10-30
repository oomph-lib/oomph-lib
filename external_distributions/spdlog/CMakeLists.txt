# ==============================================================================
# Build the spdlog library
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered spdlog subdirectory")

# Spdlog build option(s)
set(SPDLOG_BUILD_SHARED ${BUILD_SHARED_LIBS})
set(SPDLOG_INSTALL ON)

# Download and populate the spdlog target
FetchContent_Declare(
  spdlog_project
  GIT_REPOSITORY https://github.com/puneetmatharu/spdlog.git
  GIT_TAG v1.10.1
  GIT_PROGRESS TRUE
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(spdlog_project)

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and make the updated
# variable values visible to the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_SPDLOG)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_SPDLOG TRUE CACHE INTERNAL "")

message(STATUS "Leaving spdlog subdirectory")
# ------------------------------------------------------------------------------
