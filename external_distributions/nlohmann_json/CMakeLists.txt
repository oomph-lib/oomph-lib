# ==============================================================================
# Build the nlohmann JSON library
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered nlohmann_json subdirectory")

# Disable the build of tests, amalgamate the library into a single header,
# disable the CI build targets and install the library, so that we can find it
# via find_package(...)
set(JSON_BuildTests OFF)
set(JSON_MultipleHeaders OFF)
set(JSON_CI OFF)
set(JSON_Install ON)

# Download and populate the nlohmann JSON library
FetchContent_Declare(
  nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_PROGRESS TRUE
  GIT_SHALLOW TRUE
  GIT_TAG v3.11.2)
FetchContent_MakeAvailable(nlohmann_json)

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and make the updated
# variable values visible to the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_NLOHMANN_JSON)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_NLOHMANN_JSON TRUE CACHE INTERNAL "")

message(STATUS "Leaving nlohmann_json subdirectory")
# ------------------------------------------------------------------------------
