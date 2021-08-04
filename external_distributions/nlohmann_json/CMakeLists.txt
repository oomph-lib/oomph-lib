# ==============================================================================
# Build the nlohmann JSON library
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered nlohmann_json subdirectory")

# Disable the build of tests
set(JSON_BuildTests OFF)

# Download and populate the nlohmann JSON library
fetchcontent_declare(
  nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_SHALLOW TRUE
  GIT_TAG v3.9.1)
fetchcontent_makeavailable(nlohmann_json)

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and list of libraries
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_NLOHMANN_JSON)
list(APPEND EXTERNAL_DIST_LIBS nlohmann_json::nlohmann_json)

# Make the updated variable values visible to the parent scope
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)
set(EXTERNAL_DIST_LIBS ${EXTERNAL_DIST_LIBS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_NLOHMANN_JSON TRUE CACHE INTERNAL "")

message(STATUS "Leaving nlohmann_json subdirectory")
# ------------------------------------------------------------------------------
