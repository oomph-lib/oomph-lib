# ==============================================================================
# Build the nlohmann JSON library
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered oomph_superlu_dist_7.0.0 subdirectory")

# Disable the build of tests
set(JSON_BuildTests OFF)

# Information about ParMETIS required by SuperLU Dist
set(TPL_ENABLE_PARMETISLIB ON)
get_target_property(TPL_PARMETIS_INCLUDE_DIRS oomph_parmetis_4.0.3
                    INCLUDE_DIRECTORIES)
get_target_property(TPL_PARMETIS_LIBRARIES oomph_parmetis_4.0.3 LINK_LIBRARIES)

# Name of the library we're creating
set(LIBNAME oomph_superlu_dist_7.0.0)

# Download and populate the nlohmann JSON library
include(FetchContent)
FetchContent_Declare(
  ${LIBNAME}
  GIT_REPOSITORY https://github.com/xiaoyeli/superlu_dist.git
  GIT_SHALLOW TRUE
  GIT_TAG v7.0.0)
FetchContent_MakeAvailable(${LIBNAME})

# Create an alias library name
add_library(${PROJECT_NAMESPACE}::${LIBNAME} ALIAS ${LIBNAME})

# Update list of external sources preprocessor definitions and list of libraries
list(APPEND EXTERNAL_SRC_CXX_DEFINITIONS USING_OOMPH_SUPERLU_DIST)
list(APPEND EXTERNAL_SRC_LIBS ${LIBNAME})

# Make the updated variable values visible to the parent scope
set(EXTERNAL_SRC_CXX_DEFINITIONS ${EXTERNAL_SRC_CXX_DEFINITIONS} PARENT_SCOPE)
set(EXTERNAL_SRC_LIBS ${EXTERNAL_SRC_LIBS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(USING_OOMPH_SUPERLU_DIST TRUE CACHE INTERNAL "")

message(STATUS "Leaving oomph_superlu_dist_7.0.0 subdirectory")
# ------------------------------------------------------------------------------
