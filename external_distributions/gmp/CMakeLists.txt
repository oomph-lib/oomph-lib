# ==============================================================================
# Install the GMP library required by CGAL for multiprecision numbers.
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered gmp subdirectory")

# Define the name of the library
set(LIBNAME gmp)
set(ALREADY_HAVE_BUILT_GMP TRUE)

# Download the nlohmann JSON library Git repository.
FetchContent_Declare(gmp_project
                     URL "${CMAKE_CURRENT_LIST_DIR}/gmp-6.1.2.tar.xz")
FetchContent_MakeAvailable(gmp_project)

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and make the updated
# variable values visible to the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_GMP)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_GMP TRUE CACHE INTERNAL "")

message(STATUS "Leaving gmp subdirectory")
# ------------------------------------------------------------------------------
