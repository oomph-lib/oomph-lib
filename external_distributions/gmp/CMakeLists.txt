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

# Add the OOMPH_HAS_GMP symbol to the compiler flags to define and make the
# compiler flag(s) visible to the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_GMP)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Add the gmp library to the list of external distribution libraries and make
# the updated variable value visible to the parent scope
list(APPEND EXTERNAL_DIST_LIBS GMP::GMP)
set(EXTERNAL_DIST_LIBS ${EXTERNAL_DIST_LIBS} PARENT_SCOPE)

message(STATUS "Leaving gmp subdirectory")
# ------------------------------------------------------------------------------
