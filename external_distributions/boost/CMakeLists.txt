# ==============================================================================
# Build the Boost library. Be warned, it's a (really) huge library...
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered boost subdirectory")

# Download the Boost library
FetchContent_Declare(
  boost_project
  GIT_REPOSITORY https://github.com/boostorg/boost.git
  GIT_TAG boost-1.65.1
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(boost_project)

# Add the OOMPH_HAS_BOOST symbol to the compiler flags to define and make the
# compiler flag(s) visible to the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_BOOST)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Add the nlohmann JSON library to the list of external distribution libraries
# and make the updated variable value visible to the parent scope
list(APPEND EXTERNAL_DIST_LIBS Boost::Boost)
set(EXTERNAL_DIST_LIBS ${EXTERNAL_DIST_LIBS} PARENT_SCOPE)

message(STATUS "Leaving boost subdirectory")
# ------------------------------------------------------------------------------
