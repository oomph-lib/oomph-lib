# ==============================================================================
# Build/install the CGAL library
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered cgal subdirectory")

# HANDLE WHAT HAPPENS IF WE WANT FIND BOOST...
find_package(Boost)

# Specify build options
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_LIST_DIR}/cgal_default_installation")
set(BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_C_COMPILER ${CMAKE_C_COMPILER})
set(CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER})
set(CMAKE_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
# set(WITH_GMP TRUE) set(GMP_INCLUDE_DIR "${GMP_INCLUDES}") set(GMP_LIBRARIES
# "${GMP_LINK}"))

# Download the CGAL library
FetchContent_Declare(
  CGAL_project
  GIT_REPOSITORY https://github.com/CGAL/cgal.git
  GIT_TAG v5.2.2
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(CGAL_project)

# Add the OOMPH_HAS_CGAL symbol to the compiler flags to define and make the
# compiler flag(s) visible to the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_CGAL)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Add the nlohmann JSON library to the list of external distribution libraries
# and make the updated variable value visible to the parent scope
list(APPEND EXTERNAL_DIST_LIBS CGAL::CGAL)
set(EXTERNAL_DIST_LIBS ${EXTERNAL_DIST_LIBS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_CGAL TRUE CACHE INTERNAL "")

message(STATUS "Leaving cgal subdirectory")
# ------------------------------------------------------------------------------
