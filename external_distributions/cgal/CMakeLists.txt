# ==============================================================================
# Build/install the CGAL library
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered cgal subdirectory")

# Find dependencies
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
find_package(Boost 1.66.0 REQUIRED COMPONENTS Thread System Program_options)

# CGAL currently defines a variable called BUILD_TESTING which they shouldn't do
# and it messes up the use of the demo drivers. For now, we have to clone CGAL
# ourselves so that we can patch it. There is an open GitHub issue for the fix
# here: https://github.com/CGAL/cgal/issues/7121. Once the patch has been merged
# in, we'll be able to allow users to use their own pre-installed version of
# CGAL which will save a lot of time at configure time. To do this, we can just
# remove the use of the FORCE_ variable below.
set(FORCE_BUILD_OWN_VERSION_OF_CGAL TRUE)

# Configuration flags
set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE CACHE INTERNAL "")

# See if the user already has CGAL installed
find_package(CGAL 5.5.1)
if(CGAL_FOUND AND (NOT FORCE_BUILD_OWN_VERSION_OF_CGAL))
  message(
    STATUS
      "Found an existing installation of CGAL so we're not going to build it ourselves."
  )
else()
  # Path to the patch file
  set(PATCH_FILE 0001-Rename-BUILD_TESTING-CGAL_BUILD_TESTING.patch)
  set(FULL_PATCH_FILE_PATH "${CMAKE_CURRENT_LIST_DIR}/patches/${PATCH_FILE}")

  # Download the CGAL library
  FetchContent_Declare(
    cgal_project
    GIT_REPOSITORY https://github.com/CGAL/cgal.git
    GIT_TAG v5.5.1
    GIT_PROGRESS TRUE
    GIT_SHALLOW TRUE
    # If the user configures the project multiple times, the "git apply" will
    # fail on subsequent runs because the patch file will have already been
    # applied. To get around this, we reset the changes each time with "git
    # checkout". Hacky solution but it works for now...
    PATCH_COMMAND
      ${CMAKE_COMMAND} -E copy_if_different ${FULL_PATCH_FILE_PATH}
      <SOURCE_DIR> && git checkout <SOURCE_DIR> && git apply
      <SOURCE_DIR>/${PATCH_FILE})
  FetchContent_MakeAvailable(cgal_project)
endif()

find_package(CGAL REQUIRED PATHS "${CGAL_BINARY_DIR}")

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and make the updated
# variable values visible to the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_CGAL)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_CGAL TRUE CACHE INTERNAL "")

message(STATUS "Leaving cgal subdirectory")
# ------------------------------------------------------------------------------
