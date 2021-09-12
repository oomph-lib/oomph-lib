# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered tests subdirectory")

# Download Doctest
include(OomphGetDoctest)

# The subdirectories to build
set(TEST_SUBDIRS json)

# Add each subdirectory target to the build
foreach(SUBDIR ${TEST_SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving tests subdirectory")
# ------------------------------------------------------------------------------
