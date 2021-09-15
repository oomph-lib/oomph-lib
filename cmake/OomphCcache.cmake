# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Provides ccache support within the CMake build if the user has ccache.
#
# USAGE:
# ------
#  include(OomphCcache)
#
# =============================================================================
# cmake-format: on
include_guard()

# ------------------------------------------------------------------------------
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  message(STATUS "Found ccache: ${CCACHE_PROGRAM}")
  # Set up wrapper scripts
  set(C_LAUNCHER "${CCACHE_PROGRAM}")
  set(CXX_LAUNCHER "${CCACHE_PROGRAM}")
  configure_file(${CMAKE_CURRENT_LIST_DIR}/oomph-launch-c.in oomph-launch-c)
  configure_file(${CMAKE_CURRENT_LIST_DIR}/oomph-launch-cxx.in oomph-launch-cxx)
  execute_process(COMMAND chmod a+rx "${CMAKE_BINARY_DIR}/oomph-launch-c"
                          "${CMAKE_BINARY_DIR}/oomph-launch-cxx")

  if(CMAKE_GENERATOR STREQUAL "Xcode")
    # Set Xcode project attributes to route compilation and linking through our
    # scripts
    set(CMAKE_XCODE_ATTRIBUTE_CC "${CMAKE_BINARY_DIR}/oomph-launch-c")
    set(CMAKE_XCODE_ATTRIBUTE_CXX "${CMAKE_BINARY_DIR}/oomph-launch-cxx")
    set(CMAKE_XCODE_ATTRIBUTE_LD "${CMAKE_BINARY_DIR}/oomph-launch-c")
    set(CMAKE_XCODE_ATTRIBUTE_LDPLUSPLUS "${CMAKE_BINARY_DIR}/oomph-launch-cxx")
  else()
    # Support Unix Makefiles and Ninja
    set(CMAKE_C_COMPILER_LAUNCHER "${CMAKE_BINARY_DIR}/oomph-launch-c")
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CMAKE_BINARY_DIR}/oomph-launch-cxx")
  endif()
else()
  message(WARNING "Requested 'ccache', but I was unable to find it!")
endif()
# ------------------------------------------------------------------------------
