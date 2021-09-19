# ==============================================================================
# Locate/download the Google benchmarking suite, benchmark.
# ==============================================================================
# We've already found the benchmark library so we don't need to do anything
if(benchmark_FOUND)
  # Set the version of Google benchmark that we want
  set(OOMPH_GOOGLE_BENCHMARK_VERSION 1.6.0 CACHE INTERNAL "")

  # If the user passed a path to their installation of Google benchmark, add it
  # to the list of paths to search when calling find_package() further below
  if(WITH_GOOGLE_BENCHMARK)
    list(APPEND CMAKE_PREFIX_PATH ${WITH_GOOGLE_BENCHMARK})
  endif()

  # Now quietly search for a previously-installed version of the benchmark
  # library
  find_package(benchmark ${OOMPH_GOOGLE_BENCHMARK_VERSION} QUIET)

  # If we couldn't find it then we need to install it ourselves
  if(NOT benchmark_FOUND)
    # Build options
    set(BENCHMARK_ENABLE_TESTING OFF)
    set(BENCHMARK_ENABLE_LTO OFF)
    set(BENCHMARK_INSTALL_DOCS OFF)
    set(BENCHMARK_DOWNLOAD_DEPENDENCIES ON)
    set(CMAKE_BUILD_TYPE "Release")

    # Set the installation location
    set(OOMPH_GOOGLE_BENCHMARK_INSTALL_DIR
        "${CMAKE_CURRENT_BINARY_DIR}/google_benchmark_install" CACHE PATH
        "Path to the Google benchmark installation")

    # This is where the magic happens: define the rule for how to download the
    # library and where to install it
    # cmake-format: off
  FetchContent_Declare(
    google_benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_SHALLOW TRUE
    GIT_TAG v${OOMPH_GOOGLE_BENCHMARK_VERSION}
    INSTALL_DIR "${OOMPH_GOOGLE_BENCHMARK_INSTALL_DIR}"
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>)
  FetchContent_MakeAvailable(google_benchmark)
  # cmake-format: on
  endif()

  # If the benchmark library isn't available now, something went wrong
  if(NOT TARGET benchmark::benchmark)
    message(
      FATAL_ERROR
        "I tried to install the 'benchmark' library but it doesn't seem to have worked!"
    )
  endif()
endif()

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_GOOGLE_BENCHMARK TRUE CACHE INTERNAL "")
# -----------------------------------------------------------------------------
