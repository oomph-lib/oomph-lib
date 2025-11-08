# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Set up and test MPI functionality.
#
# USAGE:
# ------
#     include(OomphMPI)
#
# =============================================================================
# cmake-format: on

# ------------------------------------------------------------------------------
function(oomph_check_mpi)
  include(CheckCXXSourceRuns)
  set(CMAKE_REQUIRED_LIBRARIES MPI::MPI_CXX)
  check_cxx_source_runs(
    "#include <mpi.h>
    #include <iostream>
    int main(int argc, char** argv) {
        MPI_Init(&argc, &argv);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::cout << \"Hello from process \" << rank << std::endl;
        MPI_Finalize();
        return 0;
    }"
    OOMPH_MPI_CXX_WORKS)

  if(OOMPH_MPI_CXX_WORKS)
    message(STATUS "MPI test program compiled successfully!")
  else()
    message(FATAL_ERROR "MPI test program failed to compile!")
  endif()
endfunction()
# ------------------------------------------------------------------------------

# Look for MPI functionality if we haven't found it yet
if(NOT MPI_FOUND)
  if(OOMPH_USE_MPI_FROM)
    set(BACKUP_CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}")
    set(CMAKE_PREFIX_PATH "${OOMPH_USE_MPI_FROM}" CACHE INTERNAL "" FORCE)
    find_package(MPI REQUIRED)
    set(CMAKE_PREFIX_PATH "${BACKUP_CMAKE_PREFIX_PATH}" CACHE INTERNAL "" FORCE)
  endif()
  find_package(MPI REQUIRED)
  oomph_check_mpi()
endif()

# ------------------------------------------------------------------------------
