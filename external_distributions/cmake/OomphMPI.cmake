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
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES MPI::MPI_CXX)

  check_cxx_source_runs(
    "#include <iostream>
    #include <mpi.h>
    int main(int argc, char **argv)
    {
        int myid, numprocs;
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        MPI_Finalize();
        if (numprocs > 2)
        {
            std::cerr << \"Too many cores! MPI self-tests must be run on two cores.\"
                      << std::endl;
            return 4;
        }
        std::cout << \"This worked\" << std::endl;
        return 0;
    }"
    OOMPH_MPI_CXX_WORKS)

  if(NOT OOMPH_MPI_CXX_WORKS)
    message(FATAL_ERROR "MPI_CXX didn't work!")
  endif()
endfunction()
# ------------------------------------------------------------------------------

# Look for MPI functionality if we haven't found it yet
if(NOT MPI_FOUND)
  if(OOMPH_USE_MPI_FROM)
    set(BACKUP_CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}")
    set(CMAKE_PREFIX_PATH "${OOMPH_USE_MPI_FROM}" CACHE INTERNAL "" FORCE)
    find_package(MPI REQUIRED COMPONENTS C CXX Fortran)
    set(CMAKE_PREFIX_PATH "${BACKUP_CMAKE_PREFIX_PATH}" CACHE INTERNAL "" FORCE)
  endif()
  find_package(MPI REQUIRED COMPONENTS C CXX Fortran)
  oomph_check_mpi()
endif()

# ------------------------------------------------------------------------------
