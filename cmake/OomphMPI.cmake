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
  # Requires CMake 3.19: include(CheckSourceCompiles)
  include(CheckCXXSourceRuns)

  # TODO: I have a feeling I should be using the MPI_C library instead. Find out
  # then come back and change the code below accordingly. See TOOD-LIST.md.

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
    message(FATAL_ERROR "MPI_CXX not working.")
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

# Set the command used to run MPI-enabled self-tests
if(NOT DEFINED MPI_RUN_COMMAND)
  set(MPI_RUN_COMMAND
      "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${OOMPH_MPI_NUM_PROC}")
endif()

# Set the more complex command used to run MPI-enabled self-tests with a
# variable number of processes. The user can sed replace 'OOMPHNP' with the
# number of processes they wish to use
if(NOT DEFINED MPI_VARIABLENP_RUN_COMMAND)
  set(MPI_VARIABLENP_RUN_COMMAND
      "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} OOMPHNP ")
endif()

# Add a preprocessor definition and a CMake cache variable to indicate that MPI
# is available and works (if oomph_check_mpi() ran)
add_compile_definitions(OOMPH_HAS_MPI)
set(OOMPH_HAS_MPI TRUE CACHE INTERNAL "")
# ------------------------------------------------------------------------------
