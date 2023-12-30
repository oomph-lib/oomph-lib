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
include_guard()

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
    # NOTE: CMake provides FindMPI.cmake to find MPI via find_package(MPI).
    # However, this requires find_package(...) to operate in module mode rather
    # than config mode (where it would search for MPIConfig.cmake). As a result,
    # we can't use a simpler call of the form:
    # ~~~
    #       find_package(MPI REQUIRED
    #                    COMPONENTS C CXX Fortran
    #                    PATHS ${OOMPH_USE_CGAL_FROM} NO_DEFAULT_PATH)
    # ~~~
    # to locate the MPI package. :(
    set(BACKUP_CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}")
    set(CMAKE_PREFIX_PATH "${OOMPH_USE_MPI_FROM}" CACHE INTERNAL "" FORCE)
    find_package(MPI REQUIRED COMPONENTS C CXX Fortran)
    set(CMAKE_PREFIX_PATH "${BACKUP_CMAKE_PREFIX_PATH}" CACHE INTERNAL "" FORCE)
  endif()
  find_package(MPI REQUIRED COMPONENTS C CXX Fortran)
  oomph_check_mpi()
endif()

# Define a cache variable that can be overriden by the user from the
# command-line
set(OOMPH_MPI_NUM_PROC 2 CACHE INTERNAL
    "Number of processes to use with MPI-enabled tests")

# Sanity check
if(NOT OOMPH_MPI_NUM_PROC MATCHES "^[0-9]+$")
  message(
    FATAL_ERROR
      "The flag 'OOMPH_MPI_NUM_PROC' must be set an integer, not ${OOMPH_MPI_NUM_PROC}!"
  )
endif()

# FIXME: Override temporarily
set(MPIEXEC_NUMPROC_FLAG "-np")

# Set the command used to run MPI-enabled self-tests
if(NOT DEFINED OOMPH_MPI_RUN_COMMAND)
  set(OOMPH_MPI_RUN_COMMAND
      "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${OOMPH_MPI_NUM_PROC}")
  message(STATUS "oomph-lib MPI run command: '${OOMPH_MPI_RUN_COMMAND}'")
endif()

# Set the more complex command used to run MPI-enabled self-tests with a
# variable number of processes. The user can sed replace 'OOMPHNP' with the
# number of processes they wish to use
if(NOT DEFINED OOMPH_MPI_VARIABLENP_RUN_COMMAND)
  set(OOMPH_MPI_VARIABLENP_RUN_COMMAND
      "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} OOMPHNP ")
  message(
    STATUS
      "oomph-lib MPI run command (variable NP): '${OOMPH_MPI_VARIABLENP_RUN_COMMAND}'"
  )
endif()

# Add a preprocessor definition and a CMake cache variable to indicate that MPI
# is available and works (if oomph_check_mpi() ran)
oomph_add_c_compile_definitions(OOMPH_HAS_MPI)
oomph_add_cxx_compile_definitions(OOMPH_HAS_MPI)
set(OOMPH_HAS_MPI TRUE CACHE INTERNAL "")
# ------------------------------------------------------------------------------
