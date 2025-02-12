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
    find_package(MPI REQUIRED)
    set(CMAKE_PREFIX_PATH "${BACKUP_CMAKE_PREFIX_PATH}" CACHE INTERNAL "" FORCE)
  endif()
  find_package(MPI REQUIRED)
  oomph_check_mpi()
endif()

# Set MPIEXEC_NUMPROC_FLAG flag if it hasn't been set. The user can override
# this by passing -DMPIEXEC_NUMPROC_FLAG="..." when configuring the project
set(MPIEXEC_NUMPROC_FLAG "-np")

# Inform the user
message(STATUS "MPIEXEC: ${MPIEXEC_EXECUTABLE}")
message(STATUS "MPIEXEC_NUMPROC_FLAG: ${MPIEXEC_NUMPROC_FLAG}")

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
