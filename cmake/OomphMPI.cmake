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
  include(CheckSourceRuns)

  # TODO: I have a feeling I should be using the MPI_C library instead. Find out
  # then come back and change the code below accordingly. See TOOD-LIST.md.

  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES MPI::MPI_CXX)

  check_source_runs(
    CXX
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
# Only look for MPI functionality if we haven't already found it
if(NOT MPI_FOUND)
  find_package(MPI REQUIRED)
  oomph_check_mpi()
endif()

# if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.19.0") else() message( STATUS "I
# can only check if MPI works using the\n\ check_<LANG>_source_compile()\n\
# function which was made available in CMake 3.19. As you don't have a\n\ recent
# enough version of CMake, I'm not going to perform these tests so\n\ you should
# be aware that MPI *may* not work.\n") endif()

# Make MPI libraries available project-wide. Directory-wide assignments are
# typically a bad design decision with modern CMake, but it makes sense to use
# it here to add MPI libraries to the entire build if we want it
include_directories(SYSTEM ${MPI_C_INCLUDE_PATH} ${MPI_CXX_INCLUDE_PATH})
# include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
link_libraries(${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})

# Set the command used to run MPI-enabled self-tests if the user didn't specify
# the commands to use
if(NOT DEFINED MPI_RUN_COMMAND)
  set(MPI_RUN_COMMAND "mpirun -np 2")
endif()
if(NOT DEFINED MPI_VARIABLENP_RUN_COMMAND)
  set(MPI_VARIABLENP_RUN_COMMAND "mpirun -np OOMPHNP ")
endif()

# Add a preprocessor definition and a CMake cache variable to indicate that MPI
# is available and works (if oomph_check_mpi() ran)
add_compile_definitions(OOMPH_HAS_MPI)
set(OOMPH_HAS_MPI TRUE CACHE INTERNAL "")
# ------------------------------------------------------------------------------
