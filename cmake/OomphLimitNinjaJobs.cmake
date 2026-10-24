# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Limit number of parallel ninja jobs. Does nothing if
# - called without arguments
# - if Ninja isn't used
# - if required number of jobs is less than 1
#
# USAGE:
# ------
#     include(OomphLimitNinjaJobs)
#
#      oomph_limit_ninja_jobs(<number of jobs>)
#
#        or
#
#      oomph_limit_ninja_jobs(${oomph_NJOBS}) # uses what was used the building oomph-lib
#                                             # bails if variable doesn't exist, i.e. if
#                                             # number of ninja jobs wasn't limited.
#
# =============================================================================
include_guard()

function(oomph_limit_ninja_jobs)
    message(VERBOSE "Limiting number of ninja jobs")

    # Safety check
    get_property(_targets GLOBAL PROPERTY TARGETS)
    if (_targets)
       message(FATAL_ERROR
        "oomph_limit_ninja_jobs() must be called before any targets are created")
    endif()


    # No argument → intentional no-op
    if (ARGC EQUAL 0)
        message(STATUS "oomph_limit_ninja_jobs(): no job count provided, skipping")
        return()
    endif()
    set(NJOBS "${ARGV0}")

    # Empty or zero → intentional no-op
    if (NOT NJOBS OR NJOBS LESS 1)
        message(STATUS
            "oomph_limit_ninja_jobs(): invalid job count ('${NJOBS}'), skipping")
        return()
    endif()

    # Only do it if we're actually using ninja
    if (NOT CMAKE_GENERATOR STREQUAL "Ninja")
        return()
    endif()

    # Check if the number of parallel jobs has already been set in the environment;
    # if so, respect this and bail.
    if (DEFINED ENV{CMAKE_BUILD_PARALLEL_LEVEL})
       set(_tmp "$ENV{CMAKE_BUILD_PARALLEL_LEVEL}")
       message(STATUS
         "oomph_limit_ninja_jobs(): CMAKE_BUILD_PARALLEL_LEVEL is set to "
         "(${_tmp}); skipping Ninja job limiting")
       return()
    endif()

    # Define a global job pool (but labeled, to avoid clashes)
    set(_oomph_pool oomph_limit_pool)

    # Define the pool exactly once (defensive against leakage / reuse)
    get_property(_pools GLOBAL PROPERTY JOB_POOLS)
    if (NOT _pools MATCHES "(^|;)${_oomph_pool}=")
      set_property(GLOBAL PROPERTY JOB_POOLS ${_oomph_pool}=${NJOBS})
    endif()

    # Apply to all job types (COMPILE, LINK, CUSTOM) unless
    # limits are already set by user
    foreach(kind COMPILE LINK CUSTOM)
      if (NOT DEFINED CMAKE_JOB_POOL_${kind})
          set(CMAKE_JOB_POOL_${kind} ${_oomph_pool} PARENT_SCOPE)
      endif()
    endforeach()
    
    message(STATUS
    "The call to oomph_limit_ninja_jobs() from the current project's CMakeLists.txt\n"
    "    file limits Ninja's parallelism via a job pool with depth=${NJOBS}, so\n"
    "    the subsequent build will use at most ${NJOBS} parallel jobs. \n"
    "    Note: Passing '-j' to 'cmake --build' or 'ninja' will NOT increase this\n"
    "    but can be used to reduce the number of parallel jobs further.\n"
    "    To override, you can globally set the number of parallel build jobs via the\n"
    "    environment variable CMAKE_BUILD_PARALLEL_LEVEL prior to configuring, i.e. do\n"
    "       CMAKE_BUILD_PARALLEL_LEVEL=<N>\n"
    "    on the command line.\n"
    "    Alternatively, when configuring oomph-lib itself you can disable the limiter\n"
    "    by configuring oomph-lib using \n"
    "       cmake -G Ninja [...] -DOOMPH_DO_NOT_LIMIT_NINJA_PARALLEL_JOBS=ON\n"
    "    In that case Ninja uses its aggressive default of njobs = number of processors + 2\n"
    "    to build the project. Most computers can handle this but some die, so use\n"
    "    with caution."
    )   

endfunction()

message(VERBOSE "Finished limiting number of ninja jobs")
# ------------------------------------------------------------------------------
# cmake-format: on
