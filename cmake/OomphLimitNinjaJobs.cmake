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

    # Define a global job pool (but labeled, to avoid clashes)
    set(_oomph_pool oomph_limit_pool)
    set_property(GLOBAL PROPERTY JOB_POOLS ${_oomph_pool}=${NJOBS})

    # Apply to all job types (COMPILE, LINK, CUSTOM) unless
    # limits are already set by user
    foreach(kind COMPILE LINK CUSTOM)
      if (NOT DEFINED CMAKE_JOB_POOL_${kind})
          set(CMAKE_JOB_POOL_${kind} ${_oomph_pool} PARENT_SCOPE)
      endif()
    endforeach()

    message(STATUS "oomph-lib has limited Ninja parallel jobs to ${NJOBS}")
endfunction()

message(VERBOSE "Finished limiting number of ninja jobs")
# ------------------------------------------------------------------------------
# cmake-format: on
