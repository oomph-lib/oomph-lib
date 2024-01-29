# cmake-format: off
# ==============================================================================
# Locate optional third-party libraries requested by the user.
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Locating third-party libraries")

# Targets to store the accumulated compiler flags and external libs to link
set(EXTERNAL_DIST_CXX_DEFINITIONS "")

# TODO: Possibly remove this since all of these options are "required" by oomph-lib so they must
# be set by this point
if(NOT (OOMPH_USE_GKLIB_FROM AND OOMPH_USE_METIS_FROM AND OOMPH_USE_SUPERLU_FROM))
  message(
    FATAL_ERROR
      "You need to specify 'OOMPH_USE_GKLIB_FROM', 'OOMPH_USE_METIS_FROM', *AND* 'OOMPH_USE_SUPERLU_FROM'!"
  )
endif()

# SuperLU_DIST needs GKlib, METIS and ParMETIS
if(OOMPH_USE_SUPERLU_DIST_FROM)
  if(NOT (OOMPH_USE_GKLIB_FROM AND OOMPH_USE_METIS_FROM AND OOMPH_USE_PARMETIS_FROM))
    message(
      FATAL_ERROR
        "If you want SuperLU, you need to specify 'OOMPH_USE_GKLIB_FROM', 'OOMPH_USE_METIS_FROM', 'OOMPH_USE_PARMETIS_FROM', *AND* 'OOMPH_USE_SUPERLU_DIST_FROM'!"
    )
  endif()
endif()

# SuperLU/SuperLU_DIST
if(OOMPH_USE_SUPERLU_FROM OR OOMPH_USE_SUPERLU_DIST_FROM)
  if(OOMPH_USE_SUPERLU_DIST_FROM AND (NOT OOMPH_ENABLE_MPI))
    message(FATAL_ERROR "MPI must be enabled to use SuperLU_DIST!")
  endif()

  # NOTE: We wrote a custom Find*.cmake scripts to search for GKlib, METIS and ParMETIS and
  # SuperLU_DIST. They will internally read the variables OOMPH_USE_GKLIB_FROM, OOMPH_USE_METIS_FROM,
  # OOMPH_USE_PARMETIS_FROM and OOMPH_USE_SUPERLU_DIST_FROM, respectively. The CMake code for
  # SuperLU itself has been written properly so we can import it like a normal version CMake library
  find_package(GKlib REQUIRED GLOBAL)
  find_package(METIS REQUIRED GLOBAL)
  find_package(superlu 6.0.1 REQUIRED GLOBAL PATHS ${OOMPH_USE_SUPERLU_FROM} NO_DEFAULT_PATH)
  if(OOMPH_USE_SUPERLU_DIST_FROM)
    find_package(ParMETIS REQUIRED GLOBAL)
    find_package(SuperLU_DIST REQUIRED GLOBAL)
  endif()

  if(TARGET GKlib::GKlib)
    message(STATUS "Hurray! The target GKlib::GKlib is defined!")
  else()
    message(FATAL_ERROR "Target GKlib::GKlib is not defined!")
  endif()
  if(TARGET METIS::METIS)
    message(STATUS "Hurray! The target METIS::METIS is defined!")
  else()
    message(FATAL_ERROR "Target METIS::METIS is not defined!")
  endif()
  if(TARGET superlu::superlu)
    message(STATUS "Hurray! The target superlu::superlu is defined!")
  else()
    message(FATAL_ERROR "Target superlu::superlu is not defined!")
  endif()
  if(OOMPH_USE_SUPERLU_DIST_FROM)
    if(TARGET ParMETIS::ParMETIS)
      message(STATUS "Hurray! The target ParMETIS::ParMETIS is defined!")
    else()
      message(FATAL_ERROR "Target ParMETIS::ParMETIS is not defined!")
    endif()
    if(TARGET SuperLU_DIST::SuperLU_DIST)
      message(STATUS "Hurray! The target SuperLU_DIST::SuperLU_DIST is defined!")
    else()
      message(FATAL_ERROR "Target SuperLU_DIST::SuperLU_DIST is not defined!")
    endif()
  endif()

  list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_GKLIB OOMPH_HAS_METIS OOMPH_HAS_SUPERLU)
  if(OOMPH_USE_SUPERLU_DIST_FROM)
    list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_PARMETIS OOMPH_HAS_SUPERLU_DIST)
  endif()

  set(OOMPH_HAS_GKLIB TRUE CACHE INTERNAL "")
  set(OOMPH_HAS_METIS TRUE CACHE INTERNAL "")
  set(OOMPH_HAS_SUPERLU TRUE CACHE INTERNAL "")
  if(OOMPH_USE_SUPERLU_DIST_FROM)
    set(OOMPH_HAS_PARMETIS TRUE CACHE INTERNAL "")
    set(OOMPH_HAS_SUPERLU_DIST TRUE CACHE INTERNAL "")
  endif()
endif()

# CGAL
if(OOMPH_USE_GMP_FROM OR OOMPH_USE_MPFR_FROM OR OOMPH_USE_BOOST_FROM OR OOMPH_USE_CGAL_FROM)
  if(NOT (OOMPH_USE_GMP_FROM AND OOMPH_USE_MPFR_FROM AND OOMPH_USE_BOOST_FROM AND OOMPH_USE_CGAL_FROM))
    message(
      FATAL_ERROR
        "If you want CGAL, you need to specify 'OOMPH_USE_GMP_FROM', 'OOMPH_USE_MPFR_FROM', 'OOMPH_USE_BOOST_FROM' *AND* 'OOMPH_USE_CGAL_FROM'!"
    )
  endif()

  set(REQUIRED_BOOST_COMPONENTS thread system program_options)

  # NOTE: We wrote a custom FindGMP.cmake and FindMPFR.cmake script to search
  # for GMP and MPFR. They will internally read the variables OOMPH_USE_GMP_FROM
  # and OOMPH_USE_MPFR_FROM, respectively
  find_package(GMP REQUIRED GLOBAL)
  find_package(MPFR REQUIRED GLOBAL)
  find_package(Boost 1.83.0 REQUIRED COMPONENTS ${REQUIRED_BOOST_COMPONENTS} GLOBAL PATHS ${OOMPH_USE_BOOST_FROM} NO_DEFAULT_PATH)
  find_package(CGAL 5.6 REQUIRED GLOBAL PATHS ${OOMPH_USE_CGAL_FROM} NO_DEFAULT_PATH)

  if(TARGET GMP::GMP)
    message(STATUS "Hurray! The target GMP::GMP is defined!")
  else()
    message(FATAL_ERROR "Target GMP::GMP is not defined!")
  endif()
  if(TARGET MPFR::MPFR)
    message(STATUS "Hurray! The target MPFR::MPFR is defined!")
  else()
    message(FATAL_ERROR "Target MPFR::MPFR is not defined!")
  endif()
  foreach(BOOST_COMPONENT IN LISTS REQUIRED_BOOST_COMPONENTS)
    if(TARGET Boost::${BOOST_COMPONENT})
      message(STATUS "Hurray! The target Boost::${BOOST_COMPONENT} is defined!")
    else()
      message(FATAL_ERROR "Target Boost::${BOOST_COMPONENT} is not defined!")
    endif()
  endforeach()
  if(TARGET CGAL::CGAL)
    message(STATUS "Hurray! The target CGAL::CGAL is defined!")
  else()
    message(FATAL_ERROR "Target CGAL::CGAL is not defined!")
  endif()

  list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_GMP OOMPH_HAS_MPFR OOMPH_HAS_BOOST OOMPH_HAS_CGAL)
  set(OOMPH_HAS_GMP TRUE CACHE INTERNAL "")
  set(OOMPH_HAS_MPFR TRUE CACHE INTERNAL "")
  set(OOMPH_HAS_BOOST TRUE CACHE INTERNAL "")
  set(OOMPH_HAS_CGAL TRUE CACHE INTERNAL "")
endif()

# MUMPS
if(OOMPH_USE_MUMPS_FROM)
  set(SCALAPACK_ROOT "${OOMPH_USE_MUMPS_FROM}")
  find_package(MUMPS 5.6.2.4 REQUIRED GLOBAL PATHS ${OOMPH_USE_MUMPS_FROM} NO_DEFAULT_PATH)
  if(TARGET MUMPS::MUMPS)
    message(STATUS "Hurray! The target MUMPS::MUMPS is defined!")
  else()
    message(FATAL_ERROR "Target MUMPS::MUMPS is not defined!")
  endif()
  list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_MUMPS)
  set(OOMPH_HAS_MUMPS TRUE CACHE INTERNAL "")
endif()

# HYPRE
if(OOMPH_USE_HYPRE_FROM)
  find_package(HYPRE 2.29.0 REQUIRED GLOBAL PATHS ${OOMPH_USE_HYPRE_FROM} NO_DEFAULT_PATH)
  if(TARGET HYPRE::HYPRE)
    message(STATUS "Hurray! The target HYPRE::HYPRE is defined!")
  else()
    message(FATAL_ERROR "Target HYPRE::HYPRE is not defined!")
  endif()
  list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_HYPRE)
  set(OOMPH_HAS_HYPRE TRUE CACHE INTERNAL "")
endif()

# Trilinos
if(OOMPH_USE_TRILINOS_FROM)
  find_package(Trilinos 14.4.0 REQUIRED GLOBAL PATHS ${OOMPH_USE_TRILINOS_FROM} NO_DEFAULT_PATH)
  if(TARGET Trilinos::all_libs)
    message(STATUS "Hurray! The target Trilinos::all_libs is defined!")
  else()
    message(FATAL_ERROR "Target Trilinos::all_libs is not defined!")
  endif()
  list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_TRILINOS)
  set(OOMPH_HAS_TRILINOS TRUE CACHE INTERNAL "")
endif()

message(VERBOSE "Finished locating third-party libraries")
# ------------------------------------------------------------------------------
# cmake-format: on
