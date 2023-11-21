# cmake-format: off
# ==============================================================================
# Locate optional third-party libraries requested by the user.
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered external_distributions subdirectory")

# Targets to store the accumulated compiler flags and external libs to link
set(EXTERNAL_DIST_CXX_DEFINITIONS "")

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
  find_package(MPFR REQUIRED GLOBAL )
  find_package(Boost 1.83.0 REQUIRED COMPONENTS ${REQUIRED_BOOST_COMPONENTS} GLOBAL PATHS ${OOMPH_USE_BOOST_FROM} NO_DEFAULT_PATH)
  find_package(CGAL 5.6 REQUIRED GLOBAL PATHS ${OOMPH_USE_CGAL_FROM} NO_DEFAULT_PATH)

  if (TARGET GMP::GMP)
    message(STATUS "Hurray! The target GMP::GMP is defined!")
  else()
    message(FATAL_ERROR "Target GMP::GMP is not defined!")
  endif()
  if (TARGET MPFR::MPFR)
    message(STATUS "Hurray! The target MPFR::MPFR is defined!")
  else()
    message(FATAL_ERROR "Target MPFR::MPFR is not defined!")
  endif()
  foreach(BOOST_COMPONENT IN LISTS REQUIRED_BOOST_COMPONENTS)
    if (TARGET Boost::${BOOST_COMPONENT})
      message(STATUS "Hurray! The target Boost::${BOOST_COMPONENT} is defined!")
    else()
      message(FATAL_ERROR "Target Boost::${BOOST_COMPONENT} is not defined!")
    endif()
  endforeach()
  if (TARGET CGAL::CGAL)
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
  if (TARGET MUMPS::MUMPS)
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
  if (TARGET HYPRE::HYPRE)
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
  if (TARGET Trilinos::all_libs)
    message(STATUS "Hurray! The target Trilinos::all_libs is defined!")
  else()
    message(FATAL_ERROR "Target Trilinos::all_libs is not defined!")
  endif()
  list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_TRILINOS)
  set(OOMPH_HAS_TRILINOS TRUE CACHE INTERNAL "")
endif()

# Make the accumulated C++ definitions visible in the parent scope
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

message(VERBOSE "Leaving external_distributions subdirectory")
# ------------------------------------------------------------------------------
# cmake-format: on
