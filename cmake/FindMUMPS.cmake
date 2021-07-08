# cmake-format: off
# =============================================================================
# Try to find MUMPS
#
# Once done this will define:
#  MUMPS_FOUND   - system has MUMPS
#  MUMPS_INC     - MUMPS include directory
#  MUMPS_LIB     - MUMPS library
#  MUMPS_SEQ     - MUMPS MPISEQ library (if available)
#
# Usage:
#  find_package(MUMPS)
#
# Setting these changes the behavior of the search
#  PETSC_MUMPS   - If PETSc has MUMPS we use it
#  MUMPS_INC     - MUMPS include directory
#  MUMPS_LIB     - MUMPS library path
#  MUMPS_SEQ     - MUMPS MPISEQ library (if available)
# =============================================================================
# cmake-format: on

# =============================================================================
# Look for MUMPS in PETSc
# =============================================================================
if(PETSC_MUMPS)
  set(MUMPS_INC ${PETSC_MUMPS_INC})
  set(MUMPS_LIB ${PETSC_MUMPS_LIB} ${PETSC_SCALAPACK_LIB} ${PETSC_PARMETIS_LIB}
      ${PETSC_METIS_LIB})
  set(MUMPS_SEQ ${PETSC_MUMPS_SEQ})
endif()

# =============================================================================
# Look for environment variables
# =============================================================================
if(NOT PETSC_MUMPS)
  if(NOT DEFINED MUMPS_LIB)
    set(MUMPS_LIB $ENV{MUMPS_LIB})
  endif()
  if(NOT DEFINED MUMPS_INC)
    set(MUMPS_INC $ENV{MUMPS_INC})
  endif()
  if(NOT DEFINED MUMPS_SEQ)
    set(MUMPS_SEQ $ENV{MUMPS_SEQ})
  endif()
endif()

# =============================================================================
# CMake check and done
# =============================================================================
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  MUMPS
  "MUMPS could not be found: be sure to set MUMPS_LIB and MUMPS_INC in your environment variables"
  MUMPS_LIB
  MUMPS_INC)
