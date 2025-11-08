# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
#
#
# USAGE:
# ------
#  include(OomphSilenceWarnings)
# =============================================================================
# cmake-format: on
include_guard()

# ----------------------------------------------------------------------------
set(OOMPH_C_AND_CXX_FLAGS)
set(OOMPH_C_FLAGS ${OOMPH_C_FLAGS} ${OOMPH_C_AND_CXX_FLAGS})
set(OOMPH_CXX_FLAGS ${OOMPH_CXX_FLAGS} ${OOMPH_C_AND_CXX_FLAGS})

list(APPEND OOMPH_C_FLAGS -Wno-implicit-function-declaration)
list(APPEND OOMPH_CXX_FLAGS -Wno-instantiation-after-specialization)
list(APPEND OOMPH_CXX_FLAGS -Wno-undefined-var-template)

if(APPLE)
  list(APPEND OOMPH_CXX_FLAGS -Wno-deprecated-declarations)
endif()

include(CheckCCompilerFlag)
foreach(C_COMPILE_FLAG IN LISTS OOMPH_C_FLAGS)
  check_c_compiler_flag(${C_COMPILE_FLAG} C_COMPILER_HAS_FLAG)
  if(C_COMPILER_HAS_FLAG)
    add_compile_options($<$<COMPILE_LANGUAGE:C>:${C_COMPILE_FLAG}>)
  endif()
endforeach()

include(CheckCXXCompilerFlag)
foreach(CXX_COMPILE_FLAG IN LISTS OOMPH_CXX_FLAGS)
  check_cxx_compiler_flag(${CXX_COMPILE_FLAG} CXX_COMPILER_HAS_FLAG)
  if(CXX_COMPILER_HAS_FLAG)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:${CXX_COMPILE_FLAG}>)
  endif()
endforeach()
# ----------------------------------------------------------------------------
# Keep around!

# Silence warnings about arithmetic IF statements being deleted features in
# Fortran 2018
set(OOMPH_FORTRAN_FLAGS)
list(APPEND OOMPH_FORTRAN_FLAGS -std=legacy)

include(CheckFortranCompilerFlag)
foreach(FORTRAN_COMPILE_FLAG IN LISTS OOMPH_FORTRAN_FLAGS)
  check_fortran_compiler_flag(${FORTRAN_COMPILE_FLAG} FORTRAN_COMPILER_HAS_FLAG)
  if(FORTRAN_COMPILER_HAS_FLAG)
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:${FORTRAN_COMPILE_FLAG}>)
  endif()
endforeach()
# ----------------------------------------------------------------------------
