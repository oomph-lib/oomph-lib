# cmake-format: off
set(OOMPH_C_AND_CXX_FLAGS)
set(OOMPH_C_FLAGS ${OOMPH_C_FLAGS} ${OOMPH_C_AND_CXX_FLAGS})
set(OOMPH_CXX_FLAGS ${OOMPH_CXX_FLAGS} ${OOMPH_C_AND_CXX_FLAGS})

# list(APPEND OOMPH_C_FLAGS -Wno-compare-distinct-pointer-types)
list(APPEND OOMPH_C_FLAGS -Wno-implicit-function-declaration)
# list(APPEND OOMPH_C_FLAGS -Wno-incompatible-pointer-types)
# list(APPEND OOMPH_C_FLAGS -Wno-mismatched-new-delete)
# list(APPEND OOMPH_C_FLAGS -Wno-format-extra-args)
# list(APPEND OOMPH_C_FLAGS -Wno-enum-conversion)
# list(APPEND OOMPH_C_FLAGS -Wno-format-security)
# list(APPEND OOMPH_C_FLAGS -Wno-unused-variable)
# list(APPEND OOMPH_C_FLAGS -Wno-missing-braces)
# list(APPEND OOMPH_C_FLAGS -Wno-int-conversion)
# list(APPEND OOMPH_C_FLAGS -Wno-unused-value)
# list(APPEND OOMPH_C_FLAGS -Wno-implicit-int)
# list(APPEND OOMPH_C_FLAGS -Wno-parentheses)
# list(APPEND OOMPH_C_FLAGS -Wno-return-type)
# list(APPEND OOMPH_C_FLAGS -Wno-unsequenced)
# list(APPEND OOMPH_C_FLAGS "-Wno-#warnings")
# list(APPEND OOMPH_C_FLAGS -Wno-format)
# list(APPEND OOMPH_C_FLAGS -Wno-switch)
# list(APPEND OOMPH_CXX_FLAGS -Wno-c++11-compat-deprecated-writable-strings)
list(APPEND OOMPH_CXX_FLAGS -Wno-instantiation-after-specialization)
# list(APPEND OOMPH_CXX_FLAGS -Wno-potentially-evaluated-expression)
list(APPEND OOMPH_CXX_FLAGS -Wno-implicit-function-declaration)
# list(APPEND OOMPH_CXX_FLAGS -Wno-deprecated-declarations)
# list(APPEND OOMPH_CXX_FLAGS -Wno-sometimes-uninitialized)
list(APPEND OOMPH_CXX_FLAGS -Wno-undefined-var-template)
# list(APPEND OOMPH_CXX_FLAGS -Wno-mismatched-new-delete)
# list(APPEND OOMPH_CXX_FLAGS -Wno-int-to-pointer-cast)
# list(APPEND OOMPH_CXX_FLAGS -Wno-unused-variable)

include(CheckCCompilerFlag)
foreach(C_COMPILE_FLAG IN LISTS OOMPH_C_FLAGS)
  check_c_compiler_flag(C ${C_COMPILE_FLAG} C_COMPILER_HAS_FLAG)
  if (C_COMPILER_HAS_FLAG)
    add_compile_options($<$<COMPILE_LANGUAGE:C>:${C_COMPILE_FLAG}>)
  endif()
endforeach()

include(CheckCXXCompilerFlag)
foreach(CXX_COMPILE_FLAG IN LISTS OOMPH_CXX_FLAGS)
  check_cxx_compiler_flag(CXX ${CXX_COMPILE_FLAG} CXX_COMPILER_HAS_FLAG)
  if (CXX_COMPILER_HAS_FLAG)
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
  check_fortran_compiler_flag(Fortran ${FORTRAN_COMPILE_FLAG} FORTRAN_COMPILER_HAS_FLAG)
  if (FORTRAN_COMPILER_HAS_FLAG)
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:${FORTRAN_COMPILE_FLAG}>)
  endif()
# ----------------------------------------------------------------------------
endforeach()
# cmake-format: on
