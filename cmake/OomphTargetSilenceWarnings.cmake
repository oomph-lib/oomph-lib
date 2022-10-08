# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Makes it easy to provide flags to disable C, C++ or Fortran warnings. Provides
# checks to make sure that the compiler supports the specified flag then sets it.
#
# USAGE:
# ------
#    include(OomphTargetSilenceWarnings)
#    oomph_target_silence_warnings(TARGET                   <target-name>
#                                  [C_COMPILE_FLAGS         <c-flags>]
#                                  [CXX_COMPILE_FLAGS       <c++-flags>]
#                                  [FORTRAN_COMPILE_FLAGS   <fortran-flags>]
#                                  [WARN_IF_FLAG_NOT_FOUND])
#
# NOTE (1): At least one of the *_COMPILE_FLAGS arguments needs to be specified.
# NOTE (2): This module is intended for usage in cases that the warnings need
# to be disabled on a specific target. For this reason, the flags are added to
# the target under the PRIVATE keyword; this ensures that they are not forwarded
# to any targets that link to it.
# =============================================================================
# cmake-format: on
include_guard()

function(oomph_target_silence_warnings)
  # Define the supported set of keywords
  set(PREFIX ARG)
  set(FLAGS WARN_IF_FLAG_NOT_FOUND)
  set(SINGLE_VALUE_ARGS TARGET)
  set(MULTI_VALUE_ARGS C_COMPILE_FLAGS CXX_COMPILE_FLAGS FORTRAN_COMPILE_FLAGS)

  # Process the arguments passed in
  include(CMakeParseArguments)
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}"
                        "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

  # Redefine the headers, sources, and the name of the library in this scope but
  # with clearer variable names
  set(TARGET ${${PREFIX}_TARGET})
  set(WARN_IF_FLAG_NOT_FOUND ${${PREFIX}_WARN_IF_FLAG_NOT_FOUND})
  set(C_COMPILE_FLAGS ${${PREFIX}_C_COMPILE_FLAGS})
  set(CXX_COMPILE_FLAGS ${${PREFIX}_CXX_COMPILE_FLAGS})
  set(FORTRAN_COMPILE_FLAGS ${${PREFIX}_FORTRAN_COMPILE_FLAGS})

  # Don't silence these warnings if the user wants to see everything
  if(OOMPH_DONT_SILENCE_USELESS_WARNINGS)
    return()
  endif()

  # Just stop now if no flags have been provided
  if((NOT C_COMPILE_FLAGS)
     AND (NOT CXX_COMPILE_FLAGS)
     AND (NOT FORTRAN_COMPILE_FLAGS))
    message(WARNING "No compiler flags provided! Not doing anything.")
    return()
  endif()

  # Built-in CMake modules used to make sure the compiler supports the flags
  # specified by the user
  include(CheckCCompilerFlag)
  include(CheckCXXCompilerFlag)
  include(CheckFortranCompilerFlag)

  # Assign any provided C compiler flags
  foreach(C_COMPILE_FLAG IN LISTS C_COMPILE_FLAGS)
    check_c_compiler_flag(${C_COMPILE_FLAG} C_COMPILER_HAS_FLAG)
    if(C_COMPILER_HAS_FLAG)
      target_compile_options(${TARGET}
                             PRIVATE $<$<COMPILE_LANGUAGE:C>:${C_COMPILE_FLAG}>)
    elseif(WARN_IF_FLAG_NOT_FOUND)
      message(WARNING "Unrecognised C compiler flag: '${C_COMPILE_FLAG}'")
    endif()
  endforeach()

  # Assign any provided C++ compiler flags
  foreach(CXX_COMPILE_FLAG IN LISTS CXX_COMPILE_FLAGS)
    check_cxx_compiler_flag(${CXX_COMPILE_FLAG} CXX_COMPILER_HAS_FLAG)
    if(CXX_COMPILER_HAS_FLAG)
      target_compile_options(
        ${TARGET} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${CXX_COMPILE_FLAG}>)
    elseif(WARN_IF_FLAG_NOT_FOUND)
      message(WARNING "Unrecognised C++ compiler flag: '${CXX_COMPILE_FLAG}'")
    endif()
  endforeach()

  # Assign any provided Fortran compiler flags
  foreach(FORTRAN_COMPILE_FLAG IN LISTS FORTRAN_COMPILE_FLAGS)
    check_fortran_compiler_flag(${FORTRAN_COMPILE_FLAG}
                                FORTRAN_COMPILER_HAS_FLAG)
    if(FORTRAN_COMPILER_HAS_FLAG)
      target_compile_options(
        ${TARGET}
        PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${FORTRAN_COMPILE_FLAG}>)
    elseif(WARN_IF_FLAG_NOT_FOUND)
      message(
        WARNING "Unrecognised Fortran compiler flag: '${FORTRAN_COMPILE_FLAG}'")
    endif()
  endforeach()
endfunction()
