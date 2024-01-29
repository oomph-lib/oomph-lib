# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Generalise the creation of an oomph-lib executable.
#
# USAGE:
# ------
#     oomph_add_executable(NAME             <executable-name>
#                          SOURCES          <executable-sources-and-headers>
#                          LIBRARIES        <required-libraries>
#                          [CXX_DEFINITIONS <preprocessor-definitions>]
#                          [CXX_OPTIONS     <compiler-options>]
#                          [CXX_STANDARD    <14|17|20>]
#                          [LINK_OPTIONS    <link-options>]
#                          [SILENCE_NO_LIBS_SUPPLIED_WARNING])
#
# NOTE: It is slightly unusual for the user to not provide any libraries to link
# against but not insane. Just in case the user made a mistake, we'll issue a
# warning and if they're sure they don't want to link against any libraries,
# they can silence the warning with the flag SILENCE_NO_LIBS_SUPPLIED_WARNING.
#
# PM: Strictly speaking, you (the user) should not specify the CXX_STANDARD;
# it is likely to lead to issues down the road. I've kept it here for debugging
# purposes for now. Will probably delete it further down the road...
#
# Example:
#
#     oomph_add_executable(NAME one_d_poisson
#                          SOURCES one_d_poisson.cc
#                          LIBRARIES oomph::poisson
#                          CXX_DEFINITIONS REFINEABLE
#                          CXX_OPTIONS -Wall -Werror
#                          LINK_OPTIONS -Wl,-no_compact_unwind -Wl,-keep_dwarf_unwind)
#
# It is worth noting that the argument supplied to CXX_DEFINITIONS does not
# require a -D prefix (as is usually required to indicate that it is a
# preprocessor definition), but you can supply it if you wish.
# =============================================================================
# cmake-format: on
include_guard()

function(oomph_add_executable)
  # Define the supported set of keywords
  set(PREFIX ARG)
  set(FLAGS SILENCE_NO_LIBS_SUPPLIED_WARNING)
  set(SINGLE_VALUE_ARGS NAME CXX_STANDARD)
  set(MULTI_VALUE_ARGS SOURCES LIBRARIES CXX_DEFINITIONS CXX_OPTIONS
      LINK_OPTIONS)

  # Process the arguments passed in
  include(CMakeParseArguments)
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}"
                        "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

  # Redefine the variables in this scope without a prefix for clarity
  set(SILENCE_NO_LIBS_SUPPLIED_WARNING
      ${${PREFIX}_SILENCE_NO_LIBS_SUPPLIED_WARNING})
  set(NAME ${${PREFIX}_NAME})
  set(SOURCES ${${PREFIX}_SOURCES})
  set(LIBRARIES ${${PREFIX}_LIBRARIES})
  set(CXX_STANDARD ${${PREFIX}_CXX_STANDARD})
  set(CXX_DEFINITIONS ${${PREFIX}_CXX_DEFINITIONS})
  set(CXX_OPTIONS ${${PREFIX}_CXX_OPTIONS})
  set(LINK_OPTIONS ${${PREFIX}_LINK_OPTIONS})

  # Make sure the arguments are valid
  if(NOT NAME)
    message(FATAL_ERROR "No NAME argument supplied.")
  elseif(NOT SOURCES)
    message(FATAL_ERROR "No SOURCES argument supplied.")
  elseif(NOT LIBRARIES)
    if(NOT DEFINED SILENCE_NO_LIBS_SUPPLIED_WARNING)
      message(
        WARNING
          "\n\
          -------------------------------------------------------------------\n\
          No LIBRARIES argument supplied for the executable \"${NAME}\" in:\n\
          \n\
          ${CMAKE_CURRENT_LIST_DIR}\n\
          \n\
          If you definitely meant to do this then you can disable the warning\n\
          by adding the flag:\n\
          \n\
          SILENCE_NO_LIBS_SUPPLIED_WARNING\n\
          \n\
          to your call to oomph_add_executable(...).\n\
          -------------------------------------------------------------------\n\
          ")
    endif()
  elseif(CXX_STANDARD)
    set(SUPPORTED_CMAKE_STANDARDS 14 17 20)
    if(NOT CXX_STANDARD IN_LIST SUPPORTED_CMAKE_STANDARDS)
      message(
        FATAL_ERROR
          "Supplied invalid argument ${CXX_STANDARD} to CXX_STANDARD. Valid arguments: ${SUPPORTED_CMAKE_STANDARDS}."
      )
    endif()
  endif()

  # Create a build of a large project of executables, e.g. demo_drivers has the
  # possibility of creating a conflict in the target name, i.e. two executables
  # in two different folders might have the same name. To avoid this restricting
  # the possible set of names an executable can take, we'll hash the path to the
  # executable and append that to the target name, then we'll set the output
  # name (i.e. the name of the executable produced by the build) to the name
  # specified by the user. NOTE: This does assume that the user won't place two
  # targets with the same name in the same folder. However, I think this is a
  # perfectly valid assumption. If you did do that and you're reading this now,
  # what are you thinking?
  string(SHA1 PATH_HASH "${CMAKE_CURRENT_LIST_DIR}")
  string(SUBSTRING ${PATH_HASH} 0 7 PATH_HASH)

  # Add executable to the build
  add_executable(${NAME}_${PATH_HASH} ${SOURCES})

  # Specify the name of the executable that the user will see
  set_target_properties(${NAME}_${PATH_HASH} PROPERTIES OUTPUT_NAME ${NAME})

  # Link to the specified libraries
  target_link_libraries(${NAME}_${PATH_HASH} PUBLIC ${LIBRARIES})

  # ~~~
  # # FIXME: XCode 15.0 problems :(
  # if(APPLE)
  #   include(CheckLinkerFlag)
  #   foreach(FLAG IN ITEMS -no_warn_duplicate_libraries -ld_classic)
  #     check_linker_flag(CXX FLAG IS_CXX_FLAG)
  #     if(IS_CXX_FLAG)
  #       target_link_options(${NAME}_${PATH_HASH} PRIVATE LINKER:${FLAG})
  #     else()
  #       message(STATUS "C++ linker flag ${FLAG} is not supported.")
  #     endif()
  #   endforeach()
  # endif()
  # ~~~

  # Provide access to the MPI libraries if we've built oomph-lib with MPI
  if(OOMPH_HAS_MPI)
    target_link_libraries(${NAME}_${PATH_HASH} PUBLIC MPI::MPI_CXX)
  endif()

  # Set the C++ standard
  if(CXX_STANDARD)
    target_compile_features(${NAME}_${PATH_HASH} PUBLIC cxx_std_${CXX_STANDARD})
  else()
    target_compile_features(${NAME}_${PATH_HASH} PUBLIC cxx_std_14)
  endif()

  # General oomph-lib specific compiler definitions
  target_compile_definitions(${NAME}_${PATH_HASH}
                             PUBLIC ${OOMPH_COMPILE_DEFINITIONS})

  # Add any user-supplied compiler definitions
  if(CXX_DEFINITIONS)
    target_compile_definitions(${NAME}_${PATH_HASH} PUBLIC ${CXX_DEFINITIONS})
  endif()

  # Add any compiler flags
  if(CXX_OPTIONS)
    target_compile_options(${NAME}_${PATH_HASH} PUBLIC ${CXX_OPTIONS})
  endif()

  # Add any linker flags
  if(LINK_OPTIONS)
    target_link_options(${NAME}_${PATH_HASH} PUBLIC ${LINK_OPTIONS})
  endif()

  # Should we try to colourise the command output?
  option(FORCE_COLORED_OUTPUT "Produce ANSI-colored output (GNU/Clang)." ON)

  # Apply the compile options. The Apple-provided Cl
  if(FORCE_COLORED_OUTPUT)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      target_compile_options(${NAME}_${PATH_HASH}
                             PRIVATE -fdiagnostics-color=always)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
      target_compile_options(${NAME}_${PATH_HASH} PRIVATE -fcolor-diagnostics)
    endif()
  endif()
endfunction()
# ------------------------------------------------------------------------------
