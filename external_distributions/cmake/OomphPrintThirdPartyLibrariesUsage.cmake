# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Tells the user which flags to use the built third-party libraries in oomph-lib.
#
# USAGE:
# ------
#   include(OomphPrintThirdPartyLibrariesUsage)
#   oomph_print_third_party_libraries_usage()
# =============================================================================
include_guard()

# ------------------------------------------------------------------------------
function(oomph_print_third_party_libraries_usage)
  # Colourising
  if(NOT WIN32)
    string(ASCII 27 Esc)
    set(RESET "${Esc}[m")
    set(BOLD "${Esc}[1m")
    set(RED "${Esc}[31m")
    set(GREEN "${Esc}[32m")
    set(YELLOW "${Esc}[33m")
    set(BLUE "${Esc}[34m")
    set(MAGENTA "${Esc}[35m")
    set(CYAN "${Esc}[36m")
    set(WHITE "${Esc}[37m")
    set(BOLD_RED "${Esc}[1;31m")
    set(BOLD_GREEN "${Esc}[1;32m")
    set(BOLD_YELLOW "${Esc}[1;33m")
    set(BOLD_BLUE "${Esc}[1;34m")
    set(BOLD_MAGENTA "${Esc}[1;35m")
    set(BOLD_CYAN "${Esc}[1;36m")
    set(BOLD_WHITE "${Esc}[1;37m")
  endif()

  # Generate information about how to use these libraries in oomph-lib.
  set(OOMPH_CONFIGURE_FLAGS_LIST)
  set(OOMPH_JSON_PRESET_FLAGS_LIST)

  macro(_oomph_add_flag name value)
    # Adds '-Dname=value' to the shell flags and '"name": "value"' to the JSON preset list
    list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-D${name}=${value}")
    list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"${name}\": \"${value}\"")
  endmacro()

  # What's the build type?
  _oomph_add_flag(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}")

  # Do we also need to build the main project with MPI?
  if (OOMPH_ENABLE_MPI)
    _oomph_add_flag(OOMPH_ENABLE_MPI "ON")
  endif()

  # OpenBLAS
  if(OOMPH_USE_OPENBLAS_FROM)
    _oomph_add_flag(OOMPH_USE_OPENBLAS_FROM "${OOMPH_USE_OPENBLAS_FROM}")
  elseif(OOMPH_BUILD_OPENBLAS)
    _oomph_add_flag(OOMPH_USE_OPENBLAS_FROM "${OPENBLAS_INSTALL_DIR}")
  endif()

  # SuperLU/SuperLUDist and its dependencies
  if(OOMPH_BUILD_SUPERLU OR OOMPH_BUILD_SUPERLU_DIST)
    # GKlib
    if (OOMPH_USE_GKLIB_FROM)
      _oomph_add_flag(OOMPH_USE_GKLIB_FROM "${OOMPH_USE_GKLIB_FROM}")
    else()
      _oomph_add_flag(OOMPH_USE_GKLIB_FROM "${GKLIB_INSTALL_DIR}")
    endif()

    # METIS
    if (OOMPH_USE_METIS_FROM)
      _oomph_add_flag(OOMPH_USE_METIS_FROM "${OOMPH_USE_METIS_FROM}")
    else()
      _oomph_add_flag(OOMPH_USE_METIS_FROM "${METIS_INSTALL_DIR}")
    endif()

    # SuperLU
    _oomph_add_flag(OOMPH_USE_SUPERLU_FROM "${SUPERLU_INSTALL_DIR}")

    if(OOMPH_BUILD_SUPERLU_DIST)
      # ParMETIS
      if (OOMPH_USE_PARMETIS_FROM)
        _oomph_add_flag(OOMPH_USE_PARMETIS_FROM "${OOMPH_USE_PARMETIS_FROM}")
      else()
        _oomph_add_flag(OOMPH_USE_PARMETIS_FROM "${PARMETIS_INSTALL_DIR}")
      endif()

      # SuperLUDist
      _oomph_add_flag(OOMPH_USE_SUPERLU_DIST_FROM "${SUPERLU_DIST_INSTALL_DIR}")
    endif()
  endif()

  # CGAL and its dependencies
  if (OOMPH_BUILD_CGAL)
    # Boost
    if (OOMPH_USE_BOOST_FROM)
      _oomph_add_flag(OOMPH_USE_BOOST_FROM "${OOMPH_USE_BOOST_FROM}")
    else()
      _oomph_add_flag(OOMPH_USE_BOOST_FROM "${BOOST_INSTALL_DIR}")
    endif()

    # CGAL
    _oomph_add_flag(OOMPH_USE_CGAL_FROM "${CGAL_INSTALL_DIR}")
  endif()

  # MUMPS
  if (OOMPH_USE_MUMPS_FROM)
    _oomph_add_flag(OOMPH_USE_MUMPS_FROM "${OOMPH_USE_MUMPS_FROM}")
  elseif (OOMPH_BUILD_MUMPS)
    _oomph_add_flag(OOMPH_USE_MUMPS_FROM "${MUMPS_INSTALL_DIR}")
  endif()

  # HYPRE
  if (OOMPH_USE_HYPRE_FROM)
    _oomph_add_flag(OOMPH_USE_HYPRE_FROM "${OOMPH_USE_HYPRE_FROM}")
  elseif (OOMPH_BUILD_HYPRE)
    _oomph_add_flag(OOMPH_USE_HYPRE_FROM "${HYPRE_INSTALL_DIR}")
  endif()

  # TRILINOS
  if (OOMPH_USE_TRILINOS_FROM)
    _oomph_add_flag(OOMPH_USE_TRILINOS_FROM "${OOMPH_USE_TRILINOS_FROM}")
  elseif (OOMPH_BUILD_TRILINOS)
    _oomph_add_flag(OOMPH_USE_TRILINOS_FROM "${TRILINOS_INSTALL_DIR}")
  endif()

  # The text to use a single long CMake command
  string(JOIN " " OOMPH_CONFIGURE_FLAGS ${OOMPH_CONFIGURE_FLAGS_LIST})

  # The text to add to the top-level CMakeUserPreset.json
  message("${OOMPH_JSON_PRESET_TEXT}")
  list(TRANSFORM OOMPH_JSON_PRESET_FLAGS_LIST PREPEND "       ")
  string(JOIN ",\n" OOMPH_JSON_PRESET_TEXT ${OOMPH_JSON_PRESET_FLAGS_LIST})

  # Output files
  set(OOMPH_THIRD_PARTY_LIBARIES_USAGE_FILE "${CMAKE_CURRENT_BINARY_DIR}/usage.txt")
  set(CMAKE_FLAGS_FOR_OOMPH_LIB_TXT_FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake_flags_for_oomph_lib.txt")
  set(CMAKE_FLAGS_FOR_OOMPH_LIB_JSON_FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake_flags_for_oomph_lib.json")

  # The files after they have been *installed*
  set(CMAKE_FLAGS_FOR_INSTALLED_OOMPH_LIB_TXT_FILE "${CMAKE_INSTALL_PREFIX}/cmake_flags_for_oomph_lib.txt")
  set(CMAKE_FLAGS_FOR_INSTALLED_OOMPH_LIB_JSON_FILE "${CMAKE_INSTALL_PREFIX}/cmake_flags_for_oomph_lib.json")

  # Dump the text file
  file(WRITE "${CMAKE_FLAGS_FOR_OOMPH_LIB_TXT_FILE}" "${OOMPH_CONFIGURE_FLAGS}")

  # Dump the JSON file
  set(CMAKE_FLAGS_FOR_OOMPH_LIB_JSON_TEXT "{" "${OOMPH_JSON_PRESET_TEXT}" "}")
  string(JOIN "\n" CMAKE_FLAGS_FOR_OOMPH_LIB_JSON_TEXT ${CMAKE_FLAGS_FOR_OOMPH_LIB_JSON_TEXT})
  file(WRITE "${CMAKE_FLAGS_FOR_OOMPH_LIB_JSON_FILE}" "${CMAKE_FLAGS_FOR_OOMPH_LIB_JSON_TEXT}")

  # Relative paths for telling the user where to look for the text file
  cmake_path(GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH OOMPH_LIB_ROOT_DIR)
  cmake_path(
    RELATIVE_PATH CMAKE_FLAGS_FOR_INSTALLED_OOMPH_LIB_TXT_FILE
    BASE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE CURRENT_RELPATH_TO_CMAKE_FLAGS_TXT_FILE
  )
  cmake_path(
    RELATIVE_PATH CMAKE_FLAGS_FOR_INSTALLED_OOMPH_LIB_TXT_FILE
    BASE_DIRECTORY ${OOMPH_LIB_ROOT_DIR}
    OUTPUT_VARIABLE ROOT_RELPATH_TO_CMAKE_FLAGS_TXT_FILE
  )

  # Dump the usage instructions to file so we can just print the contents at the end of the build step
  set(
    OOMPH_THIRD_PARTY_LIBRARIES_USAGE_MESSAGE
    "\nNow the libraries have been built, you can use the following command to configure the main oomph-lib project:\n\n"
    " \tcmake -G Ninja ${OOMPH_CONFIGURE_FLAGS} -B build\n\n"
    "Alternatively, if you're using a preset file, place the following under the \"cacheVariables\" key:\n\n"
    "${OOMPH_JSON_PRESET_TEXT}\n\n"
    "I've also printed these values to the file:\n\n"
    "\t${CURRENT_RELPATH_TO_CMAKE_FLAGS_TXT_FILE}\n\n"
    "so you could simply jump up to the root oomph-lib directory and run:\n\n"
    "\tcmake -G Ninja $(cat ${ROOT_RELPATH_TO_CMAKE_FLAGS_TXT_FILE}) -B build\n\n"
    "Clever, I know!"
  )

  # Turn OOMPH_THIRD_PARTY_LIBRARIES_USAGE_MESSAGE from a list of strings to a single string
  list(JOIN "" OOMPH_THIRD_PARTY_LIBRARIES_USAGE_MESSAGE "${OOMPH_THIRD_PARTY_LIBRARIES_USAGE_MESSAGE}")

  # Remove the semicolon CMake adds
  string(REPLACE ";" "" OOMPH_THIRD_PARTY_LIBRARIES_USAGE_MESSAGE "${OOMPH_THIRD_PARTY_LIBRARIES_USAGE_MESSAGE}")

  # Dump the usage to a file
  file(WRITE "${OOMPH_THIRD_PARTY_LIBARIES_USAGE_FILE}" "${BOLD_MAGENTA}${OOMPH_THIRD_PARTY_LIBRARIES_USAGE_MESSAGE}${RESET}")

  # All the targets we build you might have
  set(THIRD_PARTY_LIBS
    openblas
    gklib
    metis
    parmetis
    superlu
    superlu_dist
    cgal
    mumps
    hypre
    trilinos
  )

  # Build a list of only the targets that actually exist
  set(EXISTING_TPL_DEPS)
  foreach(TPL_LIB IN LISTS THIRD_PARTY_LIBS)
    if(TARGET ${TPL_LIB})
      list(APPEND EXISTING_TPL_DEPS ${TPL_LIB})
    endif()
  endforeach()

  # Create a custom target that is only run after all of the libraries have been built (and
  # thus after they're installed since they're installed during the build step by ExternalProject)
  add_custom_target(
    print_usage_message_after_build ALL
    COMMAND ${CMAKE_COMMAND} -E cat "${OOMPH_THIRD_PARTY_LIBARIES_USAGE_FILE}"
    DEPENDS ${EXISTING_TPL_DEPS}
  )

  # Install the text file and JSON file with the build flags for the top-level project
  add_custom_command(
    TARGET print_usage_message_after_build
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_FLAGS_FOR_OOMPH_LIB_TXT_FILE}" "${CMAKE_INSTALL_PREFIX}/"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_FLAGS_FOR_OOMPH_LIB_JSON_FILE}" "${CMAKE_INSTALL_PREFIX}/"
    COMMAND ${CMAKE_COMMAND} -E rm -f "${CMAKE_FLAGS_FOR_OOMPH_LIB_TXT_FILE}"
    COMMAND ${CMAKE_COMMAND} -E rm -f "${CMAKE_FLAGS_FOR_OOMPH_LIB_JSON_FILE}")
endfunction()
# ------------------------------------------------------------------------------
# cmake-format: on
