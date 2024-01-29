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

  # Do we also need to build the main project with MPI?
  if (OOMPH_ENABLE_MPI)
    list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_ENABLE_MPI=ON")
    list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_ENABLE_MPI\": \"ON\"")
  endif()

  # OpenBLAS
  if(OOMPH_USE_OPENBLAS_FROM)
    list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_OPENBLAS_FROM=${OOMPH_USE_OPENBLAS_FROM}")
    list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_OPENBLAS_FROM\": \"${OOMPH_USE_OPENBLAS_FROM}\",")
  elseif(OOMPH_BUILD_OPENBLAS)
    list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_OPENBLAS_FROM=${OPENBLAS_INSTALL_DIR}")
    list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_OPENBLAS_FROM\": \"${OPENBLAS_INSTALL_DIR}\",")
  endif()

  # SuperLU/SuperLUDist and its dependencies
  if(OOMPH_BUILD_SUPERLU OR OOMPH_BUILD_SUPERLU_DIST)
    # GKlib
    if (OOMPH_USE_GKLIB_FROM)
      list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_GKLIB_FROM=${OOMPH_USE_GKLIB_FROM}")
      list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_GKLIB_FROM\": \"${OOMPH_USE_GKLIB_FROM}\",")
    else()
      list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_GKLIB_FROM=${GKLIB_INSTALL_DIR}")
      list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_GKLIB_FROM\": \"${GKLIB_INSTALL_DIR}\",")
    endif()

    # METIS
    if (OOMPH_USE_METIS_FROM)
      list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_METIS_FROM=${OOMPH_USE_METIS_FROM}")
      list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_METIS_FROM\": \"${OOMPH_USE_METIS_FROM}\",")
    else()
      list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_METIS_FROM=${METIS_INSTALL_DIR}")
      list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_METIS_FROM\": \"${METIS_INSTALL_DIR}\",")
    endif()

    # SuperLU
    list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_SUPERLU_FROM=${SUPERLU_INSTALL_DIR}")
    list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_SUPERLU_FROM\": \"${SUPERLU_INSTALL_DIR}\",")

    if(OOMPH_BUILD_SUPERLU_DIST)
      # ParMETIS
      if (OOMPH_USE_PARMETIS_FROM)
        list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_PARMETIS_FROM=${OOMPH_USE_PARMETIS_FROM}")
        list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_PARMETIS_FROM\": \"${OOMPH_USE_PARMETIS_FROM}\",")
      else()
        list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_PARMETIS_FROM=${PARMETIS_INSTALL_DIR}")
        list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_PARMETIS_FROM\": \"${PARMETIS_INSTALL_DIR}\",")
      endif()

      # SuperLUDist
      list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_SUPERLU_DIST_FROM=${SUPERLU_DIST_INSTALL_DIR}")
      list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_SUPERLU_DIST_FROM\": \"${SUPERLU_DIST_INSTALL_DIR}\",")
    endif()
  endif()

  # CGAL and its dependencies
  if (OOMPH_BUILD_CGAL)
    # GMP
    if (OOMPH_USE_GMP_FROM)
      list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_GMP_FROM=${OOMPH_USE_GMP_FROM}")
      list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_GMP_FROM\": \"${OOMPH_USE_GMP_FROM}\",")
    else()
      list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_GMP_FROM=${GMP_INSTALL_DIR}")
      list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_GMP_FROM\": \"${GMP_INSTALL_DIR}\",")
    endif()

    # MPFR
    if (OOMPH_USE_MPFR_FROM)
      list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_MPFR_FROM=${OOMPH_USE_MPFR_FROM}")
      list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_MPFR_FROM\": \"${OOMPH_USE_MPFR_FROM}\",")
    else()
      list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_MPFR_FROM=${MPFR_INSTALL_DIR}")
      list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_MPFR_FROM\": \"${MPFR_INSTALL_DIR}\",")
    endif()

    # Boost
    if (OOMPH_USE_BOOST_FROM)
      list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_BOOST_FROM=${OOMPH_USE_BOOST_FROM}")
      list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_BOOST_FROM\": \"${OOMPH_USE_BOOST_FROM}\",")
    else()
      list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_BOOST_FROM=${BOOST_INSTALL_DIR}")
      list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_BOOST_FROM\": \"${BOOST_INSTALL_DIR}\",")
    endif()

    # CGAL
    list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_CGAL_FROM=${CGAL_INSTALL_DIR}")
    list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_CGAL_FROM\": \"${CGAL_INSTALL_DIR}\",")
  endif()

  # MUMPS
  if (OOMPH_USE_MUMPS_FROM)
    list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_MUMPS_FROM=${OOMPH_USE_MUMPS_FROM}")
    list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_MUMPS_FROM\": \"${OOMPH_USE_MUMPS_FROM}\",")
  elseif (OOMPH_BUILD_MUMPS)
    list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_MUMPS_FROM=${MUMPS_INSTALL_DIR}")
    list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_MUMPS_FROM\": \"${MUMPS_INSTALL_DIR}\",")
  endif()

  # HYPRE
  if (OOMPH_USE_HYPRE_FROM)
    list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_HYPRE_FROM=${OOMPH_USE_HYPRE_FROM}")
    list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_HYPRE_FROM\": \"${OOMPH_USE_HYPRE_FROM}\",")
  elseif (OOMPH_BUILD_HYPRE)
    list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_HYPRE_FROM=${HYPRE_INSTALL_DIR}")
    list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_HYPRE_FROM\": \"${HYPRE_INSTALL_DIR}\",")
  endif()

  # TRILINOS
  if (OOMPH_USE_TRILINOS_FROM)
    list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_TRILINOS_FROM=${OOMPH_USE_TRILINOS_FROM}")
    list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_TRILINOS_FROM\": \"${OOMPH_USE_TRILINOS_FROM}\",")
  elseif (OOMPH_BUILD_TRILINOS)
    list(APPEND OOMPH_CONFIGURE_FLAGS_LIST "-DOOMPH_USE_TRILINOS_FROM=${TRILINOS_INSTALL_DIR}")
    list(APPEND OOMPH_JSON_PRESET_FLAGS_LIST "\"OOMPH_USE_TRILINOS_FROM\": \"${TRILINOS_INSTALL_DIR}\",")
  endif()

  string(JOIN " " OOMPH_CONFIGURE_FLAGS ${OOMPH_CONFIGURE_FLAGS_LIST})
  string(JOIN "\n       " OOMPH_JSON_PRESET_TEXT ${OOMPH_JSON_PRESET_FLAGS_LIST})

  # Log to file
  set(CMAKE_FLAGS_FOR_OOMPH_LIB_FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake_flags_for_oomph_lib.txt")
  file(WRITE "${CMAKE_FLAGS_FOR_OOMPH_LIB_FILE}" "${OOMPH_CONFIGURE_FLAGS}")

  # Dump the usage instructions to file so we can just print the contents at the end of the build step
  file(RELATIVE_PATH RELPATH_TO_CMAKE_FLAGS_TXT_FILE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_FLAGS_FOR_OOMPH_LIB_FILE})
  set(
    OOMPH_THIRD_PARTY_LIBARIES_USAGE_MESSAGE "\n${BOLD_MAGENTA}Now the libraries have been built, you can use the following command to configure the main oomph-lib project:\n\n \tcmake -G Ninja ${OOMPH_CONFIGURE_FLAGS} -B build\n\nAlternatively, if you're using a preset file, place the following under the \"cacheVariables\" key:\n
       ${OOMPH_JSON_PRESET_TEXT}\n\nI've also printed these values to the file:\n\n\t${RELPATH_TO_CMAKE_FLAGS_TXT_FILE}\n\nso you could simply jump up to the root oomph-lib directory and run:\n\n\tcmake -G Ninja $(cat external_distributions/build/cmake_flags_for_oomph_lib.txt) -B build\n\nClever, I know!${RESET}")
  list(JOIN "" OOMPH_THIRD_PARTY_LIBARIES_USAGE_MESSAGE "${OOMPH_THIRD_PARTY_LIBARIES_USAGE_MESSAGE}")

  # Create a target to print the usage instructions AFTER we've built everything
  set(OOMPH_THIRD_PARTY_LIBARIES_USAGE_FILE "${CMAKE_CURRENT_BINARY_DIR}/usage.txt")
  file(WRITE "${OOMPH_THIRD_PARTY_LIBARIES_USAGE_FILE}" "${OOMPH_THIRD_PARTY_LIBARIES_USAGE_MESSAGE}")
  add_custom_target(
    print_usage_message_after_build ALL
    ${CMAKE_COMMAND} -E cat "${OOMPH_THIRD_PARTY_LIBARIES_USAGE_FILE}")

  # Make sure to only print the message after all the targets have been built
  if(TARGET openblas)
    add_dependencies(print_usage_message_after_build openblas)
  endif()
  if(TARGET gklib)
    add_dependencies(print_usage_message_after_build gklib)
  endif()
  if(TARGET metis)
    add_dependencies(print_usage_message_after_build metis)
  endif()
  if(TARGET parmetis)
    add_dependencies(print_usage_message_after_build parmetis)
  endif()
  if(TARGET superlu)
    add_dependencies(print_usage_message_after_build superlu)
  endif()
  if(TARGET superlu_dist)
    add_dependencies(print_usage_message_after_build superlu_dist)
  endif()
  if(TARGET cgal)
    add_dependencies(print_usage_message_after_build cgal)
  endif()
  if(TARGET mumps)
    add_dependencies(print_usage_message_after_build mumps)
  endif()
  if(TARGET hypre)
    add_dependencies(print_usage_message_after_build hypre)
  endif()
  if(TARGET trilinos)
    add_dependencies(print_usage_message_after_build trilinos)
  endif()
endfunction()
# ------------------------------------------------------------------------------
# cmake-format: on
