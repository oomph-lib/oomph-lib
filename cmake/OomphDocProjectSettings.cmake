# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Docs the project configuration variables and compiler definitions.
#
# USAGE:
# ------
#   include(OomphDocProjectSettings)
#   oomph_doc_project_settings()
# =============================================================================
# cmake-format: on
include_guard()

# ------------------------------------------------------------------------------
function(get_all_targets BASE_DIR OUTPUT_VARIABLE)
  set(ALL_TARGETS)
  get_all_targets_recursive(ALL_TARGETS ${${BASE_DIR}})
  set(${OUTPUT_VARIABLE} ${ALL_TARGETS} PARENT_SCOPE)
endfunction()

macro(get_all_targets_recursive ALL_TARGETS CURRENT_DIRECTORY)
  get_property(
    subdirectories
    DIRECTORY ${CURRENT_DIRECTORY}
    PROPERTY SUBDIRECTORIES)

  foreach(subdir ${subdirectories})
    get_all_targets_recursive(${ALL_TARGETS} ${subdir})
  endforeach()

  get_property(
    CURRENT_DIRECTORY_TARGETS
    DIRECTORY ${CURRENT_DIRECTORY}
    PROPERTY BUILDSYSTEM_TARGETS)
  list(APPEND ${ALL_TARGETS} ${CURRENT_DIRECTORY_TARGETS})
endmacro()
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
function(oomph_doc_project_settings)
  # Define the supported set of keywords
  set(PREFIX ARG)
  set(FLAGS ENABLE_ALSO_PRINT_SETTINGS_AFTER_INSTALL ENABLE_SAVE_TO_FILE)
  set(SINGLE_VALUE_ARGS)
  set(MULTI_VALUE_ARGS)
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}"
                        "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

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

  # Initialise
  set(MARKER "⦿ ")
  set(OOMPH_SETTINGS_MESSAGE "\n")

  # Append configuration options
  string(
    APPEND
    OOMPH_SETTINGS_MESSAGE
    "********************************************************************************\n"
  )
  string(APPEND OOMPH_SETTINGS_MESSAGE
         "${BOLD}OOMPH-LIB THIRD-PARTY LIBRARIES OPTIONS:${RESET}\n")
  string(
    APPEND
    OOMPH_SETTINGS_MESSAGE
    "********************************************************************************\n"
  )
  foreach(OPTION IN LISTS OOMPH_CONFIG_VARS)
    if(NOT ${OPTION})
      string(
        APPEND OOMPH_SETTINGS_MESSAGE
        "  ${RED}${MARKER}${RESET} ${OPTION}: ${RED}'${${OPTION}}'${RESET}\n")
    else()
      string(
        APPEND
        OOMPH_SETTINGS_MESSAGE
        "  ${GREEN}${MARKER}${RESET} ${OPTION}: ${GREEN}'${${OPTION}}'${RESET}\n"
      )
    endif()
  endforeach()
  string(
    APPEND
    OOMPH_SETTINGS_MESSAGE
    "********************************************************************************\n"
  )

  # Doc it
  message(NOTICE "${OOMPH_SETTINGS_MESSAGE}")

  file(RELATIVE_PATH OOMPH_RELPATH_TO_BUILD_DIR ${CMAKE_CURRENT_SOURCE_DIR}
       ${CMAKE_CURRENT_BINARY_DIR})

  # Print info about the build and install steps
  message(
    NOTICE
    "${BOLD_MAGENTA}Project configured! Don't forget to run the build step with:\n\n"
    "\tcmake --build ${OOMPH_RELPATH_TO_BUILD_DIR}\n\n"
    "then the install step with:\n\n"
    "\tcmake --install ${OOMPH_RELPATH_TO_BUILD_DIR}${RESET}\n")

  # Construct a list of all of the targets we have defined to build
  set(BASE_DIR_FOR_FINDING_TARGETS "${CMAKE_CURRENT_SOURCE_DIR}/src")
  get_all_targets(BASE_DIR_FOR_FINDING_TARGETS ALL_TARGETS)

  # Write the post-build reminder to file; we'll print it from file
  file(
    WRITE "${CMAKE_CURRENT_BINARY_DIR}/oomphlib-install-step-reminder.txt"
    "\n${BOLD_MAGENTA} Project built! Don't forget to run the install step with:\n\n\tcmake --install ${OOMPH_RELPATH_TO_BUILD_DIR}${RESET}\n"
  )

  # Create a target that gets run after all other targets in the src/ directory
  add_custom_target(
    post_build_remind_user_to_install_library ALL
    ${CMAKE_COMMAND} -E cat
    "${CMAKE_CURRENT_BINARY_DIR}/oomphlib-install-step-reminder.txt"
    COMMENT "Post-build reminder to install the library")
  add_dependencies(post_build_remind_user_to_install_library ${ALL_TARGETS})

  # Log to file if needed
  set(OUTPUT_FILE
      "${CMAKE_CURRENT_BINARY_DIR}/oomphlib-third-party-libraries-info.log")
  if(${PREFIX}_ENABLE_SAVE_TO_FILE
     OR ${PREFIX}_ENABLE_ALSO_PRINT_SETTINGS_AFTER_INSTALL)
    file(WRITE
         "${CMAKE_CURRENT_BINARY_DIR}/oomphlib-third-party-libraries-info.log"
         "${OOMPH_SETTINGS_MESSAGE}")
  endif()

  # Doc during install step; oomph_doc_project_settings(...) must be called at
  # the end of the configuration step for this to get doced at the end
  if(${PREFIX}_ENABLE_ALSO_PRINT_SETTINGS_AFTER_INSTALL)
    install(CODE "execute_process(COMMAND cat \"${OUTPUT_FILE}\")")
  endif()
endfunction()
# ------------------------------------------------------------------------------
