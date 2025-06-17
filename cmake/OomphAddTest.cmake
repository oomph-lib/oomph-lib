# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Supports the creation of an oomph-lib demo driver (i.e. integration) test.
# These tests are distinct from regular unit tests which are far more
# lightweight and do not depend on data to be tested against.
#
# USAGE:
# ------
#     oomph_add_test(TEST_NAME             <name-of-test>
#                    DEPENDS_ON            <executables/targets-required-by-test>
#                    COMMAND               <command-to-run-test>
#                    [TEST_FILES           <files-required-by-test>]
#                    [SILENCE_MISSING_VALIDATA_WARNING]
#                    [ASSUME_PURE_CMAKE_TARGET_NAMES])
#
# By default we assume that a validata/ and a validate.sh script are always
# required for a test. Here, TEST_FILES is provided so that the user can
# specify any other files that are required at run-time, like triangle input
# mesh files. The COMMAND argument is used to run the test.
#
# Important: Every validate.sh script expects ${OOMPH_ROOT_DIR} as the first
# argument to COMMAND so that it can find fpdiff.py and validate_ok_count.
#
# EXAMPLE:
# --------
#
#     oomph_add_test(TEST_NAME poisson.one_d_poisson
#                    DEPENDS_ON one_d_poisson
#                    COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
#                    TEST_FILES validate.sh validata my_extra_data_file.dat)
#
#     oomph_add_test(TEST_NAME poisson.one_d_poisson
#                    DEPENDS_ON one_d_poisson
#                    COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_VARIABLENP_RUN_COMMAND}
#                    TEST_FILES validate.sh validata my_extra_data_file.dat)
#
# NOTE 1: Arguments to DEPENDS_ON must be already-defined executables
# or targets (i.e. defined via add_executable() or oomph_add_executable()).
# =============================================================================
# cmake-format: on
include_guard()

# ------------------------------------------------------------------------------
function(oomph_add_test)
  # Define the supported set of keywords
  set(PREFIX ARG)
  set(FLAGS SILENCE_MISSING_VALIDATA_WARNING ASSUME_PURE_CMAKE_TARGET_NAMES)
  set(SINGLE_VALUE_ARGS TEST_NAME)
  set(MULTI_VALUE_ARGS TEST_FILES DEPENDS_ON COMMAND)

  # Process the arguments passed in
  include(CMakeParseArguments)
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}"
                        "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

  # Redefine the variables in this scope without a prefix for clarity
  set(SILENCE_MISSING_VALIDATA_WARNING
      ${${PREFIX}_SILENCE_MISSING_VALIDATA_WARNING})
  set(ASSUME_PURE_CMAKE_TARGET_NAMES
      ${${PREFIX}_ASSUME_PURE_CMAKE_TARGET_NAMES})
  set(TEST_NAME ${${PREFIX}_TEST_NAME})
  set(VALIDATE_SH_COMMAND ${${PREFIX}_COMMAND})
  set(DEPENDS_ON ${${PREFIX}_DEPENDS_ON})
  set(TEST_FILES ${${PREFIX}_TEST_FILES})

  # Make sure the arguments are valid
  if(NOT TEST_NAME)
    message(FATAL_ERROR "No TEST_NAME argument supplied.")
  endif()

  if(NOT DEPENDS_ON)
    message(WARNING "No DEPENDS_ON argument supplied. Can't create a test.")
    return()
  endif()

  find_program(BASH_PROGRAM bash)
  if(NOT BASH_PROGRAM)
    message(STATUS "You don't have 'bash', so I can't construct any tests!")
  endif()

  # We *nearly* always need validata, so warn if we don't have it, just in case
  # the user's forgotten to provide it
  if(NOT "validata" IN_LIST TEST_FILES)
    if(NOT SILENCE_MISSING_VALIDATA_WARNING)
      message(
        WARNING
          "\n\
        Missing validata folder in project folder:\n\
        \n\
              ${CMAKE_CURRENT_LIST_DIR}\n\
        \n\
        Did you possibly forget to provide the validation data? If not then\n\
        you can silence this warning by supplying the following flag to the\n\
        call of oomph_add_test(...):\n\
        \n\
              SILENCE_MISSING_VALIDATA_WARNING\n\
        ")
    endif()
  endif()

  # ----------------------------------------------------------------------------
  # Hash the path to create a unique ID for our targets but shorten it to the
  # first 7 characters for brevity. A unique ID is required to avoid clashes
  # with targets in other directories
  string(SHA1 PATH_HASH "${CMAKE_CURRENT_LIST_DIR}")
  string(SUBSTRING ${PATH_HASH} 0 7 PATH_HASH)
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # Grab the validate.sh script if we have one
  set(REQUIREMENTS_WITH_PATHS "")
  set(TEST_BYPRODUCTS "")

  # Add on the extra requirements
  foreach(REQUIREMENT IN LISTS TEST_FILES)
    cmake_path(ABSOLUTE_PATH REQUIREMENT BASE_DIRECTORY
               "${CMAKE_CURRENT_LIST_DIR}")
    cmake_path(GET REQUIREMENT FILENAME FILENAME)
    list(APPEND REQUIREMENTS_WITH_PATHS "${REQUIREMENT}")
    list(APPEND TEST_BYPRODUCTS "${CMAKE_CURRENT_BINARY_DIR}/${FILENAME}")
  endforeach()
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # Flag used to control whether files are symlinked instead of copied
  set(SYMLINK_TEST_FILES_INSTEAD_OF_COPY TRUE)

  # Copy/symlink test files
  if(NOT TARGET copy_${PATH_HASH})
    # Declare a copy_... target to copy/symlink the required files to the build
    # directory
    add_custom_target(copy_${PATH_HASH} ALL)

    # Add each requirement to the copy target as a file-copy command or as a
    # directory-copy command. All of these commands will be executed when the
    # copy_<path-hash> target is called
    foreach(REQUIREMENT IN LISTS REQUIREMENTS_WITH_PATHS)
      cmake_path(GET REQUIREMENT FILENAME FILENAME)
      if(SYMLINK_TEST_FILES_INSTEAD_OF_COPY)
        add_custom_command(
          TARGET copy_${PATH_HASH}
          POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E create_symlink "${REQUIREMENT}"
                  "${CMAKE_CURRENT_BINARY_DIR}/${FILENAME}")
      else()
        if(IS_DIRECTORY "${REQUIREMENT}")
          add_custom_command(
            TARGET copy_${PATH_HASH}
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_directory "${REQUIREMENT}"
                    "${CMAKE_CURRENT_BINARY_DIR}/${FILENAME}")
        else()
          add_custom_command(
            TARGET copy_${PATH_HASH}
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${REQUIREMENT}"
                    "${CMAKE_CURRENT_BINARY_DIR}/${FILENAME}")
        endif()
      endif()
    endforeach()

    # To silence the warning about a missing COMMAND key in the
    # add_custom_command(...) command below. We're going to just enable the
    # policy locally (by pushing to the policy stack then popping after)
    if(POLICY CMP0175)
      cmake_policy(PUSH)
      cmake_policy(VERSION 3.31)
      cmake_policy(SET CMP0175 OLD)
    endif()

    # Identify the files that we'll copy as by-products so that they can be
    # cleaned up by running "make clean" if the user uses Makefile Generators
    add_custom_command(
      TARGET copy_${PATH_HASH}
      POST_BUILD
      BYPRODUCTS ${TEST_BYPRODUCTS})

    if(POLICY CMP0175)
      # Remove the policy enabled above
      cmake_policy(POP)
    endif()
  endif()
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # Create a target to build the targets we're going to test with validate.sh.
  # We need to build from the top-level build directory as this is where the
  # Makefile/build.ninja file (which contains the build recipes) lives
  add_custom_target(build_targets_${PATH_HASH} ALL)

  # Add on commands to build the targets we need
  foreach(TARGET_DEPENDENCY IN LISTS DEPENDS_ON)
    if(ASSUME_PURE_CMAKE_TARGET_NAMES)
      add_dependencies(build_targets_${PATH_HASH} ${TARGET_DEPENDENCY})
    else()
      add_dependencies(build_targets_${PATH_HASH}
                       ${TARGET_DEPENDENCY}_${PATH_HASH})
    endif()
  endforeach()
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # Create a target to wipe the Validation/ directory if it exists
  if(NOT TARGET clean_validation_dir_${PATH_HASH})
    add_custom_target(
      clean_validation_dir_${PATH_HASH}
      COMMAND ${CMAKE_COMMAND} -E remove_directory
              "${CMAKE_CURRENT_BINARY_DIR}/Validation")
  endif()
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # Command-line arguments for validate.sh. We can't run fpdiff.py if we don't
  # have Python
  if(NOT Python3_FOUND)
    list(APPEND VALIDATE_SH_COMMAND "no_fpdiff")
  endif()
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # Create run_test.sh to run the validate.sh script and handle the correct
  # propagation of the exit code, as well as appending to the "global"
  # validation.log file.
  #
  # cmake-format: off
  set(TEST_SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/run_test.sh")

  # The shebang
  file(WRITE "${TEST_SCRIPT}" "#!/bin/bash\n\n")

  # Jump to the directory of the script
  file(APPEND "${TEST_SCRIPT}" "# Entering directory where this bash script lives\n")
  file(APPEND "${TEST_SCRIPT}" "DIR=\"$(cd \"$(dirname \"\$\{BASH_SOURCE\[0\]\}\")\" && pwd)\"\n")
  file(APPEND "${TEST_SCRIPT}" "cd \"\$\{DIR\}\"\n\n")

  # Run the command
  file(APPEND "${TEST_SCRIPT}" "# Run the validation command\n")
  list(JOIN VALIDATE_SH_COMMAND " " VALIDATE_SH_COMMAND_STRING)
  file(APPEND "${TEST_SCRIPT}" "${VALIDATE_SH_COMMAND_STRING}\n\n")

  # Check exit code of test command
  file(APPEND "${TEST_SCRIPT}" "# Store the exit code\n")
  file(APPEND "${TEST_SCRIPT}" "EXIT_CODE=$?\n\n")

  # Jump to the original directory
  file(APPEND "${TEST_SCRIPT}" "# Jump back to the original calling directory\n")
  file(APPEND "${TEST_SCRIPT}" "cd - > /dev/null\n\n")

  # Now exit if nonzero
  file(APPEND "${TEST_SCRIPT}" "# Stop here if we exited with a non-zero exit code\n")
  file(APPEND "${TEST_SCRIPT}" "if [ $EXIT_CODE -ne 0 ]; then\n  echo \"Test stopped with exit code $EXIT_CODE!\"\n  exit $EXIT_CODE\nfi\n\n")

  # Check for the validation.log file
  file(APPEND "${TEST_SCRIPT}" "# Stop here if there's no validation log file\n")
  file(APPEND "${TEST_SCRIPT}" "if [ ! -e \"${CMAKE_CURRENT_BINARY_DIR}/Validation/validation.log\" ]; then\n")
  file(APPEND "${TEST_SCRIPT}" "  printf '\\n%s:\\n\\t%s\\n%s\\n' 'Unable to locate validation log file:' '\"${CMAKE_CURRENT_BINARY_DIR}/Validation/validation.log\"' 'Stopping here...'\n")
  file(APPEND "${TEST_SCRIPT}" "  exit 1\n")
  file(APPEND "${TEST_SCRIPT}" "fi\n\n")

  # Append validation.log to top-level validation.log
  file(APPEND "${TEST_SCRIPT}" "# Append the validation log file to the 'global' log file\n")
  file(APPEND "${TEST_SCRIPT}" "cat \"${CMAKE_CURRENT_BINARY_DIR}/Validation/validation.log\" >> \"${CMAKE_BINARY_DIR}/validation.log\"\n")

  # Make the script executable
  file(CHMOD "${TEST_SCRIPT}" PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  # cmake-format: on
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # Define the custom target that runs the script
  add_custom_target(
    check_${PATH_HASH}
    COMMAND "${TEST_SCRIPT}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    DEPENDS copy_${PATH_HASH} build_targets_${PATH_HASH}
            clean_validation_dir_${PATH_HASH})
  # ----------------------------------------------------------------------------

  # Create a test target that depends on the check_${PATH_HASH} target. When the
  # user runs "ninja <TEST-NAME>", it will cause the test dependencies to get
  # built and the validata files/validate.sh script to get placed in the build
  # directory
  add_custom_target(${TEST_NAME})
  add_dependencies(${TEST_NAME} copy_${PATH_HASH} build_targets_${PATH_HASH}
                   clean_validation_dir_${PATH_HASH})

  # Create the test to be run by CTest. When we run the test, it will call, e.g.
  # 'ninja check_...' which will call the 'check_...' target defined above. As
  # this target depends on the copy_... and build_targets_... targets, they will
  # be called first, resulting in the required test files to be symlinked or
  # copied to the build directory, and the targets the test depends on to get
  # built. Once the data is in place and the executables have been built, the
  # check_... commands will run, which will cause the validate.sh script or
  # individual executables to be run. Finally, the test validation.log file will
  # be appended to the root build directory validation.log file.
  add_test(
    NAME ${TEST_NAME}
    COMMAND ${CMAKE_MAKE_PROGRAM} check_${PATH_HASH}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
endfunction()
# ------------------------------------------------------------------------------
