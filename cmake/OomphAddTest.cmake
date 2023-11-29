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
#                    LABELS                <string-list-of-labels>
#                    COMMAND               <command-to-run-test>
#                    [TEST_FILES           <files-required-by-test>]
#                    [SILENCE_MISSING_VALIDATA_WARNING]
#                    [NO_VALIDATE_SH])
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
#                    TEST_FILES validate.sh validata my_extra_data_file.dat
#                    LABELS "poisson;one_d_poisson")
#
#     oomph_add_test(TEST_NAME poisson.one_d_poisson
#                    DEPENDS_ON one_d_poisson
#                    COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_VARIABLENP_RUN_COMMAND}
#                    TEST_FILES validate.sh validata my_extra_data_file.dat
#                    LABELS "poisson;one_d_poisson")
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
  set(FLAGS NO_VALIDATE_SH SILENCE_MISSING_VALIDATA_WARNING)
  set(SINGLE_VALUE_ARGS TEST_NAME)
  set(MULTI_VALUE_ARGS TEST_FILES LABELS DEPENDS_ON COMMAND)

  # Process the arguments passed in
  include(CMakeParseArguments)
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}"
                        "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

  # Redefine the variables in this scope without a prefix for clarity
  set(NO_VALIDATE_SH ${${PREFIX}_NO_VALIDATE_SH})
  set(SILENCE_MISSING_VALIDATA_WARNING
      ${${PREFIX}_SILENCE_MISSING_VALIDATA_WARNING})
  set(TEST_NAME ${${PREFIX}_TEST_NAME})
  set(COMMAND ${${PREFIX}_COMMAND})
  set(DEPENDS_ON ${${PREFIX}_DEPENDS_ON})
  set(TEST_FILES ${${PREFIX}_TEST_FILES})
  set(LABELS ${${PREFIX}_LABELS})

  # Make sure the arguments are valid
  if(NOT TEST_NAME)
    message(FATAL_ERROR "No TEST_NAME argument supplied.")
  elseif(NOT DEPENDS_ON)
    message(VERBOSE "No DEPENDS_ON argument supplied. Can't create a test.")
    return()
  elseif(NOT LABELS)
    message(WARNING "No LABELS supplied. These are helpful for running CTest\
    with subsets of the tests. We recommend that you set this!")
  endif()

  find_program(BASH_PROGRAM bash)
  if(NOT BASH_PROGRAM)
    message(STATUS "You don't have 'bash', so I can't construct any tests!")
  endif()

  # ----------------------------------------------------------------------------
  # Grab the validate.sh script if we have one
  set(REQUIREMENTS_WITH_PATHS)
  set(TEST_BYPRODUCTS)

  # Add on the extra requirements
  foreach(REQUIREMENT IN LISTS TEST_FILES)
    set(FILENAME ${REQUIREMENT})
    if(IS_ABSOLUTE ${REQUIREMENT})
      list(APPEND REQUIREMENTS_WITH_PATHS ${REQUIREMENT})
      cmake_path(GET REQUIREMENT FILENAME FILENAME)
    else()
      list(APPEND REQUIREMENTS_WITH_PATHS
           "${CMAKE_CURRENT_LIST_DIR}/${REQUIREMENT}")
    endif()
    list(APPEND TEST_BYPRODUCTS "${CMAKE_CURRENT_BINARY_DIR}/${FILENAME}")
  endforeach()

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

  # ----------------------------------------------------------------------------
  # Hash the path to create a unique ID for our targets but shorten it to the
  # first 7 characters for brevity. A unique ID is required to avoid clashes
  # with targets in other directories
  string(SHA1 PATH_HASH "${CMAKE_CURRENT_LIST_DIR}")
  string(SUBSTRING ${PATH_HASH} 0 7 PATH_HASH)
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # Declare a copy_... target to copy the required files to the build directory
  if(NOT TARGET copy_${PATH_HASH})
    add_custom_target(copy_${PATH_HASH} ALL
                      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
  endif()

  # Flag used to control whether files are symlinked instead of copied; keeping
  # the option to copy files around just in case we need it later on (but I
  # doubt it)
  set(SYMLINK_TEST_FILES_INSTEAD_OF_COPY TRUE)

  # Add each requirement to the copy target as a file-copy command or as a
  # directory-copy command. All of these commands will be executed when the
  # copy_<path-hash> target is called
  foreach(REQUIREMENT IN LISTS REQUIREMENTS_WITH_PATHS)
    if(SYMLINK_TEST_FILES_INSTEAD_OF_COPY)
      add_custom_command(
        TARGET copy_${PATH_HASH}
        POST_BUILD
        COMMAND ln -sf "${REQUIREMENT}" "${CMAKE_CURRENT_BINARY_DIR}")
    else()
      if(IS_DIRECTORY "${REQUIREMENT}")
        add_custom_command(
          TARGET copy_${PATH_HASH}
          POST_BUILD
          COMMAND cp -ur "${REQUIREMENT}" "${CMAKE_CURRENT_BINARY_DIR}")
      else()
        add_custom_command(
          TARGET copy_${PATH_HASH}
          POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy_if_different "${REQUIREMENT}"
                  "${CMAKE_CURRENT_BINARY_DIR}")
      endif()
    endif()
  endforeach()

  # Identify the files that we'll copy as by-products so that they can be
  # cleaned up by running "make clean" if the user uses Makefile Generators
  add_custom_command(
    TARGET copy_${PATH_HASH}
    POST_BUILD
    BYPRODUCTS ${TEST_BYPRODUCTS})
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # Create a target to build the targets we're going to test with validate.sh.
  # We need to build from the top-level build directory as this is where the
  # Makefile/build.ninja file (which contains the build recipes) lives
  add_custom_target(build_targets_${PATH_HASH}
                    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")

  # Add on commands to build the targets we need
  foreach(TARGET_DEPENDENCY IN LISTS DEPENDS_ON)
    add_custom_command(
      TARGET build_targets_${PATH_HASH}
      POST_BUILD
      COMMAND ${CMAKE_MAKE_PROGRAM} ${TARGET_DEPENDENCY}_${PATH_HASH}
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
  endforeach()
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # Create a target to wipe the Validation/ directory if it exists
  if(NOT TARGET clean_validation_dir_${PATH_HASH})
    add_custom_target(clean_validation_dir_${PATH_HASH}
                      COMMAND rm -rf "${CMAKE_CURRENT_BINARY_DIR}/Validation")
  endif()
  # ----------------------------------------------------------------------------

  # Command-line arguments for validate.sh. We can't run fpdiff.py if we don't
  # have Python
  if(NOT Python3_FOUND)
    list(APPEND COMMAND "no_fpdiff")
  endif()

  # TODO: Try rewriting with POST_BUILD actions:
  # https://cmake.org/cmake/help/latest/command/add_custom_command.html?highlight=add_custom_command#build-events
  # It will simply append the new command (with &&; check build.ninja) and not
  # overwrite the previous ones

  # TODO: Try simplifying things with
  # https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html#using-ctest-to-drive-complex-tests

  # FIXME: If we move away from validate.sh scripts, we need to run executables
  # with our own mpirun command if needed

  # ----------------------------------------------------------------------------
  # Run the dependencies to copy the test data, build the (sub)project(s)
  # targets then run the validate.sh script and pass the location of the
  # oomph-lib root directory so they have access to the scripts they require,
  # like fpdiff.py and validate_ok_count, and so it knows where to place the
  # validation.log output. The VERBATIM argument is absolutely necessary here to
  # ensure that the "mpirun ..." commands are correctly escaped.
  if(NOT NO_VALIDATE_SH)
    add_custom_target(
      check_${PATH_HASH}
      COMMAND ${BASH_PROGRAM} ${COMMAND}
      # Check for the validation.log file. Stop here if we can't
      COMMAND
        ${BASH_PROGRAM} -c
        "test -e \"${CMAKE_CURRENT_BINARY_DIR}/Validation/validation.log\" || ( \
          printf \"\\nUnable to locate validation log file:\\n\\n\\t${CMAKE_CURRENT_BINARY_DIR}/Validation/validation.log\\n\\nStopping here...\\n\\n\" && \
          exit 0 )"
      # Append the validation.log to the top-level validation.log
      COMMAND cat "${CMAKE_CURRENT_BINARY_DIR}/Validation/validation.log" >>
              "${CMAKE_BINARY_DIR}/validation.log"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
      DEPENDS copy_${PATH_HASH} build_targets_${PATH_HASH}
              clean_validation_dir_${PATH_HASH}
      VERBATIM)
  else()
    # Run all of the dependent targets and then append the validation.log output
    # to the global one.
    list(JOIN DEPENDS_ON " ./" RUN_DEPENDENCIES_STRING)
    add_custom_target(
      check_${PATH_HASH}
      # Run each executable
      COMMAND ${BASH_PROGRAM} -c ./${RUN_DEPENDENCIES_STRING}
      # Check for the validation.log file
      COMMAND
        ${BASH_PROGRAM} -c
        "test -e \"${CMAKE_CURRENT_BINARY_DIR}/Validation/validation.log\" || ( \
          printf \"\\nUnable to locate file:\\n\\n\\t${CMAKE_CURRENT_BINARY_DIR}/Validation/validation.log\\n\\nStopping here...\\n\\n\" && \
          exit 1 )"
      # Append validation.log to top-level validation.log
      COMMAND cat "${CMAKE_CURRENT_BINARY_DIR}/Validation/validation.log" >>
              "${CMAKE_BINARY_DIR}/validation.log"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
      DEPENDS copy_${PATH_HASH} build_targets_${PATH_HASH}
              clean_validation_dir_${PATH_HASH}
      VERBATIM)
  endif()
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
  # individual executables to be run (depending on whether 'NO_VALIDATE_SH' was
  # specified). Finally, the output validation.log file will be appended to the
  # global validation.log file in the oomph-lib root directory.
  add_test(
    NAME ${TEST_NAME}
    COMMAND ${CMAKE_MAKE_PROGRAM} check_${PATH_HASH}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")

  # Add the user-provided test labels to the test so that it can be used by
  # CTest to filter tests
  set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${LABELS}")
endfunction()
# ------------------------------------------------------------------------------
