# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Supports the creation of an oomph-lib demo driver (i.e. integration) test
# that does not rely on a validate.sh script.
#
# These tests are distinct from regular unit tests (which are far more
# lightweight and do not depend on data to be tested against).
#
# USAGE:
# ------
#  oomph_add_pure_cpp_test(TEST_NAME           <executable-name>
#                          DEPENDS_ON          <executable/target-required-by-test>
#                          LOG_FILE            <output-log-file>
#                          LABELS              <labels>
#                          [RUN_WITH           <command-to-run-executable>]
#                          [TEST_ARGS          <arguments-for-executable>]
#                          [SILENCE_MISSING_VALIDATA_WARNING])
#
# The argument to DEPENDS_ON must be an already-defined executable or target
# (i.e. defined via add_executable() or oomph_add_executable()). The LOG_FILE is the
# relative path to the validation log file under the Validation/ directory. We
# require a unique LOG_FILE argument for each test in the directory as multiple
# targets could write to the same validation.log file at the same time causing
# the resulting output to get mangled.
#
# EXAMPLE:
# --------
# (1) Run with --disable_melting --validate flags
#
#   # Define executable
#   oomph_add_executable(NAME melt
#                        SOURCES melt.cc heat_transfer_and_melt_elements.h
#                        LIBRARIES oomph::constitutive oomph::solid oomph::unsteady_heat
#                                  oomph::meshes oomph::generic)
#
#   # Define test
#   oomph_add_pure_cpp_test(TEST_NAME heat_transfer_and_melting.melt
#                           DEPENDS_ON melt
#                           TEST_ARGS --disable_melting --validate
#                           LABELS heat_transfer_and_melting solid unsteady_heat)
#
# (2) Run with MPI
#
#   # Define executable
#   oomph_add_executable(NAME melt
#                        SOURCES melt.cc heat_transfer_and_melt_elements.h
#                        LIBRARIES oomph::constitutive oomph::solid oomph::unsteady_heat
#                                  oomph::meshes oomph::generic)
#
#   # Define test
#   oomph_add_pure_cpp_test(TEST_NAME heat_transfer_and_melting.melt
#                           DEPENDS_ON melt
#                           RUN_WITH mpirun -np 2
#                           LABELS heat_transfer_and_melting solid unsteady_heat)
#
# NOTE: You can also specify the MPI run command by setting
#
#   # Define test
#   oomph_add_pure_cpp_test(TEST_NAME heat_transfer_and_melting.melt
#                           DEPENDS_ON melt
#                           RUN_WITH ${OOMPH_MPI_RUN_COMMAND}
#                           LABELS heat_transfer_and_melting solid unsteady_heat)
# =============================================================================
# cmake-format: on
include_guard()

# ------------------------------------------------------------------------------
function(oomph_add_pure_cpp_test)
  # Define the supported set of keywords
  set(PREFIX ARG)
  set(FLAGS SILENCE_MISSING_VALIDATA_WARNING)
  set(SINGLE_VALUE_ARGS TEST_NAME DEPENDS_ON LOG_FILE)
  set(MULTI_VALUE_ARGS LABELS RUN_WITH TEST_ARGS)

  # Process the arguments passed in
  include(CMakeParseArguments)
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}"
                        "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

  # Redefine the variables in this scope without a prefix for clarity
  set(TEST_NAME ${${PREFIX}_TEST_NAME})
  set(DEPENDS_ON ${${PREFIX}_DEPENDS_ON})
  set(LOG_FILE ${${PREFIX}_LOG_FILE})
  set(RUN_WITH ${${PREFIX}_RUN_WITH})
  set(TEST_ARGS ${${PREFIX}_TEST_ARGS})
  set(LABELS ${${PREFIX}_LABELS})
  set(SILENCE_MISSING_VALIDATA_WARNING
      ${${PREFIX}_SILENCE_MISSING_VALIDATA_WARNING})

  # Make sure the arguments are valid
  if(NOT TEST_NAME)
    message(FATAL_ERROR "No TEST_NAME argument supplied.")
  elseif(NOT DEPENDS_ON)
    message(FATAL_ERROR "No DEPENDS_ON argument supplied.")
  elseif(NOT LOG_FILE)
    message(FATAL_ERROR "No LOG_FILE argument supplied.")
  elseif(NOT LABELS)
    message(WARNING "No LABELS supplied. These are helpful for running CTest\
    with subsets of the tests. We recommend that you set this!")
  endif()

  # Check that we can run the test
  find_program(BASH_PROGRAM bash)
  if(NOT BASH_PROGRAM)
    message(
      STATUS "You don't have 'bash', so I can't construct any tests. Sorry!")
  endif()

  # Hash the path to create a unique ID for our targets but shorten it to the
  # first 7 characters for brevity. A unique ID is required to avoid clashes
  # with targets in other directories
  string(SHA1 PATH_HASH "${CMAKE_CURRENT_LIST_DIR}")
  string(SUBSTRING ${PATH_HASH} 0 7 PATH_HASH)

  # Make sure we've been given a proper target
  if(NOT TARGET ${DEPENDS_ON}_${PATH_HASH})
    message(FATAL_ERROR "Argument ${DEPENDS_ON} is not a target.")
  endif()

  # ----------------------------------------------------------------------------
  # Create a target to wipe the Validation/ directory if it exists
  if(NOT TARGET clean_validation_dir_${PATH_HASH})
    add_custom_target(clean_validation_dir_${PATH_HASH}
                      COMMAND rm -rf "${CMAKE_CURRENT_BINARY_DIR}/Validation")
  endif()
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # Indicate that if we build ${DEPENDS_ON}, we should copy over the required
  # files
  #
  # add_dependencies(${DEPENDS_ON}_${PATH_HASH} copy_${PATH_HASH})
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # Run the dependencies to copy the test data, build the (sub)project(s)
  # targets, run all of the dependent targets then append the validation.log
  # output to the global one.. The VERBATIM argument is necessary here to ensure
  # that the "mpirun ..." commands are correctly escaped
  set(RUN_COMMAND ${RUN_WITH} ./${DEPENDS_ON} ${TEST_ARGS})
  list(JOIN RUN_COMMAND " " RUN_COMMAND)
  add_custom_target(
    check_${TEST_NAME}_${PATH_HASH}
    # Run the executable
    COMMAND ${BASH_PROGRAM} -c ${RUN_COMMAND}
    # Check for the validation.log file
    COMMAND
      ${BASH_PROGRAM} -c
      "test -e \"${CMAKE_CURRENT_BINARY_DIR}/Validation/${LOG_FILE}\" || ( \
          printf \"\\nUnable to locate file:\\n\\n\\t${CMAKE_CURRENT_BINARY_DIR}/Validation/${LOG_FILE}\\n\\nStopping here...\\n\\n\" && \
          exit 1 )"
    # Append validation.log to top-level validation.log
    COMMAND cat "${CMAKE_CURRENT_BINARY_DIR}/Validation/${LOG_FILE}" >>
            "${CMAKE_BINARY_DIR}/validation.log"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    DEPENDS ${DEPENDS_ON}_${PATH_HASH} copy_${PATH_HASH}
            clean_validation_dir_${PATH_HASH}
    VERBATIM)
  # ----------------------------------------------------------------------------

  # Create the test to be run by CTest. When we run the test, it will call, e.g.
  # 'ninja check_...' which will call the 'check_...' target defined above. As
  # this target depends on the copy_... and build_targets_... targets, they will
  # be called first, resulting in the required test files to be symlinked or
  # copied to the build directory, and the targets the test depends on to get
  # built. Once the data is in place and the executables have been built, the
  # check_... commands will run, which will cause the individual executables to
  # be run. Finally, the output validation.log file will be appended to the
  # global validation.log file in the oomph-lib root directory.
  add_test(
    NAME ${TEST_NAME}
    COMMAND ${CMAKE_MAKE_PROGRAM} check_${TEST_NAME}_${PATH_HASH}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")

  # Add the user-provided test labels to the test so that it can be used by
  # CTest to filter tests
  set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${LABELS}")
endfunction()
# ------------------------------------------------------------------------------
