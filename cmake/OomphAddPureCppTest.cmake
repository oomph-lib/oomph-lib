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
#                           TEST_ARGS --disable_melting --validate)
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
#                           RUN_WITH mpirun -np 2)
#
# NOTE: You can also specify the MPI run command by setting
#
#   # Define test
#   oomph_add_pure_cpp_test(TEST_NAME heat_transfer_and_melting.melt
#                           DEPENDS_ON melt
#                           RUN_WITH ${OOMPH_MPI_RUN_COMMAND})
# =============================================================================
# cmake-format: on
include_guard()

# ------------------------------------------------------------------------------
function(oomph_add_pure_cpp_test)
  # Define the supported set of keywords
  set(PREFIX ARG)
  set(FLAGS SILENCE_MISSING_VALIDATA_WARNING)
  set(SINGLE_VALUE_ARGS TEST_NAME DEPENDS_ON LOG_FILE)
  set(MULTI_VALUE_ARGS RUN_WITH TEST_ARGS)

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
  set(SILENCE_MISSING_VALIDATA_WARNING
      ${${PREFIX}_SILENCE_MISSING_VALIDATA_WARNING})

  # Make sure the arguments are valid
  if(NOT TEST_NAME)
    message(FATAL_ERROR "No TEST_NAME argument supplied.")
  elseif(NOT DEPENDS_ON)
    message(FATAL_ERROR "No DEPENDS_ON argument supplied.")
  elseif(NOT LOG_FILE)
    message(FATAL_ERROR "No LOG_FILE argument supplied.")
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
  # output to the global one.
  set(RUN_COMMAND ${RUN_WITH} ./${DEPENDS_ON} ${TEST_ARGS})
  list(JOIN RUN_COMMAND " " RUN_COMMAND)
  set(TEST_SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/run_test_${TEST_NAME}.sh")

  # cmake-format: off
  file(WRITE "${TEST_SCRIPT}" "#!/bin/bash\n\n")
  file(APPEND "${TEST_SCRIPT}" "# auto-generated script (created by cmake/OomphAddPureCppTest.cmake)\n")
  file(APPEND "${TEST_SCRIPT}" "# Do not edit!\n\n")

  file(APPEND "${TEST_SCRIPT}" "# Entering directory where this bash script lives\n")
  file(APPEND "${TEST_SCRIPT}" "DIR=\"$(cd \"$(dirname \"\$\{BASH_SOURCE\[0\]\}\")\" && pwd)\"\n")
  file(APPEND "${TEST_SCRIPT}" "cd \"\$\{DIR\}\"\n\n")

  file(APPEND "${TEST_SCRIPT}" "# Run the validation command\n")
  file(APPEND "${TEST_SCRIPT}" "${RUN_COMMAND}\n\n")

  file(APPEND "${TEST_SCRIPT}" "# Store the exit code\n")
  file(APPEND "${TEST_SCRIPT}" "EXIT_CODE=$?\n\n")

  file(APPEND "${TEST_SCRIPT}" "# Define validation-log paths\n")
  file(APPEND "${TEST_SCRIPT}" "VALIDATION_DIR=\"${CMAKE_CURRENT_BINARY_DIR}/Validation\"\n")
  file(APPEND "${TEST_SCRIPT}" "VALIDATION_LOG=\"${CMAKE_CURRENT_BINARY_DIR}/Validation/${LOG_FILE}\"\n")
  file(APPEND "${TEST_SCRIPT}" "GLOBAL_VALIDATION_LOG=\"${CMAKE_BINARY_DIR}/validation.log\"\n\n")

  file(APPEND "${TEST_SCRIPT}" "# Append the validation log file to the 'global' log file\n")
  file(APPEND "${TEST_SCRIPT}" "if [ -e \"$VALIDATION_LOG\" ]; then\n")
  file(APPEND "${TEST_SCRIPT}" "  cat \"$VALIDATION_LOG\" >> \"$GLOBAL_VALIDATION_LOG\" || exit 1\n")
  file(APPEND "${TEST_SCRIPT}" "fi\n\n")

  file(APPEND "${TEST_SCRIPT}" "# Stop here if we exited with a non-zero exit code\n")
  file(APPEND "${TEST_SCRIPT}" "if [ $EXIT_CODE -ne 0 ]; then\n  echo \"Test stopped with exit code $EXIT_CODE!\"\n  exit $EXIT_CODE\nfi\n\n")

  file(APPEND "${TEST_SCRIPT}" "# Stop here if there's no validation log file\n")
  file(APPEND "${TEST_SCRIPT}" "if [ ! -e \"$VALIDATION_LOG\" ]; then\n")
  file(APPEND "${TEST_SCRIPT}" "  printf '\\n%s:\\n\\t%s\\n%s\\n' 'Unable to locate validation log file:' \"$VALIDATION_LOG\" 'Stopping here...'\n")
  file(APPEND "${TEST_SCRIPT}" "  exit 1\n")
  file(APPEND "${TEST_SCRIPT}" "fi\n\n")

  file(APPEND "${TEST_SCRIPT}" "# Optionally remove successful validation output to save disk space\n")
  file(APPEND "${TEST_SCRIPT}" "if [ \"$OOMPH_DELETE_VALIDATION_DIRECTORY_AFTER_SUCCESSFUL_TEST\" = \"yes\" ]; then\n")
  file(APPEND "${TEST_SCRIPT}" "  printf '\\n%s\\n\\t%s\\n' 'Deleting successful Validation directory to save disk space:' \"$VALIDATION_DIR\" >> \"$GLOBAL_VALIDATION_LOG\"\n")
  file(APPEND "${TEST_SCRIPT}" "  rm -rf \"$VALIDATION_DIR\"\n")
  file(APPEND "${TEST_SCRIPT}" "fi\n")

  file(CHMOD "${TEST_SCRIPT}" PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  # cmake-format: on

  add_custom_target(
    check_${TEST_NAME}_${PATH_HASH}
    COMMAND "${TEST_SCRIPT}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    DEPENDS ${DEPENDS_ON}_${PATH_HASH} copy_${PATH_HASH}
            clean_validation_dir_${PATH_HASH}
    VERBATIM)
  # ----------------------------------------------------------------------------

  # Create the test to be run by CTest. When we run the test, it will call
  # 'cmake --build ... --target check_...' which will call the check_... target
  # defined above. We pin the inner build to a single job because parallelism
  # belongs at the ctest level (ctest -j N runs N tests concurrently); without
  # this, each of those N tests would spawn its own full-core sub-build and
  # over-subscribe the machine. In the common case the inner build is a no-op
  # because the main build has already produced everything it needs.
  add_test(
    NAME ${TEST_NAME}
    COMMAND ${CMAKE_COMMAND} --build "${CMAKE_BINARY_DIR}" --target
            check_${TEST_NAME}_${PATH_HASH} --parallel 1
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
endfunction()
# ------------------------------------------------------------------------------
