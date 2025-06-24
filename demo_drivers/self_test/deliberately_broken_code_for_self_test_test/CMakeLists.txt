# ------------------------------------------------------------------------------
# OOMPH auto-inserted: warn once if user configures inside an existing build.
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt"
   AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(
    WARNING
      "You're re-configuring an existing build tree. For a clean configuration, "
      "delete the whole build folder and start with a fresh directory.")
endif()
# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE
        "Entered deliberately_broken_code_for_self_test_test subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(deliberately_broken_code_for_self_test_test C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

# Make sure we get warned about unused variables
oomph_add_executable(NAME code_with_warning SOURCES code_with_warning.cc)
set(CODE_WITH_WARNING_TARGET)
oomph_get_target_name(code_with_warning CODE_WITH_WARNING_TARGET)
target_compile_options(${CODE_WITH_WARNING_TARGET} PRIVATE -Wunused-variable)

# Define the executable the normal way but just make sure we avoid building the
# broken target normally
add_executable(deliberately_broken_code deliberately_broken_code.cc)
set_target_properties(
  deliberately_broken_code PROPERTIES EXCLUDE_FROM_ALL TRUE
                                      EXCLUDE_FROM_DEFAULT_BUILD TRUE)

# ---------------------------------- TESTING ----------------------------------

# Define this test with oomph_add_test() so we get help carrying the validate.sh
# and validata/ stuff over to the build directory when we run the tests
oomph_add_test(
  TEST_NAME
    self_test.deliberately_broken_code_for_self_test_test.code_with_warning
  DEPENDS_ON code_with_warning
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

# Just try to build the driver code
add_test(
  NAME self_test.deliberately_broken_code_for_self_test_test.deliberately_broken_code_for_self_test_test
  COMMAND ${CMAKE_COMMAND} --build . --target deliberately_broken_code
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")

# Expect the deliberately-broken test to fail (i.e. cmake --build should return
# a non-zero value)
set_tests_properties(
  self_test.deliberately_broken_code_for_self_test_test.deliberately_broken_code_for_self_test_test
  PROPERTIES WILL_FAIL TRUE)

message(VERBOSE
        "Leaving deliberately_broken_code_for_self_test_test subdirectory")
# ------------------------------------------------------------------------------
