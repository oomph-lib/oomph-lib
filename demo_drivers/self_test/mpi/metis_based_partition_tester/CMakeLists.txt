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
message(VERBOSE "Entered metis_based_partition_tester subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(metis_based_partition_tester C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME adaptive_driven_cavity
  SOURCES adaptive_driven_cavity.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

# hierher: Puneet: I'm confused about the test name; I thought this was set
# automatically. If not, should we explain our conventions here? Can I do
# anything I want? Same with labels. Is there a master list? Puneet's answer:
# Labels are different, they're just like "tags"; just a way of picking out a
# bunch of tests. You can tell CTest, for example, to run tests that only have a
# specific flag (like "mpi"). I don't have a master list for this (or should I
# say ahem main list)
#
# Nope; you set the test name manually. In the demo drivers, I set it to the
# path of the specific self-test, minus the /home/.../demo_driver/ prefix and
# replacing the slashes with periods. E.g. the test in
#
# /home/.../oomph-lib/demo_drivers/poisson/one_d_poisson/
#
# would have
#
# TEST_NAME poisson.one_d_poisson
#
# The same is true with self_test/ so we should remove the self_test. prefix.
# Since we never build the self_tests/ and demo_drivers/ directory together, we
# don't need to worry about clashes in test names (so we don't need to add a
# demo_driver. or self_test. prefix to the test name)

oomph_add_test(
  TEST_NAME self_test.mpi.metis_based_partition_tester
  DEPENDS_ON adaptive_driven_cavity
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata)

# Don't run this test alongside any other tests hierher Puneet: how does ensure
# that it's not run alongside other tests? Does the first arg have to match the
# test name? Puneet's answer: Quoting the Professional CMake book:
#
# By default, ctest will assume each test consumes one CPU. For test cases that
# use more than one CPU, their PROCESSORS test property can be set to indicate
# how many CPUs they are expected to use. ctest will then use that value when
# determining whether enough CPU resources are free before starting the test. If
# PROCESSORS is set to a value higher than the job limit, ctest will behave as
# though it was set to the job limit when determining whether the test can be
# started.
#
# and yes -- the first value needs to match the argument to TEST_NAME!

set_tests_properties(self_test.mpi.metis_based_partition_tester
                     PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving metis_based_partition_tester subdirectory")
# ------------------------------------------------------------------------------
