# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Adds the compiler flags needed for collecting code coverage results.
#
# USAGE:
# ------
#  include(OomphEnableCodeCoverage)
#  oomph_enable_code_coverage()
#
# =============================================================================
# cmake-format: on
include_guard()

# ------------------------------------------------------------------------------
macro(oomph_enable_code_coverage)
  find_program(GCOV NAMES gcov)
  set(CTEST_COVERAGE_COMMAND "${GCOV}")

  if(CMAKE_CXX_COMPILER_ID MATCHES "(GNU|Clang)")
    include(CheckCXXCompilerFlag)
    message(STATUS "Adding code coverage flags")

    # Assign the required code coverage flags
    if(CMAKE_COMPILER_IS_GNUCXX)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
    else()
      set(CMAKE_C_FLAGS
          "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage -fprofile-instr-generate -fcoverage-mapping"
      )
      set(CMAKE_CXX_FLAGS
          "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage -fno-elide-constructors -fprofile-instr-generate -fcoverage-mapping"
      )
    endif()
  endif()
endmacro()
# ------------------------------------------------------------------------------
