# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# This file provides CMake logic to enable sanitizer support (AddressSanitizer,
# UndefinedBehaviorSanitizer, etc.) in Debug builds when the user sets the
# OOMPH_ENABLE_SANITISER_SUPPORT option to ON. The sanitizers help detect memory
# errors such as out-of-bounds writes, use-after-free, and other undefined
# behaviors at runtime.
#
# USAGE:
# ------
# 1) Include this file in your top-level CMakeLists.txt:
#
#    include(AddSanitizerSupport)
#
# 2) Optionally enable sanitizers (only in Debug mode):
#
#    cmake -DOOMPH_ENABLE_SANITISER_SUPPORT=ON [other cmake args...]
#
# 3) Then call the provided function on any target you wish to instrument:
#
#    oomph_enable_sanitizers(<target-name>)
#
# EXAMPLE:
# --------
#
#   include(AddSanitizerSupport)
#   add_library(some_oomph_library ...)
#   oomph_enable_sanitizers(some_oomph_library)
#
# As a result, any executables linking to the oomph-lib library will also
# inherit these flags.
#
# IMPORTANT:
# ----------
# Sanitisers will only be enabled if:
#   1. OOMPH_ENABLE_SANITISER_SUPPORT is ON,
#   2. CMAKE_BUILD_TYPE is "Debug",
#   3. Compiler is GCC or Clang,
#   4. and the compiler supports the requested flags.
# =============================================================================
include_guard()

# ------------------------------------------------------------------------------
# Sanity checks

# Did they actually request sanitizer support?
if(NOT OOMPH_ENABLE_SANITISER_SUPPORT)
  return()
endif()

# Only enable sanitizers with a Debug build
if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(
    STATUS "Sanitizer support enabled, but build type is not Debug. Skipping..."
  )
  return()
endif()
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Check if the compiler supports the sanitizer flags

# Create a temporary test source file
set(OOMPH_SANITIZER_TEST_SRC "${CMAKE_BINARY_DIR}/test_sanitizers/test_sanitizer.cpp")
file(WRITE "${OOMPH_SANITIZER_TEST_SRC}" "int main() { return 0; }")

# Function to check compiler support using try_compile(...). Note, we're using
# this over check_cxx_compiler_flag(...) as try_compile(...) does better at
# avoiding false negatives
function(oomph_check_sanitizer_support FLAG RESULT_VAR)
  try_compile(
    ${RESULT_VAR}
    "${CMAKE_BINARY_DIR}/test_sanitizers/build"
    "${CMAKE_BINARY_DIR}/test_sanitizers"
    SOURCES "${OOMPH_SANITIZER_TEST_SRC}"
    CMAKE_FLAGS "-DCMAKE_CXX_FLAGS=${FLAG}")
endfunction()

# Check support for various sanitizers
oomph_check_sanitizer_support("-fsanitize=address" OOMPH_HAS_FSANITIZE_ADDRESS)
oomph_check_sanitizer_support("-fsanitize=undefined" OOMPH_HAS_FSANITIZE_UNDEFINED)

# If we are on AppleClang + Darwin + arm64, skip -fsanitize=leak; it is included in
# the address sanitiser
if (APPLE AND (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang") AND (CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64)$"))
  message(STATUS "Detected arm64 AppleClang (macOS); explicit -fsanitize=leak is not supported.")
  set(OOMPH_HAS_FSANITIZE_LEAK FALSE)
else()
  oomph_check_sanitizer_support("-fsanitize=leak" OOMPH_HAS_FSANITIZE_LEAK)
endif()

# Initialize sanitizer flags
set(OOMPH_SANITIZER_FLAGS "")

# AddressSanitizer (ASan)
if(OOMPH_HAS_FSANITIZE_ADDRESS)
  message(STATUS "${CMAKE_CXX_COMPILER_ID} supports AddressSanitizer (ASan)")
  list(APPEND OOMPH_SANITIZER_COMPILE_FLAGS "$<$<COMPILE_LANGUAGE:CXX>:-fsanitize=address>")
  list(APPEND OOMPH_SANITIZER_LINK_FLAGS "$<$<LINK_LANGUAGE:CXX>:-fsanitize=address>")
else()
  message(WARNING "${CMAKE_CXX_COMPILER_ID} does NOT support AddressSanitizer (ASan) on this platform!")
endif()

# UndefinedBehaviorSanitizer (UBSan)
if(OOMPH_HAS_FSANITIZE_UNDEFINED)
  message(STATUS "${CMAKE_CXX_COMPILER_ID} supports UndefinedBehaviorSanitizer (UBSan)")
  list(APPEND OOMPH_SANITIZER_COMPILE_FLAGS "$<$<COMPILE_LANGUAGE:CXX>:-fsanitize=undefined>")
  list(APPEND OOMPH_SANITIZER_LINK_FLAGS "$<$<LINK_LANGUAGE:CXX>:-fsanitize=undefined>")
else()
  message(WARNING "${CMAKE_CXX_COMPILER_ID} does NOT support UndefinedBehaviorSanitizer (UBSan) on this platform!")
endif()

# LeakSanitizer (LSan)
if(OOMPH_HAS_FSANITIZE_LEAK)
  message(STATUS "${CMAKE_CXX_COMPILER_ID} supports LeakSanitizer (LSan)")
  list(APPEND OOMPH_SANITIZER_COMPILE_FLAGS "$<$<COMPILE_LANGUAGE:CXX>:-fsanitize=leak>")
  list(APPEND OOMPH_SANITIZER_LINK_FLAGS "$<$<LINK_LANGUAGE:CXX>:-fsanitize=leak>")
else()
  message(STATUS "${CMAKE_CXX_COMPILER_ID} does NOT support LeakSanitizer (LSan) on this platform!")
endif()

# If no sanitizer flags are supported by this compiler do nothing
if(NOT (OOMPH_SANITIZER_COMPILE_FLAGS OR OOMPH_SANITIZER_LINK_FLAGS))
  message(
    FATAL_ERROR
      "Sanitizers requested but -fsanitize flags not supported. "
      "Disable them by reconfiguring with -DOOMPH_ENABLE_SANITISER_SUPPORT=OFF."
  )
endif()

# Always add -fno-omit-frame-pointer for better backtraces (in C++ compile/link contexts)
message(STATUS "Enabling -fno-omit-frame-pointer for more useful sanitizer backtraces.")
list(APPEND OOMPH_SANITIZER_COMPILE_FLAGS "$<$<COMPILE_LANGUAGE:CXX>:-fno-omit-frame-pointer>")
list(APPEND OOMPH_SANITIZER_LINK_FLAGS "$<$<LINK_LANGUAGE:CXX>:-fno-omit-frame-pointer>")

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
function(oomph_enable_sanitizers target_name)
  # Make sure the target actually exists
  if(NOT TARGET ${target_name})
    message(FATAL_ERROR "Target '${target_name}' does not exist!")
  endif()

  # Query the target type: EXECUTABLE, SHARED_LIBRARY, INTERFACE_LIBRARY, etc.
  get_property(
    TARGET_TYPE
    TARGET ${target_name}
    PROPERTY TYPE)

  # Decide which usage requirement to apply based on the target type
  if(TARGET_TYPE STREQUAL "INTERFACE_LIBRARY")
    # For interface libraries, use INTERFACE
    target_compile_options(${target_name} INTERFACE ${OOMPH_SANITIZER_COMPILE_FLAGS})
    target_link_options(${target_name} INTERFACE ${OOMPH_SANITIZER_LINK_FLAGS})
  else()
    # For normal libraries, executables, etc., use PUBLIC
    target_compile_options(${target_name} PUBLIC ${OOMPH_SANITIZER_COMPILE_FLAGS})
    target_link_options(${target_name} PUBLIC ${OOMPH_SANITIZER_LINK_FLAGS})
  endif()
endfunction()
# ------------------------------------------------------------------------------
# cmake-format: on
