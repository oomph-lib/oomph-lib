# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Provides linting for C++ code. Not finished off. Need to confine its use to
# oomph-lib sources only and not external sources/distributions too.
#
# WORK IN PROGRESS.
#
# For extra additions like cppcheck and IWYU, see:
# https://github.com/lefticus/cpp_starter_project/blob/master/cmake/StaticAnalyzers.cmake
#
# USAGE:
# ------
#     include(OomphStaticAnalysers)
#
# =============================================================================
# cmake-format: on

if(OOMPH_ENABLE_CLANG_TIDY)
  set(OOMPH_CLANG_TIDY_CHECKS)

  # ~~~
  # CHECKS TO RUN:
  # --------------
  # * modernize-use-nullptr
  # * modernize-use-override
  # * modernize-use-auto
  # * modernize-loop-convert
  # * modernize-use-using
  # * modernize-replace-auto-ptr
  # * modernize-use-equals-delete
  # * modernize-use-bool-literals
  # * modernize-make-shared
  # * modernize-make-unique
  # * modernize-replace-random-shuffle
  # * modernize-deprecated-headers
  # * modernize-avoid-bind
  # * modernize-raw-string-literal
  # * modernize-redundant-void-arg
  # * modernize-shrink-to-fit
  # * modernize-use-noexcept
  # * modernize-use-nodiscard
  # * modernize-use-equals-default
  # * modernize-concat-nested-namespaces
  # * modernize-replace-disallow-copy-and-assign-macro
  # * modernize-unary-static-assert
  # * modernize-use-uncaught-exceptions
  # * modernize-use-emplace                   <--- might not be exception-safe.
  # * modernize-use-transparent-functors      <--- useful but potentially dangerous...
  # * modernize-pass-by-value                 <--- CAREFUL!
  # * modernize-return-braced-init-list       <--- issue with readability
  # * modernize-use-trailing-return-type      <--- potential readability issue and dangerous
  #
  # * performance-for-range-copy
  # * performance-inefficient-algorithm
  # * performance-inefficient-vector-operation
  # * performance-move-const-arg
  # * performance-move-constructor-init
  # * performance-noexcept-move-constructor
  # * performance-trivially-destructible
  # * performance-type-promotion-in-math-fn
  # * performance-unnecessary-value-param
  #
  # * readability-braces-around-statements
  # * readability-convert-member-functions-to-static
  # * readability-delete-null-pointer
  # * readability-else-after-return
  # * readability-implicit-bool-conversion
  # * readability-inconsistent-declaration-parameter-name
  # * readability-isolate-declaration
  # * readability-make-member-function-const
  # * readability-string-compare
  # * readability-static-accessed-through-instance
  # * readability-simplify-subscript-expr
  # * readability-simplify-boolean-expr
  # * readability-redundant-member-init
  # * readability-redundant-function-ptr-dereference
  # * readability-redundant-declaration
  # * readability-redundant-control-flow
  # * readability-non-const-parameter
  # * readability-misplaced-array-index
  #
  # * misc-unused-alias-decls
  # * misc-redundant-expression
  # * misc-static-assert
  # ~~~

  # Concatenate all of the chosen options into a comma-delimited string
  list(JOIN OOMPH_CLANG_TIDY_CHECKS "," OOMPH_CLANG_TIDY_CHECKS)

  # Search the system paths for the clang-tidy tool
  find_program(CLANG_TIDY_PROGRAM "clang-tidy")

  # Run clang-tidy during the build if found
  if(CLANG_TIDY_PROGRAM)
    set(CMAKE_CXX_CLANG_TIDY
        ${CLANG_TIDY_PROGRAM}
        --use-color;
        --enable-check-profile;
        -checks=-*,${OOMPH_CLANG_TIDY_CHECKS};
        -header-filter=.*;
        -fix;
        --extra-arg=-std=c++14)
  else()
    message(
      STATUS "Requested static code analysis but you don't have clang-tidy!")
  endif()
endif()

if(OOMPH_ENABLE_INCLUDE_WHAT_YOU_USE)
  # Concatenate all of the chosen options into a comma-delimited string
  list(JOIN OOMPH_IWYU_CHECKS "," OOMPH_IWYU_CHECKS)

  # Search the system paths for the clang-tidy tool
  find_program(INCLUDE_WHAT_YOU_USE_PROGRAM "iwyu")

  # Run clang-tidy during the build if found
  if(INCLUDE_WHAT_YOU_USE_PROGRAM)
    message(
      STATUS
        "Requested static code analysis with the tool 'iwyu' but you don't have it!"
    )
  endif()
endif()
