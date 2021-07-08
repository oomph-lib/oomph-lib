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
  set(OOMPH_CLANG_TIDY_CHECKS
      modernize-avoid-bind
      modernize-deprecated-headers
      modernize-loop-convert
      modernize-make-shared
      modernize-make-unique
      modernize-pass-by-value
      modernize-raw-string-literal
      modernize-redundant-void-arg
      modernize-replace-auto-ptr
      modernize-shrink-to-fit
      modernize-use-auto
      modernize-use-bool-literals
      modernize-use-emplace
      modernize-use-noexcept
      modernize-use-nodiscard
      modernize-use-equals-default
      modernize-use-equals-delete
      modernize-use-nullptr
      modernize-use-override
      modernize-use-transparent-functors
      modernize-use-using
      modernize-concat-nested-namespaces
      modernize-replace-disallow-copy-and-assign-macro
      modernize-replace-random-shuffle
      modernize-return-braced-init-list
      modernize-unary-static-assert
      modernize-use-uncaught-exceptions)

  # performance-for-range-copy performance-inefficient-algorithm
  # performance-inefficient-vector-operation performance-move-const-arg
  # performance-move-constructor-init performance-noexcept-move-constructor
  # performance-trivially-destructible performance-type-promotion-in-math-fn
  # performance-unnecessary-value-param readability-braces-around-statements
  # readability-convert-member-functions-to-static
  # readability-delete-null-pointer readability-else-after-return
  # readability-implicit-bool-conversion
  # readability-inconsistent-declaration-parameter-name
  # readability-isolate-declaration readability-make-member-function-const
  # readability-string-compare misc-unused-alias-decls misc-redundant-expression
  # misc-static-assert readability-static-accessed-through-instance
  # readability-simplify-subscript-expr readability-simplify-boolean-expr
  # readability-redundant-member-init
  # readability-redundant-function-ptr-dereference
  # readability-redundant-declaration readability-redundant-control-flow
  # readability-non-const-parameter readability-misplaced-array-index)

  # Concatenate all of the chosen options into a comma-delimited string
  list(JOIN OOMPH_CLANG_TIDY_CHECKS "," OOMPH_CLANG_TIDY_CHECKS)

  # Search the system paths for the clang-tidy tool
  find_program(CLANG_TIDY_PROGRAM "clang-tidy")

  # Run clang-tidy during the build if found
  if(CLANG_TIDY_PROGRAM)
    set(CMAKE_CXX_CLANG_TIDY
        ${CLANG_TIDY_PROGRAM}
        -j4
        -checks=-*,${OOMPH_CLANG_TIDY_CHECKS}
        -header-filter=${CMAKE_SOURCE_DIR}/src/.*
        -fix
        --extra-arg-before=-std=c++11)
  else()
    message(
      STATUS "Requested static code analysis but you don't have clang-tidy!")
  endif()
endif()
