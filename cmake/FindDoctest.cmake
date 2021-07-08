# =============================================================================
# Download the unit-testing suite, doctest.
# =============================================================================
include(FetchContent)
FetchContent_Declare(
  doctest_project
  GIT_REPOSITORY "https://github.com/onqtam/doctest.git"
  GIT_TAG 2.4.6)
FetchContent_MakeAvailable(doctest_project)
