# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered joe_cool subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(one_d_poisson C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()

# Add a block of this type for your own executable
oomph_add_executable(
  NAME one_d_poisson # <--- name of the executable
  SOURCES one_d_poisson.cc # <--- source code(s)
  LIBRARIES oomph::poisson oomph::meshes oomph::generic) # <-- oomph-lib
                                                         # libraries required
                                                         # (order matters:
                                                         # specific ones go
                                                         # first; more general
                                                         # ones at the end)

message(VERBOSE "Leaving joe_cool subdirectory")
# ------------------------------------------------------------------------------
