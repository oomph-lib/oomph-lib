# ==============================================================================
# Build the Trilinos library
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered Trilinos subdirectory")

# Cache the Trilinos package version number
set(TRILINOS_VERSION 11.8.1 CACHE STRING "Trilinos version number" FORCE)

# cmake-format: off
# TODO: Update Trilinos generic code to use new Trilinos functionality.
# The patch file and how to apply it
# set(PATCH_FILE 0001-Patch-of-CMakeLists.txt-for-FetchContents.patch)
# Patch above is for the latest version of Trilinos.
#
# Patch below is for Trilinos 11.8.1.
set(PATCH_FILE 0001-Patch-of-CMakeLists.txt-for-oomph-lib-export.patch)
set(PATCH_COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/patches/${PATCH_FILE} <SOURCE_DIR>
    && git checkout <SOURCE_DIR>/CMakeLists.txt
    && git apply <SOURCE_DIR>/${PATCH_FILE})
# cmake-format: on

# The compiler issues an error about an implicit function declaration, which
# appears in a C file in Trilinos v11.8.1. I'm hoping this will be fixed in a
# newer version, but to get around this for now, we'll just disable the warning
set(TRILINOS_C_FLAGS ${CMAKE_C_FLAGS} -Wno-implicit-function-declaration)

# Specify the absolute paths to the BLAS/LAPACK libraries; if we're using
# OpenBLAS, specify the path to the libraries that we found. If not, and we're
# just going to build our own, then use a generator expression to define the
# path to the library (the value will be populated properly at build-time, when
# the actual library artifact is generated)
set(OOMPH_BLAS_LIBRARIES)
set(OOMPH_LAPACK_LIBRARIES)
if(OOMPH_USE_OPENBLAS AND BLAS_FOUND)
  set(OOMPH_BLAS_LIBRARIES "${BLAS_LIBRARIES}")
  set(OOMPH_LAPACK_LIBRARIES "${LAPACK_LIBRARIES}")
else()
  set(OOMPH_BLAS_LIBRARIES "$<TARGET_LINKER_FILE:oomph_blas>")
  set(OOMPH_LAPACK_LIBRARIES "$<TARGET_LINKER_FILE:oomph_lapack>")
endif()

# Build flags
set(OOMPH_TRILINOS_FLAGS
    -DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_C_FLAGS:STRING=${TRILINOS_C_FLAGS}
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
    -DTPL_ENABLE_MPI:BOOL=${OOMPH_ENABLE_MPI}
    -DTPL_ENABLE_BLAS:BOOL=ON
    -DTrilinos_ENABLE_TESTS:BOOL=OFF
    -DTrilinos_ENABLE_EXAMPLES:BOOL=OFF
    -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF
    -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF
    -DTrilinos_ENABLE_Amesos:BOOL=ON
    -DTrilinos_ENABLE_Anasazi:BOOL=ON
    -DTrilinos_ENABLE_AztecOO:BOOL=ON
    -DTrilinos_ENABLE_Epetra:BOOL=ON
    -DTrilinos_ENABLE_EpetraExt:BOOL=ON
    -DTrilinos_ENABLE_Ifpack:BOOL=ON
    -DTrilinos_ENABLE_ML:BOOL=ON
    -DTrilinos_ENABLE_Teuchos:BOOL=ON
    -DTrilinos_ENABLE_Triutils:BOOL=ON
    -DTPL_ENABLE_BLAS:BOOL=ON
    -DTPL_ENABLE_LAPACK:BOOL=ON
    -DTPL_BLAS_LIBRARIES:PATH=${OOMPH_BLAS_LIBRARIES}
    -DTPL_LAPACK_LIBRARIES:PATH=${OOMPH_LAPACK_LIBRARIES})

if(OOMPH_ENABLE_MPI AND MPI_FOUND)
  list(APPEND OOMPH_TRILINOS_FLAGS -DMPI_BASE_DIR:PATH=${MPI_CXX_INCLUDE_DIRS})
endif()

# Make OOMPH_TRILINOS_MODULE_PATH a cache variable so that we can get its value
# when the oomphlibConfig.cmake.in file is processed. It'll come in handy when
# we install Trilinos to a location which isn't in the user's PATH variable
set(OOMPH_TRILINOS_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/trilinos_install"
    CACHE PATH "Path to Trilinos installation directory")
set(OOMPH_TRILINOS_INCLUDE_DIR "${OOMPH_TRILINOS_INSTALL_DIR}/include")
set(OOMPH_TRILINOS_LIBRARY_DIR "${OOMPH_TRILINOS_INSTALL_DIR}/lib")

# Trilinos library names for linking
set(TRILINOS_LIBRARIES
    amesos
    anasazi
    aztecoo
    epetra
    epetraext
    ifpack
    ml
    teuchosremainder
    teuchosnumerics
    teuchoscomm
    teuchosparameterlist
    teuchoscore
    triutils
    CACHE
    INTERNAL
    "List of Trilinos libraries to link against"
    FORCE)

# The libraries built by Trilinos with full paths and lib decorations (e.g.
# "lib" prefix and ".a" suffix)
set(TRILINOS_LIBRARIES_WITH_PATHS)
foreach(TRILINOS_LIBRARY IN LISTS TRILINOS_LIBRARIES)
  list(
    APPEND
    TRILINOS_LIBRARIES_WITH_PATHS
    "${OOMPH_TRILINOS_LIBRARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}${TRILINOS_LIBRARY}${CMAKE_STATIC_LIBRARY_SUFFIX}"
  )
endforeach()

# Cache the library artifacts to link against (needed for oomphlibConfig.cmake)
set(OOMPH_TRILINOS_LIBRARIES_WITH_PATHS ${TRILINOS_LIBRARIES_WITH_PATHS} CACHE
    INTERNAL "" FORCE)

# If we're not using the OpenBLAS BLAS/LAPACK libraries then we're building our
# own. To ensure that they're built before Trilinos is (so that Trilinos can
# link against them), we need to list them as Trilinos dependencies
set(OOMPH_TRILINOS_DEPENDS oomph::oomph_blas oomph::oomph_lapack)
if(OOMPH_USE_OPENBLAS AND BLAS_FOUND)
  set(OOMPH_TRILINOS_DEPENDS)
endif()

# The recipe for building the actual library
include(ExternalProject)
ExternalProject_Add(
  oomph_trilinos_project
  GIT_REPOSITORY https://github.com/trilinos/Trilinos
  GIT_TAG trilinos-release-11-8-1
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE
  PATCH_COMMAND ${PATCH_COMMAND}
  INSTALL_DIR "${OOMPH_TRILINOS_INSTALL_DIR}"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> ${OOMPH_TRILINOS_FLAGS}
  BUILD_BYPRODUCTS ${OOMPH_TRILINOS_LIBRARIES_WITH_PATHS}
  DEPENDS ${OOMPH_TRILINOS_DEPENDS})

# Build it as a static/shared library depending on how oomph-lib is being built
set(LIBTYPE STATIC)
if(BUILD_SHARED_LIBS)
  set(LIBTYPE SHARED)
endif()

# Create the include directory to silence the CMake error about a non-existent
# path for the "imported target Trilinos::Trilinos"
file(MAKE_DIRECTORY ${OOMPH_TRILINOS_INCLUDE_DIR})

# Make a target that combines the Trilinos libraries into one and indicates to
# our library that this library is provided externally. We need to do this
# because the Trilinos library only gets configured, built and installed during
# our build phase, so we can't get hold of variables defined by Trilinos, e.g. a
# list of the libraries that it built, until then
add_library(Trilinos::Trilinos ${LIBTYPE} IMPORTED GLOBAL)
set_target_properties(
  Trilinos::Trilinos
  PROPERTIES IMPORTED_LOCATION "${OOMPH_TRILINOS_LIBRARIES_WITH_PATHS}"
             INTERFACE_INCLUDE_DIRECTORIES "${OOMPH_TRILINOS_INCLUDE_DIR}"
             INTERFACE_LINK_LIBRARIES "${TRILINOS_LIBRARIES}")

# Make it clear that our aliased library depends on the ExternalProject library
# being built to handle the race condition for linking
add_dependencies(Trilinos::Trilinos oomph_trilinos_project)

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and make the updated
# variable values visible to the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_TRILINOS)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_TRILINOS TRUE CACHE INTERNAL "")

message(VERBOSE "Leaving trilinos subdirectory")
# ------------------------------------------------------------------------------
