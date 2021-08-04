# ==============================================================================
# Build the Trilinos library
# ==============================================================================
include(ExternalProject)

# Specify the locations of the BLAS and LAPACK library artifacts
set(TPL_BLAS_LIBRARIES "$<TARGET_LINKER_FILE:oomph_blas>")
set(TPL_LAPACK_LIBRARIES "$<TARGET_LINKER_FILE:oomph_lapack>")

# The compiler issues an error about an implicit function declaration, which
# appears in a C file in Trilinos v11.8.1. I'm hoping this will be fixed in a
# newer version, but to get around this for now, we'll just disable the warning
set(TRILINOS_C_FLAGS ${CMAKE_C_FLAGS} -Wno-implicit-function-declaration)

# Build flags
set(OOMPH_TRILINOS_FLAGS
    -DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_Fortran_COMPILER=/usr/local/bin/gfortran # ${CMAKE_Fortran_COMPILER}
    -DCMAKE_C_FLAGS:STRING=${TRILINOS_C_FLAGS}
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
    -DTPL_ENABLE_MPI:BOOL=${OOMPH_HAS_MPI}
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
    -DTPL_BLAS_LIBRARIES:PATH="${TPL_BLAS_LIBRARIES}"
    -DTPL_LAPACK_LIBRARIES:PATH="${TPL_LAPACK_LIBRARIES}")

# TODO: Add FindTrilinos.cmake?
set(TRILINOS_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/trilinos_install")
set(TRILINOS_INCLUDE_DIR "${TRILINOS_INSTALL_DIR}/include")
set(TRILINOS_LIBRARY_DIR "${TRILINOS_INSTALL_DIR}/lib")
# file(MAKE_DIRECTORY ${TRILINOS_INCLUDE_DIR})

# Trilinos library names for linking
set(TRILINOS_LIBRARIES
    amesos
    anasazi
    aztecoo
    epetra
    epetraext
    ifpack
    ml
    teuchos
    triutils)

# The libraries built by Trilinos with full paths and lib decorations (e.g.
# "lib" prefix and ".a" suffix)
set(TRILINOS_LIBRARIES_WITH_PATHS)
foreach(TRILINOS_LIBRARY IN LISTS TRILINOS_LIBRARIES)
  list(
    APPEND
    TRILINOS_LIBRARIES_WITH_PATHS
    "${TRILINOS_LIBRARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}${TRILINOS_LIBRARY}${CMAKE_STATIC_LIBRARY_SUFFIX}"
  )
endforeach()

# The recipe for building the actual library
externalproject_add(
  oomph_trilinos_project
  INSTALL_DIR "${TRILINOS_INSTALL_DIR}"
  URL "${CMAKE_CURRENT_LIST_DIR}/trilinos-11.8.1-Source.tar.gz"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> ${OOMPH_TRILINOS_FLAGS}
  BUILD_BYPRODUCTS ${TRILINOS_LIBRARIES_WITH_PATHS}
  DEPENDS oomph::oomph_blas oomph::oomph_lapack)
externalproject_get_property(oomph_trilinos_project INSTALL_DIR)

# Build it as a static/shared library depending on how oomph-lib is being built
set(LIBTYPE STATIC)
if(BUILD_SHARED_LIBS)
  set(LIBTYPE SHARED)
endif()

# Make a target that combines the Trilinos libraries into one
add_library(Trilinos::Trilinos ${LIBTYPE} IMPORTED GLOBAL)
set_target_properties(
  Trilinos::Trilinos
  PROPERTIES IMPORTED_LOCATION "${TRILINOS_LIBRARIES_WITH_PATHS}"
             INTERFACE_INCLUDE_DIRECTORIES ${TRILINOS_INCLUDE_DIR})
target_link_libraries(Trilinos::Trilinos INTERFACE "${TRILINOS_LIBRARIES}")

# Make it clear that our aliased library depends on the ExternalProject library
# being built to handle the race condition for linking
add_dependencies(Trilinos::Trilinos oomph_trilinos_project)

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and list of libraries
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_TRILINOS)
list(APPEND EXTERNAL_DIST_LIBS Trilinos::Trilinos)

# Make the updated variable values visible to the parent scope
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)
set(EXTERNAL_DIST_LIBS ${EXTERNAL_DIST_LIBS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_TRILINOS TRUE CACHE INTERNAL "")

message(VERBOSE "Leaving trilinos subdirectory")
# ------------------------------------------------------------------------------
