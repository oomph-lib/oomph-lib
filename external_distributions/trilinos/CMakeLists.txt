# ==============================================================================
# Build the Trilinos library
# ==============================================================================
include(ExternalProject)

# Specify the locations of the BLAS and LAPACK library artifacts
set(TPL_BLAS_LIBRARIES "$<TARGET_LINKER_FILE:oomph_blas>")
set(TPL_LAPACK_LIBRARIES "$<TARGET_LINKER_FILE:oomph_lapack>")

# An implicit function declaration appears in a C file in Trilinos v11.8.1,
# which the compiler issues an error over(!). I'm hoping this will be fixed in a
# newer version. To get around this for now, we'll just disable the warning
set(TRILINOS_C_FLAGS ${CMAKE_C_FLAGS} -Wno-implicit-function-declaration)

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

ExternalProject_Add(
  oomph_trilinos_project
  INSTALL_DIR "${TRILINOS_INSTALL_DIR}"
  URL "${CMAKE_CURRENT_LIST_DIR}/trilinos-11.8.1-Source.tar.gz"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> ${OOMPH_TRILINOS_FLAGS}
  DEPENDS oomph_blas oomph_lapack)

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

# Build it as a static/shared library depending on how oomph-lib is being built
set(LIBTYPE STATIC)
if(BUILD_SHARED_LIBS)
  set(LIBTYPE SHARED)
endif()

# Make a target that combines the Trilinos libraries into one
add_library(Trilinos::Trilinos ${LIBTYPE} IMPORTED GLOBAL)
target_link_libraries(Trilinos::Trilinos INTERFACE "${TRILINOS_LIBRARIES}")
target_include_directories(
  Trilinos::Trilinos INTERFACE $<BUILD_INTERFACE:${TRILINOS_INCLUDE_DIR}>
                               $<INSTALL_INTERFACE:include/trilinos>)

message("Trilinos::Trilinos: ${TRILINOS_INCLUDE_DIR}")

# Make it clear that our aliased library depends on the ExternalProject library
# being built to handle the race condition for linking
add_dependencies(Trilinos::Trilinos oomph_trilinos_project)

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
