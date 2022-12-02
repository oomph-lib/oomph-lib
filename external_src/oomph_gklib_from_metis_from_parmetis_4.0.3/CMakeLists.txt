# ==============================================================================
# GKlib sources required by METIS which is required by ParMETIS 4.0.3
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE
        "Entered oomph_gklib_from_metis_from_parmetis_4.0.3 subdirectory")

# ------------------------------------------------------------------------------
# Boilerplate from GKlibSystem.cmake which ships with ParMETIS 4.0.3:

# Helper modules.
include(CheckFunctionExists)
include(CheckIncludeFile)

# Setup options.
option(GDB "enable use of GDB" OFF)
option(ASSERT "turn asserts on" OFF)
option(ASSERT2 "additional assertions" OFF)
option(DEBUG "add debugging support" OFF)
option(GPROF "add gprof support" OFF)
option(OPENMP "enable OpenMP support" OFF)
option(PCRE "enable PCRE support" OFF)
option(GKREGEX "enable GKREGEX support" OFF)
option(GKRAND "enable GKRAND support" OFF)

# Add compiler flags.
if(MSVC)
  set(GKLIB_C_OPTS "/Ox")
  set(GKLIB_C_OPTIONS "-DWIN32 -DMSC -D_CRT_SECURE_NO_DEPRECATE -DUSE_GKREGEX")
elseif(MINGW)
  set(GKLIB_C_OPTS "-DUSE_GKREGEX")
else()
  set(GKLIB_C_OPTS "-O3")
  set(GKLIB_C_OPTIONS "-DLINUX -D_FILE_OFFSET_BITS=64")
endif()
if(CYGWIN)
  set(GKLIB_C_OPTIONS "${GKLIB_C_OPTIONS} -DCYGWIN")
endif()
if(CMAKE_COMPILER_IS_GNUCC)
  # GCC opts.
  set(GKLIB_C_OPTIONS "${GKLIB_C_OPTIONS} -std=c99 -fno-strict-aliasing")
  if(NOT MINGW)
    set(GKLIB_C_OPTIONS "${GKLIB_C_OPTIONS} -fPIC")
  endif()
  # GCC warnings.
  set(GKLIB_C_OPTIONS
      "${GKLIB_C_OPTIONS} -Wall -pedantic -Wno-unused-but-set-variable -Wno-unused-variable -Wno-unknown-pragmas"
  )
elseif(${CMAKE_C_COMPILER_ID} MATCHES "Sun")
  # Sun insists on -xc99.
  set(GKLIB_C_OPTIONS "${GKLIB_C_OPTIONS} -xc99")
endif()

# Find OpenMP if it is requested.
if(OPENMP)
  include(FindOpenMP)
  if(OPENMP_FOUND)
    set(GKLIB_C_OPTIONS "${GKLIB_C_OPTIONS} -D__OPENMP__ ${OpenMP_C_FLAGS}")
  else()
    message(WARNING "OpenMP was requested but support was not found")
  endif()
endif()

# Add various definitions.
if(GDB)
  set(GKLIB_C_OPTS "${GKLIB_C_OPTS} -g")
  set(GKLIB_C_OPTIONS "${GKLIB_C_OPTIONS} -Werror")
endif()

if(DEBUG)
  set(GKLIB_C_OPTS "-g")
  set(GKLIB_C_OPTIONS "${GKLIB_C_OPTIONS} -DDEBUG")
endif()

if(GPROF)
  set(GKLIB_C_OPTS "-pg")
endif()

if(NOT ASSERT)
  set(GKLIB_C_OPTIONS "${GKLIB_C_OPTIONS} -DNDEBUG")
endif()

if(NOT ASSERT2)
  set(GKLIB_C_OPTIONS "${GKLIB_C_OPTIONS} -DNDEBUG2")
endif()

# Add various options
if(PCRE)
  set(GKLIB_C_OPTIONS "${GKLIB_C_OPTIONS} -D__WITHPCRE__")
endif()

if(GKREGEX)
  set(GKLIB_C_OPTIONS "${GKLIB_C_OPTIONS} -DUSE_GKREGEX")
endif()

if(GKRAND)
  set(GKLIB_C_OPTIONS "${GKLIB_C_OPTIONS} -DUSE_GKRAND")
endif()

# Check for features.
check_include_file(execinfo.h HAVE_EXECINFO_H)
if(HAVE_EXECINFO_H)
  set(GKLIB_C_OPTIONS "${GKLIB_C_OPTIONS} -DHAVE_EXECINFO_H")
endif()

check_function_exists(getline HAVE_GETLINE)
if(HAVE_GETLINE)
  set(GKLIB_C_OPTIONS "${GKLIB_C_OPTIONS} -DHAVE_GETLINE")
endif()

# Custom check for TLS.
if(MSVC)
  set(GKLIB_C_OPTIONS "${GKLIB_C_OPTIONS} -D__thread=__declspec(thread)")
else()
  # This if checks if that value is cached or not.
  if("${HAVE_THREADLOCALSTORAGE}" MATCHES "^${HAVE_THREADLOCALSTORAGE}$")
    try_compile(HAVE_THREADLOCALSTORAGE ${CMAKE_BINARY_DIR}
                ${CMAKE_CURRENT_LIST_DIR}/conf/check_thread_storage.c)
    if(HAVE_THREADLOCALSTORAGE)
      message(STATUS "checking for thread-local storage - found")
    else()
      message(STATUS "checking for thread-local storage - not found")
    endif()
  endif()
  if(NOT HAVE_THREADLOCALSTORAGE)
    set(GKLIB_C_OPTIONS "${GKLIB_C_OPTIONS} -D__thread=")
  endif()
endif()
# ------------------------------------------------------------------------------

# Find GKlib sources.
set(SOURCES
    b64.c
    blas.c
    csr.c
    error.c
    evaluate.c
    fkvkselect.c
    fs.c
    getopt.c
    gkregex.c
    graph.c
    htable.c
    io.c
    itemsets.c
    mcore.c
    memory.c
    omp.c
    pdb.c
    pqueue.c
    random.c
    rw.c
    seq.c
    sort.c
    string.c
    timers.c
    tokenizer.c
    util.c)

set(HEADERS
    GKlib.h
    gk_arch.h
    gk_defs.h
    gk_externs.h
    gk_getopt.h
    gk_macros.h
    gk_mkblas.h
    gk_mkmemory.h
    gk_mkpqueue.h
    gk_mkpqueue2.h
    gk_mkrandom.h
    gk_mksort.h
    gk_mkutils.h
    gk_proto.h
    gk_struct.h
    gk_types.h
    gkregex.h
    ms_inttypes.h
    ms_stat.h
    ms_stdint.h)

# Define the name of the library
set(LIBNAME oomph_gklib_from_metis_from_parmetis_4.0.3)

# Build it as a static library
set(LIBTYPE STATIC)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  LIBTYPE ${LIBTYPE}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

# Set the official C flags
target_compile_options(${LIBNAME} PRIVATE "${GKLIB_C_OPTIONS} ${GKLIB_C_OPTS}")

if(UNIX)
  target_link_libraries(${LIBNAME} m)
endif()

# ------------------------------------------------------------------------------
# BOILERPLATE: Update variables in parent scope

# Add the library to the list of external source libraries and make the updated
# variable value visible to the parent scope
list(APPEND EXTERNAL_SRC_LIBS ${LIBNAME})
set(EXTERNAL_SRC_LIBS ${EXTERNAL_SRC_LIBS} PARENT_SCOPE)

message(VERBOSE
        "Leaving oomph_gklib_from_metis_from_parmetis_4.0.3 subdirectory")
# ------------------------------------------------------------------------------
