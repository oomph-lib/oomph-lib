# ==============================================================================
# zlib sources for gzip decompression
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oomph_zlib subdirectory")

if(WIN32
   OR CYGWIN
   OR MINGW)
  message(WARNING "No support for build of oomph_gzip! Skipping...")
  return()
endif()

include(CheckTypeSize)
include(CheckFunctionExists)
include(CheckIncludeFile)

# Check to see if we have large file support
set(CMAKE_REQUIRED_DEFINITIONS -D_LARGEFILE64_SOURCE=1)

# We add these other definitions here because CheckTypeSize.cmake in CMake 2.4.x
# does not automatically do so and we want compatibility with CMake 2.4.x.
if(HAVE_SYS_TYPES_H)
  list(APPEND CMAKE_REQUIRED_DEFINITIONS -DHAVE_SYS_TYPES_H)
endif()
if(HAVE_STDINT_H)
  list(APPEND CMAKE_REQUIRED_DEFINITIONS -DHAVE_STDINT_H)
endif()
if(HAVE_STDDEF_H)
  list(APPEND CMAKE_REQUIRED_DEFINITIONS -DHAVE_STDDEF_H)
endif()
check_type_size(off64_t OFF64_T)
if(HAVE_OFF64_T)
  add_definitions(-D_LARGEFILE64_SOURCE=1)
endif()
set(CMAKE_REQUIRED_DEFINITIONS) # clear variable

# Check for fseeko
check_function_exists(fseeko HAVE_FSEEKO)
if(NOT HAVE_FSEEKO)
  add_definitions(-DNO_FSEEKO)
endif()

# Check for unistd.h
check_include_file(unistd.h Z_HAVE_UNISTD_H)

# Define the headers
set(HEADERS
    zconf.h
    zlib.h
    crc32.h
    deflate.h
    gzguts.h
    inffast.h
    inffixed.h
    inflate.h
    inftrees.h
    trees.h
    zutil.h)

# Define the sources
set(SOURCES
    adler32.c
    compress.c
    crc32.c
    deflate.c
    gzclose.c
    gzlib.c
    gzread.c
    gzwrite.c
    inflate.c
    infback.c
    inftrees.c
    inffast.c
    trees.c
    uncompr.c
    zutil.c)

# Define the name of the library
set(LIBNAME oomph_zlib)

# parse the full version number from zlib.h and include in ZLIB_FULL_VERSION
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/zlib.h _zlib_h_contents)
string(REGEX REPLACE ".*#define[ \t]+ZLIB_VERSION[ \t]+\"([-0-9A-Za-z.]+)\".*"
                     "\\1" ZLIB_FULL_VERSION ${_zlib_h_contents})

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

# Silence the (annoying) warnings produced by this library
include(OomphTargetSilenceWarnings)
oomph_target_silence_warnings(TARGET ${LIBNAME}
                              C_COMPILE_FLAGS -Wno-deprecated-non-prototype)

if(Z_HAVE_UNISTD_H)
  target_compile_definitions(${LIBNAME} PRIVATE Z_HAVE_UNISTD_H)
endif()

# ------------------------------------------------------------------------------
# BOILERPLATE: Update variables in parent scope

# Update list of external sources preprocessor definitions and list of libraries
list(APPEND EXTERNAL_SRC_LIBS ${LIBNAME})
list(APPEND EXTERNAL_SRC_CXX_DEFINITIONS OOMPH_HAS_ZLIB)

# Make the updated variable values visible to the parent scope
set(EXTERNAL_SRC_LIBS ${EXTERNAL_SRC_LIBS} PARENT_SCOPE)
set(EXTERNAL_SRC_CXX_DEFINITIONS ${EXTERNAL_SRC_CXX_DEFINITIONS} PARENT_SCOPE)

message(VERBOSE "Leaving oomph_zlib subdirectory")
# ------------------------------------------------------------------------------
