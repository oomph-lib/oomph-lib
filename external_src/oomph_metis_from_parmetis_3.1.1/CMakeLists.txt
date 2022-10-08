# ==============================================================================
# METIS sources required by ParMETIS. PM: What else is this used for?
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oomph_metis_from_parmetis_3.1.1 subdirectory")

# ------------------------------------------------------------------------
# TODO: Address following from Makefile.am (but commented out):
#
# liboomph_metis_from_parmetis_3_1_1_la_LDFLAGS = -static
#
# Additional files that come with parmetis but aren't actually used by us.
# Included for completeness.
#
# EXTRA_DIST = stats.c
# ------------------------------------------------------------------------

# Define the headers
set(HEADERS metis.h)

# Additional headers that need to be distributed
set(HEADERS_NO_COMBINE
    defs.h
    macros.h
    proto.h
    rename.h
    stdheaders.h
    struct.h)

# Define the sources
set(SOURCES
    balance.c
    graph.c
    mcoarsen.c
    mkwayrefine.c
    pqueue.c
    bucketsort.c
    initpart.c
    memory.c
    mmatch.c
    refine.c
    ccgraph.c
    kmetis.c
    mesh.c
    mmd.c
    separator.c
    checkgraph.c
    kvmetis.c
    meshpart.c
    mpmetis.c
    sfm.c
    coarsen.c
    kwayfm.c
    mfm2.c
    mrefine2.c
    srefine.c
    compress.c
    kwayrefine.c
    mfm.c
    mrefine.c
    stat.c
    debug.c
    kwayvolfm.c
    mincover.c
    mutil.c
    estmem.c
    kwayvolrefine.c
    minitpart2.c
    myqsort.c
    subdomains.c
    fm.c
    match.c
    minitpart.c
    ometis.c
    timing.c
    fortran.c
    mbalance2.c
    mkmetis.c
    parmetis.c
    util.c
    frename.c
    mbalance.c
    mkwayfmh.c
    pmetis.c)

# Define the name of the library
set(LIBNAME oomph_metis_from_parmetis_3.1.1)

# Build it as a static/shared library depending on how oomph-lib is being built
set(LIBTYPE STATIC)
if(BUILD_SHARED_LIBS)
  set(LIBTYPE SHARED)
endif()

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  LIBTYPE ${LIBTYPE}
  HEADERS ${HEADERS}
  HEADERS_NO_COMBINE ${HEADERS_NO_COMBINE}
  SOURCES ${SOURCES})

# Silence the annoying warnings produced by this library
include(OomphTargetSilenceWarnings)
oomph_target_silence_warnings(
  TARGET ${LIBNAME}
  C_COMPILE_FLAGS -Wno-implicit-function-declaration -Wno-shift-op-parentheses
                  -Wno-format-security
  WARN_IF_FLAG_NOT_FOUND)

# Add the library to the list of external source libraries and make the updated
# variable value visible to the parent scope
list(APPEND EXTERNAL_SRC_LIBS ${LIBNAME})
set(EXTERNAL_SRC_LIBS ${EXTERNAL_SRC_LIBS} PARENT_SCOPE)

message(VERBOSE "Leaving oomph_metis_from_parmetis_3.1.1 subdirectory")
# ------------------------------------------------------------------------------
