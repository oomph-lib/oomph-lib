# ==============================================================================
# Lapack sources required by SuperLU. We turn it into the library oomph_lapack
# to avoid clashes with any existing installations of SuperLU (likely to be
# accessible via -llapack)
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oomph_metis_from_parmetis_3.1.1 subdirectory")

# ------------------------------------------------------------------------
# TODO: Address following from Makefile.am (but commented out):
#
# liboomph_metis_from_parmetis_3_1_1_la_LDFLAGS = -static
# ------------------------------------------------------------------------

# Define the headers
set(HEADERS metislib.h)

# Additional headers that need to be distributed
set(HEADERS_NO_COMBINE
    defs.h
    gklib_defs.h
    gklib_rename.h
    macros.h
    metis.h
    proto.h
    rename.h
    stdheaders.h
    struct.h)

# Define the sources
set(SOURCES
    auxapi.c
    balance.c
    bucketsort.c
    checkgraph.c
    coarsen.c
    compress.c
    contig.c
    debug.c
    fm.c
    fortran.c
    frename.c
    gklib.c
    graph.c
    initpart.c
    kmetis.c
    kwayfm.c
    kwayrefine.c
    mcutil.c
    mesh.c
    meshpart.c
    minconn.c
    mincover.c
    mmd.c
    ometis.c
    options.c
    parmetis.c
    pmetis.c
    refine.c
    separator.c
    sfm.c
    srefine.c
    stat.c
    timing.c
    util.c
    wspace.c)

# Define the name of the library
set(LIBNAME oomph_metis_from_parmetis_4.0.3)

# Libraries to link against
set(LINKLIBS oomph_gklib_from_metis_from_parmetis_4.0.3)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  HEADERS_NO_COMBINE ${HEADERS_NO_COMBINE}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

# ------------------------------------------------------------------------------
# BOILERPLATE: Update variables in parent scope

# Add the library to the list of external source libraries and make the updated
# variable value visible to the parent scope
list(APPEND EXTERNAL_SRC_LIBS ${LIBNAME})
set(EXTERNAL_SRC_LIBS ${EXTERNAL_SRC_LIBS} PARENT_SCOPE)

message(VERBOSE "Leaving oomph_metis_from_parmetis_3.1.1 subdirectory")
# ------------------------------------------------------------------------------
