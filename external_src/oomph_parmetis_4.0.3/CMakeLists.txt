# ==============================================================================
# PARMETIS 3.1.1: Turn into oomph_parmetis_3.1.1 to avoid clashes with any
# existing installations of metis which are likely to be accessible via -lmetis
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oomph_parmetis_3.1.1 subdirectory")

# Define the headers
set(HEADERS parmetislib.h)

# Additional headers that need to be distributed
set(HEADERS_NO_COMBINE
    defs.h
    macros.h
    parmetis.h
    proto.h
    rename.h
    struct.h)

# Define the sources
set(SOURCES
    akwayfm.c
    ametis.c
    balancemylink.c
    comm.c
    csrmatch.c
    ctrl.c
    debug.c
    diffutil.c
    frename.c
    gkmetis.c
    gkmpi.c
    graph.c
    initbalance.c
    initmsection.c
    initpart.c
    kmetis.c
    kwayrefine.c
    match.c
    mdiffusion.c
    mesh.c
    mmetis.c
    move.c
    msetup.c
    node_refine.c
    ometis.c
    pspases.c
    redomylink.c
    remap.c
    renumber.c
    rmetis.c
    selectq.c
    serial.c
    stat.c
    timer.c
    util.c
    wave.c
    weird.c
    wspace.c
    xyzpart.c)

# Define the name of the library
set(LIBNAME oomph_parmetis_4.0.3)

# Libraries to link against
set(LINKLIBS oomph_metis_from_parmetis_4.0.3)

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

message(VERBOSE "Leaving oomph_parmetis_3.1.1 subdirectory")
# ------------------------------------------------------------------------------
