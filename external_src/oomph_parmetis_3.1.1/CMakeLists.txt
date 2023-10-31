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
    proto.h
    rename.h
    stdheaders.h
    struct.h)

# Define the sources
set(SOURCES
    frename.c
    kmetis.c
    mmetis.c
    selectq.c
    akwayfm.c
    gkmetis.c
    kwaybalance.c
    move.c
    serial.c
    ametis.c
    grsetup.c
    kwayfm.c
    msetup.c
    setup.c
    backcompat.c
    iidxsort.c
    kwayrefine.c
    node_refine.c
    stat.c
    balancemylink.c
    iintsort.c
    lmatch.c
    ometis.c
    timer.c
    comm.c
    ikeysort.c
    match.c
    order.c
    util.c
    csrmatch.c
    ikeyvalsort.c
    mdiffusion.c
    pspases.c
    wave.c
    debug.c
    initbalance.c
    redomylink.c
    weird.c
    diffutil.c
    initmsection.c
    memory.c
    remap.c
    xyzpart.c
    fpqueue.c
    initpart.c
    mesh.c
    rmetis.c)

# Define the name of the library
set(LIBNAME oomph_parmetis_3.1.1)

# The libraries that we depend on
set(LINKLIBS oomph_metis_from_parmetis_3.1.1)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  LINKLIBS ${LINKLIBS}
  HEADERS ${HEADERS}
  HEADERS_NO_COMBINE ${HEADERS_NO_COMBINE}
  SOURCES ${SOURCES})

# Silence the annoying warnings produced by this library
include(OomphTargetSilenceWarnings)
oomph_target_silence_warnings(
  TARGET ${LIBNAME} C_COMPILE_FLAGS -Wno-implicit-function-declaration
                                    -Wno-format-security)

# ------------------------------------------------------------------------------
# BOILERPLATE: Update variables in parent scope

# Add the library to the list of external source libraries and make the updated
# variable value visible to the parent scope
list(APPEND EXTERNAL_SRC_LIBS ${LIBNAME})
set(EXTERNAL_SRC_LIBS ${EXTERNAL_SRC_LIBS} PARENT_SCOPE)

message(VERBOSE "Leaving oomph_parmetis_3.1.1 subdirectory")
# ------------------------------------------------------------------------------
