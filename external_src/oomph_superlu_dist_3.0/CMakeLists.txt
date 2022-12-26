# ==============================================================================
# DOUBLE PRECISION SOURCES FROM SuperLU_dist_3.0: We turn it into the library
# oomph_superlu_dist_3.0 to avoid clashes with any existing installations of
# superlu (likely to be accessible via -lsuperlu)
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oomph_superlu_dist_3.0 subdirectory")

# Define the headers
set(HEADERS
    superlu_defs.h
    superlu_ddefs.h
    Cnames.h
    machines.h
    psymbfact.h
    supermatrix.h
    old_colamd.h
    util_dist.h
    superlu_enum_consts.h)

# Define the sources. N.B. This is the list of files specified by SuperLU Dist
# to build their double-precision parallel library. There are a few extra files
# included here that were excluded in the old Makefile.am because of clashes
# with the LAPACK source (which we don't encounter because modern CMake uses a
# target-based system rather than a "throw all your bloody sources in a bucket
# and hope for the best"...)
set(SOURCES
    dlangs.c
    dgsequ.c
    dlaqgs.c
    dutil.c
    dmemory.c
    dmyblas2.c
    dsp_blas2.c
    dsp_blas3.c
    dldperm.c
    ddistribute.c
    pdgstrf.c
    pdgstrs_Bglobal.c
    pdgstrs1.c
    pdgssvx_ABglobal.c
    pdgsrfs_ABXglobal.c
    pdgsmv_AXglobal.c
    pdgssvx.c
    pdgstrs.c
    pddistribute.c
    pdlangs.c
    pdutil.c
    pdgsequ.c
    pdlaqgs.c
    pdgsrfs.c
    pdgsmv.c
    pdgstrs_lsum.c
    pdsymbfact_distdata.c
    dreadhb.c
    dreadrb.c
    dreadtriple.c
    sp_ienv.c
    etree.c
    sp_colorder.c
    get_perm_c.c
    mmd.c
    comm.c
    memory.c
    util.c
    superlu_grid.c
    pxerbla.c
    superlu_timer.c
    GetDiagU.c
    mc64ad.c
    symbfact.c
    psymbfact.c
    psymbfact_util.c
    get_perm_c_parmetis.c
    lsame.c
    xerbla.c
    slamch.c
    dlamch.c)

# Define the name of the library
set(LIBNAME oomph_superlu_dist_3.0)

# Specify the libraries that we depend on
set(LINKLIBS oomph_metis_from_parmetis_3.1.1 MPI::MPI_C)
if(OOMPH_USE_BLAS_FROM)
  list(APPEND LINKLIBS BLAS::BLAS)
else()
  list(APPEND LINKLIBS oomph_blas)
endif()

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  LINKLIBS ${LINKLIBS}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

# Silence the (annoying) warnings produced by this library
include(OomphTargetSilenceWarnings)
oomph_target_silence_warnings(
  TARGET ${LIBNAME}
  C_COMPILE_FLAGS
    -Wno-compare-distinct-pointer-types
    -Wno-dangling-else
    -Wno-enum-conversion
    -Wno-format
    -Wno-format-security
    -Wno-implicit-function-declaration
    -Wno-incompatible-pointer-types
    -Wno-int-conversion
    -Wno-unsequenced
    -Wno-parentheses
    -Wno-return-type
    -Wno-switch
    -Wno-unused-value)

# ------------------------------------------------------------------------------
# BOILERPLATE: Update variables in parent scope

# Update list of external sources preprocessor definitions and list of libraries
list(APPEND EXTERNAL_SRC_CXX_DEFINITIONS USING_OOMPH_SUPERLU_DIST)
list(APPEND EXTERNAL_SRC_LIBS ${LIBNAME})

# Make the updated variable values visible to the parent scope
set(EXTERNAL_SRC_CXX_DEFINITIONS ${EXTERNAL_SRC_CXX_DEFINITIONS} PARENT_SCOPE)
set(EXTERNAL_SRC_LIBS ${EXTERNAL_SRC_LIBS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(USING_OOMPH_SUPERLU_DIST TRUE CACHE INTERNAL "")

message(VERBOSE "Leaving oomph_superlu_dist_3.0 subdirectory")
# ------------------------------------------------------------------------------
