# ==============================================================================
# DOUBLE PRECISION SOURCES FROM SuperLU_4.3: We turn it into the library
# oomph_superlu_4.3 to avoid clashes with any existing installations of superlu
# (likely to be accessible via -lsuperlu)
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oomph_superlu_5.2.2 subdirectory")

set(ENABLE_SINGLE ON)
set(ENABLE_DOUBLE ON)
set(ENABLE_COMPLEX ON)
set(ENABLE_COMPLEX16 ON)

set(HEADERS
    supermatrix.h
    slu_Cnames.h
    slu_dcomplex.h
    slu_scomplex.h
    slu_util.h
    superlu_enum_consts.h)

set(SOURCES
    superlu_timer.c
    util.c
    memory.c
    get_perm_c.c
    mmd.c
    sp_coletree.c
    sp_preorder.c
    sp_ienv.c
    relax_snode.c
    heap_relax_snode.c
    colamd.c
    ilu_relax_snode.c
    ilu_heap_relax_snode.c
    mark_relax.c
    mc64ad.c
    qselect.c
    input_error.c
    dmach.c
    smach.c)
set_source_files_properties(superlu_timer.c PROPERTIES COMPILE_FLAGS -O0)

if(ENABLE_SINGLE)
  list(APPEND HEADERS slu_sdefs.h)

  list(
    APPEND
    SOURCES
    slacon2.c
    smach.c
    sgssv.c
    sgssvx.c
    ssp_blas2.c
    ssp_blas3.c
    sgscon.c
    slangs.c
    sgsequ.c
    slaqgs.c
    spivotgrowth.c
    sgsrfs.c
    sgstrf.c
    sgstrs.c
    scopy_to_ucol.c
    ssnode_dfs.c
    ssnode_bmod.c
    spanel_dfs.c
    spanel_bmod.c
    sreadhb.c
    sreadrb.c
    sreadtriple.c
    sreadMM.c
    scolumn_dfs.c
    scolumn_bmod.c
    spivotL.c
    spruneL.c
    smemory.c
    sutil.c
    smyblas2.c
    sgsisx.c
    sgsitrf.c
    sldperm.c
    ilu_sdrop_row.c
    ilu_ssnode_dfs.c
    ilu_scolumn_dfs.c
    ilu_spanel_dfs.c
    ilu_scopy_to_ucol.c
    ilu_spivotL.c
    sdiagonal.c)
  set_source_files_properties(smach.c PROPERTIES COMPILE_FLAGS -O0)
endif()

if(ENABLE_DOUBLE)
  list(APPEND HEADERS slu_ddefs.h)

  list(
    APPEND
    SOURCES
    dlacon2.c
    dmach.c
    dgssv.c
    dgssvx.c
    dsp_blas2.c
    dsp_blas3.c
    dgscon.c
    dlangs.c
    dgsequ.c
    dlaqgs.c
    dpivotgrowth.c
    dgsrfs.c
    dgstrf.c
    dgstrs.c
    dcopy_to_ucol.c
    dsnode_dfs.c
    dsnode_bmod.c
    dpanel_dfs.c
    dpanel_bmod.c
    dreadhb.c
    dreadrb.c
    dreadtriple.c
    dreadMM.c
    dcolumn_dfs.c
    dcolumn_bmod.c
    dpivotL.c
    dpruneL.c
    dmemory.c
    dutil.c
    dmyblas2.c
    dgsisx.c
    dgsitrf.c
    dldperm.c
    ilu_ddrop_row.c
    ilu_dsnode_dfs.c
    ilu_dcolumn_dfs.c
    ilu_dpanel_dfs.c
    ilu_dcopy_to_ucol.c
    ilu_dpivotL.c
    ddiagonal.c
    dGetDiagU.c)
  set_source_files_properties(dmach.c PROPERTIES COMPILE_FLAGS -O0)
endif()

if(ENABLE_COMPLEX)
  list(APPEND HEADERS slu_cdefs.h)

  list(
    APPEND
    SOURCES
    clacon2.c
    scsum1.c
    icmax1.c
    scomplex.c
    cgssv.c
    cgssvx.c
    csp_blas2.c
    csp_blas3.c
    cgscon.c
    clangs.c
    cgsequ.c
    claqgs.c
    cpivotgrowth.c
    cgsrfs.c
    cgstrf.c
    cgstrs.c
    ccopy_to_ucol.c
    csnode_dfs.c
    csnode_bmod.c
    cpanel_dfs.c
    cpanel_bmod.c
    creadhb.c
    creadrb.c
    creadtriple.c
    creadMM.c
    ccolumn_dfs.c
    ccolumn_bmod.c
    cpivotL.c
    cpruneL.c
    cmemory.c
    cutil.c
    cmyblas2.c
    cgsisx.c
    cgsitrf.c
    cldperm.c
    ilu_cdrop_row.c
    ilu_csnode_dfs.c
    ilu_ccolumn_dfs.c
    ilu_cpanel_dfs.c
    ilu_ccopy_to_ucol.c
    ilu_cpivotL.c
    cdiagonal.c)
endif()

if(ENABLE_COMPLEX16)
  list(APPEND HEADERS slu_zdefs.h)

  list(
    APPEND
    SOURCES
    zlacon2.c
    dzsum1.c
    izmax1.c
    dcomplex.c
    zgssv.c
    zgssvx.c
    zsp_blas2.c
    zsp_blas3.c
    zgscon.c
    zlangs.c
    zgsequ.c
    zlaqgs.c
    zpivotgrowth.c
    zgsrfs.c
    zgstrf.c
    zgstrs.c
    zcopy_to_ucol.c
    zsnode_dfs.c
    zsnode_bmod.c
    zpanel_dfs.c
    zpanel_bmod.c
    zreadhb.c
    zreadrb.c
    zreadtriple.c
    zreadMM.c
    zcolumn_dfs.c
    zcolumn_bmod.c
    zpivotL.c
    zpruneL.c
    zmemory.c
    zutil.c
    zmyblas2.c
    zgsisx.c
    zgsitrf.c
    zldperm.c
    ilu_zdrop_row.c
    ilu_zsnode_dfs.c
    ilu_zcolumn_dfs.c
    ilu_zpanel_dfs.c
    ilu_zcopy_to_ucol.c
    ilu_zpivotL.c
    zdiagonal.c)
endif()

# Specify the libraries that we depend on
if(OOMPH_ENABLE_USE_OPENBLAS AND BLAS_FOUND)
  set(LINKLIBS BLAS::BLAS)
else()
  set(LINKLIBS oomph::oomph_blas)
endif()
if(NOT WIN32)
  list(APPEND LINKLIBS m)
endif()

# Define the name of the library
set(LIBNAME oomph_superlu_5.2.2)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  LINKLIBS ${LINKLIBS}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

# Silence the (annoying) warnings produced by this library
include(OomphTargetSilenceWarnings)
oomph_target_silence_warnings(TARGET ${LIBNAME}
                              C_COMPILE_FLAGS -Wno-unused-result)

# ------------------------------------------------------------------------------
# BOILERPLATE: Update variables in parent scope

# Update list of external sources preprocessor definitions and list of libraries
list(APPEND EXTERNAL_SRC_CXX_DEFINITIONS USING_OOMPH_SUPERLU)
list(APPEND EXTERNAL_SRC_LIBS ${LIBNAME})

# Make the updated variable values visible to the parent scope
set(EXTERNAL_SRC_CXX_DEFINITIONS ${EXTERNAL_SRC_CXX_DEFINITIONS} PARENT_SCOPE)
set(EXTERNAL_SRC_LIBS ${EXTERNAL_SRC_LIBS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(USING_OOMPH_SUPERLU TRUE CACHE INTERNAL "")

message(VERBOSE "Leaving oomph_superlu_5.2.2 subdirectory")
# ------------------------------------------------------------------------------
