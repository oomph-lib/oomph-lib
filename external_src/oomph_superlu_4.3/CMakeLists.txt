# ==============================================================================
# DOUBLE PRECISION SOURCES FROM SuperLU_4.3: We turn it into the library
# oomph_superlu_4.3 to avoid clashes with any existing installations of superlu
# (likely to be accessible via -lsuperlu)
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oomph_superlu_4.3 subdirectory")

# Define the headers
set(HEADERS
    slu_Cnames.h
    colamd.h
    supermatrix.h
    slu_util.h
    slu_ddefs.h
    slu_dcomplex.h
    slu_zdefs.h
    superlu_enum_consts.h)

# Define the sources
set(SOURCES
    colamd.c
    dgssvx.c
    dmyblas2.c
    dsnode_dfs.c
    mmd.c
    dcolumn_bmod.c
    dgstrf.c
    dpanel_bmod.c
    dsp_blas2.c
    relax_snode.c
    dcolumn_dfs.c
    dgstrs.c
    dpanel_dfs.c
    dsp_blas3.c
    sp_coletree.c
    dcopy_to_ucol.c
    dlacon.c
    dpivotgrowth.c
    dutil.c
    sp_ienv.c
    dgscon.c
    dpivotL.c
    get_perm_c.c
    sp_preorder.c
    dgsequ.c
    dlangs.c
    dpruneL.c
    heap_relax_snode.c
    superlu_timer.c
    dgsrfs.c
    dlaqgs.c
    dreadhb.c
    util.c
    dgssv.c
    dmemory.c
    dsnode_bmod.c
    memory.c
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
    zcolumn_dfs.c
    zcolumn_bmod.c
    zpivotL.c
    zpruneL.c
    zmemory.c
    zutil.c
    zmyblas2.c
    ilu_relax_snode.c
    ilu_heap_relax_snode.c
    mark_relax.c
    mc64ad.c
    qselect.c
    dgsisx.c
    dgsitrf.c
    dldperm.c
    ilu_ddrop_row.c
    ilu_dsnode_dfs.c
    ilu_dcolumn_dfs.c
    ilu_dpanel_dfs.c
    ilu_dcopy_to_ucol.c
    ilu_dpivotL.c
    ddiagonal.c)

# Define the name of the library
set(LIBNAME oomph_superlu_4.3)

# Define the library build type
#
# set(LIBTYPE STATIC)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  # LIBTYPE ${LIBTYPE}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

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

message(VERBOSE "Leaving oomph_superlu_4.3 subdirectory")
# ------------------------------------------------------------------------------
