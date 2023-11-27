# ==============================================================================
# DOUBLE PRECISION SOURCES FROM SuperLU_dist_8.2.1: We turn it into the library
# oomph_superlu_dist_8.2.1 to avoid clashes with any existing installations of
# superlu (likely to be accessible via -lsuperlu)
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oomph_superlu_dist_8.2.1 subdirectory")

if(BUILD_SHARED_LIBS)
  message(
    FATAL_ERROR
      "Build of SuperLUDist as a shared library is likely to break. Care has not been taken when stripping down the original sources to ensure the build will work."
  )
endif()

# NOTE: All of the CUDA/HIP code has been stripped out so it cannot currently be
# enabled.

# Variables required by superlu_dist_config.h.in
set(HAVE_CUDA OFF)
set(HAVE_HIP OFF)
set(HAVE_PARMETIS ON)
set(HAVE_COLAMD ON)
set(SLU_HAVE_LAPACK ON)
set(HAVE_COMBBLAS OFF)

# Set to 64 to enable 64-bit index mode
set(XSDK_INDEX_SIZE 32)

# Generate config header; depends on values of variables above
configure_file(${CMAKE_CURRENT_LIST_DIR}/superlu_dist_config.h.in
               ${CMAKE_CURRENT_LIST_DIR}/superlu_dist_config.h)

set(ENABLE_SINGLE "Enable single precision library" ON)
set(ENABLE_DOUBLE "Enable double precision library" ON)
set(ENABLE_COMPLEX16 "Enable complex16 precision library" ON)

# Define the headers
set(HEADERS
    superlu_FCnames.h
    dcomplex.h
    machines.h
    psymbfact.h
    superlu_defs.h
    superlu_enum_consts.h
    supermatrix.h
    util_dist.h
    gpu_api_utils.h
    gpu_wrapper.h
    superlu_dist_config.h
    superlu_FortranCInterface.h)

if(MSVC)
  list(APPEND HEADERS wingetopt.h)
endif()

# Define the sources.
set(SOURCES
    sp_ienv.c
    etree.c
    sp_colorder.c
    get_perm_c.c
    mmd.c
    comm.c
    memory.c
    util.c
    gpu_api_utils.c
    superlu_grid.c
    pxerr_dist.c
    superlu_timer.c
    symbfact.c
    psymbfact.c
    psymbfact_util.c
    get_perm_c_parmetis.c
    mc64ad_dist.c
    xerr_dist.c
    smach_dist.c
    dmach_dist.c
    superlu_dist_version.c
    comm_tree.c
    superlu_grid3d.c
    supernodal_etree.c
    supernodalForest.c
    trfAux.c
    communication_aux.c
    treeFactorization.c
    sec_structs.c)

if(MSVC)
  list(APPEND SOURCES wingetopt.c)
endif()

set_source_files_properties(superlu_timer.c PROPERTIES COMPILE_FLAGS -O0)

if(ENABLE_SINGLE)
  list(APPEND HEADERS superlu_sdefs.h slustruct_gpu.h)
  list(
    APPEND
    SOURCES
    slangs_dist.c
    sgsequ_dist.c
    slaqgs_dist.c
    sutil_dist.c
    smemory_dist.c
    smyblas2_dist.c
    ssp_blas2_dist.c
    ssp_blas3_dist.c
    psgssvx.c
    psgssvx_d2.c
    psgssvx_ABglobal.c
    sreadhb.c
    sreadrb.c
    sreadtriple.c
    sreadtriple_noheader.c
    sbinary_io.c
    sreadMM.c
    psgsequ.c
    pslaqgs.c
    sldperm_dist.c
    pslangs.c
    psutil.c
    pssymbfact_distdata.c
    sdistribute.c
    psdistribute.c
    psgstrf.c
    sstatic_schedule.c
    psgstrf2.c
    psgstrs.c
    psgstrs1.c
    psgstrs_lsum.c
    psgstrs_Bglobal.c
    psgsrfs.c
    psgsrfs_d2.c
    psgsmv.c
    psgsrfs_ABXglobal.c
    psgsmv_AXglobal.c
    psGetDiagU.c
    psgssvx3d.c # 3D code
    snrformat_loc3d.c
    psgstrf3d.c
    streeFactorization.c
    streeFactorizationGPU.c
    sgather.c
    sscatter3d.c
    ps3dcomm.c
    strfAux.c
    scommunication_aux.c
    strfCommWrapper.c
    ssuperlu_blas.c
    psgssvx_d2.c # with double-precision IR
    psgsrfs_d2.c
    psgsmv_d2.c
    psgsequb.c)
endif()

if(ENABLE_DOUBLE)
  list(APPEND HEADERS superlu_ddefs.h dlustruct_gpu.h)
  list(
    APPEND
    SOURCES
    dlangs_dist.c
    dgsequ_dist.c
    dlaqgs_dist.c
    dutil_dist.c
    dmemory_dist.c
    dmyblas2_dist.c
    dsp_blas2_dist.c
    dsp_blas3_dist.c
    pdgssvx.c
    pdgssvx_ABglobal.c
    dreadhb.c
    dreadrb.c
    dreadtriple.c
    dreadtriple_noheader.c
    dbinary_io.c
    dreadMM.c
    pdgsequ.c
    pdlaqgs.c
    dldperm_dist.c
    pdlangs.c
    pdutil.c
    pdsymbfact_distdata.c
    ddistribute.c
    pddistribute.c
    pdgstrf.c
    dstatic_schedule.c
    pdgstrf2.c
    pdgstrs.c
    pdgstrs1.c
    pdgstrs_lsum.c
    pdgstrs_Bglobal.c
    pdgsrfs.c
    pdgsmv.c
    pdgsrfs_ABXglobal.c
    pdgsmv_AXglobal.c
    pdGetDiagU.c
    pdgssvx3d.c # 3D code
    dnrformat_loc3d.c
    pdgstrf3d.c
    dtreeFactorization.c
    dtreeFactorizationGPU.c
    dgather.c
    dscatter3d.c
    pd3dcomm.c
    dtrfAux.c
    dcommunication_aux.c
    dtrfCommWrapper.c
    dsuperlu_blas.c)
endif()

if(ENABLE_COMPLEX16)
  list(APPEND HEADERS superlu_zdefs.h zlustruct_gpu.h)
  list(
    APPEND
    SOURCES
    dcomplex_dist.c
    zlangs_dist.c
    zgsequ_dist.c
    zlaqgs_dist.c
    zutil_dist.c
    zmemory_dist.c
    zmyblas2_dist.c
    zsp_blas2_dist.c
    zsp_blas3_dist.c
    pzgssvx.c
    pzgssvx_ABglobal.c
    zreadhb.c
    zreadrb.c
    zreadtriple.c
    zreadtriple_noheader.c
    zbinary_io.c
    zreadMM.c
    pzgsequ.c
    pzlaqgs.c
    zldperm_dist.c
    pzlangs.c
    pzutil.c
    pzsymbfact_distdata.c
    zdistribute.c
    pzdistribute.c
    pzgstrf.c
    zstatic_schedule.c
    pzgstrf2.c
    pzgstrs.c
    pzgstrs1.c
    pzgstrs_lsum.c
    pzgstrs_Bglobal.c
    pzgsrfs.c
    pzgsmv.c
    pzgsrfs_ABXglobal.c
    pzgsmv_AXglobal.c
    pzGetDiagU.c
    pzgssvx3d.c # 3D code
    znrformat_loc3d.c
    pzgstrf3d.c
    ztreeFactorization.c
    ztreeFactorizationGPU.c
    zscatter3d.c
    zgather.c
    pz3dcomm.c
    ztrfAux.c
    zcommunication_aux.c
    ztrfCommWrapper.c
    zsuperlu_blas.c)
endif()

# Define the name of the library
set(LIBNAME oomph_superlu_dist_8.2.1)

# Specify the libraries that we depend on
set(LINKLIBS MPI::MPI_CXX MPI::MPI_C oomph_metis_from_parmetis_4.0.3
    ${OpenBLAS_LIBRARIES})

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  LINKLIBS ${LINKLIBS}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

if(MSVC AND BUILD_SHARED_LIBS)
  set_target_properties(${LIBNAME} PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

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

message(VERBOSE "Leaving oomph_superlu_dist_8.2.1 subdirectory")
# ------------------------------------------------------------------------------
