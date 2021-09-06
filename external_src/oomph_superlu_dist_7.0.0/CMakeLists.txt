# ==============================================================================
# DOUBLE PRECISION SOURCES FROM SuperLU_Dist_7.0.0: We turn it into the library
# oomph_superlu_dist_7.0.0 to avoid clashes with any existing installations of
# SuperLU (likely to be accessible via -lsuperlu)
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oomph_superlu_dist_7.0.0 subdirectory")

set(ENABLE_DOUBLE "Enable double precision library" ON)
set(ENABLE_COMPLEX16 "Enable complex16 precision library" OFF)

# Indicate that we have ParMETIS and LAPACK (will be used in the config file)
set(HAVE_PARMETIS TRUE)
set(SLU_HAVE_LAPACK TRUE)

# Use 32-bit integers for index size. Change to 64 to use 64-bit integers
set(XSDK_INDEX_SIZE 32)

# ~~~
# CombBLAS. Currently disabled for now but we could possibly ship oomph-lib with
# it if we really want it.
# option(TPL_ENABLE_COMBBLASLIB   "Build the CombBLAS library" OFF)
# option(TPL_COMBBLAS_LIBRARIES "List of absolute paths to CombBLAS link libraries [].")
# option(TPL_COMBBLAS_INCLUDE_DIRS "List of absolute paths to CombBLAS include directories [].")
# ~~~

set(HEADERS
    superlu_FCnames.h
    dcomplex.h
    machines.h
    psymbfact.h
    superlu_defs.h
    superlu_enum_consts.h
    supermatrix.h
    util_dist.h
    colamd.h
    TreeBcast_slu.hpp
    TreeReduce_slu.hpp
    superlu_dist_config.h)

set(TEMPLATE_IMPL_HEADERS TreeBcast_slu_impl.hpp TreeReduce_slu_impl.hpp)

if(HAVE_CUDA)
  list(APPEND HEADERS cublas_utils.h cub/cub.cuh)
endif()

if(MSVC)
  list(APPEND HEADERS wingetopt.h)
endif()

# Define the sources
set(SOURCES
    sp_ienv.c
    etree.c
    sp_colorder.c
    get_perm_c.c
    mmd.c
    comm.c
    memory.c
    util.c
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
    colamd.c
    superlu_dist_version.c
    TreeInterface.cpp
    superlu_grid3d.c # 3D code
    supernodal_etree.c
    supernodalForest.c
    trfAux.c
    communication_aux.c
    treeFactorization.c
    sec_structs.c)

# Identify the C++ files amongst the C files
set(CXX_FILES TreeInterface.cpp TreeBcast_slu.hpp TreeReduce_slu.hpp
    TreeBcast_slu_impl.hpp TreeReduce_slu_impl.hpp)

if(HAVE_CUDA)
  list(APPEND SOURCES cublas_utils.c)
endif()

if(MSVC)
  list(APPEND SOURCES wingetopt.c)
endif()

set_source_files_properties(superlu_timer.c PROPERTIES COMPILE_FLAGS -O0)

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

  if(HAVE_CUDA)
    list(APPEND SOURCES dsuperlu_gpu.cu)
  endif()

  if(HAVE_COMBBLAS)
    list(APPEND SOURCES d_c2cpp_GetHWPM.cpp dHWPM_CombBLAS.hpp)

    # Identify the C++ files amongst the C files
    list(APPEND CXX_FILES d_c2cpp_GetHWPM.cpp dHWPM_CombBLAS.hpp)
  endif()
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

  if(HAVE_CUDA)
    list(APPEND SOURCES zsuperlu_gpu.cu)
  endif()

  if(HAVE_COMBBLAS)
    list(APPEND SOURCES z_c2cpp_GetHWPM.cpp zHWPM_CombBLAS.hpp)

    # Identify the C++ files amongst the C files
    list(APPEND CXX_FILES z_c2cpp_GetHWPM.cpp zHWPM_CombBLAS.hpp)
  endif()
endif()

foreach(CXX_FILE IN LISTS CXX_FILES)
  set_property(SOURCE ${CXX_FILE} PROPERTY CXX_STANDARD ${CMAKE_CXX_STANDARD})
  set_property(SOURCE ${CXX_FILE} PROPERTY CXX_STANDARD_REQUIRED TRUE)
endforeach()

# Libraries to link SuperLU Dist to
set(LINKLIBS ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} oomph_blas oomph_lapack
    oomph_parmetis_4.0.3) # ${COMBBLAS_LIB}

if(OOMPH_ENABLE_CUDA)
  list(APPEND LINKLIBS ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
endif()

# See:
# https://github.com/xiaoyeli/superlu_dist#dependent-external-libraries-blas-and-parmetis
# Might be able to link against OpenBLAS...
#
# -DTPL_BLAS_LIBRARIES="<BLAS library name>"

if(UNIX)
  list(APPEND LINKLIBS m)
endif()

# Define the name of the library
set(LIBNAME oomph_superlu_dist_7.0.0)

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
  HEADERS_NO_COMBINE ${TEMPLATE_IMPL_HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

# Update list of external sources preprocessor definitions and list of libraries
list(APPEND EXTERNAL_SRC_CXX_DEFINITIONS USING_OOMPH_SUPERLU_DIST)
list(APPEND EXTERNAL_SRC_LIBS ${LIBNAME})

# Make the updated variable values visible to the parent scope
set(EXTERNAL_SRC_CXX_DEFINITIONS ${EXTERNAL_SRC_CXX_DEFINITIONS} PARENT_SCOPE)
set(EXTERNAL_SRC_LIBS ${EXTERNAL_SRC_LIBS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(USING_OOMPH_SUPERLU_DIST TRUE CACHE INTERNAL "")

# Create *_config.h files from *_config.h.in input files. Place one copy in the
# current directory and another in the build directory
configure_file(${CMAKE_CURRENT_LIST_DIR}/superlu_dist_config.h.in
               ${CMAKE_CURRENT_LIST_DIR}/superlu_dist_config.h)
configure_file(${CMAKE_CURRENT_LIST_DIR}/superlu_dist_config.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/superlu_dist_config.h)

message(VERBOSE "Leaving oomph_superlu_dist_7.0.0 subdirectory")
# ------------------------------------------------------------------------------
