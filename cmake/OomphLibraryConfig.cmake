# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# A convenience function to automate the creation and installation of an
# oomph-lib library, typically consisting of one or more headers and sources.
# The function:
#   - Creates a (potentially) compiled library or an INTERFACE library (if no
#     sources).
#   - Generates a "combined" header named <LIBNAME>.h in the build include
#     directory that includes relevant headers and any SOURCES_NO_BUILD files.
#   - Installs the library artifacts (if built) and headers to the standard
#     locations under <INSTALL-DIR>/include/<LIBNAME>.
#   - Optionally installs additional headers (HEADERS_NO_COMBINE) without
#     including them in the combined header.
#   - Excludes certain headers (HEADERS_NO_INSTALL) from installation
#     or the combined header.
#
# Arguments
# ---------
#   - LIBNAME (Required): The name of the library to create.
#   - HEADERS (Optional): Headers to be installed and included in the combined
#       header.
#   - SOURCES (Optional): Source files to compile and install (if provided, the
#       library is built as STATIC/SHARED, unless overridden by LIBTYPE).
#   - SOURCES_NO_BUILD (Optional): Source files that should be installed and
#       added to the combined header but NOT compiled (e.g., template
#       implementations).
#   - HEADERS_NO_COMBINE (Optional): Headers that should be installed but not
#       included in the combined header (e.g., external or system headers).
#   - HEADERS_NO_INSTALL (Optional): Headers that should NOT be installed or
#       included in the combined header (e.g., Fortran headers).
#   - LIBTYPE (Optional): Overrides the default library type (STATIC/SHARED).
#       If omitted and SOURCES are provided, it respects the global
#       BUILD_SHARED_LIBS setting; if no SOURCES, an INTERFACE library is used.
#   - INCLUDE_SUBDIRECTORY (Optional): Subdirectory name under
#       <INSTALL-DIR>/include/ for the headers. Defaults to <LIBNAME> if not
#       provided.
#   - LINKLIBS (Optional): Additional libraries to link against.
#
# USAGE:
# ------
#   include(OomphLibraryConfig)
#
#   oomph_library_config(
#     LIBNAME             my_library
#     HEADERS             foo.h bar.h
#     SOURCES             foo.cc bar.cc
#     SOURCES_NO_BUILD    baz.cc
#     HEADERS_NO_COMBINE  system_header.h
#     HEADERS_NO_INSTALL  fortran-wrapper.h
#     LIBTYPE             STATIC
#     INCLUDE_SUBDIRECTORY my_lib
#     LINKLIBS            other_dependency
#   )
#
# Once configured, the library can be linked via:
#   target_link_libraries(your_target PRIVATE <PROJECT_NAMESPACE>::my_library)
#
# SUMMARY:
# --------
#   +--------------------+--------------+--------------------+-------------+
#   |     ARGUMENT       | LIBRARY DEP. | IN COMBINED HEADER |  INSTALLED  |
#   +--------------------+--------------+--------------------+-------------+
#   | HEADERS            |      x       |         x          |      x      |
#   | SOURCES            |      x       |                    |      x      |
#   | SOURCES_NO_BUILD   |      x       |         x          |      x      |
#   | HEADERS_NO_COMBINE |      x       |                    |      x      |
#   | HEADERS_NO_INSTALL |      x       |                    |             |
#   +--------------------+--------------+--------------------+-------------+
#
# -----------------------------------------------------------------------------
# IMPORTANT:
# ----------
# This function should be called from the directory containing all the
# headers/sources for your library (so that CMAKE_CURRENT_SOURCE_DIR is set
# correctly).
# =============================================================================
# cmake-format: on
include_guard()

# ------------------------------------------------------------------------------
function(oomph_library_config)
  # Define the supported set of keywords
  set(PREFIX ARG)
  set(FLAGS "")
  set(SINGLE_VALUE_ARGS LIBNAME LIBTYPE INCLUDE_SUBDIRECTORY)
  set(MULTI_VALUE_ARGS
      HEADERS
      SOURCES
      HEADERS_NO_COMBINE
      SOURCES_NO_BUILD
      HEADERS_NO_INSTALL
      LINKLIBS)

  # Process the arguments passed in
  include(CMakeParseArguments)
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}"
                        "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

  # We must be given the directory containing the sources/headers, the name of
  # the library we need to create, and header files. If any of these fields are
  # empty then issue a fatal error and stop.
  if(NOT ARG_LIBNAME)
    message(FATAL_ERROR "Missing ${ARG_LIBNAME}.")
  endif()

  # Redefine the headers, sources, and the name of the library in this scope but
  # with clearer variable names
  set(HEADERS ${ARG_HEADERS})
  set(SOURCES ${ARG_SOURCES})
  set(LIBNAME ${ARG_LIBNAME})
  set(LIBTYPE ${ARG_LIBTYPE})
  set(INCLUDE_SUBDIRECTORY ${ARG_INCLUDE_SUBDIRECTORY})
  set(LINKLIBS ${ARG_LINKLIBS})
  set(HEADERS_NO_COMBINE ${ARG_HEADERS_NO_COMBINE})
  set(HEADERS_NO_INSTALL ${ARG_HEADERS_NO_INSTALL})
  set(SOURCES_NO_BUILD ${ARG_SOURCES_NO_BUILD})
  set(LIBRARY_DEPS ${SOURCES} ${HEADERS} ${HEADERS_NO_COMBINE}
      ${HEADERS_NO_INSTALL} ${SOURCES_NO_BUILD})

  # Mark any files in SOURCES_NO_BUILD as headers, so they aren't compiled
  set_source_files_properties(${SOURCES_NO_BUILD} PROPERTIES HEADER_FILE_ONLY
                                                             TRUE)

  # Create the library and specify the library type. If no sources are provided,
  # we assume it is a header-only library, in which case it should be created as
  # an INTERFACE library. Note, in this case, it makes no sense to specify the
  # library type (i.e. static/shared), as the library itself will not be
  # built/compiled.
  #
  # NOTE: If the user does not specify LIBTYPE, it will be inferred from the
  # CMake cache variable BUILD_SHARED_LIBS.
  if(SOURCES)
    add_library(${LIBNAME} ${LIBTYPE} ${LIBRARY_DEPS})
    set(INCLUDE_TYPE PUBLIC)
    set(LINK_TYPE PUBLIC)
  else()
    add_library(${LIBNAME} INTERFACE)
    set(INCLUDE_TYPE INTERFACE)
    set(LINK_TYPE INTERFACE)
  endif()

  # Provide an alias for the library so that users importing the library can
  # link to the library they require by selecting the library they need from the
  # project "namespace"
  set(LIBRARY_ALIAS ${PROJECT_NAMESPACE}::${LIBNAME})

  # Create an alias library name, e.g. generic can be linked with oomph::generic
  add_library(${LIBRARY_ALIAS} ALIAS ${LIBNAME})

  # Storage for the list of libraries exported by oomph-lib; append to the cache
  # variable in two steps: first edit the variable locally, then update the
  # variable globally. This approach avoids adding a leading ";" in the variable
  list(APPEND OOMPHLIB_LIBRARIES ${LIBRARY_ALIAS})
  set(OOMPHLIB_LIBRARIES "${OOMPHLIB_LIBRARIES}" CACHE INTERNAL "")

  # Specify the location of the library headers at build-time/install-time. Add
  # these build-time directories for the following reasons:
  #
  # * CMAKE_CURRENT_SOURCE_DIR: To get direct access to the sources
  # * CMAKE_CURRENT_SOURCE_DIR/../: To access to the sources through the
  #   enclosing folder
  # * CMAKE_CURRENT_BINARY_DIR/../: To get access to the combined headers (which
  #   are placed in the build directory so they don't pollute the source tree).
  #
  # At install-time, we just place the files under:
  #
  # ${CMAKE_INSTALL_PREFIX}/${OOMPH_INSTALL_INCLUDE_DIR}
  #
  # which will most likely just be
  #
  # <oomph-lib-root>/install/include/oomphlib
  target_include_directories(
    ${LIBNAME}
    ${INCLUDE_TYPE}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../>
    $<INSTALL_INTERFACE:${OOMPH_INSTALL_INCLUDE_DIR}>)

  # Add the location of oomph-lib-config.h header if required
  if(OOMPH_ADD_CONFIG_H)
    target_include_directories(${LIBNAME} ${INCLUDE_TYPE}
                               $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src>)
  endif()

  # Specify library dependencies
  if(LINKLIBS)
    target_link_libraries(${LIBNAME} ${LINK_TYPE} ${LINKLIBS})
  endif()

  # Silence the annoying warnings produced by the library if it is going to be
  # compiled, i.e. it is not an INTERFACE library
  if(NOT (LINK_TYPE STREQUAL INTERFACE))
    include(OomphTargetSilenceWarnings)
    oomph_target_silence_warnings(TARGET ${LIBNAME}
                                  CXX_COMPILE_FLAGS -Wno-undefined-var-template)
  endif()

  # If the user requested that we enable sanitisers
  if(OOMPH_ENABLE_SANITISER_SUPPORT)
    oomph_enable_sanitizers(${LIBNAME})
  endif()

  # ----------------------------------------------------------------------------
  # The install rules: we want to do the following
  #
  # 1. Install ${LIBNAME} in <INSTALL-DIR>/lib/${PROJECT_NAME},
  # 2. Install ${LIBNAME}.h in <INSTALL-DIR>/include/${PROJECT_NAME}, and
  # 3. Install ${HEADERS} in <INSTALL-DIR>/include/${PROJECT_NAME}/${LIBNAME}
  #
  # The header ${LIBNAME}.h in <INSTALL-DIR>/include will include the relevant
  # headers from <INSTALL-DIR>/include/${LIBNAME}/. This allows us to group
  # associated headers together.
  # ----------------------------------------------------------------------------
  # ~~~
  # When headers get installed to "<INSTALL-PATH>/oomphlib/include/<SUBDIR>/",
  # we provide a combined header in the ".../include/" directory that can be
  # included by the user to include all files shipped with a library. The
  # include paths in the combined header have the form
  #                     #include <<SUBDIR>/<FILE 1>>
  #                                 ...
  #                     #include <<SUBDIR>/<FILE N>>
  # The "<SUBDIR>" is simply the folder that encloses it.
  # Therefore we need to the correct "<SUBDIR>" value. In most cases, this is
  # obvious -- it is just the name of the library itself, e.g. the <SUBDIR> of
  # the "generic" library, src/generic/, contains the header files.
  #
  # FIXME: We should, for good practice, include the full path from the root
  # project directory to preserve the folder structure. E.g. for the
  # 'space_time_block_preconditioner' library defined in
  #             src/space_time/space_time_block_preconditioner/
  # we install the combined header to
  #                 <INSTALL-DIR>/include/${PROJECT_NAME}/
  # and the files to
  #    <INSTALL-DIR>/include/${PROJECT_NAME}/space_time_block_preconditioner/
  # instead of
  #  <INSTALL-DIR>/include/${PROJECT_NAME}/space_time/space_time_block_preconditioner/
  # This could be problematic if we create libraries with the same name. Ideally
  # we would install the combined header to
  #             <INSTALL-DIR>/include/${PROJECT_NAME}/space_time/
  # so the user gets access to all the files with the following include:
  #         #include "space_time/space_time_block_preconditioner.h"
  # instead of
  #             #include "space_time_block_preconditioner.h"
  # as is currently done.
  # ~~~
  #
  # We want to provide a combined header above the level of a library directory
  # which includes the path to files in the library The include paths we provide
  # for #headers In general, the name of the folder that encloses the definition
  # of a library is the same as the name of the library, e.g. src/generic/
  # contains the definition of the "generic" library. However, in some cases the
  # name of the directory may be different to the name of the library, e.g. the
  # library "space_time_unsteady_heat_equal_order_galerkin" is contained in the
  # folder src/ galerkin_equal_order
  #
  # TODO: Update this to automatically compute the enclosing subdirectory!
  set(LIBRARY_INCLUDE_SUBDIR ${LIBNAME})
  if(INCLUDE_SUBDIRECTORY)
    set(LIBRARY_INCLUDE_SUBDIR ${INCLUDE_SUBDIRECTORY})
  endif()

  # Create the combined header and install it to the build directory, but one
  # folder above
  include(OomphCreateCombinedHeader)
  oomph_create_combined_header(
    TARGET "${CMAKE_CURRENT_BINARY_DIR}/../${LIBNAME}.h"
    HEADERS ${HEADERS} ${SOURCES_NO_BUILD}
    SUBDIRECTORY ${LIBRARY_INCLUDE_SUBDIR})

  # Define the installation locations and export target. Although this only
  # installs the library, it's important that the includes directory is added so
  # that anything linking to the exported library knows where the associated
  # headers live.
  install(
    TARGETS ${LIBNAME}
    EXPORT ${TARGETS_EXPORT_NAME}
    LIBRARY DESTINATION "${OOMPH_INSTALL_LIB_DIR}"
    ARCHIVE DESTINATION "${OOMPH_INSTALL_LIB_DIR}"
    RUNTIME DESTINATION "${OOMPH_INSTALL_BIN_DIR}")

  # Install the combined header
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/../${LIBNAME}.h"
          DESTINATION "${OOMPH_INSTALL_INCLUDE_DIR}")

  # The directory to install the library headers to relative(!) to the
  # installation directory
  set(INCLUDE_DIR_FOR_THIS_LIBRARY
      "${OOMPH_INSTALL_INCLUDE_DIR}/${LIBRARY_INCLUDE_SUBDIR}")

  # Create the headers install directory. Do it now so that there's a valid
  # directory to create symlinks from (if used)
  install(DIRECTORY DESTINATION "${INCLUDE_DIR_FOR_THIS_LIBRARY}")

  # Combine everything that shouldn't be built into a single variable
  set(ALL_HEADERS ${HEADERS} ${HEADERS_NO_COMBINE} ${SOURCES_NO_BUILD})

  # Instead of "install(FILES ...)", create a symlink for each header
  if(OOMPH_INSTALL_HEADERS_AS_SYMLINKS)
    foreach(HEADER IN LISTS ALL_HEADERS)
      get_filename_component(HEADER_NAME ${HEADER} NAME)
      install(
        CODE "execute_process(
        COMMAND \"${CMAKE_COMMAND}\" -E create_symlink
                \"${CMAKE_CURRENT_SOURCE_DIR}/${HEADER}\"
                \"\${CMAKE_INSTALL_PREFIX}/${INCLUDE_DIR_FOR_THIS_LIBRARY}/${HEADER_NAME}\"
      )")
    endforeach()
  else()
    # ...or just install copies of the headers
    install(FILES ${ALL_HEADERS} DESTINATION "${INCLUDE_DIR_FOR_THIS_LIBRARY}")
  endif()
  # ----------------------------------------------------------------------------
endfunction()
# ------------------------------------------------------------------------------
