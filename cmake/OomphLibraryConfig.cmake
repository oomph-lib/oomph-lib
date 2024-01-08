# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Generalise the creation of an oomph-lib library for the cases the process is
# easy to automate. This typically occurs when we have a set of headers and
# (possibly) a set of sources, and all we want to do is create a library from
# these headers/sources then install them in the standard location (i.e.
# build/include). A combined header will be created in build/include/ with the
# name ${LIBNAME}.h which includes the relevant headers/sources that will be
# installed to build/include/${LIBNAME}/.
#
# If no sources are provided then we assume that the library should be
# header-only and it will thus be created as an INTERFACE library.
#
# IMPORTANT:
# ----------
# This function must be called from the directory containing the files we need
# to build the library as the ${CMAKE_CURRENT_SOURCE_DIR} variable will hold the
# path to the calling directory, so we won't need to ask the user to provide it.
# For more details, see pg. 59 of:
#
#           "Scott Craig (2021), Professional CMake: A Practical Guide"
#
# USAGE:
# ------
#     include(OomphLibraryConfig)
#     oomph_library_config(LIBNAME             <library-name>
#                          [HEADERS            <headers-of-the-library>]
#                          [SOURCES            <sources-of-the-library>]
#                          [SOURCES_NO_BUILD   <sources-to-install-only>]
#                          [HEADERS_NO_COMBINE <headers-to-install-only>]
#                          [HEADERS_NO_INSTALL <headers-to-not-install>])
#
# The keyword HEADERS_NO_COMBINE is used to define extra headers that need to be
# installed in build/include but should not be placed in the combined header.
# These are typically headers that the system defines and not those associated
# with the library itself. The keyword SOURCES_NO_BUILD is used to refer to
# sources that are part of the library but should not be built such as template
# implementation headers (for pure template headers, e.g. meshes). These files
# will be installed to the include directory and also added in the combined
# header. The keyword HEADERS_NO_INSTALL is used to refer to headers that the
# library depends on but should not be part of the combined header or be placed
# in build/include/. This is helpful, for example, for fortran headers, like in
# src/generic.
#
# The above is summarised in the table below:
#
#   +--------------------+--------------+--------------------+-------------+
#   |           ARGUMENT | LIBRARY DEP. | IN COMBINED HEADER |  INSTALLED  |
#   +--------------------+--------------+--------------------+-------------+
#   |            HEADERS |      x       |         x          |      x      |
#   |            SOURCES |      x       |                    |      x      |
#   |   SOURCES_NO_BUILD |      x       |         x          |      x      |
#   | HEADERS_NO_COMBINE |      x       |                    |      x      |
#   | HEADERS_NO_INSTALL |      x       |                    |             |
#   +--------------------+--------------+--------------------+-------------+
# =============================================================================
# cmake-format: on
include_guard()

# cmake-format: off
# TODO: Add functionality to include .cc files that shouldn't be compiled:
# KEYWORD: SOURCES_NO_BUILD = ...
# cmake-format: on

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
  if(NOT ${PREFIX}_LIBNAME)
    message(FATAL_ERROR "Missing ${${PREFIX}_LIBNAME}.")
  endif()

  # Redefine the headers, sources, and the name of the library in this scope but
  # with clearer variable names
  set(HEADERS ${${PREFIX}_HEADERS})
  set(SOURCES ${${PREFIX}_SOURCES})
  set(LIBNAME ${${PREFIX}_LIBNAME})
  set(LIBTYPE ${${PREFIX}_LIBTYPE})
  set(INCLUDE_SUBDIRECTORY ${${PREFIX}_INCLUDE_SUBDIRECTORY})
  set(LINKLIBS ${${PREFIX}_LINKLIBS})
  set(HEADERS_NO_COMBINE ${${PREFIX}_HEADERS_NO_COMBINE})
  set(HEADERS_NO_INSTALL ${${PREFIX}_HEADERS_NO_INSTALL})
  set(SOURCES_NO_BUILD ${${PREFIX}_SOURCES_NO_BUILD})
  set(LIBRARY_DEPS ${SOURCES} ${HEADERS} ${HEADERS_NO_COMBINE}
      ${HEADERS_NO_INSTALL} ${SOURCES_NO_BUILD})

  # If sources that should not be built (SOURCES_NO_BUILD) are provided then
  # mark them as headers so that they do not get compiled in the library build
  set_source_files_properties(${SOURCES_NO_BUILD} PROPERTIES HEADER_FILE_ONLY
                                                             TRUE)

  # Create the library and specify the library type. If no sources are provided,
  # we assume it is a header-only library, in which case it should be created as
  # an INTERFACE library. Note, in this case, it makes no sense to specify the
  # library type, as the library itself will not be built/compiled.
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

  if(APPLE)
    target_link_options(${LIBNAME} ${LINK_TYPE} -Wl,-no_compact_unwind)
  endif()

  # ----------------------------------------------------------------------------
  # The install rules: we want to do the following
  #
  # 1. Install ${LIBNAME} in <INSTALL-DIR>/lib,
  # 2. Install ${LIBNAME}.h in <INSTALL-DIR>/include, and
  # 3. Install ${HEADERS} in <INSTALL-DIR>/include/${LIBNAME}
  #
  # The header ${LIBNAME}.h in <INSTALL-DIR>/include will include the relevant
  # headers from <INSTALL-DIR>/include/${LIBNAME}/. This allows us to group
  # associated headers together.
  # ----------------------------------------------------------------------------
  # ~~~
  # TODO: Add a note somewhere that the library directory should either match
  # the library name or be set
  # TODO: The <SUBDIR> value should NOT be the path from src/ as this would
  # break for the build.
  # FIXME: We should, for good practice, include the full path from the root
  # source directory...
  #
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
  # the "generic" library src/generic/ contains the
  # definition of
  # ~~~
  #
  # We want to provide a combined header above the level of a library directory
  # which includes the path to files in the livrary The include paths we provide
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
            COMPONENT ${PROJECT_NAME}_Runtime
            NAMELINK_COMPONENT ${PROJECT_NAME}_Development
    ARCHIVE DESTINATION "${OOMPH_INSTALL_LIB_DIR}"
            COMPONENT ${PROJECT_NAME}_Development
    RUNTIME DESTINATION "${OOMPH_INSTALL_BIN_DIR}"
            COMPONENT ${PROJECT_NAME}_Runtime
    INCLUDES
    DESTINATION "${OOMPH_INSTALL_INCLUDE_DIR}")

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

  # Install (or symlink) the headers
  install(FILES ${ALL_HEADERS} DESTINATION "${INCLUDE_DIR_FOR_THIS_LIBRARY}")
  # ----------------------------------------------------------------------------
endfunction()
# ------------------------------------------------------------------------------
