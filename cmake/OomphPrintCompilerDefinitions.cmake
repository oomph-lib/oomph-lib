# ==============================================================================
# RANDOM SHITE TO BE DELETED AT SOME POINT
# ==============================================================================
# PRINT COMPILER DEFINITIONS: USED FOR DEBUGGING WHILE PORTING PROJECT
#
# # TODO: DELETE
get_directory_property(OOMPH_COMPILE_DEFINITIONS DIRECTORY ${CMAKE_SOURCE_DIR}
                                                           COMPILE_DEFINITIONS)
message(STATUS "COMPILER DEFINITIONS:")
foreach(DEFN ${OOMPH_COMPILE_DEFINITIONS})
  message(STATUS "  -D" ${DEFN})
endforeach()
# ------------------------------------------------------------------------------
