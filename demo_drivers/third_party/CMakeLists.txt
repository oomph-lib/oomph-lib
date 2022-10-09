# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered third_party subdirectory")

set(SUBDIRS)

if(OOMPH_WANT_NLOHMANN_JSON)
  if(NOT OOMPH_HAS_NLOHMANN_JSON)
    message(WARNING "User requested 'nlohmann_json' but it could not be found!")
  else()
    list(APPEND SUBDIRS "${CMAKE_CURRENT_LIST_DIR}/nlohmann_json")
  endif()
endif()

# Add desired subdirectories to the build
foreach(SUBDIR ${SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving third_party subdirectory")
# ------------------------------------------------------------------------------
