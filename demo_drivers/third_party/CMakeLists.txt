# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered third_party subdirectory")

set(SUBDIRS)

if(OOMPH_HAS_NLOHMANN_JSON)
  list(APPEND SUBDIRS test_nlohmann_json)
endif()

if(OOMPH_HAS_SPDLOG)
  list(APPEND SUBDIRS test_spdlog)
endif()

# Add desired subdirectories to the build
foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving third_party subdirectory")
# ------------------------------------------------------------------------------
