# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered test subdirectory")

if(OOMPH_WANT_NLOHMANN_JSON)
  add_executable(test_nlohmann_json test_nlohmann_json.cc)
  target_link_libraries(test_nlohmann_json PRIVATE nlohmann_json::nlohmann_json)
  add_test(NAME test_nlohmann_json COMMAND test_nlohmann_json)
endif()

message(VERBOSE "Leaving test subdirectory")
# ------------------------------------------------------------------------------
