# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered reaction_diffusion subdirectory")

set(SUBDIRS one_d_act_inhibit two_d_act_inhibit two_d_act_inhibit_adapt
    unstructured_two_d_act_inhibit)

foreach(SUBDIR ${SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving reaction_diffusion subdirectory")
# ------------------------------------------------------------------------------
