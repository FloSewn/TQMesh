set( EXAMPLES run_examples )

add_executable( ${EXAMPLES}
  example_1.cpp
  run_examples.cpp
  main.cpp
)

target_link_libraries( ${EXAMPLES}
  util
  algorithm
)

install( TARGETS ${EXAMPLES} RUNTIME DESTINATION ${BIN} )
