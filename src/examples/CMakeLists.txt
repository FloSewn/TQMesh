set( EXAMPLES run_examples )

add_executable( ${EXAMPLES}
  example_8.cpp
  example_7.cpp
  example_6.cpp
  example_5.cpp
  example_4.cpp
  example_3.cpp
  example_2.cpp
  example_1.cpp
  run_examples.cpp
  main.cpp
)

target_link_libraries( ${EXAMPLES}
  util
  algorithm
)

install( TARGETS ${EXAMPLES} RUNTIME DESTINATION ${BIN} )
