set( TESTS TQMesh )

add_executable( ${TESTS}
  size_function.cpp
  main.cpp
)

target_link_libraries( ${TESTS}
  util
  algorithm
)

install( TARGETS ${TESTS} RUNTIME DESTINATION ${BIN} )
