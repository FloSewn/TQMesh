set( TESTS run_tests )

add_executable( ${TESTS}
  geometry_test.cpp
  qtree_test.cpp
  container_test.cpp
  vertex_test.cpp
  edgelist_test.cpp
  boundary_test.cpp
  sizefunction_test.cpp
  front_test.cpp
  triangle_test.cpp
  quad_test.cpp
  mesh_test.cpp
  smoother_test.cpp
  main.cpp
)

target_link_libraries( ${TESTS}
  util
  algorithm
)

install( TARGETS ${TESTS} RUNTIME DESTINATION ${BIN} )
