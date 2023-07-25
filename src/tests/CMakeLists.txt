set( TESTS run_tests )

add_executable( ${TESTS}
  tests_Vertex.cpp
  tests_Triangle.cpp
  tests_Quad.cpp
  tests_Front.cpp
  tests_EdgeList.cpp
  tests_Boundary.cpp
  tests_SizeFunction.cpp
  tests_Smoother.cpp
  tests_Mesh.cpp
  tests_MeshGenerator.cpp
  tests_Cleanup.cpp
  tests.cpp
  main.cpp
)

target_link_libraries( ${TESTS}
  util
  algorithm
)

install( TARGETS ${TESTS} RUNTIME DESTINATION ${BIN} )
