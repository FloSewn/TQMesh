set( TESTS run_tests )

add_executable( ${TESTS}
  tests_Vertex.cpp
  tests_Triangle.cpp
  tests_Quad.cpp
  tests_Front.cpp
  tests_EdgeList.cpp
  tests_Boundary.cpp
  tests_SizeFunction.cpp
  tests_Mesh.cpp
  tests_MeshGenerator.cpp
  tests_MeshCleanup.cpp
  tests_SmoothingStrategy.cpp
  tests.cpp
  main.cpp
)

target_link_libraries( ${TESTS} PRIVATE
  util
  algorithm
)

install( TARGETS ${TESTS} RUNTIME DESTINATION ${BIN} )

add_test(NAME Vertex COMMAND ${TESTS} "Vertex")
add_test(NAME Triangle COMMAND ${TESTS} "Triangle")
add_test(NAME Quad COMMAND ${TESTS} "Quad")
add_test(NAME Front COMMAND ${TESTS} "Front")
add_test(NAME EdgeList COMMAND ${TESTS} "EdgeList")
add_test(NAME Boundary COMMAND ${TESTS} "Boundary")
add_test(NAME SizeFunction COMMAND ${TESTS} "SizeFunction")
add_test(NAME Mesh COMMAND ${TESTS} "Mesh")
add_test(NAME MeshGenerator COMMAND ${TESTS} "MeshGenerator")
add_test(NAME MeshCleanup COMMAND ${TESTS} "MeshCleanup")
add_test(NAME SmoothingStrategy COMMAND ${TESTS} "SmoothingStrategy")
