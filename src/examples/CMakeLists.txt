set( EXAMPLES run_examples )

add_executable( ${EXAMPLES}
  #example_9.cpp
  #example_8.cpp
  #example_7.cpp
  tqmesh_banner.cpp
  merge_meshes.cpp
  fixed_vertices.cpp
  boundary_shapes.cpp
  square_in_channel.cpp
  simple_triangular_mesh.cpp
  run_examples.cpp
  main.cpp
)

target_link_libraries( ${EXAMPLES}
  util
  algorithm
)

install( TARGETS ${EXAMPLES} RUNTIME DESTINATION ${BIN} )
