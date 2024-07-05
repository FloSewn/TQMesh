set( EXAMPLES run_examples )

add_executable( ${EXAMPLES}
  #example_9.cpp
  #example_8.cpp
  #example_7.cpp
  10_tqmesh_banner.cpp
  09_fixed_edges.cpp
  08_thin_fracture.cpp
  07_multiple_meshes.cpp
  06_airfoil_from_csv.cpp
  05_merge_meshes.cpp
  04_fixed_vertices.cpp
  03_boundary_shapes.cpp
  02_square_in_channel.cpp
  01_simple_triangular_mesh.cpp
  run_examples.cpp
  main.cpp
)

target_link_libraries( ${EXAMPLES} PRIVATE
  util
  algorithm
)

install( TARGETS ${EXAMPLES} RUNTIME DESTINATION ${BIN} )
