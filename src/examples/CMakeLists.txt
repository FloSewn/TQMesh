set( EXAMPLES run_examples )

add_executable( ${EXAMPLES}
<<<<<<< HEAD
  example_9.cpp
  example_8.cpp
  example_7.cpp
  example_6.cpp
  example_5.cpp
  example_4.cpp
  example_3.cpp
  example_2.cpp
  example_1.cpp
=======
  #example_9.cpp
  #example_8.cpp
  #example_7.cpp
  07_tqmesh_banner.cpp
  06_airfoil_from_csv.cpp
  05_merge_meshes.cpp
  04_fixed_vertices.cpp
  03_boundary_shapes.cpp
  02_square_in_channel.cpp
  01_simple_triangular_mesh.cpp
>>>>>>> fa0899f5faedbc3de2d30dba4c8c9fc7b7288940
  run_examples.cpp
  main.cpp
)

target_link_libraries( ${EXAMPLES}
  util
  algorithm
)

install( TARGETS ${EXAMPLES} RUNTIME DESTINATION ${BIN} )
