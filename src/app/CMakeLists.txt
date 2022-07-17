set( TQMESH_APP TQMesh )

add_executable( ${TQMESH_APP}
  size_function.cpp
  main.cpp
)

target_link_libraries( ${TQMESH_APP}
  util
  algorithm
  extern_libs
)

install( TARGETS ${TQMESH_APP} RUNTIME DESTINATION ${BIN} )
