#***********************************************************
# Module: algorithm
#***********************************************************
set( MODULE_MESH algorithm )

add_library( ${MODULE_MESH} INTERFACE )

target_include_directories( ${MODULE_MESH}
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} )

target_link_libraries( ${MODULE_MESH}
  INTERFACE util )
