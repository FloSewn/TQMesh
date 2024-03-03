#***********************************************************
# Module: util
#***********************************************************
set( MODULE_UTIL util )

add_library( ${MODULE_UTIL} INTERFACE )

target_include_directories( ${MODULE_UTIL}
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} )
