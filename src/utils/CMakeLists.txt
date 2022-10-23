#***********************************************************
# Module: util
#***********************************************************
set( MODULE_UTIL util )

add_library( ${MODULE_UTIL} INTERFACE )

target_include_directories( ${MODULE_UTIL}
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} )

if (NOT WIN32)
target_link_libraries( ${MODULE_UTIL}
  INTERFACE m )
endif()

