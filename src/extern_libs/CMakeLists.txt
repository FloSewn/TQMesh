set( MODULE_EXTERN extern_libs )

add_library( ${MODULE_EXTERN} INTERFACE )

target_include_directories( ${MODULE_EXTERN} INTERFACE . )

install( TARGETS ${MODULE_EXTERN} RUNTIME DESTINATION ${BIN} )
