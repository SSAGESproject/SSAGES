
# Define variables
set (BUILDSCHEMA_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/buildschema.py.in")
set (BUILDSCHEMA_PY "${CMAKE_CURRENT_BINARY_DIR}/buildschema.py")

set (SCHEMA_H "${CMAKE_CURRENT_BINARY_DIR}/schema.h")

# Replace CMAKE variables in the .in file with their current values
configure_file (${BUILDSCHEMA_PY_IN} ${BUILDSCHEMA_PY})

# The following command will generate the schema header and source files.
message (STATUS "Generating schema header and source files")

execute_process (COMMAND ${PYTHON_INTERP_EXE} ${BUILDSCHEMA_PY})

# The following command will make sure that the schema is updated on bulid time.
add_custom_target (schema
                   COMMAND ${PYTHON_INTERP_EXE} ${BUILDSCHEMA_PY}
                   COMMENT "Updating schema files")
