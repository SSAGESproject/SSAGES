# Make sure, python is available
find_program (PYTHON_BIN "python")

if (PYTHON_BIN)
    # Define variables
    set (BUILDSCHEMA_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/buildschema.py.in")
    set (BUILDSCHEMA_PY "${CMAKE_CURRENT_BINARY_DIR}/buildschema.py")

    set (SCHEMA_H "${CMAKE_CURRENT_BINARY_DIR}/schema.h")

    # Replace CMAKE variables in the .in file with their current values
    configure_file (${BUILDSCHEMA_PY_IN} ${BUILDSCHEMA_PY})

    # The following command will generate the schema header and source files.
    message (STATUS "Generating schema header and source files")
    execute_process (COMMAND ${PYTHON_BIN} ${BUILDSCHEMA_PY})

else (PYTHON_BIN)
    message (SEND_ERROR "Unable to generate schema files. Python binary not found.")
endif (PYTHON_BIN)