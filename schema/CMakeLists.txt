# Make sure, python is available
if (CMAKE_VERSION VERSION_LESS 3.12)
    find_package (PythonInterp REQUIRED)
    set (PYTHON_INTERP_FOUND ${PYTHONINTERP_FOUND})
    set (PYTHON_INTERP_EXE ${PYTHON_EXECUTABLE})
else ()
    find_package (Python REQUIRED COMPONENTS Interpreter)
    set (PYTHON_INTERP_FOUND ${Python_Interpreter_FOUND})
    set (PYTHON_INTERP_EXE ${Python_EXECUTABLE})
endif ()

if (PYTHON_INTERP_FOUND)
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
else ()
    message (SEND_ERROR "\nERROR: Unable to generate schema files. Python "
                        "binary not found.")
endif ()
