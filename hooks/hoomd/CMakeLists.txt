# List of supported versions of HOOMD
#set(SUPPORTED_HOOMD_VERSIONS "2.x")

if (HOOMD_SRC)
    # Get absolute path in case user specified relative path
    if (IS_ABSOLUTE ${HOOMD_SRC})
        get_filename_component(HOOMD_SRC ${HOOMD_SRC} ABSOLUTE)
    else ()
        get_filename_component(HOOMD_SRC ${PROJECT_BINARY_DIR}/${HOOMD_SRC}
                                ABSOLUTE)
    endif ()

    # Check if given HOOMD directory exists.
    if (NOT EXISTS ${HOOMD_SRC})
        message (FATAL_ERROR "\nERROR: Directory ${HOOMD_SRC} does not exist.")
    endif ()

    # Check if directory does contain HOOMD source files
    if (NOT EXISTS ${HOOMD_SRC}/ExecutionConfiguration.h)
        # Check if subdirectory was not specified
        if (EXISTS ${HOOMD_SRC}/hoomd/ExecutionConfiguration.h)
            set(HOOMD_SRC ${HOOMD_SRC}/hoomd PARENT_SCOPE)
        else (EXISTS ${HOOMD_SRC}/hoomd/ExecutionConfiguration.h)
            # Directory does not contain HOOMD source
            message (FATAL_ERROR "\nERROR: Directory ${HOOMD_SRC} does not "
                                    "contain HOOMD source files (file ExecutionConfiguration.h "
                                    "not found).")
        endif (EXISTS ${HOOMD_SRC}/hoomd/ExecutionConfiguration.h)
    endif (NOT EXISTS ${HOOMD_SRC}/ExecutionConfiguration.h)
else ()
    set(HOOMD_SRC ${HOOMD_INCLUDE_DIR} PARENT_SCOPE)
endif ()

set(HOOKSRC
    hooks/hoomd/HOOMDHook.cpp
    PARENT_SCOPE
)

set(DRIVERSRC
    hooks/hoomd/Driver.cpp
    PARENT_SCOPE
)

set(MD_ENGINE "HOOMD" PARENT_SCOPE)
set(ENABLE_HOOMD true PARENT_SCOPE)

# Function to do post library configuration.
function(CONFIGURE_ENGINE)
    include_directories(hooks/hoomd)
    target_link_libraries(ssages libssages ${HOOMD_LIBRARIES} ${PYTHON_LIBRARY})
endfunction()
