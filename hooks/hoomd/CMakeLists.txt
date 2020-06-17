# List of supported versions of HOOMD
set(SUPPORTED_HOOMD_VERSIONS
    "2.4.0" "2.4.1" "2.4.2"
    "2.5.0" "2.5.1" "2.5.2"
)

if (HOOMD_ROOT)
    if (IS_ABSOLUTE ${HOOMD_ROOT})
        get_filename_component(HOOMD_ROOT ${HOOMD_ROOT} ABSOLUTE)
    else ()
        get_filename_component(HOOMD_ROOT ${PROJECT_BINARY_DIR}/${HOOMD_ROOT}
                               ABSOLUTE)
    endif ()
    set(HOOMD_ROOT ${HOOMD_ROOT} PARENT_SCOPE)
endif ()

execute_process (
    COMMAND "sed" "-n"
    "s/\#define HOOMD_VERSION_\\(MAJOR\\|MINOR\\|PATCH\\) \\([0-9]\\+\\)/\\2/p"
    "${HOOMD_ROOT}/include/HOOMDVersion.h"
    COMMAND "sed" "-e" ":a" "-e" "N" "-e" "\$!ba" "-e" "s/\\n/./g"
    OUTPUT_VARIABLE HOOMD_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if (NOT ";${SUPPORTED_HOOMD_VERSIONS};" MATCHES ";${HOOMD_VERSION};")
    message (FATAL_ERROR "\nERROR: HOOMD-blue version \"${HOOMD_VERSION}\" is "
                           "currently not supported by SSAGES.")
endif ()

set(HOOKSRC
    hooks/hoomd/HOOMDHook.cpp
    PARENT_SCOPE
)

set(DRIVERSRC
    hooks/hoomd/Driver.cpp
    PARENT_SCOPE
)

# Function to do post library configuration.
function(CONFIGURE_ENGINE)
    include_directories(hooks/hoomd)
    include(hooks/hoomd/FindHOOMD.cmake)
    target_link_libraries(ssages libssages ${HOOMD_LIBRARIES} ${PYTHON_LIBRARY})
endfunction()
