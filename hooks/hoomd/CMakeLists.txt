# List of supported versions of HOOMD
set(SUPPORTED_HOOMD_VERSIONS
    "2.4.0" "2.4.1" "2.4.2"
    "2.5.0" "2.5.1" "2.5.2"
    "2.6.0"
    "2.7.0"
    "2.8.0" "2.8.1" "2.8.2"
    "2.9.0" "2.9.1" "2.9.2" "2.9.3"
    PARENT_SCOPE
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
    if ("${HOOMD_VERSION}" VERSION_LESS "2.4.0")
        message (FATAL_ERROR "\nERROR: HOOMD-blue versions before 2.4.0 are "
                             "not supported by SSAGES.\n")
    elseif (NOT ";${SUPPORTED_HOOMD_VERSIONS};" MATCHES ";${HOOMD_VERSION};")
        message (WARNING "\nWARNING: HOOMD-blue version \"${HOOMD_VERSION}\" is "
                           "not officially supported by SSAGES.\n")
    endif ()
    if ("${HOOMD_VERSION}" VERSION_LESS "2.6.0")
        add_definitions(-DHOOMD_PRE_2_6)
    endif ()
    target_link_libraries(ssages libssages ${HOOMD_LIBRARIES} ${PYTHON_LIBRARY})
endfunction()
