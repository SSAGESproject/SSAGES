# List of supported versions of HOOMD
#set(SUPPORTED_HOOMD_VERSIONS "2.x")

if (IS_ABSOLUTE ${HOOMD_ROOT})
    get_filename_component(HOOMD_ROOT ${HOOMD_ROOT} ABSOLUTE)
else ()
    get_filename_component(HOOMD_ROOT ${PROJECT_BINARY_DIR}/${HOOMD_ROOT}
                           ABSOLUTE)
endif ()
set(HOOMD_ROOT ${HOOMD_ROOT} PARENT_SCOPE)

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
    include(hooks/hoomd/FindHOOMD.cmake)
    target_link_libraries(ssages libssages ${HOOMD_LIBRARIES} ${PYTHON_LIBRARY})
endfunction()
