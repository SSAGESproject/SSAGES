# List of supported versions of HOOMD
#set(SUPPORTED_HOOMD_VERSIONS "2.x")

set(HOOMD_SRC ${HOOMD_INCLUDE_DIR} PARENT_SCOPE)

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
