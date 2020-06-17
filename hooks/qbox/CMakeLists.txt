set(HOOKSRC
	hooks/qbox/QboxHook.cpp
    PARENT_SCOPE
)

set(DRIVERSRC
    hooks/qbox/Driver.cpp
    PARENT_SCOPE
)

function(CONFIGURE_ENGINE)
    add_dependencies(ssages libssages)
    include_directories(hooks/qbox/)
    target_link_libraries(ssages libssages)
endfunction()