set(HOOKSRC
	src/Drivers/QboxDriver.cpp
	hooks/qbox/QboxHook.cpp
    PARENT_SCOPE
)

set(MD_ENGINE "QBOX" PARENT_SCOPE)
set(ENABLE_QBOX true PARENT_SCOPE)

function(CONFIGURE_ENGINE)
    add_dependencies(ssages libssages)
    target_link_libraries(ssages libssages)
endfunction()