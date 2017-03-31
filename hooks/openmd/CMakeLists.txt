set(OPENMD_CXX_ARGS "-std=c++11 -I${PROJECT_SOURCE_DIR}/src -I${PROJECT_SOURCE_DIR}/include -I${PROJECT_SOURCE_DIR}/hooks/openmd")

# Configure patch.
configure_file (
	"${CMAKE_CURRENT_SOURCE_DIR}/apply_patch.sh.in"
	"${CMAKE_CURRENT_BINARY_DIR}/apply_patch.sh"
)

ExternalProject_Add(
	openmd
	PREFIX ${CMAKE_CURRENT_BINARY_DIR}/openmd
	DOWNLOAD_COMMAND ""
	SOURCE_DIR ${OPENMD_SRC}
	PATCH_COMMAND sh ${CMAKE_CURRENT_BINARY_DIR}/apply_patch.sh
	CMAKE_ARGS
	-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
	-DCMAKE_CXX_FLAGS=${OPENMD_CXX_ARGS}
	-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
	-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
	BUILD_COMMAND $(MAKE) openmd_core openmd_parallel
	BUILD_IN_SOURCE false
	INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/openmd
	INSTALL_COMMAND cp
	-a ${CMAKE_CURRENT_BINARY_DIR}/openmd/src/openmd-build/lib/. ${CMAKE_CURRENT_BINARY_DIR}/
)

# Always build openmd
ExternalProject_Add_Step(
	openmd forcebuild
	COMMAND echo "Force build OpenMD"
	DEPENDEES configure
	DEPENDERS build
	ALWAYS 1
)

set(HOOKSRC
	hooks/openmd/OpenMDHook.cpp
	PARENT_SCOPE
)

set(ENABLE_OPENMD true PARENT_SCOPE)
set(MD_ENGINE "OpenMD" PARENT_SCOPE)
set(EXTRA_CONFIG "${CMAKE_CURRENT_BINARY_DIR}/openmd/src/openmd-build/config.h" PARENT_SCOPE)

# Function to do post library configuration.
function(CONFIGURE_ENGINE)
	include_directories(${OPENMD_SRC}/src)
	add_dependencies(libssages openmd)
	add_dependencies(ssages libssages openmd)
	set_target_properties(ssages PROPERTIES
	COMPILE_DEFINITIONS IS_MPI
	)
	target_link_libraries(ssages 
		libssages
		${CMAKE_CURRENT_BINARY_DIR}/libopenmd_core.a
		${CMAKE_CURRENT_BINARY_DIR}/libopenmd_parallel.a
		${CMAKE_CURRENT_BINARY_DIR}/libopenmd_core.a
		)
endfunction()
