# List of supported versions of OpenMD
set(SUPPORTED_OPENMD_VERSIONS
    "2.5.0" "2.6.0"
)

# Get absolute path in case user specified relative path
if (IS_ABSOLUTE ${OPENMD_SRC})
    get_filename_component(OPENMD_SRC ${OPENMD_SRC} ABSOLUTE)
else ()
    get_filename_component(OPENMD_SRC ${PROJECT_BINARY_DIR}/${OPENMD_SRC}
                           ABSOLUTE)
endif ()

# Check if directory OPENMD_SRC exists
if (NOT EXISTS ${OPENMD_SRC})
	message (FATAL_ERROR "\nERROR: Directory ${OPENMD_SRC} does not exist.")
endif (NOT EXISTS ${OPENMD_SRC})

# Check if directory contains OpenMD source files
if (NOT EXISTS ${OPENMD_SRC}/src/applications/openmd/openmd.cpp)
	# Check if subdirectory was specified
	if (EXISTS ${OPENMD_SRC}/applications/openmd/openmd.cpp)
		get_filename_component (OPENMD_SRC ${OPENMD_SRC}/.. ABSOLUTE)
	else ()
		# Directory does not contain OpenMD source
		message (FATAL_ERROR "\nERROR: Directory ${OPENMD_SRC} does not "
		                       "contain OpenMD source files "
		                       "(file src/applications/openmd/openmd.cpp "
		                       "not found).")
	endif ()
endif ()

execute_process (
	COMMAND "sed" "-n"
	"s/set (VERSION_\\(MAJOR\\|MINOR\\|TINY\\) \"\\([0-9]\\+\\)\")/\\2/p"
	"${OPENMD_SRC}/CMakeLists.txt"
	COMMAND "sed" "-e" ":a" "-e" "N" "-e" "\$!ba" "-e" "s/\\n/./g"
	OUTPUT_VARIABLE OPENMD_VERSION
	OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Check if OpenMD version is supported
if (NOT ";${SUPPORTED_OPENMD_VERSIONS};" MATCHES ";${OPENMD_VERSION};")
	message (FATAL_ERROR "\nERROR: OpenMD version \"${OPENMD_VERSION}\" is "
	"currently not supported by SSAGES.")
endif (NOT ";${SUPPORTED_OPENMD_VERSIONS};" MATCHES ";${OPENMD_VERSION};")

set (MD_ENGINE_VERSION "OpenMD version ${OPENMD_VERSION}" PARENT_SCOPE)

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

set(DRIVERSRC
    hooks/openmd/Driver.cpp
    PARENT_SCOPE
)

# Libraries from the engine to link
set(ENGINE_LIBRARIES
    "${CMAKE_BINARY_DIR}/hooks/openmd/libopenmd_core.a"
    "${CMAKE_BINARY_DIR}/hooks/openmd/libopenmd_parallel.a"
    "${CMAKE_BINARY_DIR}/hooks/openmd/libopenmd_core.a"
    PARENT_SCOPE
)

set(EXTRA_CONFIG "${CMAKE_CURRENT_BINARY_DIR}/openmd/src/openmd-build/config.h" PARENT_SCOPE)

# Function to do post library configuration.
function(CONFIGURE_ENGINE)
	include_directories(${OPENMD_SRC}/src)
	include_directories(hooks/openmd/)
	add_dependencies(libssages openmd)
	add_dependencies(ssages openmd)
	set_target_properties(ssages PROPERTIES
		COMPILE_DEFINITIONS IS_MPI
	)
endfunction()
