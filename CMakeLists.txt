project(SSAGES)
cmake_minimum_required(VERSION 2.8.5)

# SAPHRON version number.
set(SSAGES_VERSION_MAJOR "0")
set(SSAGES_VERSION_MINOR "4")
set(SSAGES_VERSION_TINY "2") 

# Default type is release
if (NOT CMAKE_BUILD_TYPE)
    if (BUILD_TESTS)
		set(CMAKE_BUILD_TYPE "Debug")
    else (BUILD_TESTS)
		set(CMAKE_BUILD_TYPE "Release")
    endif (BUILD_TESTS)
endif()

string (TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE_NAME)

# List of supported versions of LAMMPS
set (SUPPORTED_LAMMPS_VERSIONS "10 Aug 2015"
                               "7 Dec 2015"
                               "14 May 2016"
                               "15 Jul 2016"
                               "30 Jul 2016")

# List of supported versions of Gromacs
set (SUPPORTED_GROMACS_VERSIONS "5.1.3")

# Resolve dependencies
find_package(MPI 1.6 REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Require Boost MPI.
find_package(Boost
	1.55
	REQUIRED
	COMPONENTS mpi serialization
)

# Set default values for variables
set(LAMMPS NO CACHE BOOL "Use LAMMPS simulation engine.")
set(GROMACS NO CACHE BOOL "Use GROMACS simulation engine.")
set(LAMMPS_SRC "" CACHE STRING "LAMMPS source directory.")
set(GROMACS_SRC "" CACHE STRING "GROMACS source directory.")

### Determine which simulation engine to use.
### Only one of the following values should be set on the command line:
###
### * LAMMPS=YES
### * GROMACS=YES
### * LAMMPS_SRC=/path/to/lammps/src -- Use LAMMPS, use local LAMMPS source.
### * GROMACS_SRC=/path/to/gromacs   -- Use GROMACS, use local GROMACS source.

## Check if new values have been given on the command line
if (LAMMPS_SRC AND NOT LAMMPS_SRC_LAST_RUN)
	# LAMMPS_SRC has been set
	set (LAMMPS NO CACHE BOOL "Use LAMMPS simulation engine." FORCE)
	set (GROMACS NO CACHE BOOL "Use GROMACS simulation engine." FORCE)
	set (GROMACS_SRC "" CACHE STRING "GROMACS source directory." FORCE)
endif (LAMMPS_SRC AND NOT LAMMPS_SRC_LAST_RUN)

if (GROMACS_SRC AND NOT GROMACS_SRC_LAST_RUN)
	# GROMACS_SRC has been set
	set (LAMMPS NO CACHE BOOL "Use GROMACS simulation engine." FORCE)
	set (GROMACS NO CACHE BOOL "Use GROMACS simulation engine." FORCE)
	set (LAMMPS_SRC "" CACHE BOOL "LAMMPS source directory." FORCE)
endif (GROMACS_SRC AND NOT GROMACS_SRC_LAST_RUN)

if (LAMMPS AND NOT LAMMPS_LAST_RUN)
	# LAMMPS has been set to YES
	set (GROMACS NO CACHE BOOL "Use GROMACS simulation engine." FORCE)
	set (LAMMPS_SRC "" CACHE BOOL "LAMMPS source directory." FORCE)
	set (GROMACS_SRC "" CACHE STRING "GROMACS source directory." FORCE)
endif (LAMMPS AND NOT LAMMPS_LAST_RUN)

if (GROMACS AND NOT GROMACS_LAST_RUN)
	# GROMACS has been set to YES
	set (LAMMPS NO CACHE BOOL "Use GROMACS simulation engine." FORCE)
	set (LAMMPS_SRC "" CACHE BOOL "LAMMPS source directory." FORCE)
	set (GROMACS_SRC "" CACHE STRING "GROMACS source directory." FORCE)
endif (GROMACS AND NOT GROMACS_LAST_RUN)

# Store current values of LAMMPS_SRC and GROMACS_SRC in the
# cache to be able to determine if they have been changed on the command line
set (LAMMPS_LAST_RUN ${LAMMPS} CACHE BOOL
	"Value of LAMMPS on last cmake run." FORCE)
set (GROMACS_LAST_RUN ${GROMACS} CACHE BOOL
	"Value of GROMACS on last cmake run." FORCE)
set (LAMMPS_SRC_LAST_RUN ${LAMMPS_SRC} CACHE BOOL
	"Use GROMACS simulation engine." FORCE)
set (GROMACS_SRC_LAST_RUN ${GROMACS_SRC} CACHE BOOL
	"Use GROMACS simulation engine." FORCE)

# Check that either LAMMPS_SRC or GROMACS_SRC have been set
if ( NOT LAMMPS AND NOT GROMACS AND NOT LAMMPS_SRC AND NOT GROMACS_SRC )
    message (FATAL_ERROR
             "\nERROR: SSAGES needs to be built with a simulation engine.\n"
             "The following simulation engines are supported by SSAGES:\n"
             "    - LAMMPS\n"
             "    - Gromacs\n"
             "Please call cmake with \"-DLAMMPS=YES\" or \"-DGROMACS=YES\" to "
             "build with SSAGES for LAMMPS and Gromacs, respectively.\n"
             "Alternatively, you can call with the option "
             "\"-DLAMMPS_SRC=/path/to/lammps/src\" or "
             "\"-DGROMACS_SRC=/path/to/gromacs\", where "
             "\"/path/to/lammps/src/\" and \"/path/to/gromacs/\" points to the "
             "local directory containing the source code for LAMMPS or Gromacs."
    )
endif ( NOT LAMMPS AND NOT GROMACS AND NOT LAMMPS_SRC AND NOT GROMACS_SRC )

include(ExternalProject)


include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}")

if(LAMMPS OR LAMMPS_SRC)
    if (LAMMPS_SRC)
		# Get absolute path in case user specified relative path
		if (IS_ABSOLUTE ${LAMMPS_SRC})
			get_filename_component(LAMMPS_SRC ${LAMMPS_SRC} ABSOLUTE)
		else (IS_ABSOLUTE ${LAMMPS_SRC})
			get_filename_component(LAMMPS_SRC ${CMAKE_CURRENT_BINARY_DIR}/${LAMMPS_SRC}
								   ABSOLUTE)
		endif (IS_ABSOLUTE ${LAMMPS_SRC})

		# Check if given LAMMPS directory exists.
		if (NOT EXISTS ${LAMMPS_SRC})
			message (FATAL_ERROR "\nERROR: Directory ${LAMMPS_SRC} does not exist.")
		endif (NOT EXISTS ${LAMMPS_SRC})

		# Check if directory does contain LAMMPS source files
		if (NOT EXISTS ${LAMMPS_SRC}/lammps.h)
			# Check if subdirectory was not specified
			if (EXISTS ${LAMMPS_SRC}/src/lammps.h)
				set (LAMMPS_SRC ${LAMMPS_SRC}/src)
			else (EXISTS ${LAMMPS_SRC}/src/lammps.h)
				# Directory does not contain LAMMPS source
				message (FATAL_ERROR "\nERROR: Directory ${LAMMPS_SRC} does not "
									 "contain LAMMPS source files (file lammps.h "
									 "not found).")
			endif (EXISTS ${LAMMPS_SRC}/src/lammps.h)
		endif (NOT EXISTS ${LAMMPS_SRC}/lammps.h)

		# Check LAMMPS version
		set (LAMMPS_VERSION "Version not found")
		execute_process (COMMAND "sed" "s/^[^\"]*\"\\([^\"]*\\)\".*/\\1/" "${LAMMPS_SRC}/version.h"
						 OUTPUT_VARIABLE LAMMPS_VERSION
						 OUTPUT_STRIP_TRAILING_WHITESPACE)

		# Check if LAMMPS version is supported
		if (NOT ";${SUPPORTED_LAMMPS_VERSIONS};" MATCHES ";${LAMMPS_VERSION};")
			message(WARNING "\nWARNING: Lammps version \"${LAMMPS_VERSION}\" is not "
							"officially supported by SSAGES.")
		endif (NOT ";${SUPPORTED_LAMMPS_VERSIONS};" MATCHES ";${LAMMPS_VERSION};")
	else (LAMMPS_SRC)
		set (LAMMPS_SRC "${CMAKE_CURRENT_BINARY_DIR}/lammps-download-prefix/src/lammps-download/src")
	endif (LAMMPS_SRC)

	set(ENABLE_LAMMPS true)

	if (LAMMPS)
        ExternalProject_Add(
            lammps-download
            DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}
			GIT_REPOSITORY git://git.lammps.org/lammps-ro.git
			GIT_TAG r15407
			TIMEOUT 600
			CONFIGURE_COMMAND ""
			BUILD_COMMAND ""
			INSTALL_COMMAND ""
		)
	endif (LAMMPS)

	# Add lammps as external project
	ExternalProject_Add(
		lammps
		DOWNLOAD_COMMAND ""
		SOURCE_DIR ${LAMMPS_SRC}
		CONFIGURE_COMMAND ""
		BUILD_COMMAND $(MAKE) mpi mode=shlib 
		CCFLAGS=-I${PROJECT_SOURCE_DIR}/src  
		CCFLAGS+=-I${PROJECT_SOURCE_DIR}/include  
		CCFLAGS+=-I${Boost_INCLUDE_DIRS}     
		CCFLAGS+=-std=c++11                  
		CCFLAGS+=-L${Boost_LIBRARY_DIRS}      
		LIB+=-lboost_mpi                      
		LIB+=-lboost_serialization           
		LIB+=${CMAKE_BINARY_DIR}/libssages.a
		INSTALL_COMMAND ""
		BUILD_IN_SOURCE true
	)

	ExternalProject_Add_Step(
		lammps forceconfigure
		COMMENT "Configuring LAMMPS"
		DEPENDERS build
		DEPENDEES update
		ALWAYS 1
	)

	add_custom_command(
	    TARGET lammps PRE_BUILD
	    COMMAND ln -sf ${PROJECT_SOURCE_DIR}/hooks/lammps/fix_ssages.h ${LAMMPS_SRC}/fix_ssages.h
	    COMMAND touch ${LAMMPS_SRC}/modify.cpp
	)

	include_directories("${LAMMPS_SRC}")

	# Add source for appropriate package hook.
	set(HOOKSRC
		hooks/lammps/fix_ssages.cpp
	)
elseif(GROMACS OR GROMACS_SRC)
    if (GROMACS_SRC)
		# Get absolute path in case user specified relative path
		if (IS_ABSOLUTE ${GROMACS_SRC})
			get_filename_component(GROMACS_SRC ${GROMACS_SRC} ABSOLUTE)
		else (IS_ABSOLUTE ${GROMACS_SRC})
			get_filename_component(GROMACS_SRC ${CMAKE_CURRENT_BINARY_DIR}/${GROMACS_SRC}
								   ABSOLUTE)
		endif (IS_ABSOLUTE ${GROMACS_SRC})

		# Check if directory GROMACS_SRC exists
		if (NOT EXISTS ${GROMACS_SRC})
			message (FATAL_ERROR "ERROR: Directory ${GROMACS_SRC} does not exist.")
		endif (NOT EXISTS ${GROMACS_SRC})

		# Get Gromacs versions
		set (GROMACS_VERSION "Version not found")
		set (GROMACS_VERSION_INFO ${GROMACS_SRC}/cmake/gmxVersionInfo.cmake)
		if (NOT EXISTS ${GROMACS_VERSION_INFO})
			set (GROMACS_VERSION "< 5.1")
		else (NOT EXISTS ${GROMACS_VERSION_INFO})
			execute_process (
				COMMAND "sed" "-n"
				"s/set(GMX_VERSION_\\(MAJOR\\|MINOR\\|PATCH\\) \\([0-9]\\+\\))/\\2/p"
				"${GROMACS_SRC}/cmake/gmxVersionInfo.cmake"
				COMMAND "sed" "-e" ":a" "-e" "N" "-e" "\$!ba" "-e" "s/\\n/./g"
				OUTPUT_VARIABLE GROMACS_VERSION
				OUTPUT_STRIP_TRAILING_WHITESPACE
			)
		endif (NOT EXISTS ${GROMACS_VERSION_INFO})

		# Check if Gromacs version is supported
		if (NOT ";${SUPPORTED_GROMACS_VERSIONS};" MATCHES ";${GROMACS_VERSION};")
			message (WARNING "\nWARNING: Gromacs version ${GROMACS_VERSION} is "
							 "currently not supported by SSAGES.")
		endif (NOT ";${SUPPORTED_GROMACS_VERSIONS};" MATCHES ";${GROMACS_VERSION};")

		# Check that files to be patched do exist
		execute_process (COMMAND "sed" "-n"
								 "s/^.*gromacs-ssages\\/\\([^t]*\\)\t.*/\\1/p"
								 "${PROJECT_SOURCE_DIR}/hooks/gromacs/gmx_diff.patch"
						 COMMAND "sed" "-e" ":a" "-e" "N" "-e" "\$!ba" "-e" "s/\\n/;/g"
						 OUTPUT_VARIABLE GROMACS_FILES_TO_BE_PATCHED
						 OUTPUT_STRIP_TRAILING_WHITESPACE)

		foreach (PATCH_FILE ${GROMACS_FILES_TO_BE_PATCHED})
			if (NOT EXISTS ${GROMACS_SRC}/${PATCH_FILE})
				message (FATAL_ERROR "\nERROR: File ${PATCH_FILE} does not exist.\n"
									 "This error is probably caused by an unsupported "
									 "version of Gromacs.")
			endif (NOT EXISTS ${GROMACS_SRC}/${PATCH_FILE})
		endforeach (PATCH_FILE IN LISTS ${GROMACS_FILES_TO_BE_PATCHED})
	else (GROMACS_SRC)
		set (GROMACS_SRC ${CMAKE_CURRENT_BINARY_DIR}/gromacs-download-prefix/src/gromacs-download)
	endif (GROMACS_SRC)

	set(ENABLE_GROMACS true)

	set(GRO_CXX_ARGS "-std=c++11 -I${PROJECT_SOURCE_DIR}/src -I${PROJECT_SOURCE_DIR}/include -I${PROJECT_SOURCE_DIR}/hooks/gromacs")

	# Configure patch.
	configure_file (
		"${PROJECT_SOURCE_DIR}/hooks/gromacs/apply_patch.sh.in"
		"${PROJECT_BINARY_DIR}/apply_patch.sh"
	)

	if (GROMACS)
        ExternalProject_Add(
			gromacs-download
			DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}/gromacs_source/
			GIT_REPOSITORY git://git.gromacs.org/gromacs.git
			GIT_TAG v5.1.3
			TIMEOUT 600
			CONFIGURE_COMMAND ""
			BUILD_COMMAND ""
			INSTALL_COMMAND ""
		)
	endif (GROMACS)

	ExternalProject_Add(
		gromacs
		DOWNLOAD_COMMAND ""
		BUILD_COMMAND $(MAKE)
		SOURCE_DIR ${GROMACS_SRC}
		PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gromacs
		BUILD_ALWAYS true
		CMAKE_ARGS
		-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
		-DGMX_MPI=ON
		-DCMAKE_CXX_FLAGS=${GRO_CXX_ARGS}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		-DGMX_EXTRA_LIBRARIES=${CMAKE_BINARY_DIR}/libssagesgmx.so
		-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
		-DCUDA_PROPAGATE_HOST_FLAGS=OFF
		-DCUDA_NVCC_FLAGS="-std=c++11 -I${PROJECT_SOURCE_DIR}/src -I${PROJECT_SOURCE_DIR}/include -I${PROJECT_SOURCE_DIR}/hooks/gromacs"
		BUILD_IN_SOURCE false
		INSTALL_DIR ${CMAKE_BINARY_DIR}/gromacs
		PATCH_COMMAND sh ${CMAKE_BINARY_DIR}/apply_patch.sh
	)

	set(HOOKSRC
		hooks/gromacs/GromacsHook.cpp
	)

	include_directories(${GROMACS_SRC}/src)
endif()

# Schema builder.
add_subdirectory (schema)

# General source files
set(SOURCE
	src/Methods/Meta.cpp
	src/Methods/Umbrella.cpp
	src/Methods/ElasticBand.cpp
	src/Methods/FiniteTempString.cpp
    	src/Methods/Swarm.cpp
	src/Methods/ForwardFlux.cpp
	src/Methods/ABF.cpp
    	src/Methods/BasisFunc.cpp
	src/Validator/RequirementLoader.cpp
	src/JSON/jsoncpp.cpp
	src/CVs/CollectiveVariable.cpp
	src/Methods/Method.cpp
	src/Grids/Grid.cpp
	src/Observers/JSONObserver.cpp
	src/Simulations/SimObserver.cpp
	src/Simulations/SimObservable.cpp
	src/Hook.cpp
	src/Constraints/Constraint.cpp
	src/Constraints/COPSSImage.cpp
)

# Include the directory with the newly generated schema header file
include_directories(${CMAKE_CURRENT_BINARY_DIR}/schema/)

include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories(${Boost_INCLUDE_DIRS})
add_definitions(${MPI_CXX_COMPILE_FLAGS})
link_directories(${MPI_CXX_LIBRARIES})

# SSAGES library.
add_library(libssages ${HOOKSRC} ${SOURCE})
set_target_properties(libssages PROPERTIES PREFIX "")
target_link_libraries(libssages ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set_target_properties(libssages PROPERTIES COMPILE_FLAGS "-Wall -Wextra -fPIC -g")
elseif (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set_target_properties(libssages PROPERTIES COMPILE_FLAGS "-O2 -fPIC -g")
else ()
    set_target_properties(libssages PROPERTIES COMPILE_FLAGS "-O3 -fPIC")
endif ()

# Require C++11.
set_property(TARGET libssages PROPERTY CXX_STANDARD 11)
set_property(TARGET libssages PROPERTY CXX_STANDARD_REQUIRED ON)

# SSAGES executable.
add_executable(ssages
	src/main.cpp
)

# Configure header.
configure_file (
"${PROJECT_SOURCE_DIR}/src/config.h.in"
"${PROJECT_BINARY_DIR}/config.h"
)

# Dependencies
add_dependencies(libssages schema)
if(LAMMPS OR LAMMPS_SRC)
	add_dependencies(lammps libssages)
	if (LAMMPS)
		add_dependencies(libssages lammps-download)
	endif (LAMMPS)
	target_link_libraries(ssages libssages ${LAMMPS_SRC}/liblammps_mpi.so)
	add_dependencies(ssages lammps)
elseif(GROMACS OR GROMACS_SRC)
	
	# Build custom shared library for dependencies. 
	add_library(ssagesgmx SHARED ${HOOKSRC} ${SOURCE})
	target_link_libraries(ssagesgmx ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})

	add_dependencies(gromacs ssagesgmx)

	if (GROMACS)
		add_dependencies(ssagesgmx gromacs-download)
	endif (GROMACS)

	# To determine Gromacs install dir. 
	include(GNUInstallDirs)
	target_link_libraries(ssages 
		ssagesgmx 
		${CMAKE_BINARY_DIR}/gromacs/${CMAKE_INSTALL_LIBDIR}/libgromacs_mpi.so
		${CMAKE_BINARY_DIR}/gromacs/${CMAKE_INSTALL_LIBDIR}/libmdrun.a
		)
	add_dependencies(ssages gromacs)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set_target_properties(ssages PROPERTIES COMPILE_FLAGS "-Wall -Wextra -g")
elseif (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set_target_properties(ssages PROPERTIES COMPILE_FLAGS "-O2 -g")
else ()
    set_target_properties(ssages PROPERTIES COMPILE_FLAGS "-O3")

endif ()

# Configure flags 
# - Globals
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -pthread -std=c++11")

# Add doc subdirectory for building the documentation
add_subdirectory(doc)

if (BUILD_TESTS)
	message(STATUS "")

	# Find Google Test Framework
	include(ExternalProject)

	# Set default ExternalProject root directory
	set_directory_properties(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/extern)

	# Add gtest
	# http://stackoverflow.com/questions/9689183/cmake-googletest
	ExternalProject_Add(
		googletest
        DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}
		GIT_REPOSITORY https://github.com/google/googletest.git
		GIT_TAG 800f542
		# TIMEOUT 10
		# # Force separate output paths for debug and release builds to allow easy
		# # identification of correct lib in subsequent TARGET_LINK_LIBRARIES commands
		# CMAKE_ARGS -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=DebugLibs
		# -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=ReleaseLibs
		# -Dgtest_force_shared_crt=ON
		# Disable install step
		INSTALL_COMMAND ""
		# Wrap download, configure and build steps in a script to log output
		LOG_DOWNLOAD ON
		LOG_CONFIGURE ON
		LOG_BUILD ON
	)

	# Specify include dir
	ExternalProject_Get_Property(googletest source_dir)
	set(GTEST_INCLUDE_DIR ${source_dir}/include)

	ExternalProject_Get_Property(googletest binary_dir)
	set(GTEST_LIBRARY_PATH ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a)
	set(GTEST_LIBRARY_MAIN_PATH ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main.a)
	set(GTEST_LIBRARY gtest)
	add_library(${GTEST_LIBRARY} UNKNOWN IMPORTED)
	set_property(TARGET ${GTEST_LIBRARY} PROPERTY IMPORTED_LOCATION
	${GTEST_LIBRARY_PATH} )
	add_dependencies(${GTEST_LIBRARY} googletest)

	message(STATUS "")
	message(STATUS "*** Google Test Framework will be used for unit tests")
	message(STATUS "*** GTEST_LIBRARY_PATH = ${GTEST_LIBRARY_PATH}")
	message(STATUS "*** GTEST_INCLUDE_DIR  = ${GTEST_INCLUDE_DIR}")
	message(STATUS "")
	# End gtest

	enable_testing()
	add_subdirectory(test)
endif (BUILD_TESTS)
