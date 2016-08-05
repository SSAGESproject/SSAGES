project(SSAGES)
cmake_minimum_required(VERSION 2.8.5)

# SAPHRON version number.
set(SSAGES_VERSION_MAJOR "0")
set(SSAGES_VERSION_MINOR "0")
set(SSAGES_VERSION_TINY "8") 

# Default type is release
if (NOT CMAKE_BUILD_TYPE)
	if (TEST)
		set(CMAKE_BUILD_TYPE "Debug")
	else (TEST)
		set(CMAKE_BUILD_TYPE "Release")
	endif (TEST)
endif()

string (TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE_NAME)

# list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Resolve dependencies
find_package(MPI 1.6 REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Require Boost MPI.
find_package(Boost
	1.55
	REQUIRED
	COMPONENTS mpi serialization
)

# Get lammps source directory
set(LAMMPS_SRC "" CACHE STRING "LAMMPS source directory")
set(GROMACS_SRC "" CACHE STRING "GROMACS source directory")

include(ExternalProject)

if(LAMMPS_SRC)
	set(ENABLE_LAMMPS true)
	# Add lammps as external project
	ExternalProject_Add(
		lammps
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
		LOG_BUILD ON
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
elseif(GROMACS_SRC)
		set(ENABLE_GROMACS true)

		ExternalProject_Add(
		gromacs
		SOURCE_DIR ${GROMACS_SRC}
		PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gromacs
		CMAKE_ARGS 
		-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> 
		-DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON 
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		BUILD_IN_SOURCE false
		INSTALL_DIR ${CMAKE_BINARY_DIR}/gromacs
		LOG_BUILD ON)
	#include_directories(${GROMACS_INCLUDE_DIRS})
	#add_definitions(${GROMACS_DEFINITIONS})
endif()

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}")


# Schema builder.
set (SCHEMA_CPP ${CMAKE_CURRENT_BINARY_DIR}/schema/schema.cpp)
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
	${SCHEMA_CPP}
	src/CVs/CollectiveVariable.cpp
	src/Methods/Method.cpp
	src/Grids/Grid.cpp
	src/Observers/JSONObserver.cpp
	src/Simulations/SimObservable.cpp
	src/Simulations/SimObserver.cpp
	src/Hook.cpp
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
    set_target_properties(libssages PROPERTIES COMPILE_FLAGS "-Wall -Wextra -fPIC")
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
if(LAMMPS_SRC)
	add_dependencies(lammps libssages)
	target_link_libraries(ssages libssages ${LAMMPS_SRC}/liblammps_mpi.so)
	add_dependencies(ssages lammps)
elseif(GROMACS_SRC)
	add_dependencies(gromacs libssages)
	target_link_libraries(ssages libssages ${CMAKE_BINARY_DIR}/gromacs/lib64/libgromacs.so)
	add_dependencies(ssages gromacs)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set_target_properties(ssages PROPERTIES COMPILE_FLAGS "-Wall -Wextra")
else ()
    set_target_properties(ssages PROPERTIES COMPILE_FLAGS "-O3")

endif ()

# Configure flags 
# - Globals
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -pthread -std=c++11")

# Add doc subdirectory for building the documentation
add_subdirectory(doc)

if (TEST)
	message(STATUS "")

	# Find Google Test Framework
	include(ExternalProject)

	# Set default ExternalProject root directory
	set_directory_properties(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/extern)

	# Add gtest
	# http://stackoverflow.com/questions/9689183/cmake-googletest
	ExternalProject_Add(
		googletest
		URL http://googletest.googlecode.com/files/gtest-1.6.0.zip
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
endif (TEST)
