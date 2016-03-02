project(SSAGES)
cmake_minimum_required(VERSION 2.8.5)

# Default type is release
if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release")
endif()

# Resolve dependencies
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Require Boost MPI.
find_package(Boost
	1.55
	REQUIRED
	COMPONENTS mpi serialization
)

# Get lammps source directory
set(LAMMPS_SRC "" CACHE STRING "LAMMPS source directory")

include(ExternalProject)

# Add lammps as external project
ExternalProject_Add(
	lammps
	SOURCE_DIR ${LAMMPS_SRC}
	CONFIGURE_COMMAND ""
	BUILD_COMMAND $(MAKE) mpi mode=shlib 
	CCFLAGS=-I${PROJECT_SOURCE_DIR}/src  
	CCFLAGS+=-I${Boost_INCLUDE_DIRS}     
	CCFLAGS+=-std=c++11                  
	CCFLAGS+=-L${Boost_LIBRARY_DIRS}      
	LIB+=-lboost_mpi                      
	LIB+=-lboost_serialization           
	LIB+=${PROJECT_SOURCE_DIR}/build/libssages.a
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

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/include")

# Include directory for the pacakge we're interested in. 
include_directories("${LAMMPS_SRC}")

# Add source for appropriate package hook.
set(HOOKSRC
	hooks/lammps/fix_ssages.cpp
)

# General source files
set(SOURCE
	src/Methods/Meta.cpp
	src/Methods/Umbrella.cpp
	src/Methods/ElasticBand.cpp
)

include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories(${Boost_INCLUDE_DIRS})
add_definitions(${MPI_CXX_COMPILE_FLAGS})
link_directories(${MPI_CXX_LIBRARIES})

# SSAGES library.
add_library(libssages ${HOOKSRC} ${SOURCE})
set_target_properties(libssages PROPERTIES PREFIX "")
target_link_libraries(libssages ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
set_target_properties(libssages PROPERTIES COMPILE_FLAGS "-Wall -Wextra -fPIC")

# Require C++11.
set_property(TARGET libssages PROPERTY CXX_STANDARD 11)
set_property(TARGET libssages PROPERTY CXX_STANDARD_REQUIRED ON)

# Dependency.
add_dependencies(lammps libssages)

# SSAGES executable.
add_executable(ssages
	src/main.cpp
)
target_link_libraries(ssages libssages ${LAMMPS_SRC}/liblammps_mpi.so)
add_dependencies(ssages lammps)
set_target_properties(ssages PROPERTIES COMPILE_FLAGS "-Wall -Wextra")

# Configure flags 
# - Globals
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -pthread -std=c++11")
