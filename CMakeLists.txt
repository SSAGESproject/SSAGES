project(SSAGES)
cmake_minimum_required(VERSION 2.8.5)

# Default type is release
if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release")
endif()

# Get lammps source directory
set(LAMMPS_SRC "" CACHE STRING "LAMMPS source directory")

# Build LAMMPS.
include(ExternalProject)

ExternalProject_Add(
	lammps
	SOURCE_DIR ${LAMMPS_SRC}
	CONFIGURE_COMMAND ""
	BUILD_COMMAND $(MAKE) mpi CCFLAGS=-I${PROJECT_SOURCE_DIR}/src CCFLAGS+=-std=c++11
	INSTALL_COMMAND ""
	BUILD_IN_SOURCE true
	LOG_BUILD ON
)

add_custom_command(
    TARGET lammps PRE_BUILD
    COMMAND ln -sf ${PROJECT_SOURCE_DIR}/hooks/lammps/fix_ssages.h ${LAMMPS_SRC}/fix_ssages.h
    COMMAND ln -sf ${PROJECT_SOURCE_DIR}/hooks/lammps/fix_ssages.cpp ${LAMMPS_SRC}/fix_ssages.cpp
)	

