project(SSAGES)
cmake_minimum_required(VERSION 2.8.8)

# SSAGES version number.
file(STRINGS "VERSION" SSAGES_VERSION)

# Default type is release
if (NOT CMAKE_BUILD_TYPE)
    if (BUILD_TESTS)
        set(CMAKE_BUILD_TYPE "Debug")
    else (BUILD_TESTS)
        set(CMAKE_BUILD_TYPE "Release")
    endif (BUILD_TESTS)
endif()

string (TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE_NAME)

# Set the custom modules path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Resolve dependencies
find_package(Git REQUIRED)
find_package(MPI 1.6 REQUIRED)

# Make sure Python is available
if (CMAKE_VERSION VERSION_LESS 3.12)
    find_package (PythonInterp REQUIRED)
    set (PYTHON_INTERP_FOUND ${PYTHONINTERP_FOUND})
    set (PYTHON_INTERP_EXE ${PYTHON_EXECUTABLE})
else ()
    find_package (Python REQUIRED COMPONENTS Interpreter)
    set (PYTHON_INTERP_FOUND ${Python_Interpreter_FOUND})
    set (PYTHON_INTERP_EXE ${Python_EXECUTABLE})
endif ()
if (NOT PYTHON_INTERP_FOUND)
    message (SEND_ERROR "\nERROR: Python binary not found!")
endif ()

# Fetch external dependencies
include(ExternalProject)
include(FetchEigen)
include(FetchJsonCpp)
include(FetchMxx)

# Schema builder.
add_subdirectory(schema)

include_directories(${MPI_INCLUDE_PATH})
include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}")
include_directories(${CMAKE_CURRENT_BINARY_DIR}/schema/)
add_definitions(${MPI_CXX_COMPILE_FLAGS})

# General source files
set(SOURCE
    src/Methods/ABF.cpp
    src/Methods/ANN.cpp
    src/Methods/BasisFunc.cpp
    src/Methods/CFF.cpp
    src/Methods/DirectForwardFlux.cpp
    src/Methods/ElasticBand.cpp
    src/Methods/FiniteTempString.cpp
    src/Methods/ForwardFlux.cpp
    src/Methods/Meta.cpp
    src/Methods/Method.cpp
    src/Methods/StringMethod.cpp
    src/Methods/Swarm.cpp
    src/Methods/Umbrella.cpp
#    src/Constraints/Constraint.cpp
#    src/Constraints/COPSSImage.cpp
    src/CVs/CollectiveVariable.cpp
    src/Loggers/Logger.cpp
    src/Validator/RequirementLoader.cpp
    src/JSON/jsoncpp.cpp
    src/Utility/PairwiseKernel.cpp
    src/Utility/Basis.cpp
    src/Hook.cpp
    src/ResourceHandler.cpp
    "include/nnet/nnet.cpp"
)

####################### MD ENGINES #########################

# Initialize list of engine libraries to link to SSAGES
set(ENGINE_LIBRARIES "")

# Set variables for supported engines.
set(LAMMPS NO CACHE BOOL "Use LAMMPS simulation engine.")
set(GROMACS NO CACHE BOOL "Use GROMACS simulation engine.")
set(LAMMPS_SRC "" CACHE PATH "LAMMPS source directory.")
set(GROMACS_SRC "" CACHE PATH "GROMACS source directory.")
set(QBOX NO CACHE BOOL "Use QBox simulation engine.")
set(OPENMD_SRC "" CACHE PATH "OpenMD source directory.")
set(HOOMD NO CACHE BOOL "Use HOOMD-blue simulation engine.")
set(HOOMD_ROOT "" CACHE PATH "HOOMD-blue installation directory.")

# Determine which simulation engine to use.
if(LAMMPS_SRC OR LAMMPS)
    add_subdirectory(hooks/lammps/)
elseif(GROMACS OR GROMACS_SRC)
    add_subdirectory(hooks/gromacs/)
elseif(QBOX)
    add_subdirectory(hooks/qbox/)
elseif(OPENMD_SRC)
	add_subdirectory(hooks/openmd/)
elseif(HOOMD OR HOOMD_ROOT)
	add_subdirectory(hooks/hoomd/)
else()
    message( FATAL_ERROR "You cannot build SSAGES without specifying an engine." )
endif()

################### END MD ENGINES ######################

# Configure header.
configure_file (
    "${PROJECT_SOURCE_DIR}/src/config.h.in"
    "${PROJECT_BINARY_DIR}/config.h"
)

# SSAGES library.
add_library(libssages SHARED ${HOOKSRC} ${SOURCE})
set_target_properties(libssages PROPERTIES PREFIX "")
target_link_libraries(libssages ${MPI_CXX_LIBRARIES})
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

# Dependencies
add_dependencies(libssages eigen jsoncpp mxx schema)

# SSAGES executable.
add_executable(ssages
    src/main.cpp
    ${DRIVERSRC}
)

# Call MD engine configuration.
CONFIGURE_ENGINE()

target_link_libraries(ssages
    libssages
    ${ENGINE_LIBRARIES}
)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set_target_properties(ssages PROPERTIES COMPILE_FLAGS "-Wall -Wextra -g")
elseif (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set_target_properties(ssages PROPERTIES COMPILE_FLAGS "-O2 -g")
else ()
    set_target_properties(ssages PROPERTIES COMPILE_FLAGS "-O3")
endif ()

# Configure flags
# - Globals
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++11")

if (SSAGES_PERF)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -march=native")
endif()

# Add doc subdirectory for building the documentation
add_subdirectory(doc)

if (BUILD_TESTS)
    include(FetchGoogletest)
    enable_testing()
    add_subdirectory(test)
endif (BUILD_TESTS)
