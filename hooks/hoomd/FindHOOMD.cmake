# CMake script for finding HOOMD and setting up all needed compile options to create and link a plugin library
#
# Variables taken as input to this module:
# HOOMD_ROOT :          location to look for HOOMD, if it is not in the python path
#
# Variables defined by this module:
# FOUND_HOOMD :         set to true if HOOMD is found
# HOOMD_LIBRARIES :     a list of all libraries needed to link to to access hoomd (uncached)
# HOOMD_INCLUDE_DIR :   a list of all include directories that need to be set to include HOOMD
# HOOMD_LIB :           a cached var locating the hoomd library to link to
#
# various ENABLE_ flags translated from hoomd_config.h so this plugin build can match the ABI of the installed hoomd
#
# as a convenience (for the intended purpose of this find script), all include directories and definitions needed
# to compile with all the various libs (boost, python, winsoc, etc...) are set within this script

if(HOOMD_ROOT)
  message(STATUS "Using HOOMD installation at " ${HOOMD_ROOT})
else(HOOMD_ROOT)
  message(FATAL_ERROR "HOOMD_ROOT must be set to the HOOMD installation path.")
endif(HOOMD_ROOT)

# search for the hoomd include directory
find_path(HOOMD_INCLUDE_DIR
          NAMES HOOMDVersion.h
          HINTS ${HOOMD_ROOT}/include
          )

if(HOOMD_INCLUDE_DIR)
    message(STATUS "Found HOOMD include directory: ${HOOMD_INCLUDE_DIR}")
else(HOOMD_INCLUDE_DIR)
    message(FATAL_ERROR "Could not find HOOMD include files at HOOMD_ROOT/include.")
endif(HOOMD_INCLUDE_DIR)

set(FOUND_HOOMD TRUE)

#############################################################
## Now that we've found hoomd, lets do some setup

include_directories(${HOOMD_INCLUDE_DIR})

# run all of HOOMD's generic lib setup scripts
set(CMAKE_MODULE_PATH ${HOOMD_ROOT}
                      ${HOOMD_ROOT}/CMake/hoomd
                      ${HOOMD_ROOT}/CMake/thrust
                      ${CMAKE_MODULE_PATH}
                      )

# grab previously-set hoomd configuration
include (hoomd_cache)

# Handle user build options
include (CMake_build_options)
include (CMake_preprocessor_flags)
# setup the install directories
include (CMake_install_options)

# Find the python executable and libraries
include (HOOMDPythonSetup)
# setup numpy
include (HOOMDNumpySetup)
# Find CUDA and set it up
include (HOOMDCUDASetup)
# Set default CFlags
include (HOOMDCFlagsSetup)
# include some os specific options
include (HOOMDOSSpecificSetup)
# setup common libraries used by all targets in this project
include (HOOMDCommonLibsSetup)
# setup macros
include (HOOMDMacros)
# setup MPI support
include (HOOMDMPISetup)

set(HOOMD_LIB ${HOOMD_ROOT}/_hoomd${PYTHON_MODULE_EXTENSION})
set(HOOMD_MD_LIB ${HOOMD_ROOT}/md/_md${PYTHON_MODULE_EXTENSION})
set(HOOMD_DEM_LIB ${HOOMD_ROOT}/dem/_dem${PYTHON_MODULE_EXTENSION})
set(HOOMD_HPMC_LIB ${HOOMD_ROOT}/hpmc/_hpmc${PYTHON_MODULE_EXTENSION})
set(HOOMD_CGCMM_LIB ${HOOMD_ROOT}/cgcmm/_cgcmm${PYTHON_MODULE_EXTENSION})
set(HOOMD_METAL_LIB ${HOOMD_ROOT}/metal/_metal${PYTHON_MODULE_EXTENSION})
set(HOOMD_DEPRECATED_LIB ${HOOMD_ROOT}/deprecated/_deprecated${PYTHON_MODULE_EXTENSION})

set(HOOMD_LIBRARIES ${HOOMD_LIB} ${HOOMD_COMMON_LIBS})
