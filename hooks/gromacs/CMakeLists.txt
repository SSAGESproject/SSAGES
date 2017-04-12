# List of supported versions of Gromacs
set (SUPPORTED_GROMACS_VERSIONS "5.1.2"
                                "5.1.3")

set(DGMX_BUILD_OWN_FFTW OFF CACHE BOOL "Let Gromacs build its own libfftw.")

if (GROMACS_SRC)
    # Get absolute path in case user specified relative path
    if (IS_ABSOLUTE ${GROMACS_SRC})
        get_filename_component(GROMACS_SRC ${GROMACS_SRC} ABSOLUTE)
    else (IS_ABSOLUTE ${GROMACS_SRC})
        get_filename_component(GROMACS_SRC ${PROJECT_BINARY_DIR}/${GROMACS_SRC}
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
    
    # Set the appropriate patch.
    if(GROMACS_VERSION MATCHES "^5.1.[0-9]+$")
        set(GMX_PATCH "gmx_diff_5.1.x.patch")
    else()
        message(WARNING "\nWARNING: Could not find Gromacs version specific patch. "
                           "Attempting to use generic one.")
    endif()
    
    # Check that files to be patched do exist
    execute_process (COMMAND "sed" "-n"
                                "s/^.*gromacs-ssages\\/\\([^t]*\\)\t.*/\\1/p"
                                "${PROJECT_SOURCE_DIR}/hooks/gromacs/${GMX_PATCH}"
                        COMMAND "sed" "-e" ":a" "-e" "N" "-e" "\$!ba" "-e" "s/\\n/;/g"
                        OUTPUT_VARIABLE GROMACS_FILES_TO_BE_PATCHED
                        OUTPUT_STRIP_TRAILING_WHITESPACE)

    foreach(PATCH_FILE ${GROMACS_FILES_TO_BE_PATCHED})
        if (NOT EXISTS ${GROMACS_SRC}/${PATCH_FILE})
            message (FATAL_ERROR "\nERROR: File ${PATCH_FILE} does not exist.\n"
                                    "This error is probably caused by an unsupported "
                                    "version of Gromacs.")
        endif (NOT EXISTS ${GROMACS_SRC}/${PATCH_FILE})
    endforeach()
else ()
    set (GROMACS_SRC ${CMAKE_CURRENT_BINARY_DIR}/gromacs-download-prefix/src/gromacs-download)
    set (GROMACS_SRC ${GROMACS_SRC} PARENT_SCOPE )
endif ()

set(GRO_CXX_ARGS "-std=c++11 -I${PROJECT_SOURCE_DIR}/src -I${PROJECT_SOURCE_DIR}/include -I${CMAKE_CURRENT_SOURCE_DIR}")

# Configure patch.
configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/apply_patch.sh.in"
    "${CMAKE_CURRENT_BINARY_DIR}/apply_patch.sh"
)

if (GROMACS)
    set(GMX_PATCH "gmx_diff_5.1.x.patch")
    
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
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gromacs
    DOWNLOAD_COMMAND ""
    SOURCE_DIR ${GROMACS_SRC}
    PATCH_COMMAND sh ${CMAKE_CURRENT_BINARY_DIR}/apply_patch.sh ${GMX_PATCH}
    CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    -DGMX_MPI=ON
    -DCMAKE_CXX_FLAGS=${GRO_CXX_ARGS}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DGMX_EXTRA_LIBRARIES=${CMAKE_BINARY_DIR}/libssagesgmx.so
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DGMX_BUILD_OWN_FFTW=${GMX_BUILD_OWN_FFTW}
    -DCUDA_PROPAGATE_HOST_FLAGS=OFF
    -DCUDA_NVCC_FLAGS="-std=c++11 -I${PROJECT_SOURCE_DIR}/src -I${PROJECT_SOURCE_DIR}/include -I${CMAKE_CURRENT_SOURCE_DIR}"
    BUILD_COMMAND $(MAKE)
    BUILD_IN_SOURCE false
    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/gromacs
)

# Always build gromacs
ExternalProject_Add_Step(
    gromacs forcebuild
    COMMAND echo "Force build gromacs"
    DEPENDEES configure
    DEPENDERS build
    ALWAYS 1
)

set(HOOKSRC
    hooks/gromacs/GromacsHook.cpp
    PARENT_SCOPE
)

set(ENABLE_GROMACS true PARENT_SCOPE)
set(MD_ENGINE "Gromacs" PARENT_SCOPE)

# Function to do post library configuration.
function(CONFIGURE_ENGINE)
    # Build custom shared library for dependencies. 
	add_library(ssagesgmx SHARED ${HOOKSRC} ${SOURCE})
	target_link_libraries(ssagesgmx ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
    include_directories(${GROMACS_SRC}/src)

	add_dependencies(gromacs ssagesgmx)

	if (GROMACS)
		add_dependencies(ssagesgmx gromacs-download)
	endif (GROMACS)

	# To determine Gromacs install dir. 
	include(GNUInstallDirs)
	target_link_libraries(ssages 
		ssagesgmx 
		${CMAKE_BINARY_DIR}/hooks/gromacs/gromacs/${CMAKE_INSTALL_LIBDIR}/libgromacs_mpi.so
		${CMAKE_BINARY_DIR}/hooks/gromacs/gromacs/${CMAKE_INSTALL_LIBDIR}/libmdrun.a
		)
	add_dependencies(ssages gromacs)
endfunction()