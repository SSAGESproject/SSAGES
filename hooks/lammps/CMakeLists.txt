# List of supported versions of LAMMPS
set (SUPPORTED_LAMMPS_VERSIONS "10 Aug 2015"
                               "7 Dec 2015"
                               "14 May 2016"
                               "15 Jul 2016"
                               "30 Jul 2016")
if (LAMMPS_SRC)
    # Get absolute path in case user specified relative path
    if (IS_ABSOLUTE ${LAMMPS_SRC})
        get_filename_component(LAMMPS_SRC ${LAMMPS_SRC} ABSOLUTE)
    else ()
        get_filename_component(LAMMPS_SRC ${PROJECT_BINARY_DIR}/${LAMMPS_SRC}
                                ABSOLUTE)
    endif ()

    # Check if given LAMMPS directory exists.
    if (NOT EXISTS ${LAMMPS_SRC})
        message (FATAL_ERROR "\nERROR: Directory ${LAMMPS_SRC} does not exist.")
    endif ()

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
else ()
    set (LAMMPS_SRC "${CMAKE_CURRENT_BINARY_DIR}/lammps-download-prefix/src/lammps-download/src")
endif ()

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

# List of optional LAMMPS packages (as of Oct 2016)
set(LAMMPS_PACKAGES "all"      # All packages
                    "standard" # All standard packages
                    "user"     # All user packages
                    # Standard packages:
                    "asphere"
                    "body"
                    "class2"
                    "colloid"
                    "compress"
                    "coreshell"
                    "dipole"
                    "gpu"
                    "granular"
                    "kim"
                    "kokkos"
                    "kspace"
                    "manybody"
                    "mc"
                    "meam"
                    "misc"
                    "molecule"
                    "mpiio"
                    "opt"
                    "peri"
                    "poems"
                    "python"
                    "qeq"
                    "reax"
                    "replica"
                    "rigid"
                    "shock"
                    "snap"
                    "srd"
                    "voronoi"
                    # User packages:
                    "user-atc"
                    "user-awpmd"
                    "user-cg-cmm"
                    "user-colvars"
                    "user-diffraction"
                    "user-dpd"
                    "user-drude"
                    "user-eff"
                    "user-fep"
                    "user-h5md"
                    "user-intel"
                    "user-lb"
                    "user-manifold"
                    "user-mgpt"
                    "user-misc"
                    "user-molfile"
                    "user-omp"
                    "user-phonon"
                    "user-qmmm"
                    "user-qtb"
                    "user-quip"
                    "user-reaxc"
                    "user-smd"
                    "user-smtbq"
                    "user-sph"
                    "user-tally"
                    "user-vtk")

foreach (PACKAGE ${LAMMPS_PACKAGES})
    add_custom_target (yes-${PACKAGE} $(MAKE) -C ${LAMMPS_SRC} yes-${PACKAGE})
    add_custom_target (no-${PACKAGE} $(MAKE) -C ${LAMMPS_SRC} no-${PACKAGE})
endforeach (PACKAGE ${LAMMPS_PACKAGES})

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
    CCFLAGS+=-O3
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
    COMMAND ln -sf ${CMAKE_CURRENT_SOURCE_DIR}/fix_ssages.h ${LAMMPS_SRC}/fix_ssages.h
    COMMAND touch ${LAMMPS_SRC}/modify.cpp
)

# Add source for appropriate package hook.
set(HOOKSRC
    hooks/lammps/fix_ssages.cpp
    PARENT_SCOPE
)

set(MD_ENGINE "LAMMPS" PARENT_SCOPE)
SET(LAMMPS_SRC ${LAMMPS_SRC} PARENT_SCOPE)

# Function to do post library configuration.
function(CONFIGURE_ENGINE)
    include_directories("${LAMMPS_SRC}")
    add_dependencies(lammps libssages)
	if (LAMMPS)
		add_dependencies(libssages lammps-download)
	endif (LAMMPS)
	target_link_libraries(ssages libssages ${LAMMPS_SRC}/liblammps_mpi.so)
	add_dependencies(ssages lammps)    
endfunction()
