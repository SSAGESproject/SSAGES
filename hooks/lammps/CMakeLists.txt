# List of supported versions of LAMMPS
set (SUPPORTED_LAMMPS_VERSIONS
     "10 Aug 2015"  "7 Dec 2015"
     "16 Feb 2016" "14 May 2016" "30 Jul 2016"
     "5 Nov 2016" "17 Nov 2016"
     "31 Mar 2017" "11 Aug 2017"
     "16 Mar 2018" "22 Aug 2018" "12 Dec 2018"
     "5 Jun 2019" "7 Aug 2019"
     "3 Mar 2020"
    )
# List of git tags for above versions
# Starting with version 15 Sep 2016, the tags are named by date
set (SUPPORTED_LAMMPS_TAGS
     "r13864" "r14304"
     "r14624" "r15061" "r15407"
     "stable_5Nov2016"  "stable_17Nov2016"
     "stable_31Mar2017" "stable_11Aug2017"
     "stable_16Mar2018" "stable_22Aug2018" "stable_12Dec2018"
     "stable_5Jun2019" "stable_7Aug2019"
     "stable_3Mar2020"
    )
if (LAMMPS_SRC)
    # Get absolute path in case user specified relative path
    if (IS_ABSOLUTE ${LAMMPS_SRC})
        get_filename_component(LAMMPS_SRC ${LAMMPS_SRC} ABSOLUTE)
    else ()
        get_filename_component(LAMMPS_SRC ${PROJECT_BINARY_DIR}/${LAMMPS_SRC}
                                ABSOLUTE)
    endif ()

    # Check if directory LAMMPS_SRC exists
    if (NOT EXISTS ${LAMMPS_SRC})
        message (FATAL_ERROR "\nERROR: Directory ${LAMMPS_SRC} does not exist.")
    endif ()

    # Check if directory contains LAMMPS source files
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
        message (WARNING "\nWARNING: LAMMPS version \"${LAMMPS_VERSION}\" is not "
                           "officially supported by SSAGES.")
    endif ()
elseif (LAMMPS)
    set (LAMMPS_VERSION ${LAMMPS})

    # Check if LAMMPS version is supported, otherwise download default
    if (NOT ";${SUPPORTED_LAMMPS_VERSIONS};" MATCHES ";${LAMMPS_VERSION};")
        set (LAMMPS_VERSION "30 Jul 2016")
        set (LAMMPS_TAG r15407)
        message (STATUS "\n-- Using default version of LAMMPS (${LAMMPS_VERSION}).\n-- ")
    else ()
        string (REGEX REPLACE " " "" LAMMPS_TAG "${LAMMPS_VERSION}")
        set (LAMMPS_TAG "stable_${LAMMPS_TAG}")
        if (NOT ";${SUPPORTED_LAMMPS_TAGS};" MATCHES ";${LAMMPS_TAG};")
            list(FIND SUPPORTED_LAMMPS_VERSIONS "${LAMMPS_VERSION}" IDX)
            list(GET SUPPORTED_LAMMPS_TAGS ${IDX} LAMMPS_TAG)
        endif ()
    endif ()

    # Check if LAMMPS tag is supported
    if (NOT ";${SUPPORTED_LAMMPS_TAGS};" MATCHES ";${LAMMPS_TAG};")
        message (FATAL_ERROR "\nWARNING: LAMMPS tag \"${LAMMPS_TAG}\" is malformed.")
    endif ()

    ExternalProject_Add(
        lammps-download
        DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}
        GIT_REPOSITORY git://github.com/lammps/lammps.git
        GIT_TAG ${LAMMPS_TAG}
        TIMEOUT 600
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
    )
    set (LAMMPS_SRC "${CMAKE_CURRENT_BINARY_DIR}/lammps-download-prefix/src/lammps-download/src")
endif ()

set (LAMMPS_SRC ${LAMMPS_SRC} PARENT_SCOPE)

# List of optional LAMMPS packages (as of 1 Jun 2018)
set (LAMMPS_PACKAGES "all"      # All packages
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
                     "latte"
                     "manybody"
                     "mc"
                     "meam"
                     "misc"
                     "molecule"
                     "mpiio"
                     "mscg"
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
                     "user-bocs"
                     "user-cgdna"
                     "user-cgsdk"
                     "user-cg-cmm" # renamed "user-cgsdk" on 4 May 2017
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
                     "user-meamc"
                     "user-meso"
                     "user-mgpt"
                     "user-misc"
                     "user-mofff"
                     "user-molfile"
                     "user-nc-dump" # renamed "user-netcdf" on 4 May 2017
                     "user-netcdf"
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
                     "user-uef"
                     "user-vtk"
                     )

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
    CCFLAGS+=-std=c++11
    CCFLAGS+=-O3
    LIB+=$<TARGET_FILE:libssages>
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE true
)

ExternalProject_Add_Step(
    lammps forceconfigure
    COMMENT "Configuring LAMMPS"
    DEPENDERS build
    DEPENDEES update
    COMMAND ln -sf ${CMAKE_CURRENT_SOURCE_DIR}/fix_ssages.h ${LAMMPS_SRC}/fix_ssages.h &&
            touch ${LAMMPS_SRC}/modify.cpp &&
            sh ${LAMMPS_SRC}/Make.sh style
    ALWAYS 1
)

# Add source for appropriate package hook.
set(HOOKSRC
    hooks/lammps/fix_ssages.cpp
    PARENT_SCOPE
)

set(DRIVERSRC
    hooks/lammps/Driver.cpp
    PARENT_SCOPE
)

# Libraries from the engine to link
set(ENGINE_LIBRARIES
    "${LAMMPS_SRC}/liblammps_mpi.so"
    PARENT_SCOPE
)

# Function to do post library configuration.
function(CONFIGURE_ENGINE)
    include_directories("${LAMMPS_SRC}")
    include_directories(hooks/lammps/)
    if (LAMMPS)
        add_dependencies(libssages lammps-download)
    endif (LAMMPS)
    add_dependencies(lammps libssages)
    add_dependencies(ssages lammps)
endfunction()
