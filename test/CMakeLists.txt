# Add testing 
set(TEST_DEPS  
	${GTEST_LIBRARY_PATH} 
	${GTEST_LIBRARY_MAIN_PATH}
	libssages)

set  (UNIT_TESTS "ParticlePositionCVTests"
                 "ABFTests"
                 "TorsionalCVTests"
                 "SnapshotTests"
                 "BasisFuncTest"
                 "ParticleSeparationCVTests"
                 "RMSD_Test"
                 "AngleCVTests"
                 "ReadFileTest"
                 "MockCVTests"
				 "COMTests"
                 "DirectForwardFluxTests"
				 "SwarmTest"
				 "FTSTest"
                 "StringMethodTest"
                 "GyrationTensorCVTests"
                 "RouseModeCVTests"
                 "CoordinationNumberCVTests"
				 "NearestNeighborsCVTests")

set (INTEGRATION_TESTS "Meta_Single_Atom_Test"
                       "Basis_ADP_Test")

foreach (TEST ${UNIT_TESTS})
    add_executable (${TEST} unit_tests/${TEST}.cpp)
    target_link_libraries (${TEST} ${TEST_DEPS})
    target_include_directories (${TEST} PRIVATE "${GTEST_INCLUDE_DIR}")
    add_dependencies (${TEST} googletest)
    add_test(UNIT_TEST_${TEST} ${TEST})
endforeach (TEST ${UNIT_TESTS})

find_program (PYTHON_BIN "python")

if (PYTHON_BIN)
    # Create working directory for output files to be written to
    file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/working_directory/)

    foreach (TEST ${INTEGRATION_TESTS})
        set (SCRIPT "integration_tests/${TEST}.py")
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${SCRIPT}.in
                       ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT})
        add_test (NAME INTEGRATION_TEST_${TEST}
            COMMAND ${PYTHON_BIN} ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/working_directory/)

    endforeach (TEST ${INTEGRATION_TESTS})
endif (PYTHON_BIN)
