# Add testing 
set(TEST_DEPS  
    ${GTEST_LIBRARY_PATH}
    ${GTEST_LIBRARY_MAIN_PATH}
    libssages
)

set (UNIT_TESTS
    # General Tests
    "COMTests"
    "GridTests"
    "ReadFileTests"
    "SnapshotTests"
    # CV Tests
    "MockCVTests"
    "AngleCVTests"
    "GyrationTensorCVTests"
    "PairwiseCVTests"
    "ParticlePositionCVTests"
    "ParticleSeparationCVTests"
    "RMSDCVTests"
    "RouseModeCVTests"
    "TorsionalCVTests"
    # Method Tests
    "ABFTests"
    "BasisFuncTests"
    "DirectForwardFluxTests"
    "FTSTests"
    "StringMethodTests"
    "SwarmTests"
)

set (INTEGRATION_TESTS
    "Meta_Single_Atom_Test"
    "Basis_ADP_Test")

foreach (TEST ${UNIT_TESTS})
    add_executable (${TEST} unit_tests/${TEST}.cpp)
    target_link_libraries (${TEST} ${TEST_DEPS})
    target_include_directories (${TEST} PRIVATE "${GTEST_INCLUDE_DIR}")
    add_dependencies (${TEST} googletest)
    add_test(UNIT_TEST_${TEST} ${TEST})
endforeach (TEST ${UNIT_TESTS})

if (CMAKE_VERSION VERSION_LESS 3.12)
    find_package (PythonInterp REQUIRED)
    set (PYTHON_INTERP_FOUND ${PYTHONINTERP_FOUND})
    set (PYTHON_INTERP_EXE ${PYTHON_EXECUTABLE})
else ()
    find_package (Python REQUIRED COMPONENTS Interpreter)
    set (PYTHON_INTERP_FOUND ${Python_Interpreter_FOUND})
    set (PYTHON_INTERP_EXE ${Python_EXECUTABLE})
endif ()

if (PYTHON_INTERP_FOUND)
    # Create working directory for output files to be written to
    file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/working_directory/)

    foreach (TEST ${INTEGRATION_TESTS})
        set (SCRIPT "integration_tests/${TEST}.py")
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${SCRIPT}.in
                       ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT})
        add_test (NAME INTEGRATION_TEST_${TEST}
            COMMAND ${PYTHON_INTERP_EXE} ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/working_directory/)

    endforeach (TEST ${INTEGRATION_TESTS})
else ()
    message (WARNING "\nWARNING: Unable to generate integration test files. "
                     "Python binary not found.")
endif ()
