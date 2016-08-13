# Add testing 
set(TEST_DEPS  
	${GTEST_LIBRARY_PATH} 
	${GTEST_LIBRARY_MAIN_PATH}
	libssages)

set  (UNIT_TESTS "AtomPositionCVTests"
                 "GridTests"
                 "NearestNeighborTests"
                 "TorsionalCVTests"
                 "UCVectorTests"
                 "UnitCellConversionTests"
                 "UnwrapCoordinatesTests"
                 "SnapshotTests")

foreach (TEST ${UNIT_TESTS})
    add_executable (${TEST} ${TEST}.cpp)
    target_link_libraries (${TEST} ${TEST_DEPS})
    target_include_directories (${TEST} PRIVATE "${GTEST_INCLUDE_DIR}")
    add_dependencies (${TEST} googletest)
    add_test(${TEST} ${TEST})
endforeach (TEST ${UNIT_TESTS})
