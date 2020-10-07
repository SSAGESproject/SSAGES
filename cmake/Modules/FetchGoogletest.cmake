set(GTEST_REPOSITORY "https://github.com/google/googletest.git")
set(GTEST_TAG "release-1.10.0")

# http://stackoverflow.com/questions/9689183/cmake-googletest
include(ExternalProject)
ExternalProject_Add(
    googletest
    DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}
    GIT_REPOSITORY ${GTEST_REPOSITORY}
    GIT_TAG ${GTEST_TAG}
    TIMEOUT 600
    CMAKE_ARGS
    -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DBUILD_GMOCK=OFF -DINSTALL_GTEST=OFF
    # Disable install step
    INSTALL_COMMAND ""
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)

# Specify include dir
ExternalProject_Get_Property(googletest SOURCE_DIR)
set(GTEST_INCLUDE_DIR ${SOURCE_DIR}/googletest/include)

ExternalProject_Get_Property(googletest BINARY_DIR)
message(STATUS ${BINARY_DIR})
set(GTEST_LIBRARY_PATH ${BINARY_DIR}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a)

message(STATUS "")
message(STATUS "*** Google Test Framework will be used for unit tests")
message(STATUS "*** GTEST_LIBRARY_PATH = ${GTEST_LIBRARY_PATH}")
message(STATUS "*** GTEST_INCLUDE_DIR  = ${GTEST_INCLUDE_DIR}")
message(STATUS "")
