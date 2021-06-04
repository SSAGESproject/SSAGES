set(JSONCPP_REPOSITORY "https://github.com/open-source-parsers/jsoncpp.git")
set(JSONCPP_TAG 1.9.4)

ExternalProject_Add(jsoncpp
    # Using DOWNLOAD_COMMAND instead of GIT_REPOSITORY to pass some convenient git flags.
    # For instance, GIT_SHALLOW is available from CMake 3.6
    DOWNLOAD_COMMAND ${GIT_EXECUTABLE} clone --branch=${JSONCPP_TAG} --depth=1
        -c advice.detachedHead=false ${JSONCPP_REPOSITORY}
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE ON
    BUILD_COMMAND ${PYTHON_INTERP_EXE} amalgamate.py
                  -s ${CMAKE_CURRENT_SOURCE_DIR}/src/JSON/jsoncpp.cpp
                  -i ${CMAKE_CURRENT_SOURCE_DIR}/include/json/json.h
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    LOG_UPDATE ON
)

set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/src/JSON/jsoncpp.cpp PROPERTIES GENERATED TRUE)
