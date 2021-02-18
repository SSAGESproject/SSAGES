set(MXX_REPOSITORY "https://github.com/patflick/mxx.git")
set(MXX_TAG "master")

ExternalProject_Add(mxx
    # Using DOWNLOAD_COMMAND instead of GIT_REPOSITORY to pass some convenient git flags.
    # For instance, GIT_SHALLOW is available from CMake 3.6
    DOWNLOAD_COMMAND ${GIT_EXECUTABLE} clone --branch=${MXX_TAG} --depth=1
        -c advice.detachedHead=false ${MXX_REPOSITORY}
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${CMAKE_COMMAND} -E copy_directory
        <SOURCE_DIR>/include/mxx ${CMAKE_CURRENT_SOURCE_DIR}/include/mxx
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    LOG_UPDATE ON
)
