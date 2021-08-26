set(EIGEN_REPOSITORY "https://gitlab.com/libeigen/eigen.git")
set(EIGEN_TAG 3.4.0)

ExternalProject_Add(eigen
    # Using DOWNLOAD_COMMAND instead of GIT_REPOSITORY to pass some convenient git flags.
    # For instance, GIT_SHALLOW is available from CMake 3.6
    DOWNLOAD_COMMAND ${GIT_EXECUTABLE} clone --branch=${EIGEN_TAG} --depth=1
        -c advice.detachedHead=false ${EIGEN_REPOSITORY}
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${CMAKE_COMMAND} -E copy_directory
        <SOURCE_DIR>/Eigen ${CMAKE_CURRENT_SOURCE_DIR}/include/Eigen
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    LOG_UPDATE ON
)
