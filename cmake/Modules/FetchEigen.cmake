set(EIGEN_REPOSITORY "https://gitlab.com/libeigen/eigen.git")
set(EIGEN_PATCH "${CMAKE_CURRENT_SOURCE_DIR}/include/patches/eigen-using-std-real.patch")
set(EIGEN_TAG 3.3.7)

if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include/eigen")
    execute_process(
        COMMAND ${GIT_EXECUTABLE} clone --branch=${EIGEN_TAG} --depth=1
            -c advice.detachedHead=false ${EIGEN_REPOSITORY}
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include"
        RESULT_VARIABLE EIGEN_GIT_CLONE_RESULT
    )
    if(NOT EIGEN_GIT_CLONE_RESULT EQUAL "0")
        message(FATAL_ERROR
            "`git clone ${EIGEN_REPOSITORY}` failed with ${EIGEN_GIT_CLONE_RESULT}"
        )
    endif()
endif()

execute_process(
    COMMAND ${GIT_EXECUTABLE} apply --check ${EIGEN_PATCH}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/eigen"
    RESULT_VARIABLE EIGEN_GIT_PATCH_RESULT
)
if(EIGEN_GIT_PATCH_RESULT EQUAL "0")
    execute_process(
        COMMAND ${GIT_EXECUTABLE} apply ${EIGEN_PATCH}
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/eigen"
    )
else()
    message(WARNING
        "`git apply --check ${EIGEN_PATCH}` failed with ${EIGEN_GIT_PATCH_RESULT}"
    )
endif()

if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include/Eigen")
    file(
        COPY "${CMAKE_CURRENT_SOURCE_DIR}/include/eigen/Eigen"
        DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/include"
    )
endif()

file(REMOVE_RECURSE "${CMAKE_CURRENT_SOURCE_DIR}/include/eigen")
