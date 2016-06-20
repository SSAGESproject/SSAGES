configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen.config.in
               ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config)

add_custom_target(apiref
                  COMMAND doxygen ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config
                  COMMENT "Building API Reference Documentation")

add_custom_target(manual
                  COMMAND sphinx-build -b html  ${CMAKE_CURRENT_SOURCE_DIR}/source
                                                ${CMAKE_CURRENT_BINARY_DIR}/Manual
                  COMMAND sphinx-build -b latex ${CMAKE_CURRENT_SOURCE_DIR}/source
                                                ${CMAKE_CURRENT_BINARY_DIR}/Manual
                  COMMENT "Building User Manual")

add_custom_target(doc
                  DEPENDS apiref manual)
