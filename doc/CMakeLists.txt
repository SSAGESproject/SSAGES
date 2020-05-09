set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})

find_package (Doxygen)

if (DOXYGEN_FOUND)
    if (DOXYGEN_DOT_FOUND)
        set (HAVE_DOT YES)
    else (DOXYGEN_DOT_FOUND)
        message (WARNING "\nWARNING: Please install Graphviz's dot tool for "
                         "full API documentation.")
        set (HAVE_DOT NO)
    endif (DOXYGEN_DOT_FOUND)

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen.config.in
                   ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config)

    add_custom_target(apiref
                      COMMAND doxygen ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config
                      COMMENT "Building API Reference Documentation")
else (DOXYGEN_FOUND)
    message (WARNING "\nWARNING: To generate API documentation, please install "
                     "Doxygen.")
endif (DOXYGEN_FOUND)

find_package (Sphinx)

if (SPHINX_FOUND)
    configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/source/conf.py.in"
                    "${CMAKE_CURRENT_BINARY_DIR}/conf.py")

    add_custom_target(
        manual
        COMMAND sphinx-build -b html -c ${CMAKE_CURRENT_BINARY_DIR}
                ${CMAKE_CURRENT_SOURCE_DIR}/source
                ${CMAKE_CURRENT_BINARY_DIR}/Manual
        COMMENT "Building User Manual"
    )

    add_custom_target(
        pdfmanual
        COMMAND sphinx-build -b latex -c ${CMAKE_CURRENT_BINARY_DIR}
                ${CMAKE_CURRENT_SOURCE_DIR}/source
                ${CMAKE_CURRENT_BINARY_DIR}/Manual
        COMMAND make
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Manual
        COMMENT "Building PDF version of Manual"
    )
else (SPHINX_FOUND)
    message (WARNING "\nWARNING: To generate User Manual, please install Sphinx.")
endif (SPHINX_FOUND)

if (DOXYGEN_FOUND AND SPHINX_FOUND)
    add_custom_target(doc
                      DEPENDS apiref manual)
endif (DOXYGEN_FOUND AND SPHINX_FOUND)
