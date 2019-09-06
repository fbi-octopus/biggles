# Doxygen support (optional)
set(DOC_TARGET_NAME "doc-cpp")
include(UseDoxygen)
if(DOXYGEN_FOUND)
    add_dependencies(doc doc-cpp)

    # Add an install target for the documentation. Make it optional so that a lac
    # of documentation does not, in itself, mean that the install will fail.
    install(
        DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/
        DESTINATION share/doc/biggles/cpp/
        OPTIONAL
    )
endif(DOXYGEN_FOUND)
