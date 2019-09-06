find_package(PythonInterp REQUIRED)

function(add_python_test _name _module)
    add_test(NAME BigglesPython${_name}
        COMMAND ${PYTHON_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/${_module}"
            "${CMAKE_CURRENT_BINARY_DIR}/../../test-data/"
    )
    set_tests_properties(BigglesPython${_name} PROPERTIES
        ENVIRONMENT "PYTHONPATH=${BIGGLES_INTERNAL_PYTHONPATH}"
    )
endfunction(add_python_test)

add_python_test(Import test_import.py)
add_python_test(Engine test_engine.py)
add_python_test(ConvertMSMM test_convert.py)

# Copy the features.dat test file to the test data directory
file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/../../test-data/msmm-features.dat
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../../test-data
)
