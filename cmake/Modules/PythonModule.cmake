# check_python_module(module_name)
#
#   Test to see if a Python module is installed. Set ${module_name}_FOUND if a Python interpreter was found and the
#   module could successfully be imported.

if(NOT PYTHONINTERP_FOUND)
    find_package(PythonInterp)
endif(NOT PYTHONINTERP_FOUND)

function(check_python_module module_name)
    message(STATUS "Checking for Python module: ${module_name}")
    set(_module_found NO)

    if(PYTHONINTERP_FOUND)
        execute_process(
            COMMAND ${PYTHON_EXECUTABLE} -c "import ${module_name}"
            RESULT_VARIABLE _result
            ERROR_QUIET OUTPUT_QUIET
        )
        if(${_result} EQUAL 0)
            set(_module_found YES)
        endif(${_result} EQUAL 0)
    endif(PYTHONINTERP_FOUND)

    if(_module_found)
        message(STATUS "Python module ${module_name} was found.")
    else(_module_found)
        message(STATUS "Python module ${module_name} was NOT found.")
    endif(_module_found)

    set("PYTHON_MODULE_${module_name}_FOUND" ${_module_found} PARENT_SCOPE)
endfunction(check_python_module)

