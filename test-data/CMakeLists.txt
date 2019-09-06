# Add a custom target to simulate test data
add_custom_target(simulate-test-data ALL)

if(BUILD_PYTHON_BINDINGS)
    include(PythonModule)
    check_python_module(matplotlib)

    if(PYTHON_MODULE_matplotlib_FOUND)
        # Add a target to plot ground truths if plot partition is available.
        option(PLOT_PARTITIONS "Set to true to generate partition plots" ON)
        add_custom_target(plot-test-data ALL)
        add_dependencies(plot-test-data python)
    endif(PYTHON_MODULE_matplotlib_FOUND)
endif(BUILD_PYTHON_BINDINGS)

if(PLOT_PARTITIONS)
    message(STATUS "Generated test partitions will be plotted.")
else(PLOT_PARTITIONS)
    message(STATUS "Generated test partitions will NOT be plotted.")
endif(PLOT_PARTITIONS)

# A function which will simulate data based on a configuration file. If you want the output from this command to be
# consistent across builds consider setting the 'seed' variable in the config.
function(add_test_data_config source_prefix)
    set(_config "${CMAKE_CURRENT_SOURCE_DIR}/${source_prefix}.ini")
    set(_ground_truth "${CMAKE_CURRENT_BINARY_DIR}/${source_prefix}-ground-truth.json")
    set(_input_data "${CMAKE_CURRENT_BINARY_DIR}/${source_prefix}-input.json")

    set(_ground_truth_plot "${CMAKE_CURRENT_BINARY_DIR}/${source_prefix}-ground-truth.pdf")
    set(_input_data_plot "${CMAKE_CURRENT_BINARY_DIR}/${source_prefix}-input.pdf")

    add_custom_command(
        OUTPUT "${_ground_truth}"
        COMMAND "${PYTHON_EXECUTABLE}"
            ${BIGGLES_INTERNAL_PY} simulate
            --config "${_config}" > "${_ground_truth}"
            DEPENDS "${_config}" ${BIGGLES_INTERNAL_PY}
    )

    add_custom_command(
        OUTPUT "${_input_data}"
        COMMAND "${PYTHON_EXECUTABLE}"
            ${BIGGLES_INTERNAL_PY} demote < "${_ground_truth}" > "${_input_data}"
        DEPENDS "${_ground_truth}" ${BIGGLES_INTERNAL_PY}
    )

    add_custom_target(simulate-test-data-${source_prefix}
        DEPENDS "${_ground_truth}" "${_input_data}"
    )

    # Add the generated data to the simulate target
    add_dependencies(simulate-test-data simulate-test-data-${source_prefix})

    if(PLOT_PARTITIONS)
        add_custom_command(
            OUTPUT "${_ground_truth_plot}"
            COMMAND "${PYTHON_EXECUTABLE}"
                ${BIGGLES_INTERNAL_PY} plot_partition
                --title "${source_prefix}: ground truth"
                "${_ground_truth}" "${_ground_truth_plot}"
            DEPENDS "${_ground_truth}" "${PLOT_PARTITION_SCRIPT}"
        )

        add_custom_command(
            OUTPUT "${_input_data_plot}"
            COMMAND "${PYTHON_EXECUTABLE}"
                ${BIGGLES_INTERNAL_PY} plot_partition
                --title "${source_prefix}: input to tracker"
                "${_input_data}" "${_input_data_plot}"
            DEPENDS "${_input_data}" "${PLOT_PARTITION_SCRIPT}"
        )

        add_custom_target(plot-test-data-${source_prefix}
            DEPENDS "${_ground_truth_plot}" "${_input_data_plot}"
        )

        add_dependencies(plot-test-data plot-test-data-${source_prefix})
    endif(PLOT_PARTITIONS)
endfunction(add_test_data_config)

# A simple test case: should be unchallenging
add_test_data_config(simple)

# A more crowded case
add_test_data_config(crowded)

# A constrained motion
add_test_data_config(constrained)
