# The CCAN TAP library
add_library(test_ccan ccan/tap/tap.c)
#add_library(test_utility utility.cpp)
#target_link_libraries(test_utility ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Do we have valgrind installed?
find_program(VALGRIND_EXECUTABLE valgrind)
if(VALGRIND_EXECUTABLE)
    message(STATUS "Valgrind was found at ${VALGRIND_EXECUTABLE}. It will be used to run tests.")
else(VALGRIND_EXECUTABLE)
    message(STATUS "Valgrind was NOT found. Tests will run without any memory checking.")
endif(VALGRIND_EXECUTABLE)

# Need an option to switch matplotlib backend since there may be cases where X is not
# available to run with a Window backend
# TODO: Use it more extensively
option(BIGGLES_MATPLOTLIB_BACKEND "The backend to use for matplotlib" "QT4AGG")

# Function to make a Biggles test
include(CMakeParseArguments)
function(add_biggles_test)
    set(oneValueArgs NAME)
    set(multiValueArgs SOURCES)
    cmake_parse_arguments(BIGGLES_TEST "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    SET(test_target "test_${BIGGLES_TEST_NAME}")
    add_executable(${test_target} ${BIGGLES_TEST_SOURCES})
    target_link_libraries(${test_target} biggles biggles_tools ${Boost_RANDOM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY} test_ccan)

    add_test(NAME ${BIGGLES_TEST_NAME} COMMAND ${test_target})
    # Use valgrind if it is present
    if(VALGRIND_EXECUTABLE AND USE_VALGRIND)
        get_target_property(test_target_output ${test_target} LOCATION)
        add_test(NAME "${BIGGLES_TEST_NAME}Valgrind" COMMAND ${VALGRIND_EXECUTABLE} --error-exitcode=1 ${test_target_output})
    endif(VALGRIND_EXECUTABLE AND USE_VALGRIND)
endfunction(add_biggles_test)

add_biggles_test(NAME GibbsSampler SOURCES gibbs.cpp)
add_biggles_test(NAME BigglesCore SOURCES biggles_core.cpp)
add_biggles_test(NAME BigglesEngine SOURCES biggles_engine.cpp)
add_biggles_test(NAME BigglesSimulate SOURCES biggles_simulate.cpp)
add_biggles_test(NAME BigglesObservationCollection SOURCES biggles_observation_collection.cpp)
add_biggles_test(NAME BigglesClutter SOURCES biggles_clutter.cpp)
add_biggles_test(NAME BigglesSamplers SOURCES biggles_samplers.cpp)
add_biggles_test(NAME BigglesSplitMove SOURCES biggles_split_move.cpp)
add_biggles_test(NAME BigglesMergeMove SOURCES biggles_merge_move.cpp)
add_biggles_test(NAME BigglesUpdateMove SOURCES biggles_update_move.cpp)
add_biggles_test(NAME BigglesBirthMove SOURCES biggles_birth_move.cpp)
add_biggles_test(NAME BigglesDeathMove SOURCES biggles_death_move.cpp)
add_biggles_test(NAME BigglesCrossOverMove SOURCES biggles_cross_over_move.cpp)
add_biggles_test(NAME BigglesExtendMove SOURCES biggles_extend_move.cpp)
add_biggles_test(NAME BigglesReduceMove SOURCES biggles_reduce_move.cpp)
add_biggles_test(NAME BigglesTransferMove SOURCES biggles_transfer_move.cpp)
add_biggles_test(NAME BigglesSimpleSamplers SOURCES biggles_simple_samplers.cpp)
add_biggles_test(NAME BigglesProposal SOURCES biggles_proposal.cpp)
add_biggles_test(NAME BigglesPartitionSampler SOURCES biggles_partition_sampler.cpp)
add_biggles_test(NAME BigglesTracker SOURCES biggles_tracker.cpp)
add_biggles_test(NAME BigglesPriors SOURCES biggles_priors.cpp)
add_biggles_test(NAME BigglesLikelihood SOURCES biggles_likelihood.cpp)
add_biggles_test(NAME BigglesKalmanFilter SOURCES biggles_kalman_filter.cpp)
add_biggles_test(NAME BigglesObservationRate SOURCES biggles_observation_rate.cpp)
add_biggles_test(NAME BigglesLinkProbabilities SOURCES biggles_link_probabilities.cpp)
add_biggles_test(NAME BigglesPartitionDistance SOURCES biggles_partition_ged.cpp)
add_biggles_test(NAME BigglesMaxCreation SOURCES biggles_max_creation.cpp)
add_biggles_test(NAME BigglesIO SOURCES biggles_io.cpp)
add_biggles_test(NAME BigglesFun SOURCES biggles_fun.cpp)
add_biggles_test(NAME BigglesCapabilities SOURCES biggles_capability_recorder.cpp)
add_biggles_test(NAME BigglesIndexedPartition SOURCES biggles_indexed_partition.cpp)

# Add tests for simulate command line tool
add_test(
    NAME "SimulateTool"
    COMMAND "${PYTHON_EXECUTABLE}" "${BIGGLES_INTERNAL_PY}" simulate
        sim-default-output.json
)

add_test(
    NAME "SimulateToolWithConfig"
    COMMAND "${PYTHON_EXECUTABLE}" "${BIGGLES_INTERNAL_PY}" simulate
        --config "${CMAKE_CURRENT_SOURCE_DIR}/sim-config.ini"
        sim-config-output.json
)

# Add tests for plotting command line tool
#add_test(
#    NAME "PlotPartitionTool"
#    COMMAND "${PYTHON_EXECUTABLE}" "${BIGGLES_INTERNAL_PY}" plot_partition
#        "${CMAKE_CURRENT_SOURCE_DIR}/simple-ground-truth.json"
#        simple-ground-truth.pdf
#)
#set_property(
#    TEST "PlotPartitionTool"
#    PROPERTY ENVIRONMENT "BIGGLES_MATPLOTLIB_BACKEND=${BIGGLES_MATPLOTLIB_BACKEND}"
#)

# Add tests for demote tool
add_test(
    NAME "DemoteTool"
    COMMAND "${PYTHON_EXECUTABLE}" "${BIGGLES_INTERNAL_PY}" demote
        "${CMAKE_CURRENT_SOURCE_DIR}/simple-ground-truth.json"
        simple-demoted.json
)

# Add tests for track tool
add_test(
    NAME "TrackTool"
    COMMAND "${PYTHON_EXECUTABLE}" "${TICKLES_INTERNAL_PY}" track
        "${CMAKE_CURRENT_SOURCE_DIR}/simple-input.json"
        --samples 2000
)
