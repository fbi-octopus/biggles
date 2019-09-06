# Make sure we can #include the utility library
include_directories(${CMAKE_SOURCE_DIR})

add_library(biggles SHARED
    detail/random.cpp
    detail/track_pair_data.cpp
    detail/transfer_possibility.cpp
    detail/cross_over_possibility.cpp
    detail/flip_possibility.cpp
    detail/possible_move_recorder.cpp
    detail/fun.cpp
    kalman_filter.cpp
    kalman_filter_cache.cpp
    io.cpp
    partition.cpp
    mh_moves/mh_moves.cpp
    mh_moves/move_birth.cpp
    mh_moves/move_death.cpp
    mh_moves/move_extend.cpp
    mh_moves/move_flip.cpp
    mh_moves/move_merge.cpp
    mh_moves/move_reduce.cpp
    mh_moves/move_split.cpp
    mh_moves/move_update.cpp
    mh_moves/move_cross_over.cpp
    mh_moves/move_transfer.cpp
    mh_moves/utility.cpp
    model.cpp
    samplers.cpp
    partition_sampler.cpp
    sampling/simple.cpp
    sampling/gibbs.cpp
    sampling/metropolis_hastings.cpp
    server/engine.cpp
    server/stepper.cpp
    simulate.cpp
    track.cpp
    track_collection.cpp
    tracker.cpp
    third-party/gen_beta.c
)

add_library(biggles_tools SHARED
    tools/sundries.cpp
    tools/partition_graph_edit_distance.cpp
    tools/max_creation.cpp
)

include_directories(${CMAKE_SOURCE_DIR}/include/biggles)

target_link_libraries(biggles
    utility
    ${Boost_RANDOM_LIBRARY}
    ${Boost_THREAD_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_FILESYSTEM_LIBRARY}
)

target_link_libraries(biggles_tools
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${Boost_FILESYSTEM_LIBRARY}
    biggles
)

# Set library name and version
set_target_properties(biggles PROPERTIES
        OUTPUT_NAME biggles
        SOVERSION   1
        VERSION     0.0.1
)

# Generate pkg-config metadata for library.
configure_file(biggles.in.pc biggles.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/biggles.pc DESTINATION lib/pkgconfig)

# Specify library installation destinations
install(TARGETS biggles
        RUNTIME DESTINATION         bin
        LIBRARY DESTINATION         lib
)

install(TARGETS biggles_tools
        RUNTIME DESTINATION         bin
        LIBRARY DESTINATION         lib
)

MESSAGE( STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}" )
