#include "biggles/kalman_filter.hpp"
#include "biggles/observation.hpp"
#include "biggles/partition.hpp"
#include "biggles/samplers.hpp"
#include "biggles/simulate.hpp"
#include "biggles/track.hpp"
#include "biggles/server/engine.hpp"
#include "biggles/server/stepper.hpp"
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <cmath>
#include <ctime>
#include <inttypes.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <utility>
#include <vector>

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

static boost::mt19937 random_generator(std::time(0));

using namespace biggles;


void test_engine(const model::parameters& params)
{
    partition_ptr_t partition_ptr;

    // simulate full data
    simulate::generate_partition(partition_ptr, params, 0, 128 + 128, 128.f, 128.f);

    diag("Generated %zu new tracks.", partition_ptr->tracks().size());
    diag("Generated %zu clutter observations.", partition_ptr->clutter().size());

    partition_ptr_t dataset_part_ptr;
    simulate::demote_all_tracks_from_partition(dataset_part_ptr, partition_ptr);

    server::engine worker;
    worker.set_initial_conditions(params, dataset_part_ptr);
    worker.start();
    diag("Started");
    boost::timer timer;
    while (timer.elapsed()<1) {}
    diag("Time over");
    worker.stop();
    diag("Stopped");

    ok1(true);
}

void test_stepper(const model::parameters& params) {
    /*
    partition_ptr_t partition_ptr;

    // simulate full data
    simulate::generate_partition(partition_ptr, params, 0, 128 + 128, 128.f, 128.f);

    diag("Generated %zu new tracks.", partition_ptr->tracks().size());
    diag("Generated %zu clutter observations.", partition_ptr->clutter().size());
    */

    boost::shared_ptr<track_collection> tracks(new track_collection());
    clutter_ptr clutter(new clutter_t());
    partition_ptr_t partition_ptr(new partition(tracks, clutter));

    /*
    partition_ptr_t dataset_part_ptr;
    simulate::demote_all_tracks_from_partition(dataset_part_ptr, partition_ptr);

    server::stepper dancer(params, dataset_part_ptr);
    */
    diag("OK");
    server::stepper dancer(params, partition_ptr);
    /*
    boost::shared_ptr<server::tracking_state> state = dancer.tracking_state_ptr();
    diag("current_partition %zu", state->current_partition.use_count());
    diag("best_partition %zu", state->best_partition.use_count());
    diag("last_proposed_partition %zu", state->last_proposed_partition.use_count());
    */
    /*
    dancer.step(1);
    state = dancer.tracking_state_ptr();
    diag("current_partition %zu", state->current_partition.use_count());
    diag("best_partition %zu", state->best_partition.use_count());
    diag("last_proposed_partition %zu", state->last_proposed_partition.use_count());
    */
}

int main(int argc, char** argv)
{
    plan_no_plan();

    float sigma_R = 0.1f;

    // initialise model parameters
    model::parameters params;
    model::mean_new_tracks_per_frame(params)           = 0.1f;
    model::mean_false_observations_per_frame(params)   = 0.05f;
    model::frame_to_frame_survival_probability(params) = 0.95f;
    model::generate_observation_probability(params)    = 0.8f;
    model::observation_error_covariance(params)       << sigma_R * sigma_R, 0, 0, sigma_R * sigma_R;

    //test_engine(params);
    test_stepper(params);
    return exit_status();
}
