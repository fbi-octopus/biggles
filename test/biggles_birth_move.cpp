#include "biggles/observation.hpp"
#include "biggles/detail/random.hpp"
#include "biggles/mh_moves/mh_moves.hpp"
#include "biggles/mh_moves/utility.hpp"
#include "biggles/simulate.hpp"
#include <boost/foreach.hpp>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <string>
#include "biggles/tools/sundries.hpp"

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

void print_obs(observation obs) {
    std::cout << "x=" << x(obs) << ", y= " << y(obs) << ", t= " << t(obs) << std::endl;
}

bool minimum_partition_test() {
    observation_collection cl_obs;
    cl_obs.insert(new_obs(1, 1, 0));
    cl_obs.insert(new_obs(1, 1, 1));
    boost::shared_ptr<track_collection> tracks(new track_collection());
    clutter_ptr clutter(new clutter_t(cl_obs.begin(), cl_obs.end()));
    partition_ptr_t original_part_ptr(
        new partition (tracks, clutter, expansion_from_observations(*tracks, *clutter)));

    bool done = false;
    for (int i = 0; i < 100; ++i) {
        float prob=0.f;
        partition_ptr_t target_part_ptr;
        done = mh_moves::birth(original_part_ptr, target_part_ptr, prob);
        if (done) break;
    }
    return done;
}

int test_pdr(const int total, const random_gen_paras& paras) {
    partition_ptr_t original_part_ptr(random_debug_partition(paras));

    capability_recorder_ptr cap_rec_ptr(new capability_recorder(
        original_part_ptr->first_time_stamp(), original_part_ptr->last_time_stamp(), original_part_ptr->tracks()));
    original_part_ptr->set_capability_recorder(cap_rec_ptr);

    test_t tests;
    BOOST_FOREACH (const std::string& msg,  move_pdr_test(original_part_ptr, "Birth", total, tests)) {
        diag("%s", msg.c_str());
    }
    BOOST_FOREACH (const test_t::value_type& item, tests) {
        ok(item.second, "%s", item.first.c_str());
    }
    return 0;
}

int test1(const int total, const random_gen_paras& paras) {
    observation_collection cl_obs;

    for (time_stamp i = 0; i < 3; ++i)
        cl_obs.insert(new_obs(1, 1, i));
    boost::shared_ptr<track_collection> tracks(new track_collection());
    clutter_ptr clutter(new clutter_t(cl_obs.begin(), cl_obs.end()));

    partition_ptr_t original_part_ptr(new partition(tracks, clutter, expansion_from_observations(*tracks, *clutter)));

    test_t tests;
    BOOST_FOREACH (const std::string& msg,  move_pdr_test(original_part_ptr, "Birth", total, tests)) {
        diag("%s", msg.c_str());
    }
    BOOST_FOREACH (const test_t::value_type& item, tests) {
        ok(item.second, "%s", item.first.c_str());
    }
    return 0;
}

int main(int argc, char** argv)
{
    general_paras gen_pars;
    gen_pars.total = 30000;
    gen_pars.seed = 0;
    gen_pars.rgp.lambda = 2.5;
    gen_pars.rgp.p_no = 0.15;
    gen_pars.rgp.p_yes = 0.15;
    gen_pars.rgp.p_tr = 0.5;
    gen_pars.rgp.min_tracks = 1;
    gen_pars.rgp.max_tracks = 2;
    if (not parse_args(argc, argv, gen_pars))
        return exit_status();
    biggles::detail::seed_prng(gen_pars.seed);

    plan_no_plan();

    test_pdr(gen_pars.total, gen_pars.rgp);
    //test1(gen_pars.total, gen_pars.rgp);

    diag("random seed = 0x%x", gen_pars.seed);
    return exit_status();
}
