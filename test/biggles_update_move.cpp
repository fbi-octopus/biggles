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


partition_ptr_t get_partition0() {
    observation_collection obs_coll;
    obs_coll.insert(new_obs(0,0,0));
    obs_coll.insert(new_obs(0,0,1));
    obs_coll.insert(new_obs(0,0,2));
    obs_coll.insert(new_obs(0,0,3));
    boost::shared_ptr<track_collection> tracks(new track_collection());
    tracks->insert(track(0, 4, obs_coll.begin(), obs_coll.end(), 1.f));
    clutter_ptr clutter(new clutter_t());
    return partition_ptr_t(new partition(tracks, clutter, partition::expansion_2d(0.f, 2.f, 0.f, 2.f)));
}

int test_pdr(const int total, const random_gen_paras& params ) {

    partition_ptr_t orignal_partition(random_debug_partition(params));

    capability_recorder_ptr cap_rec_ptr(new capability_recorder(
        orignal_partition->first_time_stamp(), orignal_partition->last_time_stamp(), orignal_partition->tracks()));
    orignal_partition->set_capability_recorder(cap_rec_ptr);

    test_t tests;
    BOOST_FOREACH (const std::string& msg,  move_pdr_test(orignal_partition, "Update", total, tests)) {
        diag("%s", msg.c_str());
    }
    BOOST_FOREACH (const test_t::value_type& item, tests) {
        ok(item.second, "%s", item.first.c_str());
    }
    return 0;
}

void test_pdr2(const int total, const random_gen_paras& params ) {

    test_t tests;
    BOOST_FOREACH (const std::string& msg,  move_pdr_test(get_partition0(), "Update", total, tests)) {
        diag("%s", msg.c_str());
    }
    BOOST_FOREACH (const test_t::value_type& item, tests) {
        ok(item.second, "%s", item.first.c_str());
    }
}

void test_3() {
    partition_ptr_t orig_part(get_partition0());
    biggles_move testmove("Update");
    partition_ptr_t forw_part;
    partition_ptr_t back_part;
    float pdr = 0.f;
    while (not testmove(orig_part, forw_part, pdr));
    diag("%s", partition_to_string(*orig_part).c_str());
    diag("%s", partition_to_string(*forw_part).c_str());
    diag("%s", clutter_to_str(forw_part->clutter(), forw_part->last_time_stamp()).c_str());
    int max_attempts = 100;
    int attempts = 0;
    bool success = false;
    std::set<std::string> part_set;
    while (not success and attempts++ < max_attempts) {
        bool done = testmove(forw_part, back_part, pdr);
        if (done) {
            diag("  %s, %s", partition_to_string(*back_part).c_str(),
                    clutter_to_str(back_part->clutter(), back_part->last_time_stamp()).c_str());
            part_set.insert(partition_to_string(*back_part));
        }
        success = *back_part == *orig_part;
    }
    //for (std::set<std::string>::iterator it = part_set.begin(); it != part_set.end(); ++it) { diag("  %s", it->c_str()); }
    if (success)
        diag("success after %d", attempts);
    else
        diag("failure");
}


int main(int argc, char** argv)
{
    general_paras gen_pars;
    gen_pars.total = 10000;
    gen_pars.seed = 0;
    gen_pars.rgp.lambda = 2.5;
    gen_pars.rgp.p_no = 0.2;
    gen_pars.rgp.p_yes = 0.5;
    gen_pars.rgp.p_tr = 0.5;
    gen_pars.rgp.min_tracks = 2;
    gen_pars.rgp.max_tracks = 3;
    if (not parse_args(argc, argv, gen_pars))
        return exit_status();
    biggles::detail::seed_prng(gen_pars.seed);

    plan_no_plan();


    test_pdr(gen_pars.total, gen_pars.rgp);
    //test_pdr2(gen_pars.total, gen_pars.rgp);

    diag("random seed = 0x%x", gen_pars.seed);
    return exit_status();
}
