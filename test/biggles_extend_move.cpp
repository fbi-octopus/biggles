#include "biggles/observation.hpp"
#include "biggles/detail/random.hpp"
#include "biggles/mh_moves/mh_moves.hpp"
#include "biggles/mh_moves/utility.hpp"
#include "biggles/simulate.hpp"
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <map>
#include <set>
#include <string>
#include "biggles/tools/sundries.hpp"

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

bool compare(const biggles::partition& part1, const biggles::partition& part2) {
    if (not are_equal(part1.clutter(), part2.clutter())) {
        diag("clutter different");
        diag("'%s' vs. '%s'", clutter_to_str(part1.clutter()).c_str(), clutter_to_str(part2.clutter()).c_str());
        return false;
    }
    biggles::track_collection tc1 = part1.tracks();
    biggles::track_collection tc2 = part2.tracks();
    if (tc1.size() not_eq tc2.size()) {
        diag("num of tracks different");
        return false;
    }
    biggles::track_collection::const_iterator it1, it2;
    for (it1 = tc1.begin(); it1 not_eq tc1.end(); ++it1) {
        bool match = false;
        const biggles::track& track_to_find = **it1;
        for (it2 = tc2.begin(); it2 not_eq tc2.end(); ++it2) {
            if (track_to_find == **it2) {
                match = true;
                break;
            }
        }
        if (not match) {
            const biggles::track &track2(**tc2.begin());
            diag("track not found '%s' ('%s')", track_to_str(track_to_find, part1.last_time_stamp()).c_str(),
                track_to_str(track2, part2.last_time_stamp()).c_str());
            diag(" * 1st %zu, last %zu, dd %f, from [%f, %f] to [%.10f, %.10f]",
                track_to_find.first_time_stamp(),
                track_to_find.last_time_stamp(),
                track_to_find.dynamic_drag(),
                track_to_find.min_location().get<0>(),
                track_to_find.min_location().get<1>(),
                track_to_find.max_location().get<0>(),
                track_to_find.max_location().get<1>()
                );
            diag(" * 1st %zu, last %zu, dd %f, from [%f, %f] to [%.10f, %.10f]",
                track2.first_time_stamp(),
                track2.last_time_stamp(),
                track2.dynamic_drag(),
                track2.min_location().get<0>(),
                track2.min_location().get<1>(),
                track2.max_location().get<0>(),
                track2.max_location().get<1>()
                );
            diag(" * min loc identical ? %d", int(track_to_find.min_location() == track2.min_location()));
            diag(" * max loc identical ? %d", int(track_to_find.max_location() == track2.max_location()));
            diag(" * max loc 0 identical ? %d", int(track_to_find.max_location().get<0>() == track2.max_location().get<0>()));
            diag(" * max loc 1 identical ? %d", int(track_to_find.max_location().get<1>() == track2.max_location().get<1>()));
            diag(" * diff max loc 1 %.50f", track_to_find.max_location().get<1>() - track2.max_location().get<1>());
            diag(" * diff max loc 0 %.50f", track_to_find.max_location().get<0>() - track2.max_location().get<0>());
            diag(" * obs identical ? %d", int(track_to_find.observations() == track2.observations()));
            diag(" * tracks identical ? %d", int(track_to_find == track2));
            return false;
        }
    }
    return true;
}


int test_pdr(const int total, const random_gen_paras& paras) {

    partition_ptr_t orignal_part_ptr(random_debug_partition(paras));

    capability_recorder_ptr cap_rec_ptr(new_cap_recorder(*orignal_part_ptr));
    orignal_part_ptr->set_capability_recorder(cap_rec_ptr);

    test_t tests;
    BOOST_FOREACH (const std::string& msg,  move_pdr_test(orignal_part_ptr, "Extend", total, tests)) {
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
    gen_pars.total = 10000;
    gen_pars.seed = 0;
    gen_pars.rgp.lambda = 2;
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

    diag("random seed = 0x%x", gen_pars.seed);
    return exit_status();
}
