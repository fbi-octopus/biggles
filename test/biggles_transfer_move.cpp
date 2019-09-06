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

typedef std::set<std::string> string_set;

partition_ptr_t get_partition1() {
    observation_collection obs_coll;
    obs_coll.insert(new_obs(0,0,5));
    obs_coll.insert(new_obs(0,0,6));
    obs_coll.insert(new_obs(0,0,8));
    boost::shared_ptr<track_collection> tracks(new track_collection());
    tracks->insert(track(5, 9, obs_coll.begin(), obs_coll.end(), 1.f));
    obs_coll.clear();
    obs_coll.insert(new_obs(1,0,5));
    obs_coll.insert(new_obs(1,0,6));
    obs_coll.insert(new_obs(1,0,7));
    obs_coll.insert(new_obs(1,0,9));
    tracks->insert(track(4, 10, obs_coll.begin(), obs_coll.end(), 1.f));
    clutter_ptr clutter(new clutter_t());
    return partition_ptr_t(new partition(tracks, clutter, partition::expansion_2d(0.f, 2.f, 0.f, 2.f)));
}

int basic_functionality_test() {
    diag("test if it does anything at all");
    observation_collection t1_obs;
    observation_collection t2_obs;
    t1_obs.insert(new_obs(1, 1, 0));
    t2_obs.insert(new_obs(1, 0, 1));
    t1_obs.insert(new_obs(1, 1, 2));
    t2_obs.insert(new_obs(1, 0, 3));
    t1_obs.insert(new_obs(1, 1, 4));
    t2_obs.insert(new_obs(1, 0, 5));
    boost::shared_ptr<track_collection> tracks(new track_collection());
    shared_const_track_ptr tr1(new track(0, 6, t1_obs.begin(), t1_obs.end(), 1.f));
    shared_const_track_ptr tr2(new track(0, 6, t2_obs.begin(), t2_obs.end(), 1.f));
    clutter_ptr clutter(new clutter_t());

    tracks->insert(tr1);
    tracks->insert(tr2);
    partition_ptr_t orignal_partition(new partition(tracks, clutter));
    partition_ptr_t end_partition;
    float pdr = 0;
    mh_moves::transfer(orignal_partition, end_partition, pdr);
    std::cout << boost::format("original '%s'") % partition_to_string(*orignal_partition) << std::endl;
    std::cout << boost::format("forward  '%s'") % partition_to_string(*end_partition) << std::endl;
    return 0;
}

int basic_functionality_test2() {
    diag("A special case test if the borders are treated correctly");
    observation_collection t1_obs;
    observation_collection t2_obs;
    t1_obs.insert(new_obs(1, 1, 2));
    t1_obs.insert(new_obs(1, 1, 4));
    t2_obs.insert(new_obs(1, 0, 2));
    t2_obs.insert(new_obs(1, 0, 4));
    t2_obs.insert(new_obs(1, 0, 5));
    boost::shared_ptr<track_collection> tracks(new track_collection());
    shared_const_track_ptr tr1(new track(0, 5, t1_obs.begin(), t1_obs.end(), 1.f));
    shared_const_track_ptr tr2(new track(2, 6, t2_obs.begin(), t2_obs.end(), 1.f));
    clutter_ptr clutter(new clutter_t());
    tracks->insert(tr1);
    tracks->insert(tr2);
    partition_ptr_t orignal_partition(new partition(tracks, clutter));
    std::cout << boost::format("original '%s'") % partition_to_string(*orignal_partition) << std::endl;
    partition_ptr_t end_partition;
    float pdr = 0;
    mh_moves::transfer(orignal_partition, end_partition, pdr);
    std::cout << boost::format("forward  '%s'") % partition_to_string(*end_partition) << std::endl;
    return 0;
}

int test_pdr(const int total, const random_gen_paras& paras) {
    partition_ptr_t orignal_partition(random_debug_partition(paras));

    capability_recorder_ptr cap_rec_ptr(new_cap_recorder(*orignal_partition));
    orignal_partition->set_capability_recorder(cap_rec_ptr);

    diag("%s", "original partition created");
    diag("num tracks = %zu", orignal_partition->tracks().size());

    test_t tests;
    BOOST_FOREACH (const std::string& msg,  move_pdr_test(orignal_partition, "Transfer", total, tests)) {
        diag("%s", msg.c_str());
    }
    BOOST_FOREACH (const test_t::value_type& item, tests) {
        ok(item.second, "%s", item.first.c_str());
    }
    return 0;
}

void test_3() {
    partition_ptr_t orig_part(get_partition1());
    biggles_move testmove("Transfer");
    partition_ptr_t forw_part;
    partition_ptr_t back_part;
    int max_attempts = 10000;
    int attempts = 0;
    float pdr = 0.f;
    std::string orig_string = partition_to_string(*orig_part);
    bool ok = false;
    std::map<std::string, partition_ptr_t> forw_partitions;
    typedef std::map<std::string, partition_ptr_t>::iterator part_map_iter;
    std::map<std::string, int> forw_count;
    std::map<std::string, float> forw_calc;
    while (attempts++ < max_attempts) {
        ok = testmove(orig_part, forw_part, pdr);
        if (not ok) {
            forw_count[orig_string]++;
        } else {
            std::string part_string(partition_to_string(*forw_part));
            forw_count[part_string]++;
            part_map_iter it = forw_partitions.find(part_string);
            if (it == forw_partitions.end()) {
                forw_partitions.insert(std::make_pair(part_string, partition_ptr_t(new partition(*forw_part))));
                forw_calc[part_string] = pdr;
            }
        }
    }
    diag("original %s", partition_to_string(*orig_part).c_str());
    std::map<std::string, int> back_count;
    std::map<std::string, float> back_calc;
    for (part_map_iter it = forw_partitions.begin(); it not_eq forw_partitions.end(); ++it) {
        std::string ps = it->first;
        diag("  %s,   %.2f%%, %.2f ", ps.c_str(), 100.f * float(forw_count[ps])/float(max_attempts),
            expf(forw_calc[ps]));
        back_count.clear();
        attempts = 0;
        forw_part = it->second;
        while (attempts++ < max_attempts) {
            ok = testmove(forw_part, back_part, pdr);
            if (ok) {
                std::string part_string(partition_to_string(*back_part));
                back_count[part_string]++;
                back_calc[part_string] = pdr;
            }
        }
        for (std::map<std::string, int>::iterator cit = back_count.begin(); cit != back_count.end(); ++cit) {
            diag("    %s, %.2f%%, %.2f", cit->first.c_str(), 100.f * float(cit->second)/float(max_attempts),
                expf(back_calc[cit->first]));
        }
    }
}


int main(int argc, char** argv) {
    general_paras gen_pars;
    gen_pars.total = 10000;
    gen_pars.seed = 0;
    gen_pars.rgp.lambda = 2;
    gen_pars.rgp.p_no = 0.2;
    gen_pars.rgp.p_yes = 0.2;
    gen_pars.rgp.p_tr = 0.5;
    gen_pars.rgp.min_tracks = 2;
    gen_pars.rgp.max_tracks = 3;
    if (not parse_args(argc, argv, gen_pars))
        return exit_status();
    biggles::detail::seed_prng(gen_pars.seed);

    plan_no_plan();
    diag("random seed = 0x%x", gen_pars.seed);
    test_pdr(gen_pars.total, gen_pars.rgp);
    //basic_functionality_test2();
    //test_3();

    diag("random seed = 0x%x", gen_pars.seed);
    return exit_status();
}
