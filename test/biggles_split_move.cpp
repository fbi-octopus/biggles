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
#include <string>
#include <cmath>
#include "biggles/tools/sundries.hpp"

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

int test1() {
    observation_collection obs_coll;
    obs_coll.insert(new_obs(0, 0, 0));
    obs_coll.insert(new_obs(0, 0, 1));
    obs_coll.insert(new_obs(0, 0, 2));
    obs_coll.insert(new_obs(0, 0, 3));
    obs_coll.insert(new_obs(0, 0, 4));
    obs_coll.insert(new_obs(0, 0, 5));
    observation_collection::const_iterator oter = obs_coll.begin();
    std::advance(oter, 5);
    track track1(0, 5, obs_coll.begin(), oter, 1.0);
    track track2(0, 6, obs_coll.begin(), obs_coll.end(), 1.0);
    ok1(long(track1.size()) == (track1.last_time_stamp()-track1.first_time_stamp()));
    track_collection tc;
    tc.insert(track1);
    tc.insert(track2);


    time_stamp ts_to_split_at(0);
    float log_ts_prob(1.f);
    const shared_const_track_ptr& track_to_split(*tc.begin());

    for (int i = 0; i < 10; ++i) {
        bool res=mh_moves::sample_time_stamp_to_split_track(track_to_split, ts_to_split_at, log_ts_prob);
        std::cout << "OK = " << res << "; time stamp = " << ts_to_split_at
            << "; prob = " << exp(log_ts_prob) << std::endl;
    }


    track_collection::const_iterator iter = tc.begin();
    std::advance(iter, 1);

    const shared_const_track_ptr& track_to_split2(*iter);
    for (int i = 0; i < 10; ++i) {
        bool res=mh_moves::sample_time_stamp_to_split_track(track_to_split2, ts_to_split_at, log_ts_prob);
        std::cout << "OK = " << res << "; time stamp = " << ts_to_split_at
            << "; prob = " << exp(log_ts_prob) << std::endl;
    }

    boost::shared_ptr<track_collection> tracks(new track_collection());
    clutter_ptr clutter(new clutter_t());
    partition_ptr_t orignal_part_ptr(new partition(tracks, clutter));

    tracks->insert(track2);
    int total = 1000;
    int num_tracks_after_split_eq2 = 0;
    int num_tracks_after_merge_eq1 = 0;
    int pdr_split_plus_pdr_merge_eq0 = 0;
    int merge_always_successful = 0;
    std::map<float, int> split_p_count;
    std::map<float, int> merge_p_count;
    for (int i = 0; i < total ; ++i) {
        bool done = false;
        partition_ptr_t split_part_ptr;
        float split_prob = 0.f;
        while (not done) {
            done = mh_moves::split(orignal_part_ptr, split_part_ptr, split_prob);
        }
        partition_ptr_t merge_part_ptr;
        float merge_prob = 0.f;
        done = mh_moves::merge(split_part_ptr, merge_part_ptr, merge_prob);
        num_tracks_after_split_eq2 += int(split_part_ptr->tracks().size() == 2);
        num_tracks_after_merge_eq1 += int(merge_part_ptr->tracks().size() == 1);
        pdr_split_plus_pdr_merge_eq0 += int(split_prob + merge_prob == 0);
        merge_always_successful += int(done);
        split_p_count[split_prob]++;
        merge_p_count[merge_prob]++;
    }
    ok(num_tracks_after_split_eq2 == total, "num_tracks_after_split_eq2");
    ok(num_tracks_after_merge_eq1 == total, "num_tracks_after_merge_eq1");
    ok(pdr_split_plus_pdr_merge_eq0 == total, "pdr_split_plus_pdr_merge_eq0");
    ok(merge_always_successful == total, "merge_always_successful");
    ok(merge_p_count.size() == 1, "one merge probability");
    ok(split_p_count.size() == 1, "one split probability");

    observation_collection track10_obs;
    observation_collection track11_obs;
    observation_collection track12_obs;
    observation_collection track13_obs;
    track10_obs.insert(new_obs(0, 0, 0));
    track10_obs.insert(new_obs(0, 0, 1));
    track10_obs.insert(new_obs(0, 0, 2));
    track10_obs.insert(new_obs(0, 0, 3));
    track10_obs.insert(new_obs(0, 0, 4));
    track11_obs.insert(new_obs(0, 0, 5));
    track11_obs.insert(new_obs(0, 0, 6));
    track11_obs.insert(new_obs(0, 0, 7));
    track11_obs.insert(new_obs(0, 0, 8));
    track11_obs.insert(new_obs(0, 0, 9));
    track12_obs.insert(new_obs(0, 0, 80));
    track12_obs.insert(new_obs(0, 0, 81));
    track12_obs.insert(new_obs(0, 0, 82));
    track12_obs.insert(new_obs(0, 0, 83));
    track12_obs.insert(new_obs(0, 0, 84));
    track13_obs.insert(new_obs(0, 0, 85));
    track13_obs.insert(new_obs(0, 0, 86));
    track13_obs.insert(new_obs(0, 0, 87));
    track13_obs.insert(new_obs(0, 0, 88));
    track13_obs.insert(new_obs(0, 0, 89));

    track track10( 0,  5, track10_obs.begin(), track10_obs.end(), 1.0);
    track track11( 5, 10, track11_obs.begin(), track11_obs.end(), 1.0);
    track track12(80, 85, track12_obs.begin(), track12_obs.end(), 1.0);
    track track13(85, 90, track13_obs.begin(), track13_obs.end(), 1.0);

    boost::shared_ptr<track_collection> tracks2(new track_collection());
    tracks2->insert(track10);
    tracks2->insert(track11);
    tracks2->insert(track12);
    tracks2->insert(track13);
    /*
    std::cout
        << "tracks could be merged?" << std::endl
        << " track10, track11 : " << mh_moves::tracks_could_be_merged(track10, track11) << std::endl
        << " track10, track12 : " << mh_moves::tracks_could_be_merged(track10, track12) << std::endl
        << " track10, track13 : " << mh_moves::tracks_could_be_merged(track10, track13) << std::endl
        << " track11, track12 : " << mh_moves::tracks_could_be_merged(track11, track12) << std::endl
        << " track11, track13 : " << mh_moves::tracks_could_be_merged(track11, track13) << std::endl
        << " track12, track13 : " << mh_moves::tracks_could_be_merged(track12, track13) << std::endl
    ;
    */
    ok1(clutter->size() == 0);
    partition_ptr_t  orignal_part_ptr2(new partition(tracks2, clutter));
    for (int i = 0; i < 100; ++i) {
        partition_ptr_t merge_part_ptr;
        float merge_prob = 0.f;
        mh_moves::merge(orignal_part_ptr2, merge_part_ptr, merge_prob);
        const track_collection &merge_tracks =  merge_part_ptr->tracks();
        std::cout << "num tracks = " << merge_tracks.size() << std::endl;
        for (track_collection::const_iterator iter= merge_tracks.begin(); iter != merge_tracks.end(); ++iter) {
            std::cout
                << " * duration = " << (*iter)->duration()
                << " - num obs = " << (*iter)->size()
                << std::endl;
        }
    }

    std::cout
        << "observations could be neighbours ? "
        << mh_moves::observations_could_be_neighbours(track12_obs.front(), track10_obs.back())
        << " dt = " << std::abs(t(track12_obs.front()) - t(track10_obs.back()))
        << std::endl
    ;
    return 0;
}

int test_pdr(const int total, const random_gen_paras& paras) {
    partition_ptr_t orignal_partition(random_debug_partition(paras));

    capability_recorder_ptr cap_rec_ptr(new capability_recorder(
        orignal_partition->first_time_stamp(), orignal_partition->last_time_stamp(), orignal_partition->tracks()));
    orignal_partition->set_capability_recorder(cap_rec_ptr);


    test_t tests;
    BOOST_FOREACH (const std::string& msg,  move_pdr_test(orignal_partition, "Split", total, tests)) {
        diag("%s", msg.c_str());
    }
    BOOST_FOREACH (const test_t::value_type& item, tests) {
        ok(item.second, "%s", item.first.c_str());
    }
    return 0;
}

int test_pdr2(const int total, const random_gen_paras& paras) {
    observation_collection obs_coll;
    obs_coll.insert(new_obs(0, 0, 0));
    obs_coll.insert(new_obs(0, 0, 1));
    obs_coll.insert(new_obs(0, 0, 2));
    obs_coll.insert(new_obs(0, 2, 5));
    obs_coll.insert(new_obs(0, 2, 6));
    obs_coll.insert(new_obs(0, 2, 7));
    obs_coll.insert(new_obs(0, 2, 8));
    obs_coll.insert(new_obs(0, 2, 9));
    track track2(0, 10, obs_coll.begin(), obs_coll.end(), 1.0);

    boost::shared_ptr<track_collection> tracks(new track_collection());
    clutter_ptr clutter(new clutter_t());
    partition_ptr_t orignal_partition(new partition(tracks, clutter));

    tracks->insert(track2);

    test_t tests;
    BOOST_FOREACH (const std::string& msg,  move_pdr_test(orignal_partition, "Split", total, tests)) {
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
    gen_pars.rgp.lambda = 4;
    gen_pars.rgp.p_no = 0.;
    gen_pars.rgp.p_yes = 0.;
    gen_pars.rgp.p_tr = 0.5;
    gen_pars.rgp.min_tracks = 2;
    gen_pars.rgp.max_tracks = 3;
    if (not parse_args(argc, argv, gen_pars))
        return exit_status();
    biggles::detail::seed_prng(gen_pars.seed);
    diag("random seed = 0x%x", gen_pars.seed);

    plan_no_plan();

    test_pdr(gen_pars.total, gen_pars.rgp);
    //test_pdr2(gen_pars.total, gen_pars.rgp);

    diag("random seed = 0x%x", gen_pars.seed);
    return exit_status();
}
