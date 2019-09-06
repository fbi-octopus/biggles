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

struct cond_track {
    std::string to, from;
    cond_track(const std::string &t, const std::string &f) : to(t), from(f) {}
    bool operator<=(const cond_track &o) const {
        return (this->from < o.from) or ((this->from == o.from) and (this->to <= o.to));
    }
    bool operator<(const cond_track &o) const {
        return (this->from < o.from) or ((this->from == o.from) and (this->to < o.to));
    }
    bool operator==(const cond_track &o) const {
        return (this->from == o.from) and (this->to == o.to);
    }
};

int test_pdr2(const int total, const random_gen_paras& paras) {
    partition_ptr_t original_part_ptr(random_debug_partition(paras));

    capability_recorder_ptr cap_rec_ptr(new_cap_recorder(*original_part_ptr));
    original_part_ptr->set_capability_recorder(cap_rec_ptr);

    test_t tests;
    BOOST_FOREACH (const std::string& msg,  move_pdr_test(original_part_ptr, "Reduce", total, tests)) {
        diag("%s", msg.c_str());
    }
    BOOST_FOREACH (const test_t::value_type& item, tests) {
        ok(item.second, "%s", item.first.c_str());
    }
    return 0;
}

int test1() {
    observation_collection t1_obs;
    observation_collection cl_obs;
    t1_obs.insert(new_obs(0, 0, 0));
    t1_obs.insert(new_obs(0, 0, 1));
    t1_obs.insert(new_obs(0, 0, 2));
    t1_obs.insert(new_obs(0, 0, 3));

    //cl_obs.insert(new_obs(0, 1, 0));
    //cl_obs.insert(new_obs(0, 1, 1));
    //cl_obs.insert(new_obs(0, 1, 2));
    //cl_obs.insert(new_obs(0, 1, 4));
    cl_obs.insert(new_obs(0, 0, 4));
    cl_obs.insert(new_obs(1, 1, 5));
    cl_obs.insert(new_obs(1, 1, 6));
    cl_obs.insert(new_obs(1, 1, 7));
    cl_obs.insert(new_obs(1, 1, 8));
    cl_obs.insert(new_obs(1, 1, 9));
    track track1(0, 5, t1_obs.begin(), t1_obs.end(), 1.0);
    boost::shared_ptr<track_collection> tracks(new track_collection());
    tracks->insert(track1);
    clutter_ptr clutter(new clutter_t(cl_obs.begin(), cl_obs.end()));

    partition_ptr_t orignal_part_ptr(new partition(tracks, clutter));

    //size_t lenc = clutter->size();

    //int iclutter = 0;
    //int itracks = 0;
    //int iboth = 0;
    int total = 50;
    //int failed = 0;

    int num_obs_after_reduce_eq4 = 0;
    int num_dur_after_reduce_eq4 = 0;
    int num_obs_after_extend_eq5 = 0;
    int num_dur_after_extend_ge5 = 0;
    int num_dur_after_extend_gt5 = 0;
    int num_dur_after_extend_eq5 = 0;
    int pdr_reduce_plus_pdr_extend_eq0 = 0;
    int reduce_always_successful = 0;

    for (int i = 0; i < total; ++i) {
        partition_ptr_t extend_part_ptr;
        float extend_pdr = 0.f;
        bool done = false;
        while (not done) {
            done = mh_moves::extend(orignal_part_ptr, extend_part_ptr, extend_pdr);
        }
        num_obs_after_extend_eq5 += int((*extend_part_ptr->tracks().begin())->size() == 5);
        num_dur_after_extend_ge5 += int((*extend_part_ptr->tracks().begin())->duration() >= 5);
        num_dur_after_extend_gt5 += int((*extend_part_ptr->tracks().begin())->duration() > 5);
        num_dur_after_extend_eq5 += int((*extend_part_ptr->tracks().begin())->duration() == 5);
        //
        float reduce_prob = 0.f;
        done = false;
        partition_ptr_t reduce_part_ptr;
        done = mh_moves::reduce(extend_part_ptr, reduce_part_ptr, reduce_prob);
        reduce_always_successful += int(done);
        num_obs_after_reduce_eq4 += int((*reduce_part_ptr->tracks().begin())->size() == 4);
        num_dur_after_reduce_eq4 += int((*reduce_part_ptr->tracks().begin())->duration() == 4);
        //track_collection::const_iterator track0 = extend_partition.tracks().begin();
        pdr_reduce_plus_pdr_extend_eq0 += int(reduce_prob + extend_pdr == 0);
        //diag("range %ld, %ld, pdr sum = %d", (*reduce_partition.tracks().begin())->first_time_stamp(),
        //        (*reduce_partition.tracks().begin())->last_time_stamp(), reduce_prob + extend_pdr);
    }
    ok(track1.size() == 4, "track1.size() == 4");
    ok(num_obs_after_reduce_eq4 == total, "num_obs_after_reduce_eq4 == total");
    diag("num_dur_after_reduce_eq4 = %d/%d", num_dur_after_reduce_eq4, total);
    //ok(num_dur_after_reduce_eq4 == total, "num_dur_after_reduce_eq4 == total");
    ok(num_obs_after_extend_eq5 == total, "num_obs_after_extend_eq5 == total");
    ok(num_dur_after_extend_ge5 == total, "num_dur_after_extend_ge5 == total");
    ok(num_dur_after_extend_gt5 > 0, "num_dur_after_extend_gt5 > 0");
    ok(num_dur_after_extend_eq5 > 0, "num_dur_after_extend_eq5 > 0");
    //ok(pdr_reduce_plus_pdr_extend_eq0 == total, "pdr_reduce_plus_pdr_extend_eq0 == total");
    ok(reduce_always_successful == total, "reduce_always_successful == total");
    return 0;
}


int main(int argc, char** argv)
{
    general_paras gen_pars;
    gen_pars.total = 100000;
    gen_pars.seed = 0;
    gen_pars.rgp.lambda = 4;
    gen_pars.rgp.p_no = 0.2;
    gen_pars.rgp.p_yes = 0.5;
    gen_pars.rgp.p_tr = 0.5;
    gen_pars.rgp.min_tracks = 2;
    gen_pars.rgp.max_tracks = 3;
    if (not parse_args(argc, argv, gen_pars))
        return exit_status();
    biggles::detail::seed_prng(gen_pars.seed);

    plan_no_plan();

    test_pdr2(gen_pars.total, gen_pars.rgp);

    diag("random seed = 0x%x", gen_pars.seed);
    //return exit_status();
}
