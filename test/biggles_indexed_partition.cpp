#include <boost/format.hpp>
#include <iostream>
#include <sstream>

#include "biggles/detail/observation_reservoir.hpp"
#include "biggles/sampling/simple.hpp"
#include "biggles/tools/clutter_transform.hpp"
#include "biggles/tools/track_coll_transform.hpp"
#include "biggles/tools/max_creation.hpp"
#include "biggles/clutter.hpp"
#include "biggles/track_collection.hpp"
#include "biggles/partition.hpp"
#include "biggles/mh_moves/mh_moves.hpp"

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

bool observation_in_track_collection(const observation& o, const track_collection& tc) {
    for (track_collection::const_iterator tcit = tc.begin(); tcit not_eq tc.end(); ++tcit) {
        if ((*tcit)->observations().contains(o))
            return true;
    }
    return false;
}

void test_proc1() {
    typedef boost::shared_ptr<ObservationReservoir> observation_reservoir_var_ptr;
    observation_reservoir_var_ptr reservoir(new ObservationReservoir);
    for (time_stamp t = 0; t < 10; ++t) {
        size_t num = sampling::uniform_int(0, 10);
        for (size_t i = 0; i < num; ++i) {
            reservoir->push(sampling::uniform_real(0.f, 10.f), sampling::uniform_real(0.f, 10.f), t);
        }
    }
    std::deque<size_t> idx_clutter;
    for (size_t i = 0; i < reservoir->size(); ++i) {
        idx_clutter.push_back(i);
    }
    clutter_ptr clutter(new clutter_t());
    partition::track_collection_ptr tracks(new track_collection());

    // convert index clutter to observation clutter
    clutter_index_to_observation(reservoir, idx_clutter, *clutter);

    ok(idx_clutter.size() == clutter->size(), "clutter size invariant");

    // statr partition only has clutter
    partition_ptr_t start_partition( new partition(reservoir, tracks, clutter));

    // create some tracks
    partition_ptr_t half;
    max_creation(start_partition, half);
    ok(start_partition->volume() == half->volume(), "volume invariant");
    ok(half->observation_count() == reservoir->size(), "obs count invariant");
    ok(2*half->clutter().size() <= half->observation_count(), "at least half the obs are in tracks");

    ok(reservoir == half->pool(), "reservoir pointers identical");

    // converting partition "half" into indexed stuff
    clutter_index_t indexed_clutter;
    track_coll_index_t indexed_tracks;
    clutter_observation_to_index(half->pool(), half->clutter(), indexed_clutter);
    track_coll_observation_to_index(half->pool(), half->tracks(), indexed_tracks);

    // compare indexed and non-indexed stuff
    ok(half->tracks().size() == indexed_tracks.size(), "num tracks invariant");
    ok(half->clutter().size() == indexed_clutter.size(), "clutter size invariant");


    // convert indexed "half" back to observed collections
    clutter_ptr clutter2(new clutter_t());
    boost::shared_ptr<track_collection> tracks2(new track_collection());
    clutter_index_to_observation(reservoir, indexed_clutter, *clutter2);
    track_coll_index_to_observation(reservoir, indexed_tracks, *tracks2);
    // do tests
    ok(tracks2->size() == half->tracks().size(), "num tracks recovered");
    ok(clutter2->size() == half->clutter().size(), "clutter size recovered");
    int failcount = 0;
    for (clutter_t::const_iterator it = clutter2->begin(); it != clutter2->end(); ++it) {
        for (clutter_t::observation_container::const_iterator oit = it->second.begin(); oit != it->second.end(); ++oit) {
            failcount += half->clutter().contains(*oit) ? 0 : 1;
        }
    }
    ok(failcount == 0, "clutter observations recovered");
    failcount = 0;
    BOOST_FOREACH(const shared_const_track_ptr& tp, *tracks2) {
        BOOST_FOREACH (const observation& o, tp->observations()) {
            failcount += observation_in_track_collection(o, half->tracks()) ? 0 : 1;
        }
    }
    ok(failcount == 0, "track observations recovered");

    ObservationReservoir::const_iterator oit1 = reservoir->begin();
    ObservationReservoir::const_iterator oit2 = half->pool()->begin();
    failcount = 0;
    while (oit1 != reservoir->end()) {
        failcount += oit1++ == oit2++ ? 0 : 1;
    }
    ok(failcount == 0, "observations order is preserved");
    ok(oit2 == half->pool()->end(), "second observation pool was exhausted");

    // do some MH moves
    mh_moves::weighted_move_sampler move_sampler;
    move_sampler.add(mh_moves::BIRTH).add(mh_moves::DEATH).add(mh_moves::EXTEND).add(mh_moves::REDUCE)
                .add(mh_moves::MERGE).add(mh_moves::SPLIT).add(mh_moves::UPDATE).add(mh_moves::TRANSFER)
                .add(mh_moves::CROSS_OVER);
    partition_ptr_t to_part;
    partition_ptr_t from_part(half);
    capability_recorder_ptr cap_rec_ptr(new_cap_recorder(*from_part));
    const size_t nloops = 1000;
    float pdr = 0.f;
    for (size_t i = 0 ; i < nloops; ++i) {
        from_part->set_capability_recorder(cap_rec_ptr);
        bool success = mh_moves::propose(move_sampler(), from_part, to_part, pdr);
        if (success) {
            from_part = to_part;
            cap_rec_ptr->commit_changes(from_part->tracks());
        }
    }

    oit1 = from_part->pool()->begin();
    oit2 = half->pool()->begin();
    failcount = 0;
    while (oit1 != from_part->pool()->end()) {
        failcount += oit1++ == oit2++ ? 0 : 1;
    }
    ok(failcount == 0, "observations order is preserved");
    ok(from_part->volume() == half->volume(), "volume is preserved");
}

int main(int argc, char** argv) {

    plan_tests(15);
    //plan_no_plan();

    test_proc1();

    return exit_status();
}
