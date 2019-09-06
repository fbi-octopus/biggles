#include <boost/assert.hpp>
#include <cmath>
#include <limits>
#include <utility>

#include "mh_moves/mh_moves.hpp"
#include "mh_moves/mh_move_weights.hpp"
#include "sampling/simple.hpp"
#include "mh_moves/utility.hpp"
#include "tools/debug.hpp"

namespace biggles {

namespace mh_moves
{

float get_log_splitting_prob(shared_const_track_ptr& track1, shared_const_track_ptr& track2) {
    shared_const_track_ptr tr_lo;
    shared_const_track_ptr tr_hi;
    if (track1->first_time_stamp() < track2->first_time_stamp()) {
        tr_lo = track1;
        tr_hi = track2;
    } else {
        tr_lo = track2;
        tr_hi = track1;
    }
    float log_splitting_prob = 0.f;
    if (tr_hi->observations().first_time_stamp() > tr_lo->observations().last_time_stamp()) // is there a gap?
    {
        time_stamp last_obs_first_track = tr_lo->observations().last_time_stamp();
        time_stamp split_ts = tr_hi->first_time_stamp();
        log_splitting_prob = -logf(split_ts - last_obs_first_track + 1);
    }
    return log_splitting_prob;
}

/// \brief gets the back merge partner for the front piece
bool get_merge_partner(const shared_const_track_ptr &piece1, const track_collection &tracks,
    shared_const_track_ptr &partner, float &log_prob)
{
    back_weight_merging back_weight(piece1);
    track_collection::const_iterator it =
        sampling::weighted_choice(tracks.begin(), tracks.end(), back_weight, log_prob);
    if (it == tracks.end())
        return false;
    partner = *it;
    return true;
}

bool merge(const partition_ptr_t& start_partition_ptr, partition_ptr_t& end_partition_ptr, float& proposal_density_ratio)
{
    //typedef std::pair<shared_const_track_ptr, shared_const_track_ptr> track_pair;

    capability_recorder_ptr cap_rec_ptr = start_partition_ptr->get_capability_recorder();
    BOOST_ASSERT(cap_rec_ptr.get() not_eq 0); // there actually is something assigned

    const track_collection& tracks(start_partition_ptr->tracks());
    if (tracks.empty())
        return false;

    float log_track_sample = 0.f;
    front_weight_merging fw(start_partition_ptr->first_time_stamp(), start_partition_ptr->last_time_stamp());
    track_collection::const_iterator it =
        sampling::weighted_choice(tracks.begin(), tracks.end(), fw, log_track_sample);
    if (it == tracks.end())
        return false;
    shared_const_track_ptr piece1(*it);

    shared_const_track_ptr piece2;
    float log_partner_prob = 0.f;
    bool partner_found = get_merge_partner(piece1, tracks, piece2, log_partner_prob);

    if (not partner_found)
        return false;

    // create a new track by merging old ones
    // we also re-sample a new dynamic drag as 1
    shared_track_ptr new_track(new track(*piece1, *piece2, 1.f));
    BOOST_ASSERT(new_track->size() == piece1->size() + piece2->size());
    BOOST_ASSERT(new_track->size() >= 4);
    BOOST_ASSERT(new_track->duration() >= 3);

    // create copies of the tracks. The clutter is unaffected
    shared_track_collection_ptr new_tracks(new track_collection(tracks));

    cap_rec_ptr->add_erase_track(piece1);
    cap_rec_ptr->add_erase_track(piece2);
    cap_rec_ptr->add_insert_track(new_track);
    cap_rec_ptr->set_editing_finished();

    // remove old tracks
    new_tracks->remove(piece1);
    new_tracks->remove(piece2);

    // insert new one
    new_tracks->insert(new_track);

    // update partition
    end_partition_ptr = partition_ptr_t(
        new partition(start_partition_ptr->pool(), new_tracks, start_partition_ptr->clutter_ptr(),
                        start_partition_ptr->expansion()));

    // calculate forward log prob.
    float forward_log_prob = log_track_sample + log_partner_prob;

    // calculate backward log prob. = prob of selecting new track and splitting at right place
    /*
    float backward_log_prob = 1.f;
    BOOST_VERIFY(sample_track_with_minimum_size(*new_tracks, 4, backward_log_prob) != new_tracks->end());
    BOOST_ASSERT(backward_log_prob <= 0.f);
    */
    float backward_log_prob = get_track_log_prob(split_weight_t(), *new_tracks, new_track);
    float log_ts_prob = 0.f;
    time_stamp ts_to_split = std::max(piece2->first_time_stamp(),
                                piece1->first_time_stamp());
    log_ts_prob = get_time_stamp_to_split_track_prob(new_track, ts_to_split);
    //if (not std::isfinite(log_ts_prob)) OK(log_ts_prob);
    //sample_time_stamp_to_split_track(new_track, foo, log_ts_prob); // this is non-sense
    backward_log_prob += log_ts_prob;
    backward_log_prob += get_log_splitting_prob(piece1, piece2);

    proposal_density_ratio = backward_log_prob - forward_log_prob;

    return true;
}

}

}
