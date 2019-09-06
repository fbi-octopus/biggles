#include <boost/assert.hpp>
#include <cmath>
#include <limits>
#include <utility>

#include "mh_moves/mh_moves.hpp"
#include "mh_moves/mh_move_weights.hpp"
#include "sampling/simple.hpp"
#include "mh_moves/utility.hpp"

namespace biggles {

namespace mh_moves
{

/// \brief get the probability that pieces are merged given tracks (which contains pieces)
float get_log_merge_prob(const partition_ptr_t &end_partition, const shared_const_track_ptr &front_piece,
    const shared_const_track_ptr &back_piece)
{
    const track_collection &tracks = end_partition->tracks();
    BOOST_ASSERT(tracks.contains(front_piece));
    BOOST_ASSERT(tracks.contains(back_piece));
    BOOST_ASSERT(tracks_could_be_merged2(*front_piece, *back_piece)>0.f);
    front_weight_merging fw(end_partition->first_time_stamp(), end_partition->last_time_stamp());
    float front_piece_weight = fw(front_piece);
    float total_weight = fun_sum(tracks.begin(), tracks.end(), fw, 0.f);
    float log_front_piece_prob = std::log(front_piece_weight) - std::log(total_weight);
    back_weight_merging bw(front_piece);
    float back_piece_weight = bw(back_piece);
    BOOST_ASSERT(back_piece_weight > 0);
    total_weight = fun_sum(tracks.begin(), tracks.end(), bw, 0.f);
    float log_back_piece_prob = std::log(back_piece_weight) - std::log(total_weight);
    return log_front_piece_prob + log_back_piece_prob;
}

bool split(const partition_ptr_t& start_partition_ptr, partition_ptr_t& end_partition, float& proposal_density_ratio)
{
    const track_collection& tracks(start_partition_ptr->tracks());

    capability_recorder_ptr cap_rec_ptr = start_partition_ptr->get_capability_recorder();
    BOOST_ASSERT(cap_rec_ptr.get() not_eq 0); // there actually is something assigned

    // early-out: there must be at least one track
    if(tracks.empty())
        return false;

    // choose a track to split
    float log_track_prob(1.f);
    /*
    track_collection::const_iterator track_it(sample_track_with_minimum_size(tracks, 4, log_track_prob));
    if(track_it == tracks.end())
        return false;

    // check we found one
    BOOST_ASSERT(track_it != tracks.end());
    BOOST_ASSERT((*track_it)->size() >= 4);
    BOOST_ASSERT(log_track_prob <= 0.f);
    const shared_const_track_ptr& track_to_split(*track_it);
        */
    shared_const_track_ptr track_to_split;
    if (not select_track(split_weight_t(), tracks, track_to_split, log_track_prob)) {
        return false;
    }

    // choose a splitting point
    time_stamp ts_to_split_at(0);
    float log_ts_prob(1.f);
    // IMPORTANT: BOOST_VERIFY will still evaluate its argument even when not in debug mode
    BOOST_VERIFY(sample_time_stamp_to_split_track(track_to_split, ts_to_split_at, log_ts_prob));
    BOOST_ASSERT(log_ts_prob <= 0.f);

    // create new tracks
    float log_splitting_prob = 0.f; // due to noobs at the splitting end of the new tracks
    shared_const_track_ptr_pair new_track_pair(
        split_track_at_time_stamp(track_to_split, ts_to_split_at, log_splitting_prob));

    // sanity check tracks
    BOOST_ASSERT(new_track_pair.first->size() + new_track_pair.second->size() == track_to_split->size());
    BOOST_ASSERT(new_track_pair.first->size() >= 2);
    BOOST_ASSERT(new_track_pair.second->size() >= 2);
    BOOST_ASSERT(new_track_pair.first->last_time_stamp() <= new_track_pair.second->first_time_stamp());
    BOOST_ASSERT(tracks_could_be_merged(*new_track_pair.first, *new_track_pair.second));

    // create copies of the tracks. The clutter is unaffected
    shared_track_collection_ptr new_tracks(new track_collection(tracks));


    cap_rec_ptr->add_erase_track(track_to_split);
    cap_rec_ptr->add_insert_track(new_track_pair.first);
    cap_rec_ptr->add_insert_track(new_track_pair.second);
    cap_rec_ptr->set_editing_finished();

    // remove the track from the collection
    new_tracks->remove(track_to_split);
    new_tracks->insert(new_track_pair.first);
    new_tracks->insert(new_track_pair.second);

    // update partition
    end_partition = partition_ptr_t(
        new partition(start_partition_ptr->pool(), new_tracks, start_partition_ptr->clutter_ptr(), start_partition_ptr->expansion())
        );
    end_partition->set_first_time_stamp(start_partition_ptr->first_time_stamp());
    end_partition->set_last_time_stamp(start_partition_ptr->last_time_stamp());


    // calculate forward (split) log prob.
    float forward_log_prob = log_ts_prob + log_track_prob + log_splitting_prob;

    float backward_log_prob = get_log_merge_prob(end_partition, new_track_pair.first, new_track_pair.second);

    proposal_density_ratio = backward_log_prob - forward_log_prob;

    return true;
}

}

}
