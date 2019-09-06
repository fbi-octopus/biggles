#include <boost/assert.hpp>
#include <cmath>
#include <limits>
#include <utility>

#include "mh_moves/mh_moves.hpp"
#include "mh_moves/mh_move_weights.hpp"
#include "sampling/simple.hpp"
#include "mh_moves/utility.hpp"
#include "detail/fun.hpp"

namespace biggles {

namespace mh_moves
{

/// \brief TODO debug only
int encode3_obs(const observation& o) {
    return (int(x(o)) % 3 + 3 * (int(y(o)) % 3));
}

/// \brief TODO debug only
std::string o2s(const observation& o) {
    std::stringstream os;
    os << "[" << x(o) << ", " << y(o) << ", " << t(o) << "]";
    return os.str();
}

float log_prob_at_time_stamp(const clutter_t& clutter,
    const boost::shared_ptr<const track>& track_to_reduce, const time_stamp& ts, observation& ref_obs)
{
    float radius = std::abs(float(t(ref_obs)) - float(ts)) * detail::speed_of_light_;
    std::deque<observation> candidates = clutter.locate_near(ts, ref_obs, radius);
    std_normal_weight obs_weight(ref_obs);
    float noobs_weight = obs_weight.sigmas(2, ts);
    float total_weight = fun_sum(candidates.begin(), candidates.end(), obs_weight, noobs_weight);
    float sample_weight = noobs_weight;

    observation_collection::const_range range = track_to_reduce->observations().at_time_stamp(ts);
    if (range.first != range.second) { // obs present in original track
        ref_obs = *range.first;
        sample_weight = obs_weight(ref_obs);
    }
    return logf(sample_weight) - logf(total_weight);
}

float get_extend_log_prob(const partition& end_partition, const boost::shared_ptr<const track>& track_to_reduce,
    const boost::shared_ptr<const track>& new_track, const bool reduce_at_front)
{
    //float extend_log_prob_ = -logf(end_partition.tracks().size());
    const track_collection& tracks(end_partition.tracks());
    extend_weight_t extend_weight(end_partition.duration());

    float extend_log_prob_ = get_track_log_prob(extend_weight, tracks, new_track);
    time_stamp extension_first_time_stamp = end_partition.first_time_stamp();
    time_stamp extension_last_time_stamp = end_partition.last_time_stamp();
    if (reduce_at_front) {
        extension_last_time_stamp = new_track->first_time_stamp();
        observation ref_obs = new_track->observations().front();
        for (time_stamp ts = extension_last_time_stamp - 1; ts >= track_to_reduce->first_time_stamp(); --ts)
            extend_log_prob_ += log_prob_at_time_stamp(end_partition.clutter(), track_to_reduce, ts, ref_obs);
    } else {
        extension_first_time_stamp = new_track->last_time_stamp();
        observation ref_obs = new_track->observations().back();
        for (time_stamp ts = extension_first_time_stamp; ts < track_to_reduce->last_time_stamp(); ++ts)
            extend_log_prob_ += log_prob_at_time_stamp(end_partition.clutter(), track_to_reduce, ts, ref_obs);
    }
    extend_log_prob_ += -logf(extension_last_time_stamp - extension_first_time_stamp);
    float num_poss = new_track->num_possible_extensions(end_partition.first_time_stamp(), end_partition.last_time_stamp());
    extend_log_prob_ += -logf(   num_poss  );
    return extend_log_prob_;
}

bool reduce(const partition_ptr_t& original_ptr, partition_ptr_t& end_partition_ptr, float& proposal_density_ratio) {
    const clutter_t& clutter(original_ptr->clutter());
    const track_collection& tracks(original_ptr->tracks());

    capability_recorder_ptr cap_rec_ptr = original_ptr->get_capability_recorder();
    BOOST_ASSERT(cap_rec_ptr.get() not_eq 0); // there actually is something assigned

    // early-out: there must be at least one track
    if(tracks.empty())
        return false;

    // choose a track to reduce
    float log_track_prob(1.f);
    /*
    track_collection::const_iterator track_it(sampling::from_range(tracks.begin(), tracks.end(), log_track_prob));
    if(track_it == tracks.end())
        return false;

    // check we found one
    BOOST_ASSERT(track_it != tracks.end());
    BOOST_ASSERT(log_track_prob <= 0.f);
    const boost::shared_ptr<const track>& track_to_reduce(*track_it);
    */
    shared_const_track_ptr track_to_reduce;
    if (not select_track(reduce_weight_t(), tracks, track_to_reduce, log_track_prob))
        return false;

    // front or back? p=0.5 doesn't need to be considered if(!) it is also ignored in extend
    const bool reduce_at_front = sampling::uniform_int(0, 2) == 1;

    // create copies of the partition clutter and tracks
    clutter_ptr new_clutter(new clutter_t(clutter));
    boost::shared_ptr<track_collection> new_tracks(new track_collection(tracks));

    const observation_collection& obs_to_reduce(track_to_reduce->observations());
    observation_collection::const_iterator obs_iter = obs_to_reduce.begin();
    time_stamp new_last_ts = track_to_reduce->last_time_stamp();
    time_stamp new_first_ts = track_to_reduce->first_time_stamp();
    float log_time_stamp_prob(1.f);
    if (reduce_at_front) {
        size_t num_obs = track_to_reduce->size();
        std::advance(obs_iter, num_obs -2); // the second last observation (we want to keep that)
        time_stamp final_ts = t(*obs_iter);
        if (final_ts == new_first_ts)
            return false;
        new_first_ts = sampling::uniform_int(new_first_ts, final_ts) + 1;
        log_time_stamp_prob = -logf(final_ts - track_to_reduce->first_time_stamp());
    } else {
        std::advance(obs_iter, 1); // the second observation (we want to keep that)
        time_stamp initial_ts = t(*obs_iter)+1;
        if (initial_ts == new_last_ts)
            return false;
        new_last_ts = sampling::uniform_int(initial_ts, new_last_ts);
        log_time_stamp_prob = -logf(track_to_reduce->last_time_stamp() - initial_ts);
    }
    BOOST_FOREACH(const observation &o, obs_to_reduce) {
        if (t(o) < new_first_ts or t(o) >= new_last_ts)
            new_clutter->insert(o);
    }
    observation_collection::const_iterator first_obs_iter = obs_to_reduce.lower_bound_for_time_stamp(new_first_ts);
    observation_collection::const_iterator last_obs_iter = obs_to_reduce.lower_bound_for_time_stamp(new_last_ts);

    const float new_dynamic_drag = 1.0f;
    boost::shared_ptr<track> new_track(new track(
            new_first_ts, new_last_ts, first_obs_iter, last_obs_iter, new_dynamic_drag));
    BOOST_ASSERT(new_track->size() >= 2);

    cap_rec_ptr->add_erase_track(track_to_reduce);
    cap_rec_ptr->add_insert_track(new_track);
    cap_rec_ptr->set_editing_finished();
    // replace old track
    new_tracks->replace(track_to_reduce, new_track);

    // update partition
    end_partition_ptr = partition_ptr_t(
        new partition(original_ptr->pool(), new_tracks, new_clutter, original_ptr->expansion())
        );
    end_partition_ptr->set_first_time_stamp(original_ptr->first_time_stamp());
    end_partition_ptr->set_last_time_stamp(original_ptr->last_time_stamp());

    float reduce_log_prob = log_track_prob + log_time_stamp_prob - std::log(2.f);
    float extend_log_prob = get_extend_log_prob(*end_partition_ptr, track_to_reduce, new_track, reduce_at_front);
    proposal_density_ratio = extend_log_prob - reduce_log_prob;

    return true;

}


}

}
