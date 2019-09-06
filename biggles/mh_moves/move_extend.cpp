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


float extend_at_beginning_(const boost::shared_ptr<const track>& track_to_extend,
    const time_stamp& final_time_stamp, const clutter_t& clutter,
    std::deque<observation>& new_track_obs)
{
    float log_prob = 0.f;
    observation last_obs = track_to_extend->observations().front();
    new_track_obs.clear();
    weights_t weights;
    for (time_stamp ts = track_to_extend->first_time_stamp() -1 ; ts >= final_time_stamp; --ts) {
        float radius = float(t(last_obs) - ts) * detail::speed_of_light_;
        std::deque<observation> candidates = clutter.locate_near(ts, last_obs, radius);

        weights.assign(candidates.size(), 0.f);
        std_normal_weight obs_weight(last_obs);
        float noobs_weight = obs_weight.sigmas(2, ts);
        std::transform(candidates.begin(), candidates.end(), weights.begin(), obs_weight);
        candidates.push_back(observation());
        weights.push_back(noobs_weight);
        float obs_prob = 1.f;
        observation new_obs =
            sampling::weighted_choice(candidates.begin(), weights.begin(), weights.end(), obs_prob);
        if (not noobs(new_obs)) {
            new_track_obs.push_back(new_obs);
            last_obs = new_obs;
        }
        log_prob += obs_prob;
    }
    return log_prob;
}
float extend_at_end_(const boost::shared_ptr<const track>& track_to_extend,
    const time_stamp& final_time_stamp, const clutter_t& clutter,
    std::deque<observation>& new_track_obs)
{
    float log_prob = 0.f;
    observation last_obs = track_to_extend->observations().back();
    new_track_obs.clear();
    weights_t weights;
    for (time_stamp ts = track_to_extend->last_time_stamp(); ts < final_time_stamp; ++ts) {
        float radius = float(ts - t(last_obs)) * detail::speed_of_light_;
        std::deque<observation> candidates = clutter.locate_near(ts, last_obs, radius);
        weights.assign(candidates.size(), 0.f);
        std_normal_weight obs_weight(last_obs);
        float noobs_weight = obs_weight.sigmas(2, ts);
        std::transform(candidates.begin(), candidates.end(), weights.begin(), obs_weight);
        candidates.push_back(observation());
        weights.push_back(noobs_weight);
        float obs_prob = 1.f;
        observation new_obs =
            sampling::weighted_choice(candidates.begin(), weights.begin(), weights.end(), obs_prob);
        if (not noobs(new_obs)) {
            new_track_obs.push_back(new_obs);
            last_obs = new_obs;
        }
        log_prob += obs_prob;
    }
    return log_prob;
}

std::pair<bool, float> new_last_time_stamp_(const boost::shared_ptr<const track>& track_to_extend,
    const partition_ptr_t& original, time_stamp& new_last_time_stamp)
{
    new_last_time_stamp = track_to_extend->last_time_stamp();
    time_stamp max_time_step = original->last_time_stamp() - new_last_time_stamp;
    if (max_time_step == 0) {
        return std::make_pair(false, 0.f); // nothing to extend at the end
    }
    new_last_time_stamp += sampling::uniform_int(0, max_time_step) + 1;
    BOOST_ASSERT(new_last_time_stamp <= original->last_time_stamp());
    BOOST_ASSERT(new_last_time_stamp > track_to_extend->last_time_stamp());
    return std::make_pair(true, -logf(max_time_step));
}

bool extend(const partition_ptr_t& start_partition_ptr, partition_ptr_t& end_partition_ptr, float& proposal_density_ratio) {

    //partition original(start_partition);
    const clutter_t& clutter(start_partition_ptr->clutter());
    const track_collection& tracks(start_partition_ptr->tracks());
    if (tracks.size() == 0)
        return false;
    capability_recorder_ptr cap_rec_ptr = start_partition_ptr->get_capability_recorder();
    BOOST_ASSERT(cap_rec_ptr.get() not_eq 0); // there actually is something assigned
    //const track_collection::track_const_ptr_set& extendible_tracks(cap_rec_ptr->extendible_tracks());

    float log_track_prob = 1.f;
    shared_const_track_ptr track_to_extend;
    extend_weight_t extend_weight(start_partition_ptr->duration());

    if (not select_track(extend_weight, tracks, track_to_extend, log_track_prob)) {
        return false;
    }

    // front or back? p=0.5 doesn't need to be considered if(!) it is also ignored in reduce
    // if the following is false, then the track must be extendible at the end
    float log_front_or_end_prob = 0.f;
    bool extend_at_front = start_partition_ptr->first_time_stamp() < track_to_extend->first_time_stamp();
    BOOST_ASSERT(extend_at_front or start_partition_ptr->last_time_stamp() > track_to_extend->last_time_stamp());
    if (extend_at_front and start_partition_ptr->last_time_stamp() > track_to_extend->last_time_stamp()) {
        extend_at_front = sampling::uniform_int(0, 2) == 1;
        log_front_or_end_prob = -std::log(2.f);
    }

    time_stamp extension_last_time_stamp = track_to_extend->first_time_stamp();
    time_stamp extension_first_time_stamp = track_to_extend->last_time_stamp();
    float log_time_stamp_prob(1.f);
    float log_new_obs_prob(1.f);
    std::deque<observation> new_track_obs;
    if (not extend_at_front) {
        std::pair<bool, float> result = new_last_time_stamp_(track_to_extend, start_partition_ptr, extension_last_time_stamp);
        if (not result.first)
            return false; // nothing to extend at the end
        log_time_stamp_prob = result.second;
        log_new_obs_prob = extend_at_end_(track_to_extend, extension_last_time_stamp, clutter, new_track_obs);
    }
    else {
        if (start_partition_ptr->first_time_stamp() == track_to_extend->first_time_stamp())
            return false; // nothing to extend at the beginning
        extension_first_time_stamp = sampling::uniform_int(start_partition_ptr->first_time_stamp(), track_to_extend->first_time_stamp());
        log_time_stamp_prob = -logf(track_to_extend->first_time_stamp() - start_partition_ptr->first_time_stamp());
        log_new_obs_prob = extend_at_beginning_(track_to_extend, extension_first_time_stamp, clutter, new_track_obs);
    }
    const float new_dynamic_drag = 1.0f;
    boost::shared_ptr<track> track_extension(new track(extension_first_time_stamp, extension_last_time_stamp,
        new_track_obs.begin(), new_track_obs.end(), new_dynamic_drag));
    boost::shared_ptr<track> new_track(new track(*track_extension, *track_to_extend, track_to_extend->dynamic_drag()));


    // create copies of the partition clutter and tracks
    clutter_ptr new_clutter(new clutter_t(clutter));
    boost::shared_ptr<track_collection> new_tracks(new track_collection(tracks));

    cap_rec_ptr->add_insert_track(new_track);
    cap_rec_ptr->add_erase_track(track_to_extend);
    cap_rec_ptr->set_editing_finished();

    // replace old track
    new_tracks->replace(track_to_extend, new_track);

    // remove obs from clutter
    BOOST_FOREACH(const observation& o, new_track_obs) {
        new_clutter->erase(o);
    }

    // update partition
    end_partition_ptr = partition_ptr_t(
        new partition(start_partition_ptr->pool(), new_tracks, new_clutter, start_partition_ptr->expansion())
        );

    // reduce prob = prob(selecting correct track) * prob(sample the correct time point)
    //float reduce_log_prob = -logf(new_tracks->size());
    float reduce_log_prob = get_track_log_prob(reduce_weight_t(), *new_tracks, new_track);
    observation_collection::const_iterator obs_iter = new_track->begin();
    if (extend_at_front) {
        size_t num_obs = new_track->size();
        std::advance(obs_iter, num_obs -2); // the second last observation (we want to keep that)
        reduce_log_prob += -logf(t(*obs_iter) - new_track->first_time_stamp());
    } else {
        std::advance(obs_iter, 1); // the second observation (we want to keep that)
        reduce_log_prob += -logf(new_track->last_time_stamp() - t(*obs_iter)-1);
    }
    reduce_log_prob +=  -std::log(2.f);

    // extend prob == prob of sampling correct observation
    float extend_log_prob = log_track_prob + log_time_stamp_prob + log_new_obs_prob + log_front_or_end_prob;

    proposal_density_ratio = reduce_log_prob - extend_log_prob;
    return true;

}


}

}
