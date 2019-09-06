#include <boost/assert.hpp>
#include <cmath>
#include <limits>
#include <utility>
#include <iostream>
#include "mh_moves/mh_moves.hpp"
#include "sampling/simple.hpp"
#include "mh_moves/utility.hpp"
#include "mh_moves/mh_move_weights.hpp"
#include "detail/physics.hpp"

namespace biggles {

namespace mh_moves
{

namespace birth_fun {
    float get_new_track_obs(const clutter_t& clutter, time_stamp first_time_stamp,
        time_stamp last_time_stamp, std::deque<observation>& new_track_obs)
    {
        new_track_obs.clear(); // this should be not necessary
        float new_obs_prob(0.0f);
        for (time_stamp new_obs_t = first_time_stamp; new_obs_t < last_time_stamp; ++new_obs_t)
        {
            std::deque<observation> candidates;
            weights_t weights;
            if (new_track_obs.empty()) {
                candidates = clutter.observations(new_obs_t);
                weights.assign(candidates.size(), 1.f);
                float noobs_weight = 0.1f;
                weights.push_back(noobs_weight);
                candidates.push_back(observation());
            } else {
                observation last_obs = new_track_obs.back();
                candidates = clutter.locate_near(new_obs_t ,last_obs, light_years(new_obs_t, t(last_obs)));
                std_normal_weight obs_weight(last_obs);
                weights.assign(candidates.size(), 1.f);
                std::transform(candidates.begin(), candidates.end(), weights.begin(), obs_weight);
                weights.push_back(obs_weight.sigmas(2, new_obs_t));
                candidates.push_back(observation());
            }

            float obs_prob = 1.f;
            observation new_obs =
                sampling::weighted_choice(candidates.begin(), weights.begin(), weights.end(), obs_prob);

            if (not noobs(new_obs)) {
                new_track_obs.push_back(new_obs);
            }
            new_obs_prob += obs_prob;

        }
        return new_obs_prob;
    }
}

/** \brief Birth move
 *
 * 1) Sample a first time stamp
 *    The time stamp is sampled so that it is guaranteed that there is at least one later observation
 * 2) Sample a last time stamp
 *    This is done so that it is guaranteed that the last time stamp is at least two time stamps away from the first
 * 3) Loop over [first time stamp, last time stamp)
 * 3.1) Find all observations, o, at the current form clutter that fulfil:
 *      - if the track so far is empty all observations are accepted
 *      - if the track so far already has observations, o is in the light cone of the last observation
 * 3.2)
 *
 */
bool birth(const partition_ptr_t& start_partition_ptr, partition_ptr_t& end_partition_ptr, float& proposal_density_ratio)
{
    // copy the start partition so that start_partition can equal end_partition and we can modify end_partition
    //partition original(start_partition);
    const clutter_t& clutter(start_partition_ptr->clutter());
    const track_collection& tracks(start_partition_ptr->tracks());

    capability_recorder_ptr cap_rec_ptr = start_partition_ptr->get_capability_recorder();
    BOOST_ASSERT(cap_rec_ptr.get() not_eq 0); // there actually is something assigned

    // early-out: there must be at least at two time stamps clutter observations
    if(clutter.last_time_stamp()-clutter.first_time_stamp() < 2) {
        return false;
    }

    time_stamp last_clutter_time_stamp = clutter.last_time_stamp();

    time_stamp last_for_first_time_stamp = last_clutter_time_stamp - 1;
    time_stamp first_time_stamp = sampling::uniform_int(start_partition_ptr->first_time_stamp(), last_for_first_time_stamp);

    float first_time_stamp_prob = -logf(last_for_first_time_stamp - start_partition_ptr->first_time_stamp());

    // the last time stamp for the new track is at least two time stamps away
    time_stamp last_time_stamp = sampling::uniform_int(first_time_stamp + 2, start_partition_ptr->last_time_stamp() + 1);
    BOOST_ASSERT(last_time_stamp - first_time_stamp > 1);

    float last_time_stamp_prob = -logf(start_partition_ptr->last_time_stamp() - first_time_stamp - 1);

    std::deque<observation> new_track_obs;

    float new_obs_prob = birth_fun::get_new_track_obs(clutter, first_time_stamp, last_time_stamp, new_track_obs);

    if(new_track_obs.size() < 2) {
        return false;
    }


    // create copies of the partition clutter and tracks
    clutter_ptr new_clutter(new clutter_t(clutter));
    boost::shared_ptr<track_collection> new_tracks(new track_collection(tracks));

    // create new track with a dynamic drag of 1
    boost::shared_ptr<track> new_track(new track(first_time_stamp, last_time_stamp,
        new_track_obs.begin(), new_track_obs.end(), 1.f));

    cap_rec_ptr->add_insert_track(new_track);
    cap_rec_ptr->set_editing_finished();
    // insert new track
    new_tracks->insert(new_track);

    // remove obs from clutter
    BOOST_FOREACH(const observation& o, new_track_obs)
    {
        new_clutter->erase(o);
    }

    // update partition
    end_partition_ptr = partition_ptr_t(
        new partition(start_partition_ptr->pool(), new_tracks, new_clutter, start_partition_ptr->expansion()));

    // death prob == prob selecting correct track
    float death_log_prob = -logf(new_tracks->size());

    // birth prob == prob of sampling correct observation
    float birth_log_prob = first_time_stamp_prob + last_time_stamp_prob + new_obs_prob;

    proposal_density_ratio = death_log_prob - birth_log_prob;
    return true;
}


}

}
