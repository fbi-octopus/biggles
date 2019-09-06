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

namespace death_fun {
    float get_birth_obs_prob(const clutter_ptr& clutter, time_stamp first_time_stamp, time_stamp last_time_stamp,
        const shared_const_track_ptr& track_to_remove)
    {
        float noobs_weight(0.1f);
        float sample_weight(0.f);
        float total_weight(0.f);

        float new_obs_prob(0.f);
        observation last_obs;
        const observation_collection &track_to_remove_obs(track_to_remove->observations());
        observation_collection::const_iterator obs_it = track_to_remove_obs.begin();
        for (time_stamp new_obs_t = first_time_stamp; new_obs_t < last_time_stamp; ++new_obs_t)
        {

            std::deque<observation> candidates;
            if (noobs(last_obs)) {
                candidates = clutter->observations(new_obs_t);
                noobs_weight = 0.1f;
                total_weight = float(candidates.size()) + noobs_weight;
                if (new_obs_t == t(*obs_it)) {
                    last_obs = *obs_it++;
                    sample_weight = 1.f;
                } else {
                    sample_weight = noobs_weight;
                }
            } else {
                float radius = float(new_obs_t - t(last_obs)) * detail::speed_of_light_;
                candidates = clutter->locate_near(new_obs_t, last_obs, radius);
                std_normal_weight obs_weight(last_obs);
                noobs_weight = obs_weight.sigmas(2, new_obs_t);
                total_weight = fun_sum(candidates.begin(), candidates.end(), obs_weight, noobs_weight);
                if (obs_it != track_to_remove_obs.end() and new_obs_t == t(*obs_it)) {
                    last_obs = *obs_it++;
                    sample_weight = obs_weight(last_obs);
                } else {
                    sample_weight = noobs_weight;
                }
            }

            new_obs_prob += logf(sample_weight) - logf(total_weight);
        }

        return new_obs_prob;
    }
}

bool death(const partition_ptr_t& start_partition_ptr, partition_ptr_t& end_partition_ptr, float& proposal_density_ratio)
{
    // copy the start partition so that start_partition can equal end_partition and we can modify end_partition
    //partition original(start_partition);
    const clutter_t& clutter(start_partition_ptr->clutter());
    const track_collection& tracks(start_partition_ptr->tracks());

    capability_recorder_ptr cap_rec_ptr = start_partition_ptr->get_capability_recorder();
    BOOST_ASSERT(cap_rec_ptr.get() not_eq 0); // there actually is something assigned


    // early-out: there must be at least one track
    if(tracks.empty())
        return false;

    // choose a track to remove
    size_t n_to_sample = sampling::uniform_int(0, tracks.size());
    track_collection::const_iterator it(tracks.begin()); std::advance(it, n_to_sample);
    shared_const_track_ptr track_to_remove(*it);

    // create copies of the partition clutter and tracks
    clutter_ptr new_clutter(new clutter_t(clutter));
    boost::shared_ptr<track_collection> new_tracks(new track_collection(tracks));

    // put track obs into clutter
    BOOST_FOREACH(const observation& o, *track_to_remove)
    {
        new_clutter->insert(o);
    }

    cap_rec_ptr->add_erase_track(track_to_remove);
    cap_rec_ptr->set_editing_finished();

    // remove old track
    new_tracks->remove(track_to_remove);

    // update partition
    end_partition_ptr = partition_ptr_t(
        new partition(start_partition_ptr->pool(), new_tracks, new_clutter, start_partition_ptr->expansion())
    );
    end_partition_ptr->set_first_time_stamp(start_partition_ptr->first_time_stamp());
    end_partition_ptr->set_last_time_stamp(start_partition_ptr->last_time_stamp());

    // calculating the probs to sample the first and the last time stamp
    // see birth move for more details
    time_stamp last_clutter_time_stamp = new_clutter->last_time_stamp();
    time_stamp last_for_first_time_stamp = last_clutter_time_stamp - 1;
    time_stamp first_time_stamp = track_to_remove->first_time_stamp();
    float first_time_stamp_prob = -logf(last_for_first_time_stamp - end_partition_ptr->first_time_stamp());
    BOOST_ASSERT(std::isfinite(first_time_stamp_prob));
    time_stamp last_time_stamp = track_to_remove->last_time_stamp();
    float last_time_stamp_prob = -logf(end_partition_ptr->last_time_stamp() - first_time_stamp - 1);
    BOOST_ASSERT(std::isfinite(last_time_stamp_prob));

    float new_obs_prob =
        death_fun::get_birth_obs_prob(new_clutter, first_time_stamp, last_time_stamp, track_to_remove);

    // death prob == prob selecting correct track
    float death_log_prob = -logf(start_partition_ptr->tracks().size());

    // birth prob == prob of sampling correct observation
    float birth_log_prob = first_time_stamp_prob + last_time_stamp_prob + new_obs_prob;

    proposal_density_ratio = birth_log_prob - death_log_prob;
    return true;
}

}

}
