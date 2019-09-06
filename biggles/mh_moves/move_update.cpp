#include <boost/assert.hpp>
#include <cmath>
#include <limits>
#include <utility>
#include <algorithm>

#include "mh_moves/mh_moves.hpp"
#include "mh_moves/mh_move_weights.hpp"
#include "sampling/simple.hpp"
#include "mh_moves/utility.hpp"
#include "tools/debug.hpp"

namespace biggles {

namespace mh_moves
{

namespace update_fun {
    void common_observations(
            std::deque<observation>& clutter_candidates_bef,
            std::deque<observation>& clutter_candidates_aft,
            std::deque<observation>& clutter_candidates )
    {
        std::sort(clutter_candidates_bef.begin(), clutter_candidates_bef.end());
        std::sort(clutter_candidates_aft.begin(), clutter_candidates_aft.end());
        clutter_candidates.resize(std::min(clutter_candidates_bef.size(), clutter_candidates_aft.size()));
        std::deque<observation>::iterator end_iter = std::set_intersection(
            clutter_candidates_bef.begin(), clutter_candidates_bef.end(),
            clutter_candidates_aft.begin(), clutter_candidates_aft.end(),
            clutter_candidates.begin());
        clutter_candidates.resize(end_iter-clutter_candidates.begin());
    }

    void get_candidates_at(
                const observation_collection& track_to_update_obs,
                const time_stamp ts_to_update,
                observation_collection::const_iterator before,
                observation_collection::const_iterator after,
                const clutter_t& clutter,
                std::deque<observation>& clutter_candidates
            )
    {
        std::deque<observation> clutter_candidates_bef;
        std::deque<observation> clutter_candidates_aft;
        bool have_bef = before != track_to_update_obs.end();
        bool have_aft = after  != track_to_update_obs.end();

        if (have_bef) {
            float rad_bef = std::abs(float(t(*before)) - float(ts_to_update)) * detail::speed_of_light_;
            clutter_candidates_bef = clutter.locate_near(ts_to_update, *before, rad_bef);
        }
        if (have_aft) {
            float rad_aft = std::abs(float(t(*after )) - float(ts_to_update)) * detail::speed_of_light_;
            clutter_candidates_aft = clutter.locate_near(ts_to_update, *after , rad_aft);
        }
        if (have_bef and have_aft) {
            std::sort(clutter_candidates_bef.begin(), clutter_candidates_bef.end());
            std::sort(clutter_candidates_aft.begin(), clutter_candidates_aft.end());
            clutter_candidates.resize(std::min(clutter_candidates_bef.size(), clutter_candidates_bef.size()));
            std::deque<observation>::iterator end_iter = std::set_intersection(
                clutter_candidates_bef.begin(), clutter_candidates_bef.end(),
                clutter_candidates_aft.begin(), clutter_candidates_aft.end(),
                clutter_candidates.begin());
            clutter_candidates.resize(end_iter-clutter_candidates.begin());
        } else if (have_bef) {
            clutter_candidates = clutter_candidates_bef;
        } else {
            clutter_candidates = clutter_candidates_aft;
        }

    }

    observation get_centre(const track& target_track, const time_stamp ts) {
        const observation_collection& track_obs(target_track.observations());
        if (ts < t(track_obs.front())) // "==" must be treated seperately
            return track_obs.front();
        if (ts == t(track_obs.front())) {
            observation_collection::const_iterator it = track_obs.begin();
            it++;
            return *it;
        }
        if (ts > t(track_obs.back()))
            return track_obs.back();
        if (ts == t(track_obs.back())) {
            observation_collection::const_iterator it = track_obs.end();
            it--; it--;
            return *it;
        }
        observation_collection::const_iterator later(track_obs.begin());
        observation_collection::const_iterator sooner(track_obs.begin());
        later++; // second obs
        while (t(*later) < ts) {
            later++;
            sooner++;
        }
        // now :: t(*sooner) < ts ; t(later) >= ts
        if (t(*later) == ts)
            later++; // not end() since ts < tmax
        // now t(*later) > ts
        time_stamp bef = t(*sooner);
        time_stamp aft = t(*later);
        BOOST_ASSERT(bef < ts);
        BOOST_ASSERT(aft > ts);
        if ((aft - ts) >= (ts - bef))
            return convex_combination(*sooner, *later, bef);
        else
            return convex_combination(*sooner, *later, aft);
    }

    void get_candidates_for_track(
        const track& target_track,
        const clutter_t& clutter,
        std::map<time_stamp, std::deque<observation> >& candidates_map,
        std::map<time_stamp, observation>& drops
        )
    {
        typedef observation_collection::observation_multimap observations_map;
        const observations_map& track_obs(target_track.observations().get());
        observations_map::const_iterator leading(track_obs.begin());
        observations_map::const_iterator current(track_obs.begin());
        observations_map::const_iterator following(leading);
        following++;
        BOOST_ASSERT(following not_eq track_obs.end()); // there should be at least two observations
        bool initial(true);
        //OK(track_obs.size());
        for (time_stamp ts = target_track.first_time_stamp(); ts != target_track.last_time_stamp(); ++ts) {
            std::deque<observation> candidates;
            while (current != track_obs.end() and time_of(current) < ts)
                ++current;
            bool drop = current != track_obs.end() and ts == time_of(current);
            if (drop)
                drops.insert(std::make_pair(ts, observation_of(current)));
            while (following != track_obs.end() and time_of(following) < ts) {
                ++following;
                ++leading;
            }
            //BOOST_ASSERT(drop and time_of(following) == ts);
            observations_map::const_iterator trailing(following);
            if (following != track_obs.end() and time_of(following) == ts) ++trailing; // make sure that trailing points after ts
            BOOST_ASSERT(trailing == track_obs.end() or time_of(trailing) > ts);
            if (leading not_eq track_obs.end() and ts < time_of(leading)) { //
                candidates = clutter.locate_near(ts, observation_of(leading), light_years(time_of(leading), ts));
                if (drop or not candidates.empty())
                    candidates_map.insert(std::make_pair(ts, candidates));
                continue;
            }
            if (initial and ts == time_of(leading)) {
                initial = false;
                candidates = clutter.locate_near(ts, observation_of(following), light_years(time_of(following), ts));
                if (drop or not candidates.empty())
                    candidates_map.insert(std::make_pair(ts, candidates));
                continue;
            }
            BOOST_ASSERT(time_of(leading) < ts);
            if (trailing == track_obs.end()) {
                candidates = clutter.locate_near(ts, observation_of(leading), light_years(ts, time_of(leading)));
                if (drop or not candidates.empty())
                    candidates_map.insert(std::make_pair(ts, candidates));
                continue;
            }
            std::deque<observation> clutter_candidates_bef =
                clutter.locate_near(ts, observation_of(leading), light_years(ts, time_of(leading)));
            std::deque<observation> clutter_candidates_aft =
                clutter.locate_near(ts, observation_of(trailing), light_years(ts, time_of(trailing)));
            common_observations(clutter_candidates_bef, clutter_candidates_aft, candidates);
            if (drop or not candidates.empty())
                candidates_map.insert(std::make_pair(ts, candidates));
        }
    }
} // namespace biggles::mh_moves::update_fun

bool update(const partition_ptr_t& start_partition_ptr, partition_ptr_t& end_partition_ptr, float& proposal_density_ratio)
{
    // early-out: there must be at least one track
    if(start_partition_ptr->tracks().empty()) {
        return false;
    }
    proposal_density_ratio = 0.f;

    capability_recorder_ptr cap_rec_ptr = start_partition_ptr->get_capability_recorder();
    BOOST_ASSERT(cap_rec_ptr.get() not_eq 0); // there actually is something assigned

    const clutter_t& clutter(start_partition_ptr->clutter());
    const track_collection& tracks(start_partition_ptr->tracks());

    // choose a track to update
    track_collection::const_iterator track_it(sampling::from_range(tracks.begin(), tracks.end()));

    // check we found one
    BOOST_ASSERT(track_it != tracks.end());
    BOOST_ASSERT((*track_it)->size() >= 2);
    const shared_const_track_ptr& track_to_update(*track_it);

    std::map<time_stamp, std::deque<observation> > candidates_map;
    std::map<time_stamp, observation > possible_drops;
    update_fun::get_candidates_for_track(*track_to_update, clutter, candidates_map, possible_drops);
    if (candidates_map.empty())
        return false;
    std::map<time_stamp, std::deque<observation> >::const_iterator
        ts_and_cand_to_update(sampling::from_range(candidates_map.begin(), candidates_map.end()));
    time_stamp ts_to_update = ts_and_cand_to_update->first;
    std::map<time_stamp, observation >::const_iterator droppings = possible_drops.find(ts_to_update);
    bool observation_dropped = droppings not_eq possible_drops.end();
    /* ********** */
    std::deque<observation>& clutter_candidates(candidates_map.at(ts_to_update));
    size_t n_options = clutter_candidates.size();
    //BOOST_ASSERT(n_options > 0);
    observation_collection new_track_obs(track_to_update->observations());
    clutter_ptr new_clutter(new clutter_t(clutter));
    if (observation_dropped) {
        new_track_obs.erase(droppings->second);
        new_clutter->insert(droppings->second);
    }
    weights_t weights(n_options, 1.f);
    observation centre = update_fun::get_centre(*track_to_update, ts_to_update);
    //OK(ts_to_update);
    //OK(t(centre));
    std_normal_weight obs_weight(centre);
    std::transform(clutter_candidates.begin(), clutter_candidates.end(), weights.begin(), obs_weight);

    float forw_obs_weight = 0.f;
    float back_obs_weight = 0.f;
    float noobs_weight = obs_weight.sigmas(1, ts_to_update);
    //OK(noobs_weight);

    if (observation_dropped) {
        forw_obs_weight += obs_weight(droppings->second);
        clutter_candidates.push_back(observation());
        weights.push_back(noobs_weight); // leaving a gap has the weight of 1 sigma
    } else {
        forw_obs_weight += noobs_weight;
    }
    if (weights.empty())
        return false;
    float forward_log_prob = 1.f;
    const observation& sample(sampling::weighted_choice(
        clutter_candidates.begin(), weights.begin(), weights.end(), forward_log_prob));
    if (not noobs(sample)) {
        new_clutter->erase(sample);
        new_track_obs.insert(sample);
        back_obs_weight += obs_weight(sample);
    } else {
        if(new_track_obs.size() < 2)
            return false;
        back_obs_weight += noobs_weight;
    }
    float total_weight = std::accumulate(weights.begin(), weights.end(), 0.f);
    float total_back_weight = total_weight - back_obs_weight + forw_obs_weight;

    float back_log_prob = 1.f;
    if (observation_dropped) {
        back_log_prob = obs_weight(droppings->second);
    } else {
        back_log_prob = noobs_weight;
    }
    back_log_prob = logf(back_log_prob/total_back_weight);


    BOOST_ASSERT(new_track_obs != track_to_update->observations()); // something must have happend

    const float new_dynamic_drag = 1.f;
    shared_track_ptr new_track(new track(track_to_update->first_time_stamp(), track_to_update->last_time_stamp(),
                                new_track_obs.begin(), new_track_obs.end(), new_dynamic_drag));

    // sanity check the new track
    BOOST_ASSERT(new_track->size() >= 2);
    BOOST_ASSERT(verify_track(*new_track));

    cap_rec_ptr->add_erase_track(track_to_update);
    cap_rec_ptr->add_insert_track(new_track);
    cap_rec_ptr->set_editing_finished();
    shared_track_collection_ptr new_tracks(new track_collection(tracks));
    new_tracks->remove(track_to_update);
    new_tracks->insert(new_track);

    proposal_density_ratio = back_log_prob - forward_log_prob;

    end_partition_ptr = partition_ptr_t(
        new partition(start_partition_ptr->pool(), new_tracks, new_clutter, start_partition_ptr->expansion())
        );

    return true;
} // bool update(...)

} // namespace mh_moves

} // namespace biggles
