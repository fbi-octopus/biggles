#include <boost/assert.hpp>
#include <cmath>
#include <limits>
#include <utility>
#include <algorithm>

#include "mh_moves/mh_moves.hpp"
#include "sampling/simple.hpp"
#include "mh_moves/utility.hpp"
#include "tools/debug.hpp"

namespace biggles {

namespace creation_fun {

    const size_t max_gap = 3;

    typedef std::deque<observation> obs_container;

    bool get_track_obs(const observation& first_obs, const clutter_t& clutter, const size_t max_gap,
        obs_container& track_obs)
    {
        BOOST_ASSERT(max_gap > 1);
        bool done = false;
        track_obs.clear();
        track_obs.push_back(first_obs);
        observation last_obs = first_obs;
        time_stamp final_ts = clutter.last_time_stamp();
        while (not done) {
            bool added_obs = false;
            for (size_t gap = 1; gap <= max_gap; ++ gap) {
                time_stamp current_ts = t(last_obs) + gap;
                if (current_ts == final_ts) { break; }
                float radius = mh_moves::light_years(t(last_obs), current_ts);
                obs_container candidate_obs = clutter.locate_near(current_ts, last_obs, radius);
                if (candidate_obs.empty()) continue;
                last_obs = *sampling::from_range(candidate_obs.begin(), candidate_obs.end());
                track_obs.push_back(last_obs);
                added_obs = true;
            }
            done = not added_obs;
        }
        return track_obs.size() > 1;
    }

    bool get_track_at_ts(const time_stamp ts, const clutter_t& clutter, const size_t max_gap, obs_container& track_obs) {
        obs_container firsts = clutter.observations(ts);
        bool done = false;
        while (firsts.size() > 0 and not done) {
            obs_container::iterator obs_it = sampling::from_range(firsts.begin(), firsts.end());
            firsts.erase(obs_it);
            done = get_track_obs(*obs_it, clutter, max_gap, track_obs);
        }
        return done;
    }


    bool create_track(const partition_ptr_t& start_partition_ptr, partition_ptr_t& end_partition_ptr) {
        const clutter_t& clutter(start_partition_ptr->clutter());
        const track_collection& tracks(start_partition_ptr->tracks());

        time_stamp first_ts = clutter.first_time_stamp();
        time_stamp last_ts = clutter.last_time_stamp();
        obs_container track_obs;
        bool success = false;
        for (time_stamp ts = first_ts; ts < last_ts - 1; ts++) {
            success = get_track_at_ts(ts, clutter, max_gap, track_obs);
            if (success) break;
        }
        if (not success) {
            end_partition_ptr = start_partition_ptr;
            return success;
        }

        // make new track
        shared_track_ptr new_track(new track(t(track_obs.front()), t(track_obs.back())+1,
                                track_obs.begin(), track_obs.end(), 1.f));
        // make new clutter
        clutter_ptr new_clutter(new clutter_t(clutter));
        for (obs_container::iterator obs_it = track_obs.begin(); obs_it != track_obs.end(); obs_it++) {
            new_clutter->erase(*obs_it);
        }
        // make new partition
        mh_moves::shared_track_collection_ptr new_tracks(new track_collection(tracks));
        new_tracks->insert(new_track);
        end_partition_ptr = partition_ptr_t(
            new partition(start_partition_ptr->pool(), new_tracks, new_clutter, start_partition_ptr->expansion())
        );
        return success;
    }

    bool create_track_reverse(const partition_ptr_t& start_partition_ptr, partition_ptr_t& end_partition_ptr) {
        const clutter_t& clutter(start_partition_ptr->clutter());
        const track_collection& tracks(start_partition_ptr->tracks());

        time_stamp first_ts = clutter.first_time_stamp();
        time_stamp last_ts = clutter.last_time_stamp();
        obs_container track_obs;
        bool success = false;
        for (time_stamp ts = first_ts; ts < last_ts - 1; ts++) {
            success = get_track_at_ts(ts, clutter, max_gap, track_obs);
            if (success) break;
        }
        if (not success) {
            end_partition_ptr = start_partition_ptr;
            return success;
        }

        // make new track
        shared_track_ptr new_track(new track(t(track_obs.front()), t(track_obs.back())+1,
                                track_obs.begin(), track_obs.end(), 1.f));
        // make new clutter
        clutter_ptr new_clutter(new clutter_t(clutter));
        for (obs_container::iterator obs_it = track_obs.begin(); obs_it != track_obs.end(); obs_it++) {
            new_clutter->erase(*obs_it);
        }
        // make new partition
        mh_moves::shared_track_collection_ptr new_tracks(new track_collection(tracks));
        new_tracks->insert(new_track);
        end_partition_ptr = partition_ptr_t(
            new partition(start_partition_ptr->pool(), new_tracks, new_clutter, start_partition_ptr->expansion())
        );
        return success;
    }

} // creation_fun


bool max_creation(const partition_ptr_t& start_partition_ptr, partition_ptr_t& end_partition_ptr) {

    partition_ptr_t input_partition = start_partition_ptr;

    float total_obs = input_partition->observation_count();

    bool done = false;
    while (not done) {
        done = not creation_fun::create_track(input_partition, end_partition_ptr)
               or float(end_partition_ptr->clutter().size())/total_obs <= 0.5;
        input_partition = end_partition_ptr;
    }
    return true;
} // max_creation

bool rev_creation(const partition_ptr_t& start_partition_ptr, partition_ptr_t& end_partition_ptr) {

    partition_ptr_t input_partition = start_partition_ptr;

    float total_obs = input_partition->observation_count();

    bool done = false;
    while (not done) {
        done = not creation_fun::create_track_reverse(input_partition, end_partition_ptr)
               or float(end_partition_ptr->clutter().size())/total_obs <= 0.5;
        input_partition = end_partition_ptr;
    }
    return true;
} // rev_creation



} // biggles
