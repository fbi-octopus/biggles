#ifndef BIGGLES_TRACK_COLL_TRANSFORM_HPP__
#define BIGGLES_TRACK_COLL_TRANSFORM_HPP__

#include "../track_collection.hpp"
#include "../detail/observation_reservoir.hpp"
#include <deque>

namespace biggles {

    struct indexed_track {
        time_stamp first_time_stamp;
        time_stamp last_time_stamp;
        std::deque<size_t> observations;
    };

    typedef std::deque<indexed_track> track_coll_index_t;

    void track_coll_observation_to_index(const observation_reservoir_ptr& obs_res_p, const track_collection& obs_track_coll,
        track_coll_index_t& idx_track_coll)
    {
        BOOST_FOREACH(const shared_const_track_ptr& t_ptr, obs_track_coll) {
            idx_track_coll.push_back(indexed_track());
            idx_track_coll.back().observations.resize(t_ptr->size());
            obs_res_p->idx(t_ptr->observations().begin(), t_ptr->observations().end(),
                idx_track_coll.back().observations.begin());
            idx_track_coll.back().first_time_stamp = t_ptr->first_time_stamp();
            idx_track_coll.back().last_time_stamp = t_ptr->last_time_stamp();
        }
    }

    void track_coll_index_to_observation(const observation_reservoir_ptr& obs_res_p, const track_coll_index_t& idx_track_coll,
        track_collection& obs_track_coll)
    {
        BOOST_FOREACH(const indexed_track& idx_t, idx_track_coll) {
            std::deque<observation> obs_list(idx_t.observations.size());
            obs_res_p->obs(idx_t.observations.begin(), idx_t.observations.end(), obs_list.begin());
            shared_track_ptr tp(new track(idx_t.first_time_stamp, idx_t.last_time_stamp, obs_list.begin(), obs_list.end()));
            obs_track_coll.insert(tp);
        }
    }

} //namespace biggles

#endif // BIGGLES_TRACK_COLL_TRANSFORM_HPP__
