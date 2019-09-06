#ifndef BIGGLES_CLUTTER_TRANSFORM_HPP__
#define BIGGLES_CLUTTER_TRANSFORM_HPP__

#include "../clutter.hpp"
#include "../detail/observation_reservoir.hpp"
#include <deque>

namespace biggles {

    typedef std::deque<size_t> clutter_index_t;

    void clutter_observation_to_index(const observation_reservoir_ptr& obs_res_p, const clutter_t& obs_clutter,
        clutter_index_t& idx_clutter)
    {
        idx_clutter.resize(obs_clutter.size());
        clutter_index_t::iterator idx_clutter_it = idx_clutter.begin();
        for (clutter_t::const_iterator obs_clutter_it = obs_clutter.begin(); obs_clutter_it != obs_clutter.end(); ++obs_clutter_it) {
            idx_clutter_it = obs_res_p->idx(obs_clutter_it->second.begin(), obs_clutter_it->second.end(), idx_clutter_it);
        }
    }

    void clutter_index_to_observation(const observation_reservoir_ptr& obs_res_p, const clutter_index_t& idx_clutter,
        clutter_t& obs_clutter)
    {
        std::deque<observation> obs_list(idx_clutter.size());
        obs_res_p->obs(idx_clutter.begin(), idx_clutter.end(), obs_list.begin());
        obs_clutter.insert(obs_list.begin(), obs_list.end());
    }


}

#endif // BIGGLES_CLUTTER_TRANSFORM_HPP__
