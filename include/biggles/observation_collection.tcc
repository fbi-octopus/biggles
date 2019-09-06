#ifndef WITHIN_BIGGLES_OBSERVATION_COLLECTION_HPP__
#error "this file must only be included by observation_collection.hpp"
#endif // WITHIN_BIGGLES_OBSERVATION_COLLECTION_HPP__

#include <boost/assert.hpp>
#include <boost/foreach.hpp>
#include <cmath>
#include <Eigen/Dense>
#include <utility>

namespace biggles
{

inline observation_collection::const_range observation_collection::at_time_stamp(time_stamp time_stamp) const
{
    std::pair<observation_multimap::const_iterator, observation_multimap::const_iterator>
        obs_its = observations_.equal_range(time_stamp);
    return std::make_pair(const_iterator(obs_its.first), const_iterator(obs_its.second));
}

inline observation_collection::observation_multimap::const_iterator observation_collection::real_find(const observation& o) const
{
    // UGH!
    return const_cast<observation_collection*>(this)->real_find(o);
}

inline observation_collection::observation_multimap::iterator observation_collection::real_find(const observation& o)
{
    // look for observations with the right timestamp
    std::pair<observation_multimap::iterator, observation_multimap::iterator> obs_it_range(observations_.equal_range(t(o)));

    // no observations? if so, return end()
    if(obs_it_range.first == obs_it_range.second)
        return observations_.end();

    // return first observation which matches
    for(observation_multimap::iterator i(obs_it_range.first);
        i != obs_it_range.second;
        ++i)
    {
        if(i->second == o)
            return i;
    }

    // no observation was found
    return observations_.end();
}

inline bool observation_collection::erase(const observation& o)
{
    observation_multimap::iterator it(real_find(o));
    if(it == observations_.end())
        return false;
    observations_.erase(it);
    return true;
}

inline bool observation_is_within_light_cone(
    const observation& o,
    const observation& centre,
    float speed_of_light,
    time_stamp first_time_stamp,
    time_stamp last_time_stamp)
{
    // observation is outside of temporal bounds?
    if((t(o) < first_time_stamp) || (t(o) >= last_time_stamp))
        return false;

    // observation is outside of spatial bounds?
    return observations_see_each_other(o, centre, speed_of_light);
    /*
    float dx(x(o) - x(centre)), dy(y(o) - y(centre));
    float ds_sq(dx*dx + dy*dy);
    time_stamp dt(abs(t(o) - t(centre)));
    if(ds_sq >= dt*dt*speed_of_light*speed_of_light)
        return false;

    return true;
    */
}

template<typename OutputIterator>
void observations_within_light_cone(const observation_collection& observations,
                                    const observation& centre,
                                    float speed_of_light,
                                    time_stamp first_time_stamp,
                                    time_stamp last_time_stamp,
                                    std::set<time_stamp>& time_stamps_with_obs,
                                    OutputIterator output)
{
    // sanity check input
    BOOST_ASSERT(last_time_stamp >= first_time_stamp);

    // early out if we've been asked for a time cone with zero extent
    if(first_time_stamp == last_time_stamp)
        return;

    // the first observation whose time stamp is >= first_time_stamp
    observation_collection::const_iterator first(
        observations.lower_bound_for_time_stamp(first_time_stamp));

    // the end of the observation collection
    observation_collection::const_iterator last(
        observations.lower_bound_for_time_stamp(last_time_stamp));

    // remove any elements from the set of time stamps where candidates have been found
    time_stamps_with_obs.clear();

    // loop over all observations in range temporally
    for(; first != last; ++first)
    {
        const observation& cand(*first);

        if((t(cand) < first_time_stamp) || (t(cand) >= last_time_stamp))
            continue;

        BOOST_ASSERT(t(cand) >= first_time_stamp);
        BOOST_ASSERT(t(cand) < last_time_stamp);
        if(observation_is_within_light_cone(cand, centre, speed_of_light, first_time_stamp, last_time_stamp))
        {
            *output = cand;
            time_stamps_with_obs.insert(t(cand));
            ++output;
        }
    }
}

inline bool observation_collection_does_allow_tracks(const observation_collection& observations, float speed_of_light) {
    observation_collection::const_iterator current_iter, loop_iter;
    bool success = false;
    for (current_iter = observations.begin(); current_iter != observations.end(); ++current_iter) {
        loop_iter = current_iter;
        while (++loop_iter != observations.end()) {
            // observations should be at different time points
            if (t(*loop_iter) == t(*current_iter)) {
                continue;
            }
            // if the observation at loop_iter is within the light cone of the observation at current_iter
            // declare success
            if (observations_see_each_other(*loop_iter, *current_iter, speed_of_light))
            {
                success = true;
                break;
            }
        }
        // if we are already successful we can stop searching
        if (success) break;
    }
    return success;
}


}
