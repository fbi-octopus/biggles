#include "kalman_filter_cache.hpp"

#include <utility>

namespace biggles
{

const kalman_filter& kalman_filter_cache::get(
    boost::shared_ptr<const track> track_ptr,
    const matrix2f &R, const matrix4f &Q)
{
    boost::weak_ptr<const track> track_weak_ptr(track_ptr);

    // Do we have an entry for this ptr? If not, make one and return it.
    cache_map::const_iterator entry_it(cache_.find(track_weak_ptr));
    if(entry_it == cache_.end())
        return make_entry(track_weak_ptr, R, Q);

    // Is the ptr still valid?
    if(entry_it->first.expired())
        return make_entry(track_weak_ptr, R, Q);

    // Does R11, R12 and R22 match?
    if((entry_it->second.get<1>() != R(0, 0)) ||
       (entry_it->second.get<2>() != R(0, 1)) ||
       (entry_it->second.get<3>() != R(1, 1)))
        return make_entry(track_weak_ptr, R, Q);

    // If we get here, the entry is valid
    return entry_it->second.get<0>();
}

const kalman_filter& kalman_filter_cache::make_entry(
    boost::weak_ptr<const track> track_weak_ptr,
    const matrix2f &R, const matrix4f &Q)
{
    // Remove the old entry if there was one
    cache_.erase(track_weak_ptr);

    // Create a new entry
    std::pair<cache_map::iterator, bool> rv = cache_.insert(
        std::make_pair(track_weak_ptr,
                       boost::make_tuple( kalman_filter(), R(0, 0), R(0, 1), R(1, 1) )
                       )
        );
    BOOST_ASSERT(rv.second);

    // Initialise Kalman filter from track (avoids a copy)
    boost::shared_ptr<const track> track_ptr(track_weak_ptr.lock());
    rv.first->second.get<0>().reinitialise(
        track_ptr->first_time_stamp(),
        track_ptr->last_time_stamp(),
        track_ptr->observations().begin(),
        track_ptr->observations().end(),
        R, Q, track_ptr->dynamic_drag());

    ++purge_timer_;
    if(purge_timer_ > 1024)
    {
        std::deque<boost::weak_ptr<const track> > to_purge;
        BOOST_FOREACH(const cache_map::value_type& kv_pair, cache_)
        {
            if(kv_pair.first.expired())
                to_purge.push_back(kv_pair.first);
        }

        BOOST_FOREACH(const boost::weak_ptr<const track>& t, to_purge)
        {
            cache_.erase(t);
        }

        purge_timer_ = 0;
    }

    // Return the filter we made
    return rv.first->second.get<0>();
}

}
