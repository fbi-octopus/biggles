#ifndef WITHIN_BIGGLES_TRACK_HPP__
#error "This file should only be included by biggles/track.hpp"
#endif

#include <boost/assert.hpp>
#include <boost/foreach.hpp>

namespace biggles
{

template<typename ObservationIndexInputIterator, typename ObervationRandomAccessIterator>
track::track(time_stamp first_time_stamp,
             time_stamp last_time_stamp,
             ObservationIndexInputIterator observation_indices_first,
             ObservationIndexInputIterator observation_indices_last,
             ObervationRandomAccessIterator observation_first,
             ObervationRandomAccessIterator observation_last,
             float dynamic_drag)
    : first_time_stamp_(first_time_stamp)
    , last_time_stamp_(last_time_stamp)
    , dynamic_drag_(dynamic_drag)
{
    BOOST_ASSERT(last_time_stamp >= first_time_stamp);

    // copy observations into track
    size_t n_obs(observation_last - observation_first);
    BOOST_FOREACH(const size_t& i, std::make_pair(observation_indices_first, observation_indices_last))
    {
        BOOST_ASSERT(i < n_obs);
        observations_.insert(observation_first[i]);
    }

    // check for observations with common timestamps
    check_common_timestamps();

    // update minimum and maximum locations
    update_boundary_cache();
}

inline bool within_bounds(const track& t, const point_2d& p)
{
    if(t.empty())
        return false;
    return (x(p) >= x(t.min_location())) && (x(p) <= x(t.max_location())) &&
           (y(p) >= y(t.min_location())) && (y(p) <= y(t.max_location()));
}

inline bool within_bounds(const track& t, const point_2d& p, time_stamp ts)
{
    if(t.empty())
        return false;
    return within_bounds(t,p) && (ts >= t.first_time_stamp()) && (ts <= t.last_time_stamp());
}

inline bool within_bounds(const track& track, const observation& o)
{
    return within_bounds(track, point_2d(x(o), y(o)), t(o));
}

inline bool track::operator != (const track& t) const
{
    // ordered in increasing cost to compute
    return (t.first_time_stamp_ != first_time_stamp_) ||
           (t.last_time_stamp_ != last_time_stamp_) ||
           (t.dynamic_drag_ != dynamic_drag_) ||
           //(t.boundary_min_location_ != boundary_min_location_) ||
           //(t.boundary_max_location_ != boundary_max_location_) ||
           (t.observations_ != observations_);
}

inline void track::insert(const observation& o)
{
    observations_.insert(o);

    // update time stamp boundaries
    first_time_stamp_ = std::min(first_time_stamp_, t(o));
    last_time_stamp_ = std::max(last_time_stamp_, t(o)+1);

    // update spatial boundaries
    x(boundary_min_location_) = std::min(x(boundary_min_location_), x(o));
    x(boundary_max_location_) = std::max(x(boundary_max_location_), x(o));

    y(boundary_min_location_) = std::min(y(boundary_min_location_), y(o));
    y(boundary_max_location_) = std::max(y(boundary_max_location_), y(o));

    // check we didn't add an observation at an existing time stamp
    BOOST_ASSERT(observations_.count_at_time_stamp(t(o)) == 1);
}

inline void track::time_stamps_from_observations() {
    first_time_stamp_ = observations_.first_time_stamp();
    last_time_stamp_  = observations_.last_time_stamp();
}

}
