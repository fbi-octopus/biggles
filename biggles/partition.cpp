#include <algorithm>
#include <boost/assert.hpp>
#include <boost/foreach.hpp>
#include <limits>

#include "partition.hpp"
#include "tools/debug.hpp"
namespace biggles
{

void partition::update_boundary_cache()
{
    // special case for no tracks and no clutter
    if(clutter().empty() && tracks().empty())
    {
        first_time_stamp_ = last_time_stamp_ = 0;
        return;
    }

    first_time_stamp_ = std::numeric_limits<time_stamp>::max();
    last_time_stamp_ = std::numeric_limits<time_stamp>::min();
    /*

    // scan the clutter

    BOOST_FOREACH(const observation& o, clutter())
    {
        first_time_stamp_ = std::min(first_time_stamp_, t(o));
        last_time_stamp_ = std::max(last_time_stamp_, t(o) + 1);
    }
    */

    if (not clutter().empty()) {
        first_time_stamp_ = clutter().first_time_stamp();
        last_time_stamp_ = clutter().last_time_stamp();
    }

    // scan the tracks
    BOOST_FOREACH(const boost::shared_ptr<const track>& tp, tracks())
    {
        if(tp->duration() > 0)
        {
            first_time_stamp_ = std::min(first_time_stamp_, tp->first_time_stamp());
            last_time_stamp_ = std::max(last_time_stamp_, tp->last_time_stamp());
        }
    }

    BOOST_ASSERT(last_time_stamp_ >= first_time_stamp_);
}

void partition::update_observation_count() {
    observation_count_ = clutter_->size();
    BOOST_FOREACH(const boost::shared_ptr<const track>& tp, tracks()) {
        observation_count_ += tp->observations().size();
    }
    BOOST_ASSERT(observation_count_ == observation_pool_->size());
}

partition::partition(const track_collection_ptr& tracks, const const_clutter_ptr& clutter)
        : tracks_(tracks)
        , clutter_(clutter)
        , expansion_(expansion_from_observations(*tracks, *clutter))
{
    set_observation_pool_from_tracks_clutter();
    update_boundary_cache();
    update_observation_count();
}

partition::partition(const observation_pool_ptr& pool, const track_collection_ptr& tracks,
    const const_clutter_ptr& clutter)
    : tracks_(tracks)
    , clutter_(clutter)
    , observation_pool_(pool)
    , expansion_(expansion_from_observations(*tracks, *clutter))
{
    update_boundary_cache();
    update_observation_count();
}



partition::partition(const track_collection_ptr& tracks,
            const const_clutter_ptr& clutter, const expansion_2d& expansion)
    : tracks_(tracks)
    , clutter_(clutter)
    , expansion_(expansion)
{
    set_observation_pool_from_tracks_clutter();
    update_boundary_cache();
    update_observation_count();
}



partition::expansion_2d expansion_from_observations(const track_collection& tracks,
                        const clutter_t& clutter)
{
    partition::expansion_2d exp(detail::myriad, -detail::myriad, detail::myriad, -detail::myriad);

    // scan the clutter
    if (clutter.size() >0) {
        exp.grow(clutter.lower_bound());
        exp.grow(clutter.upper_bound());
    }
    /*
    BOOST_FOREACH(const observation& o, clutter) {
        exp.grow(point_2d(x(o), y(o)));
    }
    */

    BOOST_FOREACH(const boost::shared_ptr<const track>& tp, tracks)
    {
        if(tp->duration() > 0) {
            exp.grow(tp->min_location(), tp->max_location());
        }
    }
    return exp;
}
/*
std::string partition::string_from_track(const shared_const_track_ptr& track_p) const {
    std::stringstream trackstr;
    trackstr << track_p->first_time_stamp() << "," << track_p->last_time_stamp() << "%";
    observation_collection::const_iterator it = track_p->observations().begin();
    observation_collection::const_iterator end_it = track_p->observations().end();
    trackstr << pool()->idx(*it++);
    while (it != end_it)
        trackstr << "," << pool()->idx(*it++);
    return trackstr.str();
}
*/

std::string partition::string_from_track(const shared_const_track_ptr& track_p) const {
    std::set<size_t> obs_idx;
    BOOST_FOREACH(const observation& obs, track_p->observations())
        obs_idx.insert(pool()->idx(obs));
    std::stringstream trackstr;
    trackstr << track_p->first_time_stamp() << "," << track_p->last_time_stamp() << "%";
    std::set<size_t>::iterator it = obs_idx.begin();
    size_t first = *it;
    size_t current = *it;
    it++;
    while (it != obs_idx.end()) {
        if ( current + 1 == *it)
            current++;
        else {
            if (current == first)
                trackstr << current << ",";
            else
                trackstr << first << "-" << current << ",";
            first = *it;
            current = *it;
        }
        it++;
    }
    if (current == first)
        trackstr << current;
    else
        trackstr << first << "-" << current;
    return trackstr.str();
}

std::string partition::as_string() const {
    std::set<std::string> track_strings;
    BOOST_FOREACH(const boost::shared_ptr<const track>& track_p, tracks()) {
        track_strings.insert(string_from_track(track_p));
    }
    if (track_strings.size() == 0)
        return std::string();
    std::set<std::string>::iterator it = track_strings.begin();
    std::string result(*it++);
    while (it != track_strings.end())
        result += "/" + *it++;
    return result;
}

bool verify_partition(partition& p) {
    BOOST_FOREACH(const boost::shared_ptr<const track>& track_p, p.tracks())
    {
        if(!verify_track(*track_p))
            return false;
    }

    return true;
}



}
