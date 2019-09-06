#include <algorithm>
#include <boost/assert.hpp>
#include <boost/foreach.hpp>
#include <deque>
#include <iterator>
#include <set>
#include <stdexcept>

#include "detail/physics.hpp"
#include "mh_moves/utility.hpp"
#include "track.hpp"
#include "track_collection.hpp"

namespace biggles
{

shared_const_track_ptr track_collection::insert(const track& t) {
    shared_const_track_ptr new_track_p(new track(t));
    tracks_.insert(new_track_p);
    return new_track_p;
}

shared_const_track_ptr track_collection::insert(const shared_const_track_ptr& t_p) {
    std::pair<std::set<shared_const_track_ptr>::iterator, bool> rv(tracks_.insert(t_p));
    if(!rv.second)
        throw std::runtime_error("track_collection::insert: track was already present in the collection.");
    return t_p;
}

size_t track_collection::unsafe_remove(const shared_const_track_ptr& h)
{
    return tracks_.erase(h);
}

void track_collection::remove(const shared_const_track_ptr& h)
{
    if(0 == unsafe_remove(h))
        throw std::runtime_error("track_collection::check_remove: track was not present in collection.");
}

const track& track_collection::find(shared_const_track_ptr t) const
{
    track_const_ptr_set::const_iterator ti(tracks_.find(t));
    if(ti == tracks_.end())
        throw std::runtime_error("track_collection::find: track handle is not in collection.");
    return **ti;
}

shared_const_track_ptr track_collection::merge(const shared_const_track_ptr& track_1,
                                                          const shared_const_track_ptr& track_2,
                                                          float dynamic_drag)
{
    // remove tracks 1 and 2 from the collection keeping a copy of their pointer.
    const shared_const_track_ptr& t1(pull_out(track_1));
    const shared_const_track_ptr& t2(pull_out(track_1));

    // insert a new merged version
    return insert(shared_const_track_ptr(new track(*t1, *t2, dynamic_drag)));
}

shared_const_track_ptr track_collection::pull_out(const shared_const_track_ptr& h)
{
    track_const_ptr_set::const_iterator ti(tracks_.find(h));
    if(ti == tracks_.end())
    {
        throw std::runtime_error("track_collection::pull_out: track is not in the collection.");
    }

    // keep a copy of the pointer
    shared_const_track_ptr tp(*ti);

    // erase it from the collection
    tracks_.erase(ti);

    return tp;
}


} // biggles
