#ifndef WITHIN_BIGGLES_TRACK_COLLECTION_HPP__
#error "This file should only be included by biggles/track_collection.hpp"
#endif

#include <stdexcept>

namespace biggles
{

template<typename Track>
shared_const_track_ptr track_collection::replace(const shared_const_track_ptr& old_track, Track new_track)
{
    remove(old_track);
    return insert(new_track);
}

template<typename Track>
shared_const_track_ptr track_collection::check_replace(const shared_const_track_ptr& old_track, Track new_track)
{
    remove(old_track);
    return insert(new_track);
}

}
