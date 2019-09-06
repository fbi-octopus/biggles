/// @file track_collection.hpp A mutable collection of tracks with implicit sharing

#ifndef BIGGLES_TRACK_COLLECTION_HPP__
#define BIGGLES_TRACK_COLLECTION_HPP__

#include <boost/shared_ptr.hpp>
#include <set>

#include "detail/pair_collection.hpp"
#include "detail/track_pair_data.hpp"
#include "detail/transfer_possibility.hpp"
#include "detail/cross_over_possibility.hpp"
#include "track.hpp"

namespace biggles
{

/// @brief WIP: A collection of tracks.
///
/// Tracks are implicitly shared between collections.
///
/// FIXME: Is this a bottleneck for multi-threading?

class track_collection
{
public:
    /// @brief The internal container for the track collection.
    typedef std::set<shared_const_track_ptr> track_const_ptr_set;

    /// @brief Locate a track by its handle, throwing an exception if it isn't present.
    ///
    /// @param t A handle to the track.
    ///
    /// @return a const reference to the track associated with handle \p t.
    ///
    /// @throw std::runtime_error if \p t is not in the collection.
    const track& find(shared_const_track_ptr t) const;


    /// \brief checks if a track is in the collection
    bool contains(shared_const_track_ptr t) const { return tracks_.find(t) != tracks_.end(); }

    /// @brief Insert a new track into the collection.
    ///
    /// The track \p t is copied into the collection.
    ///
    /// @param t A reference to the track to copy into the collection.
    ///
    /// @return A handle for the newly inserted track and a flag indicating if the track was inserted.
    shared_const_track_ptr insert(const track& t);

    /// @brief Insert a new track into the collection from a \c boost::shared_ptr<>.
    ///
    /// This is a zero-copy version of insert() which inserts the track and increments the internal reference count of
    /// \p t_p.
    ///
    /// @param t_p A pointer to the track to insert.
    ///
    /// @return A handle to the newly inserted track.
    ///
    /// @throw std::runtime_error if the track pointed to by \p t_p is already in the collection.
    shared_const_track_ptr insert(const shared_const_track_ptr& t_p);

    /// @brief Remove a track from the collection, checking if the track is present.
    ///
    /// This is like remove() except that an exception is thrown if \p is not in the collection.
    ///
    /// @param track The handle of the track to remove.
    ///
    /// @throw std::runtime_error if \p track is not in the collection.
    void remove(const shared_const_track_ptr& track);

    /// @brief Replace a track in the collection with a new one.
    ///
    /// This is semantically identical to passing \p old_track to remove() and then passing \p new_track to insert()
    /// except that there may be a more efficient specialisation for particular types.
    ///
    /// @tparam Track Any type capable of being passed as the first argument to insert().
    /// @param old_track A handle to the track to be replaced.
    /// @param new_track A new track to take the place of \p old_track.
    ///
    /// @return A handle to the track which replaced \p old_track.
    template<typename Track>
    shared_const_track_ptr replace(const shared_const_track_ptr& old_track, Track new_track);

    /// @brief Replace a track in the collection with a new one, checking the old track is present.
    ///
    /// This is semantically identical to passing \p old_track to check_remove() and then passing \p new_track to
    /// insert() except that there may be a more efficient specialisation for particular types.
    ///
    /// @tparam Track Any type capable of being passed as the first argument to insert().
    /// @param old_track A handle to the track to be replaced.
    /// @param new_track A new track to take the place of \p old_track.
    ///
    /// @return A handle to the track which replaced \p old_track.
    ///
    /// @throw std::runtime_error if \p old_track is not in the collection.
    template<typename Track>
    shared_const_track_ptr check_replace(const shared_const_track_ptr& old_track, Track new_track);

    /// @brief Merge two tracks into one.
    ///
    /// @param track_1
    /// @param track_2
    /// @param dynamic_drag
    ///
    /// @return A handle to the new track formed from merging \p track_1 and \p track_2.
    shared_const_track_ptr merge(const shared_const_track_ptr& track_1, const shared_const_track_ptr& track_2, float dynamic_drag);

    /// @brief FIXME: UNIMPLEMENTED
    ///
    /// @param track
    /// @param at
    ///
    /// @return
    boost::tuple<shared_const_track_ptr, shared_const_track_ptr> split(shared_const_track_ptr track, size_t at);

    /// @brief Remove a track from the collection but return the shared pointer to it.
    ///
    /// This is like remove() but also returns a shared pointer to the removed track.
    ///
    /// @param track
    ///
    /// @throw std::runtime_error if \p track is not in the collection.
    shared_const_track_ptr pull_out(const shared_const_track_ptr& track);

    /// @name Iterating over tracks
    ///
    /// These member functions are provided to ease interaction with STL algorithms such as \p std::copy.
    ///
    /// @note There is no non-const iterator type because all mutation of the track set should be done by a member
    /// function. These iterators are provided for convenience only and should not be used to mutate the collection.
    /// @{

    /// @brief An const iterator type for iterating over all tracks.
    typedef track_const_ptr_set::const_iterator const_iterator;

    /// @brief An const iterator type for iterating over all tracks.
    typedef track_const_ptr_set::const_iterator iterator;

    /// @brief Query if there are no tracks in the collection.
    bool empty() const { return tracks_.empty(); }

    /// @brief Query the number of tracks in the collection.
    size_t size() const { return tracks_.size(); }

    /// @brief Return an iterator pointing to the first track.
    const_iterator begin() const { return tracks_.begin(); }

    /// @brief Return an iterator pointing just beyond the last track.
    const_iterator end() const { return tracks_.end(); }

    /// @}

    track_collection() { }

    /// @brief Copy constructor.
    ///
    /// @param c
    track_collection(const track_collection& c)
        : tracks_(c.tracks_)
    { }

protected:
    /// @brief The internal collection of track pointers.
    track_const_ptr_set tracks_;

    /// @brief Remove a track from the collection.
    ///
    /// @param track The handle of the track to remove.
    ///
    /// @return The number of tracks removed. This is 1 if \p track was removed and 0 if \p track could not be found
    /// within the collection.
    ///
    /// @sa remove()
    size_t unsafe_remove(const shared_const_track_ptr& track);

    /// @brief Locate a track by its handle.
    ///
    /// This function may be significantly faster than check_find() since it doesn't need to assert that the track is in
    /// the collection.
    ///
    /// @note This member function is unsafe in that no checks are made that \p t is in the collection.
    ///
    /// @sa find()
    ///
    /// @param t A handle to the track to find.
    ///
    /// @return a const reference to the track associated with handle \p t.
    ///
    /// @throw std::runtime_error if \p t is not in the collection.
    const track& unsafe_find(shared_const_track_ptr t) const { return *t; }

private:

    // To make a track collection non-assignable.
    const track_collection& operator = (const track_collection&);
};

template <class OBS_CONTAINER>
void ingest_observations(OBS_CONTAINER& container, const observation_collection& oc) {
    container.insert(container.end(), oc.begin(), oc.end());
}

template <class OBS_CONTAINER>
void ingest_observations(OBS_CONTAINER& container, const track_collection& tc) {
    BOOST_FOREACH(const shared_const_track_ptr& track_p, tc) {
        ingest_observations(container, track_p->observations());
    }
}


}

#define WITHIN_BIGGLES_TRACK_COLLECTION_HPP__
#include "track_collection.tcc"
#undef WITHIN_BIGGLES_TRACK_COLLECTION_HPP__

#endif // BIGGLES_TRACK_COLLECTION_HPP__
