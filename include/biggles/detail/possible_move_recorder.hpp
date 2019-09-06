/// @file track_collection.hpp A mutable collection of tracks with implicit sharing

#ifndef BIGGLES_POSSIBLE_MOVE_RECORDER_HPP__
#define BIGGLES_POSSIBLE_MOVE_RECORDER_HPP__

#include <boost/shared_ptr.hpp>
#include <set>

#include "pair_collection.hpp"
#include "track_pair_data.hpp"
#include "transfer_possibility.hpp"
#include "cross_over_possibility.hpp"
#include "flip_possibility.hpp"
#include "../track_collection.hpp"

namespace biggles
{

/// \brief Shared pointer to a map-like container with all pair on which a transfer move could be conducted
typedef detail::new_track_pair_data_t< detail::track_pair_item_t<detail::transfer_possibility_ptr> > transferables_t;
/// \brief Shared pointer to a map-like container with all pair on which a cross_over move could be conducted
typedef detail::new_track_pair_data_t< detail::sym_track_pair_item_t<detail::cross_over_possibility_ptr> > crossables_t;
/// \brief Shared pointer to a map-like container with all pair on which a merge move could be conducted
typedef detail::new_track_pair_data_t< detail::sym_track_pair_item_t<int> > mergeables_t;
/// \brief Shared pointer to a map-like container with all pair on which a flip move could be conducted
typedef detail::new_track_pair_data_t< detail::sym_track_pair_item_t<detail::flip_possibility_ptr> > flippables_t;


typedef std::set<shared_const_track_ptr> extendibles_t;

class capability_recorder {
    /// @name Records
    /// @{
    /// @brief A mapping between tracks and tracks which are candidates for merging with them.
    mergeables_t mergeable_pairs_;

    /// @brief A mapping between tracks and tracks which are candidates for observation *transfer*.
    transferables_t transfer_pairs_;

    /// @brief A mapping between tracks and tracks which are candidates for cross_over.
    crossables_t cross_over_pairs_;

    /// @brief A mapping between tracks and tracks which are candidates for flip move.
    flippables_t flip_pairs_;

    /// @brief the set of pointers to the extedible tracks
    extendibles_t extensible_tracks_;
    /// @}


    /// @name Track buffers
    ///
    /// tracks will kept here until it is decided if the records need updating
    /// @{
    std::deque<shared_const_track_ptr> erase_buffer_;
    std::deque<shared_const_track_ptr> insert_buffer_;
    bool editing_finished_;
    //typedef std::deque<shared_const_track_ptr>::iterator buffer_iterator;
    /// @}

    // non-copyable
    capability_recorder(const capability_recorder&);
    void operator=(const capability_recorder&);

    void insert_track_pair_into_transferables(const shared_const_track_ptr& t1, const shared_const_track_ptr& t2);
    void insert_track_pair_into_merge_pairs(const shared_const_track_ptr& t1, const shared_const_track_ptr& t2);
    void insert_track_pair_into_crossables(const shared_const_track_ptr& t1, const shared_const_track_ptr& t2);
    void insert_track_pair_into_flippables(const shared_const_track_ptr& t1, const shared_const_track_ptr& t2);

    /// @brief Remove the pair of tracks from the merge-able pair list.
    ///
    /// @param t
    void remove_track_from_merge_pairs(const shared_const_track_ptr& t);

    /// @brief Insert the pair of tracks from the merge-able pair list.
    ///
    /// @param t
    void insert_track_into_merge_pairs(const shared_const_track_ptr& t,
        const track_collection& partition_tracks);

    void remove_track_from_transferables(const shared_const_track_ptr& track_p);
    void insert_track_into_transferables(const shared_const_track_ptr& track_p,
        const track_collection& partition_tracks);

    void remove_track_from_crossables(const shared_const_track_ptr& track_p);
    void remove_track_from_flippables(const shared_const_track_ptr& track_p);
    void insert_track_into_crossables(const shared_const_track_ptr& track_p,
        const track_collection& partition_tracks);
    /// \brief counts how many inserts would be made
    size_t count_inserts_into_crossables(const shared_const_track_ptr& new_p,
        const track_collection &tracks);
    void insert_track_into_flippables(const shared_const_track_ptr& track_p,
        const track_collection& partition_tracks);

    void remove_from_records(const shared_const_track_ptr& t) {
        //remove_track_from_merge_pairs(t);
        remove_track_from_transferables(t);
        remove_track_from_crossables(t);
        remove_track_from_flippables(t);
        //extensible_tracks_.erase(t);
    }
    void insert_into_records(const shared_const_track_ptr& t, const track_collection& tracks) {
        //insert_track_into_merge_pairs(t, tracks);
        insert_track_into_transferables(t, tracks);
        insert_track_into_crossables(t, tracks);
        insert_track_into_flippables(t, tracks);
        //if (t->first_time_stamp() > begin_ts_ or t->last_time_stamp() < end_ts_)
        //    extensible_tracks_.insert(t);
    }

    void initial_setup(const track_collection& partition_tracks);

    /// @name Partition size
    /// @{
    time_stamp begin_ts_;
    time_stamp end_ts_;
    /// @}

    /// \brief clear the buffers of possible changes
    void clear_buffers() {
        erase_buffer_.clear();
        insert_buffer_.clear();
        editing_finished_ = false;
    }

public:
    capability_recorder(time_stamp begin_ts, time_stamp end_ts, const track_collection& tracks) :
        editing_finished_(false), begin_ts_(begin_ts), end_ts_(end_ts)
    {
        initial_setup(tracks);
    }

    /// @name buffer handling
    ///
    /// The buffer keeps the tracks that would change if the proposal is accepted.
    /// @{

    void set_editing_finished() { editing_finished_ = true; }

    const std::deque<shared_const_track_ptr>& erase_buffer() const { return erase_buffer_; }
    const std::deque<shared_const_track_ptr>& insert_buffer() const { return insert_buffer_; }

    void add_insert_track(const shared_const_track_ptr &track_ptr) {
        if (editing_finished_)
            clear_buffers();
        insert_buffer_.push_back(track_ptr);
    }
    void add_erase_track(const shared_const_track_ptr &track_ptr) {
        if (editing_finished_)
            clear_buffers();
        erase_buffer_.push_back(track_ptr);
    }
    /// \brief the proposal has been accepted commit the changes to the records
    void commit_changes(const track_collection& tracks) {
        BOOST_FOREACH(const shared_const_track_ptr &track_ptr, erase_buffer_)  {
            remove_from_records(track_ptr);
        }
        BOOST_FOREACH(const shared_const_track_ptr &track_ptr, insert_buffer_) {
            insert_into_records(track_ptr, tracks);
        }
        clear_buffers();
    }
    int num_crossables_if_commited(const track_collection& tracks) {
        typedef detail::sym_track_pair_item_t<detail::cross_over_possibility_ptr> crossable_track_pair_item;
        std::list< crossable_track_pair_item > crossables(cross_over_pairs_.begin(), cross_over_pairs_.end());
        BOOST_FOREACH(const shared_const_track_ptr &track_ptr, erase_buffer_)  {
            crossables.remove_if(detail::ContainsTrack<crossable_track_pair_item>(track_ptr));
        }
        BOOST_FOREACH(const shared_const_track_ptr &tp1, insert_buffer_) {
            if(tp1->empty()) continue;

            BOOST_FOREACH(const shared_const_track_ptr& tp2, tracks) {
                if(tp2 == tp1)
                    continue;
                detail::cross_over_possibility_ptr possibility_test(new detail::cross_over_possibility_test(tp1, tp2));

                if(possibility_test->passed()) {
                    crossable_track_pair_item item(tp1, tp2, possibility_test);
                    std::list< crossable_track_pair_item >::iterator low =
                        std::lower_bound (crossables.begin(), crossables.end(), item);
                    if (low == crossables.end() or not(*low == item))
                            crossables.insert(low, item);
                }
            }
        }
        return crossables.size();
    }
    ///@}


    /// @brief Return a const reference to a collection of merge-able pairs within this collection.
    ///
    /// A pair of tracks is deemed merge-able if their terminal observations are within each other's light cones and
    /// within the maximum temporal separation distance.
    const mergeables_t& mergeable_pairs() const { return mergeable_pairs_; }

    /// @brief Return a const reference to a collection of track pairs that can transfer observations
    /// within this collection.
    const transferables_t& transfer_pairs() const { return transfer_pairs_; }

    /// @brief Return a const reference to a collection of pairs that can cross-over.
    const crossables_t& cross_over_pairs() const { return cross_over_pairs_; }

    /// @brief Return a const reference to a collection of track pairs that can be flipped.
    const flippables_t& flip_pairs() const { return flip_pairs_; }

    const extendibles_t& extendible_tracks() const { return extensible_tracks_; }

    /// \brief checks if all tracks in the records are in the track collection
    bool verify_against_track(const track_collection& tracks) const;


    size_t size() const {
        return _size_mergeable_pairs_() + _size_cross_over_pairs_() + _size_extendible_tracks_()
            + _size_transfer_pairs_();
    }

    /// @name debug functions
    /// size of internal record keeping containers.
    /// @{
    size_t _size_mergeable_pairs_() const {
        //return 0;
        return mergeable_pairs_.size() ;
    }
    size_t _size_transfer_pairs_() const { return transfer_pairs_.size() ; }
    size_t _size_cross_over_pairs_() const { return cross_over_pairs_.size() ; }
    size_t _size_extendible_tracks_() const { return extensible_tracks_.size(); }
    /// @}

};

typedef boost::shared_ptr<capability_recorder> capability_recorder_ptr;

}

#endif
