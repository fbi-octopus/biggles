#ifndef TRACK_PAIR_DATA_HPP__
#define TRACK_PAIR_DATA_HPP__

#include <list>
#include "../track.hpp"

namespace biggles { namespace detail {

    /// \brief list item to record moveable (Trf) track pairs. A pair can appear twice: (t1, t2) != (t2, t1)
    template<class VALUE_TYPE>
    struct track_pair_item_t {
        shared_const_track_ptr first; ///< \brief first track
        shared_const_track_ptr second; ///< \brief second track
        typedef VALUE_TYPE value_type;
        VALUE_TYPE value; ///< move information
        /// \brief the constructor allows any relation between t1 and t2
        track_pair_item_t(const shared_const_track_ptr &track1, const shared_const_track_ptr &track2,
            const VALUE_TYPE &val) : first(track1), second(track2), value(val) {}
        bool operator<(const track_pair_item_t<VALUE_TYPE>& o) const {
            return first<o.first or (first==o.first and second<o.second);
        }
        bool operator==(const track_pair_item_t<VALUE_TYPE>& o) const {
            return first==o.first and second==o.second;
        }
    };

    /// \brief list item to record moveable (Mrg, Crs) track pairs. A pair can only appear once: (t1, t2) == (t2, t1).
    template<class VALUE_TYPE>
    struct sym_track_pair_item_t {
        shared_const_track_ptr first; ///< \brief first track
        shared_const_track_ptr second; ///< \brief second track
        typedef VALUE_TYPE value_type;
        VALUE_TYPE value; ///< move information
        /// \brief The constructor ensures t1 <= t2
        sym_track_pair_item_t(const shared_const_track_ptr &track1, const shared_const_track_ptr &track2,
            const VALUE_TYPE &val) : first(track1< track2 ? track1 : track2),
                                    second(track1< track2 ? track2 : track1), value(val) {}
        bool operator<(const sym_track_pair_item_t<VALUE_TYPE>& o) const {
            return first<o.first or (first==o.first and second<o.second);
        }
        bool operator==(const sym_track_pair_item_t<VALUE_TYPE>& o) const {
            return first==o.first and second==o.second;
        }
    };

    template<class ITEM_TYPE>
    class ContainsTrack : public std::unary_function<ITEM_TYPE, bool>{
        const shared_const_track_ptr& track_;
    public:
        ContainsTrack(const shared_const_track_ptr& track_p) : track_(track_p) {}
        bool operator()(const ITEM_TYPE &item) {
            return item.first == track_ or item.second == track_;
        }
    };

    template<class ITEM_TYPE>
    class new_track_pair_data_t {
    public:
        typedef ITEM_TYPE item_type;
        typedef typename std::list<item_type> track_pair_data_list;
        typedef typename track_pair_data_list::const_iterator const_iterator;
        typedef typename track_pair_data_list::iterator iterator;
        typedef typename item_type::value_type value_type;
    private:
        track_pair_data_list track_pair_data_;
        new_track_pair_data_t(const new_track_pair_data_t&); ///< \brief forbidden
        void operator=(const new_track_pair_data_t&); ///< \brief forbidden
    public:
        new_track_pair_data_t() {}  ///< \brief the only constructor
        size_t size() const { return track_pair_data_.size(); }
        bool empty() const { return track_pair_data_.empty(); }
        const_iterator begin() const { return track_pair_data_.begin(); }
        const_iterator end() const { return track_pair_data_.end(); }
        iterator end() { return track_pair_data_.end(); }
        iterator begin() { return track_pair_data_.begin(); }
        /// \brief insert a new track pair and its data
        ///
        /// do we need to test if the track pair is already in the list?
        void insert(
            const shared_const_track_ptr& track_1, const shared_const_track_ptr& track_2, const value_type& data_ptr)
        {
            item_type item(track_1, track_2, data_ptr);
            insert(item);
        }

        void insert(const item_type& item) {
            iterator low = std::lower_bound (begin(), end(), item);
            if (low == end() or not(*low == item))
                    track_pair_data_.insert(low, item);
        }

        /// \brief remove all entries that contain a particular track
        void erase_items_with_track(const shared_const_track_ptr& track_p) {
            track_pair_data_.remove_if(ContainsTrack<item_type>(track_p));
        }

        size_t count_items_with_track(const shared_const_track_ptr& track_p) {
            return std::count_if(track_pair_data_.begin(), track_pair_data_.end(), ContainsTrack<item_type>(track_p));
        }


    };




    /// \brief A mimic of a map pair<track, track> -> DATA. Used to keep record of tracks that can be used in MH moves
    ///
    /// This class stores track pairs and related data that was gathered when the pair was tested for the
    /// possibility if a MH move could be successfully executed.
    ///
    /// Could be used with the moves: cross_over, transfer, merge.
    /// DATA is intended to be a class that tests if the move can be executed on the track pair.
    /// DATA contains some information useful for the respective MH move, which is gathered during the test.
    template <class DATA>
    class track_pair_data_t {
    public:
        typedef typename std::map<shared_const_track_ptr_pair, DATA> track_pair_data_map;
        typedef typename track_pair_data_map::const_iterator const_iterator;
        typedef typename track_pair_data_map::iterator iterator;
        typedef typename track_pair_data_map::value_type value_type;
        typedef typename track_pair_data_map::key_type key_type;
        typedef typename track_pair_data_map::mapped_type mapped_type;
    private:
        track_pair_data_map track_pair_data_;
    public:
        track_pair_data_t() { } ///< \brief default constructor
        /// \brief copy constructor
        track_pair_data_t(const track_pair_data_t& o) : track_pair_data_(o.track_pair_data_) { }
        /// \brief range constructor
        track_pair_data_t(const_iterator b, const_iterator e) : track_pair_data_(b, e) { }

        size_t size() const { return track_pair_data_.size(); }
        bool empty() const { return track_pair_data_.empty(); }
        const_iterator begin() const { return track_pair_data_.begin(); }
        const_iterator end() const { return track_pair_data_.end(); }

        /// \brief insert a new track pair and its data
        std::pair<iterator, bool>
        insert(const shared_const_track_ptr_pair& track_pair, const DATA& data_ptr) {
            return track_pair_data_.insert(std::make_pair(track_pair, data_ptr));
        }

        /// \brief remove all entries that contain a particular track
        void erase_items_with_track(const shared_const_track_ptr& track_p) {
            for(;;) {
                // This erases a item pair with track_p as first or second element in the key
                iterator i(track_pair_data_.begin());
                for(; (i != track_pair_data_.end()) and (i->first.first != track_p) and (i->first.second != track_p); ++i);
                if(i != track_pair_data_.end()) {
                    track_pair_data_.erase(i);
                } else {
                    break;
                }
            }
        }

    };

} } // biggles::detail

#endif // TRACK_PAIR_DATA_HPP__
