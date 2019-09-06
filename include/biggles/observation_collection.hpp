/// @file observation_collection.hpp Collection type for observations

#ifndef BIGGLES_OBSERVATION_COLLECTION_HPP__
#define BIGGLES_OBSERVATION_COLLECTION_HPP__

#include <algorithm>
#include <boost/assert.hpp>
#include <iterator>
#include <map>
#include <set>

#include "detail/pair_iterator_wrapper.hpp"
#include "observation.hpp"

namespace biggles
{

/// @brief A collection of observations.
///
/// The member functions begin() and end() may be used to get objects conforming to the \c ForwardIterator concept and
/// which are guaranteed to return the observations in ascending time stamp order.
///
/// This class differs from a usual STL container in that an effort is made to make the compiler catch if you try to
/// directly modify an observation within it. This allows the collection class to 'under the covers' provide a number of
/// optimisations which allow for spatio-temporal queries to be performed in logarithmic time. In addition, this allows
/// the class to maintain a weak ordering of observations by time stamp. This allows for convenient access to
/// observations in a narrow time stamp range. This was done in an attempt to make Biggles' underlying data structures
/// scale at worst logarithmically with number of frames.
///
/// @note Observations which share a common time stamp are allowed by this class. If this is not desired, you will have
/// to roll your own consistency checking code.
class observation_collection
{
public:
    /// @brief A multimap which stores observations sorted by time stamp.
    typedef std::map<time_stamp, observation> observation_multimap;

    /// @brief An iterator over const observation values.
    typedef detail::pair_iterator_wrapper<observation_multimap::const_iterator, const observation> const_iterator;

    /// @brief An iterator over observation values.
    ///
    /// @note unlike most STL containers, the underlying observation in this iterator is still const.
    typedef detail::pair_iterator_wrapper<observation_multimap::const_iterator, const observation> iterator;

    /// @brief The value type of this collection.
    typedef observation value_type;

    /// @brief A reference to a value stored in this collection.
    typedef observation& reference;

    /// @brief A const reference to a value stored in this collection.
    typedef const observation& const_reference;

    /// @brief A pair of const_iterator-s which define a range of observations.
    typedef std::pair<const_iterator, const_iterator> const_range;

    /// @brief Query if the collection is empty.
    ///
    /// @return \c true iff there are no observations in the collection.
    bool empty() const { return observations_.empty(); }

    /// @brief The number of observations in the collection.
    size_t size() const { return observations_.size(); }

    const observation_multimap& get() const { return observations_; }

    /// @brief Insert an observation into the collection.
    ///
    /// @param o The observation to copy and insert into the collection.
    ///
    /// @return an iterator pointing to the newly inserted observation
    iterator insert(const observation& o)
        { return iterator(observations_.insert(std::make_pair(t(o), o)).first); }

    /// @brief Insert an observation into the collection.
    ///
    /// @note The \p position parameter is currently ignored. This member function exists to keep interface equivalence
    /// with set-like STL containers.
    ///
    /// @param position An iterator hinting at the appropriate position for the observation.
    /// @param o The observation to copy and insert into the collection.
    ///
    /// @return an iterator pointing to the newly inserted observation
    iterator insert(iterator position, const observation& o)
        { return iterator(observations_.insert(std::make_pair(t(o), o)).first); }

    /// @brief Insert multiple observations into the collection.
    ///
    /// Observations are copied from the interval [\p first, \p last) into the collection.
    ///
    /// @tparam InputIterator
    /// @param first
    /// @param last
    template<typename InputIterator>
    void insert(InputIterator first, InputIterator last)
        { std::copy(first, last, std::inserter(*this, begin())); }

    /// @brief An ForwardIterator pointing to the first observation in the collection.
    const_iterator begin() const { return const_iterator(observations_.begin()); }

    /// @brief An ForwardIterator pointing just beyond the last observation in the collection.
    const_iterator end() const { return const_iterator(observations_.end()); }

    /// @brief An ForwardIterator pointing to the first observation in the collection.
    iterator begin() { return iterator(observations_.begin()); }

    /// @brief An ForwardIterator pointing just beyond the last observation in the collection.
    iterator end() { return iterator(observations_.end()); }

    /// @brief The number of observations at a given time stamp.
    ///
    /// @param time_stamp The time stamp to query.
    size_t count_at_time_stamp(size_t time_stamp) const { return observations_.count(time_stamp); }

    /// @brief Get a pair of iterators delimiting a range of observations with a given time stamp.
    ///
    /// @param time_stamp
    ///
    /// @return A pair of iterators, (first, last), so that the observations on the interval (first, last] all have time
    /// stamp \p time_stamp.
    const_range at_time_stamp(time_stamp time_stamp) const;

    /// @brief Erase all observations from the collection.
    void clear() { observations_.clear(); }

    /// @brief Create an empty collection.
    observation_collection() { }

    /// @brief Copy an existing collection.
    ///
    /// @param c
    observation_collection(const observation_collection& c) : observations_(c.observations_) { }

    /// @brief Construct from a range of observation iterators.
    ///
    /// Observations are copied from the interval [\p first, \p last).
    ///
    /// @sa insert<>(InputIterator, InputIterator)
    ///
    /// @tparam InputIterator
    /// @param first
    /// @param last
    template<typename InputIterator>
    observation_collection(InputIterator first, InputIterator last) { insert(first, last); }

    /// @brief Return an iterator to the lower bound for a time stamp
    ///
    /// The iterator returned from this member function points to the first observations whose time stamp is not less
    /// than \p t. That is to say, the first observation with a time stamp equal <em>or greater</em> than \p t.
    ///
    /// @param t A time stamp
    const_iterator lower_bound_for_time_stamp(time_stamp t) const
        { return const_iterator(observations_.lower_bound(t)); }

    /// @brief Return an iterator to the upper bound for a time stamp
    ///
    /// The iterator returned from this member function points to the first observations whose time stamp is greater
    /// than \p t. Unlike lower_bound(), this member function does not look of observations whose time stamp
    /// <em>equals</em> \p t.
    ///
    /// @param t A time stamp
    const_iterator upper_bound_for_time_stamp(time_stamp t) const
        { return const_iterator(observations_.upper_bound(t)); }

    /// @brief The first time stamp in the collection.
    ///
    /// Combined with first_time_stamp(), this may be used to construct the interval [\c first_time_stamp, \c
    /// last_time_stamp) on which time stamps of observations in this collection lie. If \c first_time_stamp = \c
    /// last_time_stamp then the collection contains no observations.
    ///
    /// @note If the collection is empty, this returns 0.
    time_stamp first_time_stamp() const { return empty() ? 0 : t(*begin()); }

    /// @brief Just beyond last time stamp in the collection.
    ///
    /// Combined with first_time_stamp(), this may be used to construct the interval [\c first_time_stamp, \c
    /// last_time_stamp) on which time stamps of observations in this collection lie. If \c first_time_stamp = \c
    /// last_time_stamp then the collection contains no observations.
    ///
    /// @note If the collection is empty, this returns 0.
    time_stamp last_time_stamp() const { return empty() ? 0 : (t(*(--end())) + 1); }

    /// @brief the difference between last time stamp and first time stamp
    time_stamp duration() const { return last_time_stamp() - first_time_stamp(); }

    /// @brief Return a const reference to the earliest observation in the collection.
    ///
    /// @note Calling this on an empty collection is undefined behaviour.
    const observation& front() const { BOOST_ASSERT(!empty()); return *begin(); }

    /// @brief Return a const reference to the latest observation in the collection.
    ///
    /// @note Calling this on an empty collection is undefined behaviour.
    const observation& back() const { BOOST_ASSERT(!empty()); return *(--end()); }

    /// @brief Find an observation within the collection.
    ///
    /// Find the first observation within the collection which compares equal to \p o and return an iterator pointing to
    /// it. If no observation was found, return end().
    ///
    /// @param o
    const_iterator find(const observation& o) const
        { observation_multimap::const_iterator it(real_find(o)); return (it == observations_.end()) ? end() : const_iterator(it); }

    /// @brief Query if an observation is within the collection.
    ///
    /// This is a convenience member function equivalent to <code>find(o) != end()</code>.
    ///
    /// @param o
    bool contains(const observation& o) const { return find(o) != end(); }

    /// @brief Erase an observation from the collection.
    ///
    /// The first observation comparing equal to \p o is removed from the collection and \c true is returned. If no
    /// observation in the collection compared equal to \p o, no change is made to the collection and \c false is
    /// returned.
    ///
    /// @param o
    bool erase(const observation& o);

    /// @brief Assign from an existing collection.
    ///
    /// @param c
    ///
    /// @return A reference to \c *this.
    const observation_collection& operator = (const observation_collection& c) { observations_ = c.observations_; return *this; }

    /// @name Comparison operators
    /// @{

    /// @brief Test for equality with another collection.
    ///
    /// @param c An observation_collection to compare \c this with.
    ///
    /// @return \c true iff all elements of \p c compare equal with corresponding elements of \c this and \p c has the
    /// same number of elements and \c this.
    bool operator == (const observation_collection& c) const { return observations_ == c.observations_; }

    /// @brief Test for inequality with another collection.
    ///
    /// @param c An observation_collection to compare \c this with.
    ///
    /// @return \c true iff an element of \p c compares not equal with a corresponding elements of \c this or if \p c
    /// has a differing number of elements to \c this.
    bool operator != (const observation_collection& c) const { return observations_ != c.observations_; }

    /// @}

protected:
    /// @brief The underlying collection of observations.
    observation_multimap observations_;

    /// @brief Real implementation of find() on underlying multimap.
    ///
    /// @param o
    ///
    /// @return <code>observations_.end()</code> if \p o is not in the collection.
    observation_multimap::const_iterator real_find(const observation& o) const;

    /// @brief Real implementation of find() on underlying multimap.
    ///
    /// @param o
    ///
    /// @return <code>observations_.end()</code> if \p o is not in the collection.
    observation_multimap::iterator real_find(const observation& o);
};

/// @brief Determine if an observation is within the light cone of a central observation.
///
/// The light cone is defined as a spatio-temporal cone centred on \p centre. Observations \c t away in time
/// from this observation are within the light cone if their distance in space is smaller than <code>abs(t) *
/// speed_of_light</code>.
///
/// Only observations whose time stamps are on the interval [\p first_time_stamp, \p last_time_stamp) will be
/// considered.
///
/// @note This is not part of biggles::observation_collection because it does not really count as 'core functionality'
/// to that data structure.
///
/// @param o the observation to determine if it is within the light cone
/// @param centre the centre of the light cone
/// @param speed_of_light the speed of light in pixels-per-frame
/// @param first_time_stamp the minimum time stamp for observations to consider
/// @param last_time_stamp the time stamp just after that of the observations to consider
bool observation_is_within_light_cone(
    const observation& o,
    const observation& centre,
    float speed_of_light,
    time_stamp first_time_stamp,
    time_stamp last_time_stamp);

/// @brief Obtain a set of observations within a light cone of a central observation.
///
/// The light cone is defined as a spatio-temporal cone centred on \p centre. Observations \c t away in time
/// from this observation are within the light cone if their distance in space is smaller than <code>abs(t) *
/// speed_of_light</code>.
///
/// Only observations whose time stamps are on the interval [\p first_time_stamp, \p last_time_stamp) will be
/// considered.
///
/// @note This is not part of biggles::observation_collection because it does not really count as 'core functionality'
/// to that data structure.
///
/// @tparam OutputIterator an iterator into which observations within the light cone will be written
/// @param observations a collection of observations
/// @param centre the centre of the light cone
/// @param speed_of_light the speed of light in pixels-per-frame
/// @param first_time_stamp the minimum time stamp for observations to consider
/// @param last_time_stamp the time stamp just after that of the observations to consider
/// @param output the iterator to write matching observations to
template<typename OutputIterator>
void observations_within_light_cone(const observation_collection& observations,
                                    const observation& centre,
                                    float speed_of_light,
                                    time_stamp first_time_stamp,
                                    time_stamp last_time_stamp,
                                    std::set<time_stamp>& time_stamps_with_obs,
                                    OutputIterator output);

/// \brief Determine if the observation collection does allow the formation of at least one valid track
///
/// A valid track consists of at least two observations that see each other (one is in the light cone of another)
/// at different time points
///
/// @param observations a collection of observations
/// @param speed_of_light the speed of light in pixels-per-frame
bool observation_collection_does_allow_tracks(const observation_collection& observations, float speed_of_light);

inline time_stamp time_of(observation_collection::observation_multimap::const_iterator& it) { return it->first; }
inline const observation& observation_of(observation_collection::observation_multimap::const_iterator& it) { return it->second; }


}

#define WITHIN_BIGGLES_OBSERVATION_COLLECTION_HPP__
#include "observation_collection.tcc"
#undef WITHIN_BIGGLES_OBSERVATION_COLLECTION_HPP__

#endif // BIGGLES_OBSERVATION_COLLECTION_HPP__
