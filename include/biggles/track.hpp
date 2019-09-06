/// @file track.hpp A collection of observations associated with a single molecule

#ifndef BIGGLES_TRACK_HPP__
#define BIGGLES_TRACK_HPP__

#include <boost/assert.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <limits>
#include <stdexcept>
#include <set>

#include "point_2d.hpp"
#include "observation.hpp"
#include "observation_collection.hpp"
#include "types.hpp"

namespace biggles
{

// forward-declare kalman_filter since kalman_filter.hpp requires track be concrete
class kalman_filter;

/// @brief An inferred track associated with a set of observations.
///
/// Biggles represents a track as a series observations and a birth and death time.
///
/// Tracks are initialised using a shared collection of observations and indices into that collection. Observations
/// associated with the track are copied internally.
///
/// In addition to the locations, each track is associated with some observations from the underlying input data-set.
class track
{
public:
    /// @brief Construct an empty track with a birth time of 0 and a dynamic drag factor of 1.
    track()
        : first_time_stamp_(0)
        , last_time_stamp_(0)
        , boundary_min_location_(std::numeric_limits<float>::max())
        , boundary_max_location_(std::numeric_limits<float>::max())
        , dynamic_drag_(1.f)
    {
        cached_posterior_.second = 1.f;
    }

    /// @brief Fully construct a track from observation indices.
    ///
    /// The track is constructed with a birth time and death time as specified. Locations for the track are copied from [\p
    /// locations_first, \p locations_last) indexed by [\p observation_indices_first, \p observation_indices_last).
    ///
    /// @throw std::runtime_error if any observations share a common time stamp.
    ///
    /// @tparam ObservationIndexInputIterator An InputIterator concept over a collection of \c size_t.
    /// @tparam ObervationRandomAccessIterator A RandomAccessIterator over a collection of biggles::observation.
    /// @param first_time_stamp The time stamp associated with the first location in the track.
    /// @param last_time_stamp The first time stamp at which the track has no observations.
    /// @param observation_indices_first An iterator pointing to the first observation index associated with the track.
    /// @param observation_indices_last An iterator pointing just after the last observation associated with the track.
    /// @param observation_first An iterator pointing to the first observation whose collection is indexed by the track.
    /// @param observation_last An iterator pointing just after the last observation whose collection is indexed by the track.
    /// @param dynamic_drag The dynamic drag factor for this track (proportion of velocity added on each time step).
    template<typename ObservationIndexInputIterator, typename ObervationRandomAccessIterator>
    track(time_stamp first_time_stamp,
          time_stamp last_time_stamp,
          ObservationIndexInputIterator observation_indices_first,
          ObservationIndexInputIterator observation_indices_last,
          ObervationRandomAccessIterator observation_first,
          ObervationRandomAccessIterator observation_last,
          float dynamic_drag);

    /// @brief Fully construct a track from a set of observations.
    ///
    /// The track is constructed with a birth time and duration as specified. Locations for the track are copied from [\p
    /// observations_first, \p observations_last).
    ///
    /// @throw std::runtime_error if any observations share a common time stamp.
    ///
    /// @tparam ObervationInputIterator An InputIterator over a collection of biggles::observation.
    /// @param first_time_stamp The time stamp associated with the first location in the track.
    /// @param last_time_stamp The first time stamp at which the track has no observations.
    /// @param observation_first An iterator pointing to the first observation whose collection is indexed by the track.
    /// @param observation_last An iterator pointing just after the last observation whose collection is indexed by the track.
    /// @param dyanamic_drag The dynamic drag factor for this track (proportion of velocity added on each time step).
    template<typename ObervationInputIterator>
    track(time_stamp first_time_stamp,
          time_stamp last_time_stamp,
          ObervationInputIterator observation_first,
          ObervationInputIterator observation_last,
          float dynamic_drag = 1.f)
        : first_time_stamp_(first_time_stamp)
        , last_time_stamp_(last_time_stamp)
        , observations_(observation_first, observation_last)
        , dynamic_drag_(dynamic_drag)
    {
        BOOST_ASSERT(last_time_stamp_ >= first_time_stamp_);
        update_boundary_cache();
        check_common_timestamps();
    }

    /// @brief Merging constructor
    ///
    /// Initialise a track by merging two existing tracks.
    ///
    /// @throw std::runtime_error if any observations share a common time stamp.
    ///
    /// @param t1 The first track to merge.
    /// @param t2 The second track to merge.
    /// @param dyanamic_drag The dynamic drag factor for this track (proportion of velocity added on each time step).
    track(const track& t1, const track& t2, float dynamic_drag)
        : first_time_stamp_(std::min(t1.first_time_stamp_, t2.first_time_stamp_))
        , last_time_stamp_(std::max(t1.last_time_stamp_, t2.last_time_stamp_))
        , observations_(t1.observations_)
        , boundary_min_location_(
            std::min(x(t1.boundary_min_location_), x(t2.boundary_min_location_)),
            std::min(y(t1.boundary_min_location_), y(t2.boundary_min_location_)))
        , boundary_max_location_(
            std::max(x(t1.boundary_max_location_), x(t2.boundary_max_location_)),
            std::max(y(t1.boundary_max_location_), y(t2.boundary_max_location_)))
        , dynamic_drag_(dynamic_drag)
    {
        BOOST_ASSERT(last_time_stamp_ >= first_time_stamp_);
        observations_.insert(t2.observations_.begin(), t2.observations_.end());
        if(size() != t1.size() + t2.size())
            throw std::runtime_error("track::track: cannot merge tracks which share observations.");
        check_common_timestamps();
    }

    /// @brief Copy constructor.
    track(const track& t)
        : first_time_stamp_(t.first_time_stamp_)
        , last_time_stamp_(t.last_time_stamp_)
        , observations_(t.observations_)
        , boundary_min_location_(t.boundary_min_location_)
        , boundary_max_location_(t.boundary_max_location_)
        , dynamic_drag_(t.dynamic_drag_)
    { }

    /// @brief Copy constructor with dynamic drag.
    /// @param dynamic_drag
    track(const track& t, float dynamic_drag)
        : first_time_stamp_(t.first_time_stamp_)
        , last_time_stamp_(t.last_time_stamp_)
        , observations_(t.observations_)
        , boundary_min_location_(t.boundary_min_location_)
        , boundary_max_location_(t.boundary_max_location_)
        , dynamic_drag_(dynamic_drag)
    { }

    /// @brief The duration in ticks of this track.
    time_stamp duration() const { return last_time_stamp_ - first_time_stamp_; }

    /// @brief Obtain a reference to the underlying observation collection.
    ///
    /// This reference is const by design. Any mutation should be done via the track member functions themselves so that
    /// cached bounding box information can be preserved.
    const observation_collection& observations() const { return observations_; }

    /// @name STL-like container access
    ///
    /// The track can be accessed like an immutable STL container (one with only a \p const_iterator) via the following
    /// member functions.
    /// @{

    /// @brief A const iterator implementing the \c BidirectionalIterator concept.
    typedef observation_collection::const_iterator const_iterator;

    /// @brief An iterator implementing the \c BidirectionalIterator concept.
    ///
    /// @note This iterator is effectively a \c const_iterator since the value it points to cannot be modified.
    typedef observation_collection::const_iterator iterator;

    /// @brief Determine if the track is empty.
    ///
    /// An empty track is one with no observations defined.
    bool empty() const { return observations_.empty(); }

    /// @brief The number of observations in this track.
    size_t size() const { return observations_.size(); }

    /// @brief Obtain an iterator pointing to the first observation.
    const_iterator begin() const { return observations_.begin(); }

    /// @brief Obtain an iterator pointing just beyond the last observation.
    const_iterator end() const { return observations_.end(); }

    /// @brief Insert an observation into the track.
    ///
    /// The first and last timestamps of the track will be extended if necessary and the spatio-temporal bounding box is
    /// updated.
    ///
    /// @note This is relatively expensive if you repeatedly insert observations.
    ///
    /// @param o
    void insert(const observation& o);

    /// @brief sets first and last time stamps from the observations
    ///
    /// this isn't needed by the main algorithm but it is useful for certain python applications
    void time_stamps_from_observations();

    /// @brief Erase the given observation from the collection if present.
    ///
    /// @param o
    void erase(const observation& o)
        { observations_.erase(o); update_boundary_cache(); }

    /// @}

    /// @name Spatio-temporal bounding box.
    /// @{

    /// @brief The birth time of the track.
    time_stamp first_time_stamp() const { return first_time_stamp_; }

    /// @brief The death time. This is the birth time plus the duration.
    time_stamp last_time_stamp() const { return last_time_stamp_; }

    /// @brief The minimum co-ordinate of the spatial bounding box.
    const point_2d& min_location() const { return boundary_min_location_; }

    /// @brief The maximum co-ordinate of the spatial bounding box.
    const point_2d& max_location() const { return boundary_max_location_; }

    float radius() const;

    /// @}

    /// @brief The dynamic drag associated with this track.
    float dynamic_drag() const { return dynamic_drag_; }

    /// @brief Equality operator
    ///
    /// @param t A track to compare with this one.
    ///
    /// @return false iff the birth and death times of \p t and this track differ or if any non-common observations are
    /// found. In addition, the bounding boxes are checked.
    bool operator == (const track& t) const { return !(t!=*this); }

    /// @brief Inequality operator
    ///
    /// @param t A track to compare with this one.
    ///
    /// @return true iff the birth and death times of \p t and this track differ or if any non-common observations are
    /// found. In addition, the bounding boxes are checked.
    bool operator != (const track& t) const;

    /// @brief Calculate the log-pdf for the track <em>observations</em> given the parameters.
    ///
    /// This function works by using a biggles::kalman_filter to sample missing states for a track and then to calculate
    /// the likelihood of the observations we've seen given the parameters.
    ///
    /// Specifically, suppose we have an observation \f$ y \f$, a predicted state, \f$ \hat{x} \f$ and state estimation
    /// error estimate, \f$ \hat{P} \f$. Then the total error in the predicted observation, \f$ \hat{y} = B \hat{x} \f$ is
    /// given by \f$ \hat{\Sigma} = B \hat{P} B^T + R \f$. We therefore calculate the likelihood of \f$ y \f$ assuming a
    /// Gaussian model:
    ///
    /// \f[
    /// P(y | \hat{x}, \hat{P}) = \mathcal{N}(y ; \hat{y}, \hat{\Sigma}).
    /// \f]
    ///
    /// We combine all these likelihoods for each observation in the track.
    ///
    /// @sa biggles::kalman_filter
    ///
    /// @param parameters
    ///
    /// @return The value of \f$ \ell(d_i|t_i, \theta) \f$.
    float log_posterior(const model::parameters& parameters) const;

    /// @brief Construct a Kalman filter instance from a set of model parameters.
    ///
    /// This may be used to access the interpolated states for the track. The track may cache this value for a given set
    /// of parameters but don't rely on that behaviour for efficiency.
    ///
    /// @param p
    boost::shared_ptr<const kalman_filter> make_kalman_filter(const model::parameters& p) const;

    /// @brief returns the the number of sides where the track is extendible (returns either 0, 1 or 2)
    size_t num_possible_extensions(time_stamp min_t, time_stamp max_t) const {
        return size_t(min_t < first_time_stamp()) + size_t(last_time_stamp() < max_t);
    }

    /// \brief clear the cached posterior
    void clear_chached_posterior() const { cached_posterior_.second = 1.f; }

protected:
    /// @brief The time stamp of the first location within the track.
    time_stamp first_time_stamp_;

    /// @brief The number of time ticks the track lasts.
    time_stamp last_time_stamp_;

    /// @brief The observations associated with the track.
    observation_collection observations_;

    /// @brief The minimum co-ordinate of the spatial bounding box.
    point_2d boundary_min_location_;

    /// @brief The maximum co-ordinate of the spatial bounding box.
    point_2d boundary_max_location_;

    /// @brief The cached posterior for a particular set of model parameters.
    mutable std::pair<matrix2f, float> cached_posterior_;

    /// @brief The dyanamic drag factor associated with this track.
    float dynamic_drag_;

    /// @brief Update cached location boundaries.
    void update_boundary_cache();

    /// @brief Check for common timestamps.
    ///
    /// @throw std::runtime_error if any observations with common timestamps are found.
    void check_common_timestamps();

private:
    /// @brief Assignment constructor. Private to disallow mutating assignment.
    const track& operator = (const track& t);
};

/// @brief A reference-counted shared const pointer to a biggles::track.
typedef boost::shared_ptr<const track> shared_const_track_ptr;

/// @brief A reference-counted shared non-const pointer to a biggles::track.
typedef boost::shared_ptr<track> shared_track_ptr;

/// @brief A pair of reference-counted shared const pointers to a biggles::track.
typedef std::pair<shared_const_track_ptr, shared_const_track_ptr> shared_const_track_ptr_pair;


/// @brief Test if point \p p is within the spatial bounding box of the track.
///
/// @param t The track whose bounds should be queried.
/// @param p The point to test.
inline bool within_bounds(const track& t, const point_2d& p);

/// @brief Test if point \p p is within the spatial bounding box and \p t is within the temporal bounding box of the
/// track
///
/// @param t The track whose bounds should be queried.
/// @param p The point to test.
/// @param ts The time stamp to test
inline bool within_bounds(const track& t, const point_2d& p, time_stamp ts);

/// @brief Test if observation \p o is within the spatio-temporal bounding box of the track.
///
/// @param t The track whose bounds should be queried.
/// @param o The observation to test.
inline bool within_bounds(const track& t, const observation& o);

/// @brief the "distance" of two tracks
///
/// @return the minimal distance between the boundary boxes of the two tracks
float dist(const track& t1, const track& t2);

/// @brief Verify that a track does not violate the laws of physics.
///
/// @param track
bool verify_track(const track& track);

} // biggles

#define WITHIN_BIGGLES_TRACK_HPP__
#include "track.tcc"
#undef WITHIN_BIGGLES_TRACK_HPP__

#endif // BIGGLES_TRACK_HPP__
