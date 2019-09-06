/// @file partition.hpp A set of tracks and clutter observations which partition all observed features

#ifndef BIGGLES_PARTITION_HPP__
#define BIGGLES_PARTITION_HPP__

#include <boost/shared_ptr.hpp>
#include <boost/any.hpp>
#include <deque>
#include <Eigen/Dense>
#include <iterator>

#include "observation_collection.hpp"
#include "track_collection.hpp"
#include "detail/possible_move_recorder.hpp"
#include "detail/observation_reservoir.hpp"
#include "clutter.hpp"

namespace biggles
{


typedef boost::shared_ptr<const clutter_t> const_clutter_ptr;
typedef boost::shared_ptr<clutter_t> clutter_ptr;


/// @brief A configuration of tracks and spurious observations.
///
/// A full configuration within Biggles consists of a collection of tracks for observations assigned to molecules and a
/// collection of observations deemed to be spurious detections. These spurious detections are called the 'clutter'.
///
/// A partition instance maintains a \c boost::shared_ptr pointing to a const biggles::track_collection and const
/// biggles::observation_collection. In this way, partitions are semi-immutable. They are designed to be lightweight and
/// used essentially to refer to a particular pair of track collections and clutter.
///
/// They are thus immutable by design except that they may be assigned to. The track collections and clutter within them
/// cannot be changed.
///
/// A partition's immutability allows it to maintain a few cached parameters, such as minimum and maximum time stamps.
class partition
{
public:

    /// \brief The spatial expansion of a partition
    ///
    /// This is required to calculate the clutter likelihood
    struct expansion_2d {
        float x_lo, x_hi, y_lo, y_hi;
        expansion_2d(float xlo, float xhi, float ylo, float yhi) : x_lo(xlo), x_hi(xhi), y_lo(ylo), y_hi(yhi) {}
        expansion_2d(const expansion_2d& e) : x_lo(e.x_lo), x_hi(e.x_hi), y_lo(e.y_lo), y_hi(e.y_hi) {}
        /// \brief Increases the expansion, if \p min_location or \p max_location are not included in it
        void grow(const point_2d &min_location, const point_2d &max_location) {
            x_lo = std::min(x_lo, min_location.get<0>());
            y_lo = std::min(y_lo, min_location.get<1>());
            x_hi = std::max(x_hi, max_location.get<0>());
            y_hi = std::max(y_hi, max_location.get<1>());
        }
        /// \brief Increases the expansion, if \p location is not included in it
        void grow(const point_2d &location){
            x_lo = std::min(x_lo, location.get<0>());
            y_lo = std::min(y_lo, location.get<1>());
            x_hi = std::max(x_hi, location.get<0>());
            y_hi = std::max(y_hi, location.get<1>());
        }
    private:
        /// \brief default construction is not allowed
        expansion_2d();
    };

    /// @brief A shared pointer to a biggles::track_collection.
    typedef boost::shared_ptr<const track_collection> track_collection_ptr;

    /// @brief A shared pointer to a biggles::observation_collection.
    typedef boost::shared_ptr<const observation_collection> observation_collection_ptr;

    typedef boost::shared_ptr<const ObservationReservoir> observation_pool_ptr;

    /// @brief Default constructor.
    partition()
        : tracks_(new track_collection())
        , clutter_(new clutter_t())
        , observation_pool_(new ObservationReservoir())
        , first_time_stamp_(0)
        , last_time_stamp_(0)
        , expansion_(0.f, 0.f, 0.f, 0.f)
        , observation_count_(0)
    { }

    /// @brief Fully initialise a partition.
    ///
    /// @param tracks
    /// @param clutter
    partition(const track_collection_ptr& tracks,
        const const_clutter_ptr& clutter);

    partition(const track_collection_ptr& tracks, const const_clutter_ptr& clutter, const expansion_2d& expansion);
    partition(const observation_pool_ptr& pool, const track_collection_ptr& tracks, const const_clutter_ptr& clutter);

    partition(const observation_pool_ptr& pool, const track_collection_ptr& tracks,
        const const_clutter_ptr& clutter, const expansion_2d& expansion)
        : tracks_(tracks), clutter_(clutter), observation_pool_(pool), expansion_(expansion)
    {
        update_boundary_cache();
        update_observation_count();
    }
    /// @brief Copy constructor.
    ///
    /// @param p A reference to a partition to copy from.
    partition(const partition& p)
        : tracks_(p.tracks_)
        , clutter_(p.clutter_)
        , observation_pool_(p.observation_pool_)
        , first_time_stamp_(p.first_time_stamp_)
        , last_time_stamp_(p.last_time_stamp_)
        , expansion_(p.expansion_)
        , observation_count_(p.observation_count_)
    { }

    /// @brief Obtain a reference to the track collection for this partition.
    const track_collection& tracks() const { return *tracks_; }

    /// @brief Obtain a reference to the clutter observations for this partition.
    const clutter_t& clutter() const { return *clutter_; }

    /// @brief Obtain a reference to the track collection shared pointer.
    const track_collection_ptr& tracks_ptr() const { return tracks_; }

    /// @brief Obtain a reference to the clutter shared pointer.
    const const_clutter_ptr& clutter_ptr() const { return clutter_; }

    /// @brief The first time stamp covered by this partition.
    time_stamp first_time_stamp() const { return first_time_stamp_; }
    void set_first_time_stamp(const time_stamp& new_time_stamp) { first_time_stamp_ = new_time_stamp; }

    /// @brief The first time stamp <em>not</em> covered by this partition after first_time_stamp().
    time_stamp last_time_stamp() const { return last_time_stamp_; }
    void set_last_time_stamp(const time_stamp& new_time_stamp) { last_time_stamp_ = new_time_stamp; }

    /// @brief The number of frames covered by this partition.
    time_stamp duration() const { return last_time_stamp() - first_time_stamp(); }

    /// @brief Assignment operator.
    ///
    /// @param p A reference to a partition to assign from.
    ///
    /// @return A const reference to \c this.
    const partition& operator = (const partition& p) {
        if (this == &p) return *this;
        tracks_ = p.tracks_;
        clutter_ = p.clutter_;
        observation_pool_ = p.observation_pool_;
        first_time_stamp_ = p.first_time_stamp_;
        last_time_stamp_ = p.last_time_stamp_;
        expansion_ = p.expansion_;
        observation_count_ = p.observation_count_;
        return *this;
    }

    /// \brief Constant access to the expansion
    const expansion_2d& expansion() const { return expansion_; }
    /// \brief Set the expansion from the minimal and maximal coordinates
    void set_expansion(const float xlo, const float xhi, const float ylo, const float yhi) {
        expansion_ = expansion_2d(xlo, xhi, ylo, yhi);
    }
    /// \brief The area of the rectangle which is the expansion
    float volume() const {
        return (expansion_.x_hi - expansion_.x_lo) * (expansion_.y_hi - expansion_.y_lo);
    }

    std::string string_from_track(const shared_const_track_ptr& track_p) const;
    std::string as_string() const;

    const std::string str() const { return "track_partition"; }

    /// \brief the total number of observations
    const size_t observation_count() const { return observation_count_; }

    /// \brief the observation pool
    const observation_pool_ptr& pool() const { return observation_pool_; }

    void set_capability_recorder(const capability_recorder_ptr &cr_ptr) { cr_ptr_ = cr_ptr; }
    capability_recorder_ptr get_capability_recorder() const { return cr_ptr_; }

    template<class VALUE_TYPE> void set_misc(const VALUE_TYPE& val) { misc_ = val; }
    boost::any misc() const { return misc_; }

protected:
    /// @brief The collection of tracks associated with this partition.
    track_collection_ptr tracks_;

    /// @brief The clutter observations associated with this partition.
    const_clutter_ptr clutter_;

    /// @brief A reference to the pool of all observations
    observation_pool_ptr observation_pool_;

    /// @brief The minimum time stamp observed in the partition.
    time_stamp first_time_stamp_;

    /// @brief A time stamp just beyond the last observed in the partition.
    time_stamp last_time_stamp_;

    /// @brief The spatial expantion
    expansion_2d expansion_;

    /// @brief The total number of observations
    size_t observation_count_;

    /// \brief This is ugly: pointer to capability recorder
    capability_recorder_ptr cr_ptr_;

    /// \brief Anything else
    boost::any misc_;

    /// @brief Update the partition boundary information from the tracks and clutter observations.
    void update_boundary_cache();

    /// @brief Update the total observation count from the tracks and clutter observations.
    void update_observation_count();

    void set_observation_pool_from_tracks_clutter() {
        std::deque<observation> all_obs;
        ingest_observations(all_obs, *clutter_);
        ingest_observations(all_obs, *tracks_);
        observation_pool_ = observation_pool_ptr(new ObservationReservoir(all_obs.begin(), all_obs.end()));
    }

};

struct internal_containers_observer {
    size_t size_cross_over_pairs;
    size_t size_extendible_tracks;
    size_t size_mergeable_pairs;
    size_t size_transfer_pairs;

    const internal_containers_observer& update(const capability_recorder& p) {
        size_cross_over_pairs = p._size_cross_over_pairs_();
        size_extendible_tracks = p._size_extendible_tracks_();
        size_mergeable_pairs = p._size_mergeable_pairs_();
        size_transfer_pairs = p._size_transfer_pairs_();
        return *this;
    }

    std::deque<uint64_t> get() const {
        std::deque<uint64_t> values(4);
        values[0] = size_mergeable_pairs;
        values[1] = size_cross_over_pairs;
        values[2] = size_transfer_pairs;
        values[3] = size_extendible_tracks;
        return values;
    }
};

typedef boost::shared_ptr<partition> partition_ptr_t;

/// \brief creates a smallest expansion that contains all \p tracks and all \p clutter observations
partition::expansion_2d expansion_from_observations(const track_collection& tracks,
                        const clutter_t& clutter);

partition value_of(const partition_ptr_t& partition_p);

inline capability_recorder_ptr new_cap_recorder(const partition& part) {
    return capability_recorder_ptr(
        new capability_recorder(part.first_time_stamp(), part.last_time_stamp(), part.tracks()));
}

/// \brief checks if a partition is valid
bool verify_partition(partition& p);
} // namespace biggles


#endif // BIGGLES_PARTITION_HPP__
