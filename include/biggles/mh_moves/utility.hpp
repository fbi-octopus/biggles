/// @file utility.hpp Common utility functions for the Metropolis-Hastings moves.

#ifndef BIGGLES_MH_MOVES_UTILITY_HPP__
#define BIGGLES_MH_MOVES_UTILITY_HPP__

#include <boost/shared_ptr.hpp>
#include <utility>

#include "../sampling/simple.hpp"
#include "../observation.hpp"
#include "../observation_collection.hpp"
#include "../track.hpp"
#include "../track_collection.hpp"

namespace biggles
{

namespace mh_moves
{

/// @brief A pair of reference-counted shared non-const pointers to a biggles::track.
typedef std::pair<shared_track_ptr, shared_track_ptr> shared_track_ptr_pair;

/// @brief A reference-counted shared non-const pointer to a biggles::track_collection.
typedef boost::shared_ptr<track_collection> shared_track_collection_ptr;

/// @brief A reference-counted shared non-const pointer to a biggles::observation_collection.
typedef boost::shared_ptr<observation_collection> shared_observation_collection_ptr;

inline bool operator<(const shared_const_track_ptr_pair& a, const shared_const_track_ptr_pair& b) {
    return (a.first < b.first) or (a.first == b.first and a.second < b.second);
}

/// \brief Calc the weights for choosing a split point
std::deque<float> get_split_weights(observation_collection::const_iterator fst,
    const observation_collection::const_iterator& lst);

/// \brief get the log probability to sample the split time stamp \e ts_to_split_at
float get_time_stamp_to_split_track_prob(const shared_const_track_ptr& track_to_split,
        const time_stamp& ts_to_split_at);

/// @brief Sample a track from the set of all those with more than \p minimum_size observations.
///
/// @param tc
/// @param minimum_size
/// @param output_log_prob
track_collection::const_iterator sample_track_with_minimum_size(
    const track_collection& tc,
    size_t minimum_size,
    float& output_log_prob);

/// @brief Split a track at a specified time stamp.
///
/// The new tracks will have last and first time stamps modified to reflect observations either side of the splitting
/// point. It may or may not leave noobs (time stamps without observation) at the splitting end
/// if there is a gap in the observations at the splitting point
///
/// @param track
/// @param split_ts
/// @param log_splitting_prob
///
/// @return A pair of shared pointers pointing to the new tracks.
shared_const_track_ptr_pair split_track_at_time_stamp(
        const shared_const_track_ptr& track,
        time_stamp split_ts, float& log_splitting_prob);

/// @brief Sample a time stamp within a track where a split could occurr.
///
/// @note This will ensure that the split tracks would always have at least two observations in them. Consequently the
/// result of calling this function with a track with fewer than four observations is undefined.
///
/// @param track
/// @param ts_to_split_at
/// @param log_prob_ts
bool sample_time_stamp_to_split_track(
        const shared_const_track_ptr& track,
        time_stamp& ts_to_split_at,
        float& log_prob_ts);

/// @brief Return \c true iff \p o1 and \p o2 could be neighbours in a track.
///
/// There are a few laws of physics built into Biggles w.r.t. maximum temporal and spatial distances. These have been
/// set to be very inclusive and probably shouldn't be changed unless you understand the implications.
///
/// @param o1
/// @param o2
bool observations_could_be_neighbours(
    const observation& o1,
    const observation& o2);

/// @brief Return \c true iff the nearest terminal observations of \p t1 and \p t2 could be neighbours in a merged track.
///
/// There are a few laws of physics built into Biggles w.r.t. maximum temporal and spatial distances. These have been
/// set to be very inclusive and probably shouldn't be changed unless you understand the implications.
///
/// @param t1
/// @param t2
bool tracks_could_be_merged(const track& t1,  const track& t2);
float tracks_could_be_merged2(const track& t1,  const track& t2);

bool track_could_transfer_obs_to_other(const track& t1,  const track& t2);

bool tracks_could_cross(const track& t1,  const track& t2);

shared_const_track_ptr_pair choose_two_tracks(const track_collection& tracks);

inline float light_years(time_stamp t1, time_stamp t2) {
    return std::abs(static_cast<float>(t1) - static_cast<float>(t2)) * detail::speed_of_light_;
}

namespace cross_over_fun {

    /// \brief the cross over data
    struct cross_over_data_t {
        time_stamp t_beg, t_end; // time limits (inclusive) at which the cut is possible
    };

    /// \brief given a track pair get the over-lapping time range
    bool no_cross_over_time(shared_const_track_ptr_pair& tracks_to_cross, cross_over_data_t& data);

    /// \brief given a track pair and the cross-over data get the time points where a cross over is possible
    /// spatial analysis
    void get_possible_times(shared_const_track_ptr_pair& tracks_to_cross, std::deque<time_stamp>& possible_times,
        const cross_over_data_t& data);
}

void print_track(const track& t1);

} // namespace mh_moves

} // namespace biggles

#define WITHIN_BIGGLES_MH_MOVES_UTILITY_HPP__
#include "utility.tcc"
#undef WITHIN_BIGGLES_MH_MOVES_UTILITY_HPP___

#endif // BIGGLES_MH_MOVES_UTILITY_HPP__
