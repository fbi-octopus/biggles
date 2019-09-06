/// @file biggles.hpp Include all other header files for the Biggles library

#ifndef BIGGLES_HPP__
#define BIGGLES_HPP__

/// @page biggles_data_model Data model
///
/// Here we briefly discuss the data model Biggles uses to represent its input and the tracking results it generates.
///
/// @section biggles_data_model_observation Observations
///
/// The input to the Biggles tracker is a set of <em>observations</em>. These are typically the output from a feature
/// detection stage. An observation is a pair of floating point values giving the (floating-point) x- and y-co-ordinate
/// of the feature along with an integral time stamp.
///
/// An observation is represented via the biggles::observation type which is internally a wrapper around \c
/// boost::tuple. Convenience functions x(), y() and t() are provided to get access to the underlying observation
/// co-ordinate values.
///
/// Observations are often collected together into a biggles::observation_collection. This class provides STL-container
/// like access to observations but also allows for specialised queries to be made on the collection such as obtaining a
/// pair of iterators delimiting those observations with a specified time stamp. Internally
/// biggles::observation_collection maintains a sorted list of observations ordered by time stamp making a) these
/// operations efficient and b) the guarantee that iterating over the collection will return observations in ascending
/// time stamp order.
///
/// @section biggles_data_model_track Tracks
///
/// A track is a light-weight wrapper around a biggles::observation_collection and is represented by the biggles::track
/// class. In the algorithms, a track is comprised of two time stamps delimiting the life-time of the track and a set of
/// observations which are associated with the track. Internally a track also maintains a spatio-temporal bounding box
/// allowing quick queries about whether a given observation could be in the track or not. The convenience function
/// within_bounds() is provided to query an observation against a track's bounding box.
///
/// In addition to the above a biggles::track also enforces the rule that no track may contain two or more observations
/// which share a common time stamp.
///
/// A biggles::track_collection provides a STL-container like collection of tracks. Each track is stored as a \c
/// boost::shared_ptr which allows for implicit sharing of tracks between collections. Since the tracks form the bulk of
/// Biggles' state and the state is often duplicated as part of the algorithm, implicit sharing like this makes sense.
/// To facilitate sharing of tracks, a biggles::track_collection provides only \c const access to the tracks. It is
/// envisages that one never modifies a track once it is placed into a collection; tracks should be modified by cloning,
/// modification and re-insertion.
///
/// @section biggles_data_model_partition Partitions
///
/// The output from the Biggles algorithm includes a <em>partition</em> of the input observations. Each observation is
/// partitioned into zero or more tracks or a special track&mdash;termed the <em>clutter</em>&mdash;which contains all
/// observations deemed to be the result of spurious feature detection rather than due to an underlying molecule. Unlike
/// normal tracks the clutter is allowed to have observations which share common time stamps.
///
/// A partition is represented via the biggles::partition class and is in essence a pointer to a
/// biggles::track_collection of tracks and a biggles::observation_collection for the clutter. Again the actual data is
/// stored via a \c boost::shared_ptr allowing partitions to be light-weight objects for passing around and returning
/// from functions. It also allows implicit sharing of clutters between partitions who differ only in their tracks.
///
/// To make this implicit sharing work biggles::partition provides only \c const access to the clutter and tracks.
/// Should you wish to modify the partition, you must first clone the clutter and/or track collection.

/// @brief The Biggles tracker API.
namespace biggles { }

#include "observation_collection.hpp"
#include "kalman_filter.hpp"
#include "mh_moves/mh_moves.hpp"
#include "model.hpp"
#include "observation.hpp"
#include "partition.hpp"
#include "point_2d.hpp"
#include "samplers.hpp"
#include "simulate.hpp"
#include "track.hpp"
#include "track_collection.hpp"
#include "tracker.hpp"

#endif // BIGGLES_HPP__
