/// @file simulate.hpp Simulate the dynamic model to generate synthetic datasets

#ifndef BIGGLES_SIMULATE_HPP__
#define BIGGLES_SIMULATE_HPP__

#include <boost/shared_ptr.hpp>

#include "model.hpp"
#include "point_2d.hpp"
#include "track.hpp"

namespace biggles {

/// @brief Simulate ground truth data and synthetic input for the Biggles tracker
namespace simulate
{

/// @brief Simulate a single track.
///
/// Given the model parameters specified in \p parameters, simulate a single track with birth time \p first_time_stamp. The
/// state evolution covariance matrix, \f$ Q \f$, is set to
///
/// \f[
/// Q = \left[
/// \begin{array}{cccc}
/// s_p^2 & 0 & 0 & 0 \\ 0 & s_v^2 & 0 & 0 \\ 0 & 0 & s_p^2 & 0 \\ 0 & 0 & 0 & s_v^2
/// \end{array}
/// \right]
/// \f]
///
/// where \f$ s_p \f$ is the value of \p position_sigma and \f$ s_v \f$ is the value of \p velocity_sigma.
///
/// @param parameters The model parameters to generate this track with.
/// @param first_time_stamp The birth time of the new track.
/// @param initial_position The initial position of the track.
/// @param initial_velocity The initial velocity of the track.
/// @param position_sigma The standard deviation of the position random walk.
/// @param velocity_sigma The standard deviation of the velocity random walk.
///
/// @return A shared pointer to the new track.
boost::shared_ptr<track> generate_track(std::deque<state_t> &states,
                                        const model::parameters& parameters,
                                        const point_2d& initial_position = point_2d(),
                                        const point_2d& initial_velocity = point_2d(),
                                        time_stamp first_time_stamp = 0,
                                        float position_sigma = 0.25f,
                                        float velocity_sigma = 0.01f);

boost::shared_ptr<track> generate_track(const model::parameters& parameters,
                                        const point_2d& initial_position = point_2d(),
                                        const point_2d& initial_velocity = point_2d(),
                                        time_stamp first_time_stamp = 0,
                                        float position_sigma = 0.25f,
                                        float velocity_sigma = 0.01f);

/// @brief Simulate an entire partition configuration.
///
/// @note All generated observations are included, even those outside of the bounds specified in the parameters. If you
/// want observations bounded to a region, use biggles::crop_partition() after calling this function.
///
/// @sa biggles::simulate_track()
/// @sa biggles::crop_partition()
///
/// @param[out] output_partition Write the generated partition to this reference.
/// @param params The model parameters.
/// @param start_time_stamp The first frame of data to generate.
/// @param end_time_stamp The last frame of data to birth tracks on. This may not be the last time stamp generated, see
/// the function documentation.
/// @param image_width The width of the output image. Used for generating clutter and choosing track initial position.
/// @param image_height The height of the output image. Used for generating clutter and choosing track initial position.
/// @param position_sigma The standard deviation of the position random walk.
/// @param velocity_sigma The standard deviation of the velocity random walk.
/// @param initial_velocity_sigma The standard deviation of the initial velocity for the track.
void generate_partition(partition_ptr_t& output_part_ptr,
                        const model::parameters& params,
                        time_stamp start_time_stamp = 0,
                        time_stamp end_time_stamp = 128,
                        float image_width = 128.f,
                        float image_height = 128.f,
                        float position_sigma = 0.25f,
                        float velocity_sigma = 0.01f,
                        float initial_velocity_sigma = 1.f);

/// @brief Crop out observations from a partition within a spatio-temporal window.
///
/// Remove observations from a partition whose x co-ordinate does not lie on [ \p first_x, \p last_x ) or whose y
/// co-ordinate does not lie on [ \p first_y, \p last_y ) or whose time stamp does not lie on [ \p first_time_stamp, \p
/// last_time_stamp ).
///
/// Tracks which end up with no observations are removed from the partition. Track birth and death times are modified to
/// be cropped to [ \p first_time_stamp, \p last_time_stamp ).
///
/// It is safe for \p output_partition and \p input_partition to refer to the same partition.
///
/// @param output_partition
/// @param input_partition
/// @param first_x
/// @param last_x
/// @param first_y
/// @param last_y
/// @param first_time_stamp
/// @param last_time_stamp
void crop_partition(partition_ptr_t& output_part_ptr,
                    const partition_ptr_t& input_part_ptr,
                    float first_x, float last_x,
                    float first_y, float last_y,
                    time_stamp first_time_stamp,
                    time_stamp last_time_stamp);

/// @brief Demote all tracks in a partition to the clutter.
///
/// Take all the tracks in \p input_partition and demote all their observations to the clutter. This function is useful
/// for generating synthetic input data sets from ground truths.
///
/// It is safe for \p output_partition and \p input_partition to refer to the same partition.
///
/// @param[out] output_partition overwritten with new partition
/// @param input_partition
void demote_all_tracks_from_partition(partition_ptr_t& output_part_ptr,
                                      const partition_ptr_t& input_part_ptr);

}

}

#endif // BIGGLES_SIMULATE_HPP__
