#ifndef BIGGLES_DETAIL_PHYSICS_HPP___
#define BIGGLES_DETAIL_PHYSICS_HPP___

#include <Eigen/Dense>

#include "../observation.hpp"

namespace biggles { namespace detail
{

/// @brief This is the default maximum separation of observations to be considered for observations to be neighbours.
static const time_stamp last_delta_t_ = 31;

/// @brief The default maximum separation of observations to be considered for extending a track.
static const time_stamp extend_last_delta_t_ = 2;

/// @brief The default maximum separation of observations to be considered for extending a track.
static const float extend_speed_of_light_ = 1.2f;

/// @brief This is the default speed of light.
///
/// @note If this were any other algorithm, specifying a 'default' speed of light would be a bit of a WTF.
static const float speed_of_light_ = 3.f;

inline Eigen::Matrix4f initQ() {
    Eigen::Matrix4f Q(Eigen::Matrix4f::Zero());
    Q(0,0) = Q(2,2) = 1e-1f * 1e-1f * speed_of_light_ * speed_of_light_;
    Q(1,1) = Q(3,3) = 1e-2f * 1e-2f * speed_of_light_ * speed_of_light_;
    return Q;
}


} }

#endif // BIGGLES_DETAIL_PHYSICS_HPP___
