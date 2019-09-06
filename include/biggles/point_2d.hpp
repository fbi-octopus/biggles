/// @file point_2d.hpp A pair of floats representing a 2D point location

#ifndef BIGGLES_POINT_2D_HPP__
#define BIGGLES_POINT_2D_HPP__

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

namespace biggles
{

/// @brief An alias to \c boost::tuple representing a 2D floating point point
typedef boost::tuples::tuple<float, float> point_2d;

/// @brief Extract a reference to a point's x-co-ordinate (non-const).
inline float& x(point_2d& p) { return boost::tuples::get<0>(p); }

/// @brief Extract a reference to a point's x-co-ordinate (const).
inline const float& x(const point_2d& p) { return boost::tuples::get<0>(p); }

/// @brief Extract a reference to a point's y-co-ordinate (non-const).
inline float& y(point_2d& p) { return boost::tuples::get<1>(p); }

/// @brief Extract a reference to a point's y-co-ordinate (const).
inline const float& y(const point_2d& p) { return boost::tuples::get<1>(p); }

}

#endif // BIGGLES_POINT_2D_HPP__
