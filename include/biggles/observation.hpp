/// @file observation.hpp A tuple type representing an observed feature

#ifndef BIGGLES_OBSERVATION_HPP__
#define BIGGLES_OBSERVATION_HPP__

#include <utility>
#include <functional>
#include <cmath>

#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "detail/fun.hpp"

namespace biggles
{

/// @brief A type used to represent an observation's time stamp or frame index.
typedef off_t time_stamp; // shall be signed integer type

/*
/// @brief An alias to \c boost::tuple representing an observation.
typedef boost::tuple<float, float, time_stamp> observation;
*/

class observation_t;
typedef boost::shared_ptr<observation_t> observation;

class observation_t {
    float x_, y_;
    time_stamp t_;
public:
    observation_t() : x_(0), y_(0), t_(0) {}
    observation_t(const float x, const float y, const float t) : x_(x), y_(y), t_(t) {}
    observation_t(const observation_t& o) : x_(o.x_), y_(o.y_), t_(o.t_) {}
    observation_t& operator=(const observation_t& o) {
        if (&o == this) return *this;
        x_ = o.x_;
        y_ = o.y_;
        t_ = o.t_;
        return *this;
    }
    //
    bool operator==(const observation_t& o) const {
        return x_ == o.x_ and y_ == o.y_ and t_ == o.t_;
    }
    bool operator<(const observation_t& o) const {
        return t_ < o.t_ or (t_ == o.t_ and x_ < o.x_) or (t_ == o.t_ and x_ == o.x_ and y_ < o.y_);
    }
    bool operator!=(const observation_t& o) const { return not (*this == o); }
    //
    friend const float x(const observation& o);
    friend const float y(const observation& o);
    friend const time_stamp t(const observation& o);
    void set_x(const float x) { x_ = x; }
    void set_y(const float y) { y_ = y; }
    void set_t(const float t) { t_ = t; }
};


/// @brief Extract a reference to a observation's x-co-ordinate (const).
inline const float x(const observation& o) { return o->x_; }

/// @brief Extract a reference to a observation's y-co-ordinate (const).
inline const float y(const observation& o) { return o->y_; }

/// @brief Extract a reference to a observation's time stamp (const).
inline const time_stamp t(const observation& o) { return o->t_; }

/*
/// @brief Extract a reference to a observation's x-co-ordinate (non-const).
inline float& x(observation& o) { return boost::tuples::get<0>(o); }

/// @brief Extract a reference to a observation's y-co-ordinate (non-const).
inline float& y(observation& o) { return boost::tuples::get<1>(o); }

/// @brief Extract a reference to a observation's time stamp (non-const).
inline time_stamp& t(observation& o) { return boost::tuples::get<2>(o); }
*/

/// \brief creates a new observation from cooridantes
inline observation new_obs(float x, float y, time_stamp t) {
    return observation(new observation_t(x, y, t));
}

/// \brief is the observation defined? If no, then it has no coordinates
inline bool noobs(const observation& o) { return o.use_count() == 0; }

/// \brief a compare functor for sets of observations
struct earlier { bool operator()(const observation& a, const observation& b) { return t(a) < t(b); } };

/// @brief Test if observations are in each others light cone
inline bool observations_see_each_other(const observation& o1, const observation& o2, const float& speed_of_light);
inline float observations_see_each_other2(const observation& o1, const observation& o2, const float& speed_of_light);

/// \brief Euclidean distance disregarding the time.
inline float dist(const observation& p, const observation& q) {
    return sqrtf(square(x(p) - x(q)) + square(y(p) - y(q)));
}

inline observation convex_combination(const observation& o1, const observation& o2, time_stamp t0) {
    const float a1 = float(t(o2) - t0) / float(t(o2) - t(o1));
    return new_obs(a1*x(o1) + (1.f-a1)*x(o2), a1*y(o1) + (1.f-a1)*y(o2), t0);
}


} // namespace biggles

#define WITHIN_BIGGLES_OBSERVATION_HPP__
#include "observation.tcc"
#undef WITHIN_BIGGLES_OBSERVATION_HPP__


#endif // BIGGLES_OBSERVATION_HPP__
