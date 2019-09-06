#ifndef WITHIN_BIGGLES_OBSERVATION_HPP__
#error "this file must only be included by observation.hpp"
#endif // WITHIN_BIGGLES_OBSERVATION_HPP__

#include <cmath>
namespace biggles
{

inline bool observations_see_each_other(const observation& o1, const observation& o2, const float& speed_of_light)
{

    time_stamp dt(std::abs(t(o1) - t(o2)));
    float dx(x(o1) - x(o2));
    float dy(y(o1) - y(o2));
    float ds_sq(dx*dx + dy*dy);
    float max_ds(dt * speed_of_light);
    if(ds_sq > max_ds*max_ds)
        return false;
    return true;
}

inline float observations_see_each_other2(const observation& o1, const observation& o2, const float& speed_of_light)
{
    float ds(dist(o1, o2));
    float max_ds(std::abs(t(o1) - t(o2)) * speed_of_light);
    if(ds > max_ds)
        return -1.f;
    return sqrtf(square(x(o1)-x(o2)) + square(y(o1)-y(o2)) + float(square(t(o2)-t(o1))));
    //return ds;
}

}
