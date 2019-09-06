#ifndef BIGGLES_TOOL_DEBUG_HPP___
#define BIGGLES_TOOL_DEBUG_HPP___

#include <iostream>
#include <boost/timer.hpp>
#include <boost/shared_ptr.hpp>

#define OK(a) std::cerr << "OK " << #a << " = " << (a) << std::endl
#define OK1(a) std::cerr << "OK " << #a << " = " << std::endl << (a) << std::endl

struct TimeObserver {
    double t0, t1, t2, t3, t4;
    boost::timer timer;
    TimeObserver() : t0(0.0), t1(0.0), t2(0.0), t3(0.0), t4(0.0) {}
};

typedef boost::shared_ptr<TimeObserver> time_observer_ptr;

#endif //BIGGLES_TOOL_DEBUG_HPP___
