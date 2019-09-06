#include "biggles/observation.hpp"
#include "biggles/observation_collection.hpp"
#include "biggles/samplers.hpp"
#include "biggles/sampling/simple.hpp"
#include <boost/tuple/tuple_io.hpp>
#include <sstream>

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

std::ostream& operator<<(std::ostream& os, const observation& o) {
    os << "(" << x(o) << ", " << y(o) << ", " << t(o) << ")";
    return os;
}

int main(int argc, char** argv)
{
    plan_tests(6 - 1 );

    // a set of dummy observations
    observation obs_array[] = {
        new_obs(0.f, 0.f, 3),
        new_obs(0.f, 2.f, 3),
        new_obs(1.f, 1.f, 3),
        new_obs(-1.f, 3.f, 3),
        new_obs(0.f, 0.f, 3), // intentionally duplicate observation

        // intentionally miss time stamp 4

        new_obs(3.f, 3.f, 5),
        new_obs(2.f, 3.f, 5),
        new_obs(3.f, 2.f, 5),
        new_obs(3.f, 4.f, 5),
        new_obs(7.f, 4.f, 5),
        new_obs(7.f, 7.f, 5),

        // out of order 6 & 7

        new_obs(4.f, 1.5f, 7),
        new_obs(2.f, 1.5f, 7),
        new_obs(1.f, 1.5f, 7),

        new_obs(1.f, 1.5f, 6),
        new_obs(1.5f, 0.5f, 6),
    };
    size_t obs_array_size = sizeof(obs_array) / sizeof(observation);
    ok1(obs_array_size > 0);
    ok1(obs_array_size < 32);

    // construct via iterator
    observation_collection oc(obs_array, obs_array + obs_array_size);
    //ok1(oc.size() == obs_array_size);
    ok1(!oc.empty());
    ok1(oc.first_time_stamp() == 3);
    ok1(oc.last_time_stamp() == 8);

    return exit_status();
}
