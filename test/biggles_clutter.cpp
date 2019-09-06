#include <boost/format.hpp>
#include "biggles/clutter.hpp"
#include "biggles/sampling/simple.hpp"
#include <iostream>
#include <sstream>
//#include <stdlib.h>

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

void test_proc1() {
    clutter_t clutter;
    clutter.insert(new_obs(1, 1, 0));
    clutter.insert(new_obs(1, 1, 3));
    clutter.insert(new_obs(3, 1, 4));
    clutter.insert(new_obs(0, 5, 6));
    clutter.insert(new_obs(1, 1, 7));
    clutter.insert(new_obs(1, 0, 7));
    clutter.insert(new_obs(1, 1, 9));
    ok1(clutter.size() == 7);
    ok1(clutter.count_at_time_stamp(7) == 2);
    ok1(clutter.count_at_time_stamp(5) == 0);
    ok1(clutter.locate_near(5, new_obs(0,0, 5), 3).size() == 0);
    ok1(clutter.locate_near(7, new_obs(0,0, 7), 3).size() == 2);
    diag("first ts = %zu, last ts = %zu", clutter.first_time_stamp(), clutter.last_time_stamp());
    ok1(clutter.duration() == 10);
    ok1(clutter.empty() == false);
    ok1(clutter.count_at_time_stamp(0) == 1);
    std::deque<observation> neigbours = clutter.locate_near(0, new_obs(0, 0, 0), 2);
    clutter.erase(neigbours[0]);
    ok1(clutter.count_at_time_stamp(0) == 0);
    clutter_t new_clutter(clutter);
    diag("new_clutter.size() = %zu, clutter.size() = %zu", new_clutter.size(), clutter.size());
    ok1(new_clutter.size() == clutter.size());
    neigbours = clutter.locate_near(7, new_obs(0, 0, 7), 2);
    new_clutter.erase(neigbours[0]);
    ok1(new_clutter.size() + 1 == clutter.size());
    ok1(x(clutter.lower_bound()) == 0);
    ok1(y(clutter.lower_bound()) == 0);
    ok1(x(clutter.upper_bound()) == 3);
    ok1(y(clutter.upper_bound()) == 5);
}

int main(int argc, char** argv)
{
    plan_tests(15);

    test_proc1();

    return exit_status();
}
