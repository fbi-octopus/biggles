#include <stdlib.h>
#include <vector>

#include "biggles/detail/pair_collection.hpp"
#include "biggles/observation.hpp"
#include "biggles/point_2d.hpp"
#include "biggles/track.hpp"

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

void test_point_2d()
{
    point_2d p1, p2(3.f, 4.f), p3(p2);
    const point_2d p4(5,6);

    ok1(x(p1) == 0.f);
    ok1(y(p1) == 0.f);
    ok1(x(p2) == 3.f);
    ok1(y(p2) == 4.f);
    ok1(x(p3) == x(p2));
    ok1(y(p3) == y(p2));
    ok1(x(p4) == 5.f);
    ok1(y(p4) == 6.f);

    x(p2) = 4.5f;
    y(p2) = x(p4);

    ok1(x(p2) == 4.5f);
    ok1(y(p2) == 5.f);
}

void test_observation()
{
    diag("test observation");
    observation o1(new observation_t());
    observation o2(new_obs(3.f, 4.f, 100));
    observation o3(o2);
    const observation o4(new_obs(5,6,7));

    ok1(x(o1) == 0.f);
    ok1(y(o1) == 0.f);
    ok1(t(o1) == 0);
    ok1(x(o2) == 3.f);
    ok1(y(o2) == 4.f);
    ok1(t(o2) == 100);
    ok1(x(o3) == x(o2));
    ok1(y(o3) == y(o2));
    ok1(t(o3) == t(o2));
    ok1(x(o4) == 5.f);
    ok1(y(o4) == 6.f);
    ok1(t(o4) == 7);

    o2->set_x(4.5f);
    o2->set_y(x(o4));

    o1->set_t(1);
    o2->set_t(t(o1));

    ok1(x(o2) == 4.5f);
    ok1(y(o2) == 5.f);
    ok1(t(o2) == 1);
}

void test_track()
{
    observation obs[] = {
        new_obs(1.0f, 2.0f, 1), // 0
        new_obs(1.0f, 3.0f, 1), // 1
        new_obs(3.0f, 4.0f, 1), // 2
        new_obs(1.0f, 0.0f, 1), // 3

        new_obs(1.3f, 2.2f, 2), // 4
        new_obs(2.1f, 1.6f, 2), // 5
        new_obs(1.3f, 0.1f, 2), // 6
        new_obs(3.2f, 4.5f, 2), // 7

        new_obs(0.2f, 4.5f, 3), // 8

        new_obs(3.2f, 7.5f, 4), // 9
    };
    const size_t n_obs = sizeof(obs) / sizeof(observation);

    size_t indices[] = { 1, 5 };
    const size_t n_indices = sizeof(indices) / sizeof(size_t);

    track t1(1, 4, indices, indices + n_indices, obs, obs + n_obs, 1.f), t2(t1), t3;

    ok1(t1.first_time_stamp() == 1);
    ok1(t1.last_time_stamp() == 4);
    ok1(t1.duration() == 3);
    ok1(t1.dynamic_drag() == 1.f);

    ok1(t1.size() == 2);

    ok1(*t1.begin() == obs[1]);
    ok1(*(--t1.end()) == obs[5]);

    ok1(x(t1.min_location()) == 1.f);
    ok1(y(t1.min_location()) == 1.6f);
    ok1(x(t1.max_location()) == 2.1f);
    ok1(y(t1.max_location()) == 3.0f);

    ok1(within_bounds(t1,point_2d(1.5f, 2.5f)));
    ok1(!within_bounds(t1,point_2d(8.5f, 2.5f)));
    ok1(!within_bounds(t1,point_2d(0.5f, 2.5f)));
    ok1(!within_bounds(t1,point_2d(1.5f, 1.5f)));
    ok1(!within_bounds(t1,point_2d(1.5f, 9.5f)));

    ok1(within_bounds(t1,point_2d(1.5f, 2.5f), 1));
    ok1(!within_bounds(t1,point_2d(8.5f, 2.5f), 1));
    ok1(!within_bounds(t1,point_2d(0.5f, 2.5f), 1));
    ok1(!within_bounds(t1,point_2d(1.5f, 1.5f), 1));
    ok1(!within_bounds(t1,point_2d(1.5f, 9.5f), 1));

    ok1(within_bounds(t1,point_2d(1.5f, 2.5f), 2));
    ok1(!within_bounds(t1,point_2d(8.5f, 2.5f), 2));
    ok1(!within_bounds(t1,point_2d(0.5f, 2.5f), 2));
    ok1(!within_bounds(t1,point_2d(1.5f, 1.5f), 2));
    ok1(!within_bounds(t1,point_2d(1.5f, 9.5f), 2));

    ok1(!within_bounds(t1,point_2d(1.5f, 2.5f), 0));
    ok1(!within_bounds(t1,point_2d(1.5f, 2.5f), 41));

    ok1(within_bounds(t1,new_obs(1.5f, 2.5f, 1)));
    ok1(!within_bounds(t1,new_obs(8.5f, 2.5f, 1)));
    ok1(!within_bounds(t1,new_obs(0.5f, 2.5f, 1)));
    ok1(!within_bounds(t1,new_obs(1.5f, 1.5f, 1)));
    ok1(!within_bounds(t1,new_obs(1.5f, 9.5f, 1)));

    ok1(within_bounds(t1,new_obs(1.5f, 2.5f, 2)));
    ok1(!within_bounds(t1,new_obs(8.5f, 2.5f, 2)));
    ok1(!within_bounds(t1,new_obs(0.5f, 2.5f, 2)));
    ok1(!within_bounds(t1,new_obs(1.5f, 1.5f, 2)));
    ok1(!within_bounds(t1,new_obs(1.5f, 9.5f, 2)));

    ok1(!within_bounds(t1,new_obs(1.5f, 2.5f, 0)));
    ok1(!within_bounds(t1,new_obs(1.5f, 2.5f, 41)));

    ok1(t1.size() == 2);

    ok1(t2.first_time_stamp() == 1);
    ok1(t2.last_time_stamp() == 4);
    ok1(t2.duration() == 3);

    ok1(t2.size() == 2);

    ok1(*t2.begin() == obs[1]);
    ok1(*(--t2.end()) == obs[5]);

    ok1(x(t2.min_location()) == 1.f);
    ok1(y(t2.min_location()) == 1.6f);
    ok1(x(t2.max_location()) == 2.1f);
    ok1(y(t2.max_location()) == 3.0f);

    ok1(within_bounds(t2,point_2d(1.5f, 2.5f)));
    ok1(!within_bounds(t2,point_2d(8.5f, 2.5f)));
    ok1(!within_bounds(t2,point_2d(0.5f, 2.5f)));
    ok1(!within_bounds(t2,point_2d(1.5f, 1.5f)));
    ok1(!within_bounds(t2,point_2d(1.5f, 9.5f)));

    ok1(within_bounds(t2,point_2d(1.5f, 2.5f), 1));
    ok1(!within_bounds(t2,point_2d(8.5f, 2.5f), 1));
    ok1(!within_bounds(t2,point_2d(0.5f, 2.5f), 1));
    ok1(!within_bounds(t2,point_2d(1.5f, 1.5f), 1));
    ok1(!within_bounds(t2,point_2d(1.5f, 9.5f), 1));

    ok1(within_bounds(t2,point_2d(1.5f, 2.5f), 2));
    ok1(!within_bounds(t2,point_2d(8.5f, 2.5f), 2));
    ok1(!within_bounds(t2,point_2d(0.5f, 2.5f), 2));
    ok1(!within_bounds(t2,point_2d(1.5f, 1.5f), 2));
    ok1(!within_bounds(t2,point_2d(1.5f, 9.5f), 2));

    ok1(!within_bounds(t2,point_2d(1.5f, 2.5f), 0));
    ok1(!within_bounds(t2,point_2d(1.5f, 2.5f), 41));

    ok1(within_bounds(t2,new_obs(1.5f, 2.5f, 1)));
    ok1(!within_bounds(t2,new_obs(8.5f, 2.5f, 1)));
    ok1(!within_bounds(t2,new_obs(0.5f, 2.5f, 1)));
    ok1(!within_bounds(t2,new_obs(1.5f, 1.5f, 1)));
    ok1(!within_bounds(t2,new_obs(1.5f, 9.5f, 1)));

    ok1(within_bounds(t2,new_obs(1.5f, 2.5f, 2)));
    ok1(!within_bounds(t2,new_obs(8.5f, 2.5f, 2)));
    ok1(!within_bounds(t2,new_obs(0.5f, 2.5f, 2)));
    ok1(!within_bounds(t2,new_obs(1.5f, 1.5f, 2)));
    ok1(!within_bounds(t2,new_obs(1.5f, 9.5f, 2)));

    ok1(!within_bounds(t2,new_obs(1.5f, 2.5f, 0)));
    ok1(!within_bounds(t2,new_obs(1.5f, 2.5f, 41)));

    ok1(t2.size() == 2);

    ok1(t3.first_time_stamp() == 0);
    ok1(t3.last_time_stamp() == 0);
    ok1(t3.duration() == 0);

    ok1(t3.size() == 0);

    ok1(!within_bounds(t3,point_2d(0.f, 0.f)));
    ok1(!within_bounds(t3,point_2d(0.f, 0.f), 0));
    ok1(!within_bounds(t3,new_obs(0.f, 0.f, 0)));

    // test merging

    size_t indices_2[] = { 8, 9 };
    const size_t n_indices_2 = sizeof(indices_2) / sizeof(size_t);

    track t4(3, 10, indices_2, indices_2 + n_indices_2, obs, obs + n_obs, 1.f);

    track t5(t4,t1,1.f), t6(t1,t4,1.f), t6a(t1,t4,2.f);
    ok1(t5 == t6);
    ok1(t5 != t4);
    ok1(t5 != t1);
    ok1(t6 != t6a);

    ok1(t5.first_time_stamp() == 1);
    ok1(t5.last_time_stamp() == 10);
    ok1(t5.min_location() == point_2d(0.2f, 1.6f));
    ok1(t5.max_location() == point_2d(3.2f, 7.5f));

    track t7(t1,t3,1.f), t8(t3,t1,1.f);
    ok1(t7 == t8);
    ok1(t7 != t1);
    ok1(t7 != t3);

    bool caught(false);
    try
    {
        track t8(t1,t1,1.f);
    }
    catch(std::runtime_error& e)
    {
        caught = true;
    }
    ok(caught, "merging overlapping tracks fails");

    track t9(t8, 5.f);
    ok1(t9 != t8);
    ok1(t9.dynamic_drag() == 5.f);
}


void test_pair_collection()
{
    typedef std::pair<int, int> int_pair;

    int_pair input_pairs[] =
    {
        std::make_pair( 1, 1 ),
        std::make_pair( 2, 1 ),
        std::make_pair( 2, 2 ),
        std::make_pair( 1, 2 ), // duplicate of 1
        std::make_pair( 3, 4 ),
        std::make_pair( 5, 4 ),
        std::make_pair( 1, 4 ),
        std::make_pair( 4, 7 ),
        std::make_pair( 7, 2 ),
    };
    size_t n_pairs = sizeof(input_pairs) / sizeof(int_pair);
    ok1(n_pairs > 0);

    detail::pair_collection<int> pairs(input_pairs, input_pairs + n_pairs);
    diag("Pair count: %zu", pairs.size());
    ok1(pairs.size() == n_pairs - 1); // -1 to account for duplicate

    BOOST_FOREACH(const detail::pair_collection<int>::value_pair& p, pairs)
    {
        diag("pair collection contains (%i, %i)", p.first, p.second);
    }

    ok1(!pairs.contains_pair_with_element(8));
    pairs.insert(std::make_pair( 4, 8 ));
    ok1(pairs.size() == n_pairs);
    ok1(pairs.contains_pair_with_element(8));

    ok1(!pairs.contains_pair_with_element(0));
    ok1(pairs.contains_pair_with_element(1));
    ok1(pairs.contains_pair_with_element(2));
    ok1(pairs.contains_pair_with_element(3));
    ok1(pairs.contains_pair_with_element(4));
    ok1(pairs.contains_pair_with_element(5));
    ok1(!pairs.contains_pair_with_element(6));
    ok1(pairs.contains_pair_with_element(7));
    ok1(pairs.contains_pair_with_element(8));

    pairs.erase_pairs_with_element(4);
    ok1(pairs.size() == n_pairs - 5);

    ok1(!pairs.contains_pair_with_element(0));
    ok1(pairs.contains_pair_with_element(1));
    ok1(pairs.contains_pair_with_element(2));
    ok1(!pairs.contains_pair_with_element(3));
    ok1(!pairs.contains_pair_with_element(4));
    ok1(!pairs.contains_pair_with_element(5));
    ok1(!pairs.contains_pair_with_element(6));
    ok1(pairs.contains_pair_with_element(7));
    ok1(!pairs.contains_pair_with_element(8));

    BOOST_FOREACH(const detail::pair_collection<int>::value_pair& p, pairs)
    {
        diag("pair collection contains (%i, %i)", p.first, p.second);
    }
}

int main(int argc, char** argv)
{
    plan_tests(10 + 15 + 102 + 24);

    test_point_2d();
    test_observation();
    test_track();
    test_pair_collection();

    return exit_status();
}
