#include "biggles/observation.hpp"
#include "biggles/observation_collection.hpp"
#include "biggles/detail/physics.hpp"

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

bool test_min_ok_collection() {
    observation_collection oc;
    oc.insert(new_obs(0.f, 0.f, 3));
    oc.insert(new_obs(0.f, 5.f, 1));
    diag("obs see each other ? %d", int(observations_see_each_other(oc.front(), oc.back(), detail::speed_of_light_)));
    diag("obs see each other ? %d", int(observations_see_each_other(oc.back(), oc.front(), detail::speed_of_light_)));
    diag("obs in lc ? %d", int(observation_is_within_light_cone(oc.front(), oc.back(), detail::speed_of_light_,
        std::min(t(oc.front()), t(oc.back())), std::max(t(oc.front()), t(oc.back())) )));
    diag("obs in lc ? %d", int(observation_is_within_light_cone(oc.back(), oc.front(), detail::speed_of_light_,
        std::min(t(oc.front()), t(oc.back())), std::max(t(oc.front()), t(oc.back())) )));
    diag("obs in lc ? %d", int(observation_is_within_light_cone(oc.back(), oc.front(), detail::speed_of_light_,
        std::min(t(oc.front()), t(oc.back())), std::max(t(oc.front()), t(oc.back()))+1 )));
    return observation_collection_does_allow_tracks(oc, detail::speed_of_light_);
}

bool test_two_obs_at_same_time() {
    observation_collection oc;
    oc.insert(new_obs(0.f, 0.f, 1));
    oc.insert(new_obs(0.f, 1.f, 1));
    return not observation_collection_does_allow_tracks(oc, detail::speed_of_light_);
}

bool test_obs_too_far_away() {
    observation_collection oc;
    oc.insert(new_obs(0.f, 0.f, 1));
    oc.insert(new_obs(2.f, 2.5f, 2));
    return not observation_collection_does_allow_tracks(oc, detail::speed_of_light_);
}

int main(int argc, char** argv)
{
    plan_tests(26 + 24 + 5 + 3 - 5 -5 - 6  - 4);

    // create empty obs. collection
    observation_collection empty_oc, another_empty_oc;

    ok1(empty_oc.empty());
    ok1(empty_oc.size() == 0);
    ok1(empty_oc.begin() == empty_oc.end());
    ok1(empty_oc.count_at_time_stamp(0) == 0);
    ok1(empty_oc.count_at_time_stamp(1) == 0);
    ok1(empty_oc == another_empty_oc);
    ok1(!(empty_oc != another_empty_oc));

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
    observation_collection oc1(obs_array, obs_array + obs_array_size);

    // construct via insert
    observation_collection oc2;
    for(size_t i=0; i<obs_array_size; ++i)
    {
        oc2.insert(obs_array[i]);
    }

    // XXX: Is this actually guaranteed?
    ok1(oc1 == oc2);
    ok1(!(oc1 != oc2));

    oc2.insert(new_obs(0.f, 0.f, 5));

    // XXX: Is this actually guaranteed?
    // observation_collection is now map
    // ok1(oc1 != oc2);
    // ok1(!(oc1 == oc2));

    /* the following doesn't work as such anymore (tests - 5)
    ok1(oc2.contains(new_obs(1.f, 1.5f, 7)));
    ok1(!oc2.erase(new_obs(1.f, 1.6f, 7)));
    ok1(oc2.erase(new_obs(1.f, 1.5f, 7)));
    ok1(!oc2.contains(new_obs(1.f, 1.5f, 7)));
    ok1(!oc2.erase(new_obs(1.f, 1.5f, 7)));
    */

    oc2.clear();
    ok1(oc2.empty());
    ok1(oc2.size() == 0);
    ok1(oc2 == empty_oc);

    // check assignment
    another_empty_oc = oc1;
    ok1(another_empty_oc == oc1);
    ok1(oc1 == another_empty_oc);
    ok1(!(another_empty_oc != oc1));
    ok1(!(oc1 != another_empty_oc));

    // check size is as expected
    //ok1(oc1.size() == obs_array_size);
    ok1(!oc1.empty());

    ok1(oc1.count_at_time_stamp(0) == 0);
    ok1(oc1.count_at_time_stamp(1) == 0);
    ok1(oc1.count_at_time_stamp(2) == 0);
    ok1(oc1.count_at_time_stamp(3) == 1);
    ok1(oc1.count_at_time_stamp(4) == 0);
    ok1(oc1.count_at_time_stamp(5) == 1);
    ok1(oc1.count_at_time_stamp(6) == 1);
    ok1(oc1.count_at_time_stamp(7) == 1);

    /* comparision with a new obs does not work (tests - 5)
    ok1(oc1.contains(new_obs(-1.f, 3.f, 3)));
    ok1(!oc1.contains(new_obs(-1.1f, 3.f, 3)));
    ok1(!oc1.contains(new_obs(-1.f, 3.1f, 3)));
    ok1(!oc1.contains(new_obs(-1.f, 3.f, 4)));
    ok1(!oc1.contains(new_obs(-1.f, 3.f, 5)));
    */

    /* comparision with a new obs does not work (tests - 6)
    ok1(*oc1.find(new_obs(4.f, 1.5f, 7)) == new_obs(4.f, 1.5f, 7));
    ok1(*oc1.find(new_obs(4.f, 1.5f, 7)) != new_obs(4.1f, 1.5f, 7));
    ok1(*oc1.find(new_obs(4.f, 1.5f, 7)) != new_obs(4.f, 1.6f, 7));
    ok1(*oc1.find(new_obs(4.f, 1.5f, 7)) != new_obs(4.f, 1.5f, 8));
    ok1(*oc1.find(new_obs(4.f, 1.5f, 7)) != new_obs(2.f, 1.5f, 7));

    ok1(oc1.find(new_obs(4.f, 1.6f, 7)) == oc1.end());
    */

    // check ordering and iterator decrement/increment
    bool is_ordered = true;
    bool decrement_it_worked = true;
    size_t n_looked_at = 0;
    for(observation_collection::iterator it(++oc1.begin()), prev_it(oc1.begin()); it != oc1.end(); ++it, ++prev_it)
    {
        observation_collection::iterator prev_it_2(it);
        --prev_it_2;

        decrement_it_worked = decrement_it_worked && (*prev_it_2 == *prev_it) && (prev_it_2 == prev_it);

        is_ordered = is_ordered && (t(*it) >= t(*prev_it));

        ++n_looked_at;
    }
    ok1(n_looked_at > 0);
    //ok1(n_looked_at + 1 == obs_array_size);
    ok1(is_ordered);
    ok1(decrement_it_worked);

    // light cones
    observation_collection light_cone_oc1, light_cone_oc2;
    std::set<time_stamp> time_stamp_with_obs;

    // check corner-case of empty collection
    observations_within_light_cone(empty_oc, new_obs(2.5f, 2.5f, 5), 1.f, 3, 8,
                                   time_stamp_with_obs,
                                   std::inserter(light_cone_oc2, light_cone_oc2.begin()));
    ok1(light_cone_oc2.empty());

    observations_within_light_cone(oc1, new_obs(2.5f, 2.5f, 5), 1.f, 3, 8,
                                   time_stamp_with_obs,
                                   std::inserter(light_cone_oc1, light_cone_oc1.begin()));

    // should have some results, but not all of oc1
    ok1(!light_cone_oc1.empty());
    ok1(light_cone_oc1.size() != oc1.size());

    // check temporal bound
    ok1(light_cone_oc1.first_time_stamp() >= 3);
    ok1(light_cone_oc1.last_time_stamp() <= 8);

    ok1(test_min_ok_collection());
    ok1(test_obs_too_far_away());
    ok1(test_two_obs_at_same_time());

    return exit_status();
}
