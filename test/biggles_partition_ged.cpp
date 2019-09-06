#include "biggles/observation.hpp"
#include "biggles/detail/random.hpp"
#include "biggles/mh_moves/mh_moves.hpp"
#include "biggles/mh_moves/utility.hpp"
#include "biggles/simulate.hpp"
#include <boost/foreach.hpp>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <string>
#include "biggles/tools/sundries.hpp"
#include "biggles/tools/partition_graph_edit_distance.hpp"

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

void test1() {
    std::deque<observation> obs_list;
    obs_list.push_back(new_obs(1, 1, 0));
    obs_list.push_back(new_obs(1, 1, 1));
    obs_list.push_back(new_obs(1, 1, 2));
    obs_list.push_back(new_obs(1, 1, 3));
    obs_list.push_back(new_obs(1, 1, 4));
    clutter_ptr clutter1_p(new clutter_t());
    clutter1_p->insert(obs_list.begin(), obs_list.end());
    boost::shared_ptr<track_collection> tracks1_p(new track_collection());
    clutter_ptr clutter2_p(new clutter_t());
    boost::shared_ptr<track_collection> tracks2_p(new track_collection());
    tracks2_p->insert(shared_track_ptr(new track(0, 5, obs_list.begin(), obs_list.end())));
    partition part1(tracks1_p, clutter1_p);
    partition part2(tracks2_p, clutter2_p);
    ok(partition_ged(part1, part2) == 4.f, "dist(part1, part2) == 4");
    ok(partition_ged(part1, part1) == 0.f, "Partition 1 is equal to itself");
    ok(partition_ged(part2, part2) == 0.f, "Partition 2 is equal to itself");
}

void test2() {
    std::deque<observation> obs_list1;
    obs_list1.push_back(new_obs(1, 1, 0));
    obs_list1.push_back(new_obs(1, 1, 1));
    obs_list1.push_back(new_obs(1, 1, 2));
    obs_list1.push_back(new_obs(1, 1, 3));
    obs_list1.push_back(new_obs(1, 1, 4));
    std::deque<observation> obs_list2;
    obs_list2.push_back(new_obs(0, 1, 0));
    obs_list2.push_back(new_obs(0, 1, 1));
    obs_list2.push_back(new_obs(0, 1, 2));
    obs_list2.push_back(new_obs(0, 1, 3));
    obs_list2.push_back(new_obs(0, 1, 4));

    clutter_ptr clutter1_p(new clutter_t());
    boost::shared_ptr<track_collection> tracks1_p(new track_collection());
    clutter_ptr clutter2_p(new clutter_t());
    boost::shared_ptr<track_collection> tracks2_p(new track_collection());
    track t11;
    track t12;
    track t21;
    track t22;
    t11.insert(obs_list1[0]);
    t11.insert(obs_list1[1]);
    t11.insert(obs_list1[2]);
    t11.insert(obs_list2[3]);
    t11.insert(obs_list2[4]);

    t12.insert(obs_list2[0]);
    t12.insert(obs_list2[1]);
    t12.insert(obs_list2[2]);
    t12.insert(obs_list1[3]);
    t12.insert(obs_list1[4]);

    t21.insert(obs_list1[0]);
    t21.insert(obs_list1[1]);
    t21.insert(obs_list1[2]);
    t21.insert(obs_list1[3]);
    t21.insert(obs_list1[4]);

    t22.insert(obs_list2[0]);
    t22.insert(obs_list2[1]);
    t22.insert(obs_list2[2]);
    t22.insert(obs_list2[3]);
    t22.insert(obs_list2[4]);

    tracks1_p->insert(t11);
    tracks1_p->insert(t12);

    tracks2_p->insert(t21);
    tracks2_p->insert(t22);

    partition part1(tracks1_p, clutter1_p);
    partition part2(tracks2_p, clutter2_p);

    diag("distance = %.1f", partition_ged(part1, part2));
    ok(partition_ged(part1, part2) == 4.f, "dist(part1, part2) == 4");
}

int main(int argc, char** argv) {

    general_paras gen_pars;
    gen_pars.total = 200000;
    gen_pars.seed = 0;
    gen_pars.rgp.lambda = 2.5;
    gen_pars.rgp.p_no = 0.2;
    gen_pars.rgp.p_yes = 0.2;
    gen_pars.rgp.p_tr = 0.5;
    gen_pars.rgp.min_tracks = 1;
    gen_pars.rgp.max_tracks = 2;
    if (not parse_args(argc, argv, gen_pars))
        return exit_status();
    biggles::detail::seed_prng(gen_pars.seed);

    plan_no_plan();

    //test1();
    test2();

    diag("random seed = 0x%x", gen_pars.seed);
    return exit_status();
}
