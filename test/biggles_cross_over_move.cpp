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

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

int test_pdr(const int total, const random_gen_paras& params ) {

    partition_ptr_t original_part_ptr(random_debug_partition_adv(params));

    capability_recorder_ptr cap_rec_ptr(new_cap_recorder(*original_part_ptr));
    original_part_ptr->set_capability_recorder(cap_rec_ptr);

    const crossables_t& crossables =  cap_rec_ptr->cross_over_pairs();

    test_t tests;
    BOOST_FOREACH (const std::string& msg,  move_pdr_test_adv(original_part_ptr, "Cross-over", total, tests)) {
        diag("%s", msg.c_str());
    }
    BOOST_FOREACH (const test_t::value_type& item, tests) {
        ok(item.second, "%s", item.first.c_str());
    }
    return 0;
}


int main(int argc, char** argv)
{
    general_paras gen_pars;
    gen_pars.total = 100000;
    gen_pars.seed = 0;
    gen_pars.rgp.lambda = 4.5;
    gen_pars.rgp.p_no = 0.2;
    gen_pars.rgp.p_yes = 0.5;
    gen_pars.rgp.p_tr = 0.5;
    gen_pars.rgp.min_tracks = 5;
    gen_pars.rgp.max_tracks = 6;
    if (not parse_args(argc, argv, gen_pars))
        return exit_status();
    biggles::detail::seed_prng(gen_pars.seed);

    plan_no_plan();

    test_pdr(gen_pars.total, gen_pars.rgp);

    diag("random seed = 0x%x", gen_pars.seed);
    return exit_status();
}
