#include "biggles/observation.hpp"
#include "biggles/detail/random.hpp"
#include "biggles/mh_moves/mh_moves.hpp"
#include "biggles/mh_moves/utility.hpp"
#include "biggles/simulate.hpp"
#include "biggles/partition.hpp"
#include "biggles/samplers.hpp"
#include "biggles/partition_sampler.hpp"
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <Eigen/Dense>

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;
using boost::format;

int main(int argc, char** argv)
{
    // seed the PRNG to ensure consistent results across runs.
    biggles::detail::seed_prng(0xfacedead);

    //bool verbose = (argc > 1) && (0 == strcmp(argv[1], "-v"));

    plan_tests(2);

    boost::shared_ptr<observation_collection> obs_col(new observation_collection());
    obs_col->insert(new_obs(0, 0, 0));
    obs_col->insert(new_obs(0, 1, 0));
    obs_col->insert(new_obs(1, 0, 0));
    obs_col->insert(new_obs(0, 0, 1));
    obs_col->insert(new_obs(0, 0.5, 1));
    obs_col->insert(new_obs(0.5, 0, 1));
    obs_col->insert(new_obs(0, 0, 2));
    obs_col->insert(new_obs(0, 1, 2));
    obs_col->insert(new_obs(1, 0, 2));
    obs_col->insert(new_obs(0, 0, 3));
    obs_col->insert(new_obs(0, 0.5, 3));
    obs_col->insert(new_obs(0.5, 0, 3));
    obs_col->insert(new_obs(0, 0, 4));
    obs_col->insert(new_obs(0, 1, 4));
    obs_col->insert(new_obs(1, 0, 4));
    obs_col->insert(new_obs(0, 0.3, 5));
    obs_col->insert(new_obs(0.3, 0, 5));

    clutter_ptr clutter(new clutter_t());

    boost::shared_ptr<track_collection> tracks(new track_collection());
    clutter->insert(obs_col->begin(), obs_col->end());

    Eigen::Matrix2f obs_var;
    obs_var << 0.01, 0.0, 0.0, 0.01;
    model::parameters input_paras;
    model::birth_rate(input_paras) = 1.0f;
    model::clutter_rate(input_paras) = 0.5f;
    model::survival_probability(input_paras) = 0.6f;
    model::observation_probability(input_paras) = 0.9f;
    model::observation_error_covariance(input_paras) = obs_var;
    model::constraint_radius(input_paras) = 5.0f;
    model::process_noise_covariance(input_paras) = detail::initQ();
    partition_ptr_t original_part_ptr(
        new partition(tracks, clutter, expansion_from_observations(*tracks, *clutter)));
    original_part_ptr->set_expansion(0.f, 2.f, 0.f, 2.f);
    partition_sampler test_sampler(original_part_ptr, input_paras);
    partition_ptr_t current_part_ptr;
    diag("original partition volume = %.1f", original_part_ptr->volume());

    std::map<int, std::string> move_name;
    std::map<int, int> move_count;
    std::map<int, int> error_move_count;
    std::map<int, int> accepted_count;
    typedef std::map<int, int>::const_iterator move_map_iter;
    move_name[0] = "none";
    move_name[1] = "birth";
    move_name[2] = "death";
    move_name[3] = "extend";
    move_name[4] = "reduce";
    move_name[5] = "split";
    move_name[6] = "merge";
    move_name[7] = "update";
    move_name[8] = "identity";
    move_name[9] = "what? 1";

    int new_inf(0);
    int old_inf(0);
    int qq_inf(0);
    int total(10000);
    float best = test_sampler.current_sample_log_density();

    diag("last log density = %.2f", test_sampler.current_sample_log_density());


    for (int i = 0; i < total; ++i) {
        float P_old = test_sampler.current_sample_log_density();
        if (best < P_old) best = P_old;
        current_part_ptr = test_sampler.draw().partition_sample_ptr;
        int proposed_move = test_sampler.proposed_sample().executed_move;
        int last_move = test_sampler.last_move();
        move_count[proposed_move]++;
        accepted_count[last_move]++;
        float P_new = test_sampler.last_log_density();
        float qq = test_sampler.last_pdr();
        /*
        std::cout
            << format("%7.3f < %10.3f - %10.3f %+10.3f = %10.3f (%6s) ? %2s #%2d")
            % alpha % P_new % P_old % qq % mh_ratio
            % move_name[proposed_move] % (alpha < mh_ratio ? "OK" : "-")
            % test_sampler.last_sample().template get<0>().tracks().size()
            << std::endl
        ;
        */
        if (isinff(qq)   ) qq_inf++;
        if (isinff(P_new)) new_inf++;
        if (isinff(P_old)) old_inf++;
        if (isinff(qq) or isinff(P_new) or isinff(P_old)) error_move_count[proposed_move]++;
    }
    model::parameters paras = test_sampler.parameters();
    diag("last volume = %.2f", test_sampler.last_partition()->volume());
    diag("best log-PDF = %.2f", best);
    diag("--------------------");
    diag("birth rate        (b) = %.4f", model::birth_rate(paras));
    diag("clutter rate      (c) = %.4f", model::clutter_rate(paras));
    diag("survival prob.    (s) = %.4f", model::survival_probability(paras));
    diag("observation prob. (o) = %.4f", model::observation_probability(paras));
    diag("--------------------");
    diag("proposed inf  = %d, (%6.2f%%)", new_inf, (float(new_inf)/float(total)*100.0));
    diag("current inf   = %d", old_inf);
    diag("QQ infinity   = %d", qq_inf);
    diag("total samples = %d", total);
    for (move_map_iter iter = error_move_count.begin(); iter != error_move_count.end(); ++iter) {
        diag("%s = %d", move_name.at(iter->first).c_str(), iter->second);
    }
    diag("--------------------");
    diag("proposed moves: ");
    for (int i = 0; i < 8; ++i) {
        diag("%s = %d", move_name.at(i).c_str(), move_count[i]);
    }
    diag("--------------------");
    diag("accepted moves: ");
    for (int i = 0; i < 8; ++i) {
        diag("%s = %d", move_name.at(i).c_str(), accepted_count[i]);
    }

    ok(mh_moves::observations_could_be_neighbours(new_obs(59.3354, 93.3354, 27), new_obs(59.3354, 90.3354, 26)),
        "Obs could be neighbours? ");

    ok1(true);

    return exit_status();
}
