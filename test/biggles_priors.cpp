#include "biggles/model.hpp"
#include "biggles/partition.hpp"
#include "biggles/tracker.hpp"
#include "biggles/simulate.hpp"
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <fstream>
#include <string>
#include <cmath>

#include "biggles/tools/sundries.hpp"

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

const float image_size(32.f);
const time_stamp time_extent(128);

struct part_para {
    float b, c, o, s;
    part_para (float bb, float cc, float oo, float ss) : b(bb), c(cc), o(oo), s(ss) {}
    part_para (const model::parameters& params) :
        b(model::birth_rate(params)),
        c(model::clutter_rate(params)),
        o(model::observation_probability(params)),
        s(model::survival_probability(params)) {}
};

bool is_in_01_incl(float f) { return 0.0f <= f and f <= 1.0f; }
bool is_positive(float f) { return f > 0.0f; }

float quartiles(data_iter b, data_iter e, float &lo, float &hi) {
    std::deque<float> data_copy(b, e);
    std::sort(data_copy.begin(), data_copy.end());
    data_iter median_ptr = data_copy.begin();
    int step = data_copy.size()/4;
    std::advance(median_ptr, step);
    lo = *median_ptr;
    std::advance(median_ptr, step);
    float result = *median_ptr;
    std::advance(median_ptr, step);
    hi = *median_ptr;
    return result;
}

bool dump_data(data_iter b, data_iter e, const std::string& fname) {
    try {
        std::ofstream dumpsite(fname.c_str());
        while (b != e) dumpsite << *b++ << std::endl;
    }
    catch (const std::exception &e) {
        return false;
    }
    return true;
}

int test_inv_wishart() {
    diag("Comparing the implemented (\"here\") Wishart distribution with the R implementation (\"there\")");
    model::parameters params;
    model::mean_new_tracks_per_frame(params)           = 2.f;
    model::mean_false_observations_per_frame(params)   = 1.f;
    model::frame_to_frame_survival_probability(params) = 0.95f;
    model::generate_observation_probability(params)    = 0.9f;
    model::constraint_radius(params) = 100.0f;

    const float eps(0.001f);

    int total = 0;
    int passed = 0;
    model::observation_error_covariance(params) << 1.430304, -1.191502, -1.191502, 1.389336;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (-3.559525)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), -3.559525);
    model::observation_error_covariance(params) << -0.985086068104712, 0.672046045379445, 0.672046045379445, -1.29904989329519;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (2.656294)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), 2.656294);
    model::observation_error_covariance(params) << -0.455529316010436, 0.548656022278254, 0.548656022278254, -1.00701139266364;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (15.80541)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), 15.80541);
    model::observation_error_covariance(params) << 1.70035796238832, -1.78545975647887, -1.78545975647887, 0.244721201367077;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (-4.233217)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), -4.233217);
    model::observation_error_covariance(params) << 1.09814576536207, -0.448974502788745, -0.448974502788745, 0.0357872453122741;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (13.40429)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), 13.40429);
    model::observation_error_covariance(params) << -0.118243745472366, 1.09596329801294, 1.09596329801294, -0.98904929569856;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (-2.20168)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), -2.20168);
    model::observation_error_covariance(params) << -1.95335035648088, -0.436234307572509, -0.436234307572509, 0.235099378404611;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (-1.776402)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), -1.776402);
    model::observation_error_covariance(params) << 0.199379431961018, 1.4705812699427, 1.4705812699427, 0.40286955401663;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (-3.501687)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), -3.501687);
    model::observation_error_covariance(params) << -0.0779657833300163, -1.837210992881, -1.837210992881, 0.532279115899551;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (-5.638954)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), -5.638954);
    model::observation_error_covariance(params) << -0.0732336965241064, -0.524441669960241, -0.524441669960241, -1.5908819356183;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (-3.986849)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), -3.986849);
    model::observation_error_covariance(params) << -0.227371066465159, 1.52811540974702, 1.52811540974702, 0.0398600298840976;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (-4.344813)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), -4.344813);
    model::observation_error_covariance(params) << 1.89778321707126, 0.250303891010732, 0.250303891010732, 1.59230258422261;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (-6.376119)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), -6.376119);
    model::observation_error_covariance(params) << -0.700789642971323, 0.470250080529942, 0.470250080529942, -1.90636169182093;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (1.046791)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), 1.046791);
    model::observation_error_covariance(params) << -0.181573785054159, -0.425369242870216, -0.425369242870216, 0.583612355113527;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (5.538613)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), 5.538613);
    model::observation_error_covariance(params) << 1.11705717417555, 2.10740753511284, 2.10740753511284, 0.0487264847222338;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (-6.505638)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), -6.505638);
    model::observation_error_covariance(params) << 1.3789250772606, 1.61765964610381, 1.61765964610381, -1.40028160559548;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (-6.920236)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), -6.920236);
    model::observation_error_covariance(params) << 1.01052383657972, -1.41049421968031, -1.41049421968031, -0.109930915230516;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (-3.39717)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), -3.39717);
    model::observation_error_covariance(params) << -0.283212686024601, -0.571187607380736, -0.571187607380736, 0.527547374581063;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (2.628803)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), 2.628803);
    model::observation_error_covariance(params) << -1.45309132768864, 0.96601199049902, 0.96601199049902, -0.911329803921612;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (8.944553)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), 8.944553);
    model::observation_error_covariance(params) << 0.457908536414127, 0.619392170299201, 0.619392170299201, -0.441046398340659;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (1.312182)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), 1.312182);
    model::observation_error_covariance(params) << 0.38793441552551, 0.97225238415363, 0.97225238415363, -1.64687763597453;
    passed += int(std::abs(model::log_parameters_prior_density(params) - (-3.491962)) < eps); ++total;
    diag("here %f, there %f", model::log_parameters_prior_density(params), -3.491962);

    //ok(passed == total, "inverse Wishart");
    diag("passed %d, total %d", passed, total);

    return 0;
}


bool generate_partition(const model::parameters& params, partition_ptr_t& partition_ptr, bool shutup=false) {

    // simulate full data
    simulate::generate_partition(partition_ptr, params, 0, (time_extent>>1) + time_extent, image_size * 1.5f, image_size * 1.5f);

    if (not shutup) {
        diag("Generated %zu new tracks.", partition_ptr->tracks().size());
        diag("Generated %zu clutter observations.", partition_ptr->clutter().size());
    }

    return true;
}

bool generate_parameters2(model::parameters& params) {
    model::mean_new_tracks_per_frame(params)           = 1.2f;
    model::mean_false_observations_per_frame(params)   = 0.001f;
    model::frame_to_frame_survival_probability(params) = 0.9;
    model::generate_observation_probability(params)    = 0.9f;
    float sigma_R = 0.3f;
    model::observation_error_covariance(params)       << sigma_R * sigma_R, 0, 0, sigma_R * sigma_R;
    model::constraint_radius(params) = sampling::uniform_real(0, 10);
    return true;
}

void print_mat_for_r(const Eigen::Matrix2f& mat) {
    std::cout << "c(" << mat(0) << ", " << mat(1) << ", " << mat(2) << ", " << mat(3) << ")";
}

int test4(int n_it, const model::parameters& params) {
    std::deque<float> log_priors;
    std::deque<float> log_part;
    std::deque<float> log_clutter;
    std::deque<float> log_tracks;
    for (int i=0 ; i < n_it; ++i) {
        partition_ptr_t original_part_ptr;
        generate_partition(params, original_part_ptr, true);
        log_priors.push_back(model::log_parameters_prior_density(params));
        log_part.push_back(model::log_partition_given_parameters_density(*original_part_ptr, params));
        log_clutter.push_back(model::log_clutter_given_parameters_density(*original_part_ptr, params));
        float log_sum = 0.0f;
        BOOST_FOREACH(const boost::shared_ptr<const track>& t_ptr, original_part_ptr->tracks())
        {
            log_sum += model::log_track_given_parameters_density(t_ptr, params);
        }
        log_tracks.push_back(log_sum);
    }

    float lo, med, hi;
    med = quartiles(log_priors.begin(), log_priors.end(), lo, hi);
    diag("log-Prior.   [%.4f, %.4f]. Median %.4f (%.4f) ",
        *std::min_element(log_priors.begin(), log_priors.end()),
        *std::max_element(log_priors.begin(), log_priors.end()),
        med, hi-lo);
    med = quartiles(log_part.begin(), log_part.end(), lo, hi);
    diag("log-Part.    [%.4f, %.4f]. Median %.4f (%.4f)",
        *std::min_element(log_part.begin(), log_part.end()),
        *std::max_element(log_part.begin(), log_part.end()),
        med, hi-lo);
    med = quartiles(log_clutter.begin(), log_clutter.end(), lo, hi);
    diag("log-Clutter. [%.4f, %.4f]. Median %.4f (%.4f)",
        *std::min_element(log_clutter.begin(), log_clutter.end()),
        *std::max_element(log_clutter.begin(), log_clutter.end()),
        med, hi-lo);
    med = quartiles(log_tracks.begin(), log_tracks.end(), lo, hi);
    diag("log-Tracks.  [%.4f, %.4f]. Median %.4f (%.4f)",
        *std::min_element(log_tracks.begin(), log_tracks.end()),
        *std::max_element(log_tracks.begin(), log_tracks.end()),
        med, hi-lo);


    return 0;
}

int test3(int n_it, const partition_ptr_t& partition_ptr, bool dump) {
    model::parameters target_parameters;
    std::deque<float> birth_rates;
    std::deque<float> clutter_rates;
    std::deque<float> obs_prob;
    std::deque<float> surv_prob;
    std::deque<float> constraint_radius;
    std::deque<float> log_priors;
    std::deque<float> log_part;
    std::deque<float> log_clutter;
    std::deque<float> log_tracks;
    std::deque<float> obs_error;
    generate_parameters2(target_parameters);
    model::process_noise_covariance(target_parameters) = detail::initQ();

    int op_duration(partition_ptr->duration());
    diag("clutter rate = %.2f", float(partition_ptr->clutter().size())/float(op_duration));
    diag("birth rate = %.2f", float(partition_ptr->tracks().size())/float(op_duration));

    int count_symmetry = 0;

    for(int i=0; i<n_it; ++i) {
        sample_model_parameters_given_partition(*partition_ptr, target_parameters);
        part_para ppars(target_parameters);
        birth_rates.push_back(ppars.b);
        clutter_rates.push_back(ppars.c);
        obs_prob.push_back(ppars.o);
        surv_prob.push_back(ppars.s);
        Eigen::Matrix2f& oe_cov = model::observation_error_covariance(target_parameters);
        count_symmetry += int(oe_cov == oe_cov.transpose());
        obs_error.push_back(sqrtf(oe_cov.determinant()));
        constraint_radius.push_back(sqrtf(model::constraint_radius(target_parameters)));
        log_priors.push_back(model::log_parameters_prior_density(target_parameters));
        //if (log_priors.back() > 0.f) {
        //    print_mat_for_r(oe_cov);
        //    std::cout << " p=" << log_priors.back() << std::endl;
        //}
        log_part.push_back(model::log_partition_given_parameters_density(*partition_ptr, target_parameters));
        log_clutter.push_back(model::log_clutter_given_parameters_density(*partition_ptr, target_parameters));
        float log_sum = 0.0f;
        BOOST_FOREACH(const boost::shared_ptr<const track>& t_ptr, partition_ptr->tracks())
        {
            log_sum += model::log_track_given_parameters_density(t_ptr, target_parameters);
        }
        log_tracks.push_back(log_sum);


    }
    ok(std::count_if(obs_prob.begin(), obs_prob.end(), is_in_01_incl) == n_it, "(o) in [0, 1]");
    ok(std::count_if(surv_prob.begin(), surv_prob.end(), is_in_01_incl) == n_it, "(s) in [0, 1]");
    ok(std::count_if(birth_rates.begin(), birth_rates.end(), is_positive) == n_it, "(b) > 0");
    ok(std::count_if(clutter_rates.begin(), clutter_rates.end(), is_positive) == n_it, "(c) > 0");
    ok(std::count_if(obs_error.begin(), obs_error.end(), is_positive) == n_it, "R > 0");
    ok(count_symmetry == n_it, "R always symmetric");

    typedef std::pair<float, float> fpair;

    fpair b_stats = mean_stdev(birth_rates.begin(), birth_rates.end());
    fpair c_stats = mean_stdev(clutter_rates.begin(), clutter_rates.end());
    fpair o_stats = mean_stdev(obs_prob.begin(), obs_prob.end());
    fpair s_stats = mean_stdev(surv_prob.begin(), surv_prob.end());
    fpair r_stats = mean_stdev(obs_error.begin(), obs_error.end());
    fpair t_stats = mean_stdev(constraint_radius.begin(), constraint_radius.end());

    diag("b %.4f +/- %.4f", b_stats.first, b_stats.second);
    diag("c %.4f +/- %.4f", c_stats.first, c_stats.second);
    diag("o %.4f +/- %.4f", o_stats.first, o_stats.second);
    diag("s %.4f +/- %.4f", s_stats.first, s_stats.second);
    diag("R %.4f +/- %.4f", r_stats.first, r_stats.second);
    diag("T %.4f +/- %.4f", t_stats.first, t_stats.second);

    float lo, med, hi;
    med = quartiles(log_priors.begin(), log_priors.end(), lo, hi);
    diag("log-Prior.   [%.4f, %.4f]. Median %.4f (%.4f) ",
        *std::min_element(log_priors.begin(), log_priors.end()),
        *std::max_element(log_priors.begin(), log_priors.end()),
        med, hi-lo);
    med = quartiles(log_part.begin(), log_part.end(), lo, hi);
    diag("log-Part.    [%.4f, %.4f]. Median %.4f (%.4f)",
        *std::min_element(log_part.begin(), log_part.end()),
        *std::max_element(log_part.begin(), log_part.end()),
        med, hi-lo);
    med = quartiles(log_clutter.begin(), log_clutter.end(), lo, hi);
    diag("log-Clutter. [%.4f, %.4f]. Median %.4f (%.4f)",
        *std::min_element(log_clutter.begin(), log_clutter.end()),
        *std::max_element(log_clutter.begin(), log_clutter.end()),
        med, hi-lo);
    med = quartiles(log_tracks.begin(), log_tracks.end(), lo, hi);
    diag("log-Tracks.  [%.4f, %.4f]. Median %.4f (%.4f)",
        *std::min_element(log_tracks.begin(), log_tracks.end()),
        *std::max_element(log_tracks.begin(), log_tracks.end()),
        med, hi-lo);

    if (dump) {
        bool ok0 = true;
        ok0 = ok0 and dump_data(birth_rates.begin(), birth_rates.end(), "/tmp/test_BigglesBirthRates.txt");
        ok0 = ok0 and dump_data(clutter_rates.begin(), clutter_rates.end(), "/tmp/test_BigglesClutterRates.txt");
        ok0 = ok0 and dump_data(surv_prob.begin(), surv_prob.end(), "/tmp/test_BigglesSurvivalProb.txt");
        ok0 = ok0 and dump_data(obs_prob.begin(), obs_prob.end(), "/tmp/test_BigglesObservationProb.txt");
        ok0 = ok0 and dump_data(obs_error.begin(), obs_error.end(), "/tmp/test_BigglesObservationError.txt");
        diag("Dumping data successful ? %d", ok0);
    }

    return 0;

}

partition_ptr_t get_test_partition() {
    observation_collection t1_obs;
    observation_collection t2_obs;
    observation_collection cl_obs;
    t1_obs.insert(new_obs(0, 0, 0));
    //t1_obs.insert(new_obs(0, 0, 1));
    t1_obs.insert(new_obs(0, 0.1, 2));
    t1_obs.insert(new_obs(0.1, 0, 3));
    t1_obs.insert(new_obs(0, 0.2, 4));

    //t2_obs.insert(new_obs(0, 1, 0));
    t2_obs.insert(new_obs(0, 1.1, 1));
    t2_obs.insert(new_obs(0.1, 1, 2));
    //t2_obs.insert(new_obs(0, 1, 3));
    t2_obs.insert(new_obs(0, 1.2, 4));

    cl_obs.insert(new_obs(1.1, 1, 0));
    cl_obs.insert(new_obs(1, 1.1, 1));
    cl_obs.insert(new_obs(1.1, 1.1, 2));
    cl_obs.insert(new_obs(1, 1.2, 3));
    cl_obs.insert(new_obs(1.2, 1, 4));
    track track1(0, 5, t1_obs.begin(), t1_obs.end(), 1.0);
    track track2(1, 5, t2_obs.begin(), t2_obs.end(), 1.0);
    boost::shared_ptr<track_collection> tracks(new track_collection());
    tracks->insert(track1);
    tracks->insert(track2);
    clutter_ptr clutter(new clutter_t(cl_obs.begin(), cl_obs.end()));

    return partition_ptr_t(new partition(tracks, clutter));
}

int main(int argc, char** argv)
{
    plan_no_plan();

    const float sigma_R = 0.1f;

    bool dump = (argc > 1) && (0 == strcmp(argv[1], "-d"));

    //biggles::detail::seed_prng(0xfacedead);

    // initialise model parameters
    model::parameters params;
    model::mean_new_tracks_per_frame(params)           = 2.f;
    model::mean_false_observations_per_frame(params)   = 1.f;
    model::frame_to_frame_survival_probability(params) = 0.95f;
    model::generate_observation_probability(params)    = 0.9f;
    model::observation_error_covariance(params)       << sigma_R * sigma_R, 0, 0, sigma_R * sigma_R;
    model::constraint_radius(params) = 100.0f;
    model::process_noise_covariance(params) = detail::initQ();


    int n_it(50);
    if(argc > 1 and 0 != strcmp(argv[1], "-d"))
        n_it = atoi(argv[1]);

    diag("number of iterations = %d", n_it);

    //test1(n_it, params);
    partition_ptr_t original_part_ptr;
    generate_partition(params, original_part_ptr);
    test3(n_it, original_part_ptr, dump);
    test3(n_it, get_test_partition(), false);
    test_inv_wishart();
    test4(n_it, params);

    return exit_status();
}
