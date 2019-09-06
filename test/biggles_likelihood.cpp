/// \file biggles_likelihood.cpp Test for p(D|T, theta)

#include "biggles/observation.hpp"
#include "biggles/detail/random.hpp"
#include "biggles/kalman_filter.hpp"
#include "biggles/mh_moves/mh_moves.hpp"
#include "biggles/mh_moves/utility.hpp"
#include "biggles/simulate.hpp"
#include "biggles/partition.hpp"
#include "biggles/samplers.hpp"
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include "biggles/tools/sundries.hpp"

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;
using boost::format;

const float INV_EXP(1.f/expf(1.f));
const float EXP_ONE(expf(1.f));
const float EPSILON(std::numeric_limits<float>::epsilon());

int clutter_likelihood_test() {
    observation_collection clutter_obs;
    clutter_obs.insert(new_obs(0.f, 0.f, 0));
    clutter_obs.insert(new_obs(0.f, 0.f, 1));
    clutter_obs.insert(new_obs(0.f, 0.f, 2));
    clutter_obs.insert(new_obs(0.f, 0.f, 3));
    clutter_obs.insert(new_obs(0.f, 0.f, 4));
    boost::shared_ptr<track_collection> tracks(new track_collection());
    clutter_ptr clutter(new clutter_t(clutter_obs.begin(), clutter_obs.end()));
    partition part(tracks, clutter, partition::expansion_2d(0.f, EXP_ONE, 0.f, 1.f));
    model::parameters paras;
    float log_pdf = model::log_clutter_given_parameters_density(part, paras);
    diag("log PDF = %f, anticipated %f", log_pdf, -5.f);
    ok(std::abs(float(clutter_obs.size()) + log_pdf) < EPSILON, "clutter test 1");
    return 0;
}

int clutter_likelihood_test2() {
    observation_collection clutter_obs;
    clutter_obs.insert(new_obs(1.f, 0.f, 0));
    clutter_obs.insert(new_obs(0.f, 0.f, 0));
    clutter_obs.insert(new_obs(0.f, 0.f, 1));
    clutter_obs.insert(new_obs(0.f, 0.f, 2));
    clutter_obs.insert(new_obs(1.f, 0.f, 2));
    clutter_obs.insert(new_obs(0.f, 1.f, 2));
    clutter_obs.insert(new_obs(0.f, 0.f, 4));
    boost::shared_ptr<track_collection> tracks(new track_collection());
    clutter_ptr clutter(new clutter_t(clutter_obs.begin(), clutter_obs.end()));
    partition part(tracks, clutter, partition::expansion_2d(0.f, 2., 0.f, 2.f));
    model::parameters paras;
    float log_pdf = model::log_clutter_given_parameters_density(part, paras);
    diag("log PDF = %f, anticipated %f", log_pdf, -float(clutter_obs.size()) * logf(4.f));
    ok(std::abs(float(clutter_obs.size()) * logf(4.f) + log_pdf) < EPSILON, "clutter test 2");
    return 0;
}

int partition_likelihood_test1() {
    observation_collection clutter_obs;
    clutter_obs.insert(new_obs(0.f, 0.f, 0));
    clutter_obs.insert(new_obs(0.f, 0.f, 1));
    clutter_obs.insert(new_obs(0.f, 0.f, 2));
    clutter_obs.insert(new_obs(0.f, 0.f, 3));
    clutter_obs.insert(new_obs(0.f, 0.f, 4));
    boost::shared_ptr<track_collection> tracks(new track_collection());
    clutter_ptr clutter(new clutter_t(clutter_obs.begin(), clutter_obs.end()));
    partition part(tracks, clutter, partition::expansion_2d(0.f, 2.718121828f, 0.f, 1.f));
    model::parameters paras;
    model::birth_rate(paras) = 1.f;
    model::clutter_rate(paras) = 1.f;
    model::observation_probability(paras) = 0.5f;
    model::survival_probability(paras) = 0.5f;
    float log_pdf = model::log_partition_given_parameters_density(part, paras);
    diag("log pdf %f, anticipated %f", log_pdf, logf(powf(INV_EXP, 10.f)));
    ok(std::abs(log_pdf + 10.f) < EPSILON, "partition test 1");
    return 0;
}

int partition_likelihood_test2() {
    observation_collection clutter_obs;
    clutter_obs.insert(new_obs(0.f, 0.f, 0));
    clutter_obs.insert(new_obs(0.f, 0.f, 1));
    clutter_obs.insert(new_obs(0.f, 0.f, 2));
    clutter_obs.insert(new_obs(0.f, 0.f, 3));
    clutter_obs.insert(new_obs(0.f, 0.f, 4));
    boost::shared_ptr<track_collection> tracks(new track_collection());
    clutter_ptr clutter(new clutter_t(clutter_obs.begin(), clutter_obs.end()));
    observation_collection obs_col;
    tracks->insert(track(0, 2, obs_col.begin(), obs_col.end(), 1.f));
    tracks->insert(track(1, 3, obs_col.begin(), obs_col.end(), 1.f));
    tracks->insert(track(2, 4, obs_col.begin(), obs_col.end(), 1.f));
    tracks->insert(track(3, 5, obs_col.begin(), obs_col.end(), 1.f));
    partition part(tracks, clutter, partition::expansion_2d(0.f, EXP_ONE, 0.f, 1.f));
    model::parameters paras;
    model::birth_rate(paras) = 1.f;
    model::clutter_rate(paras) = 1.f;
    model::observation_probability(paras) = 0.5f;
    model::survival_probability(paras) = 0.5f;
    float log_pdf = model::log_partition_given_parameters_density(part, paras);
    diag("log pdf %f, anticipated %f", log_pdf, -10.f + 6.f*logf(0.5) + 3.f*logf(0.25));
    ok(std::abs(log_pdf + 10.f - 6.f*logf(0.5) - 3.f*logf(0.25)) < EPSILON, "partition test 2");
    return 0;
}

int partition_likelihood_test3() {
    boost::shared_ptr<track_collection> tracks(new track_collection());
    clutter_ptr clutter(new clutter_t());
    observation_collection obs_col1;
    observation_collection obs_col2;
    observation_collection obs_col3;
    observation_collection obs_col4;
    obs_col1.insert(new_obs(0.f, 0.f, 0));
    obs_col1.insert(new_obs(0.f, 0.f, 2));
    obs_col2.insert(new_obs(0.f, 0.f, 1));
    obs_col2.insert(new_obs(0.f, 0.f, 2));
    obs_col2.insert(new_obs(0.f, 0.f, 3));
    obs_col3.insert(new_obs(0.f, 0.f, 2));
    obs_col3.insert(new_obs(0.f, 0.f, 3));
    obs_col4.insert(new_obs(0.f, 0.f, 2));
    obs_col4.insert(new_obs(0.f, 0.f, 4));
    tracks->insert(track(0, 3, obs_col1.begin(), obs_col1.end(), 1.f));
    tracks->insert(track(1, 5, obs_col2.begin(), obs_col2.end(), 1.f));
    tracks->insert(track(2, 4, obs_col3.begin(), obs_col3.end(), 1.f));
    tracks->insert(track(2, 5, obs_col4.begin(), obs_col4.end(), 1.f));
    partition part(tracks, clutter, partition::expansion_2d(0.f, EXP_ONE, 0.f, 1.f));
    model::parameters paras;
    model::birth_rate(paras) = 1.5f;
    model::clutter_rate(paras) = .8f;
    model::observation_probability(paras) = 0.6f;
    model::survival_probability(paras) = 0.8f;
    float expected = -5*0.8f; // clutter
    expected += - 2.f * 1.0945348918918356f - 1.3822169643436166f - 2.f * 1.5; //birth
    expected += logf(0.8f) + log(0.64f) + logf(0.4096f) + logf(0.384f); // surv
    expected += logf(0.6f) + 2.f*logf(0.48f) + logf(0.1296f) + logf(0.432f); // obs
    float log_pdf = model::log_partition_given_parameters_density(part, paras);
    diag("log pdf %f, anticipated %f", log_pdf, expected);
    ok(std::abs(log_pdf - expected) < EPSILON, "partition test 3");
    return 0;
}

void print_kalman_filter_states(const boost::shared_ptr<const track>& t_ptr, const model::parameters paras) {
    const Eigen::Matrix2f& R(model::observation_error_covariance(paras));
    boost::shared_ptr<const kalman_filter> p_kf(t_ptr->make_kalman_filter(paras));
    Eigen::Matrix<float, 2, 4> B;
    B << 1, 0, 0, 0,
         0, 0, 1, 0;
    BOOST_FOREACH(const kalman_filter::state_covariance_pair& state_and_cov, p_kf->predictions()) {
        Eigen::Matrix2f covariance = B * state_and_cov.second * B.transpose() + R;
        Eigen::Vector2f predicted_obs = B * state_and_cov.first;
        diag("cov [[%f, %f], [%f, %f]]", covariance(0, 0), covariance(0, 1), covariance(1, 0), covariance(1, 1));
        diag("loc [%f, %f]", predicted_obs(0), predicted_obs(1));
    }
};

int observation_likelihood1() {
    model::parameters paras;
    model::observation_error_covariance(paras) = Eigen::Matrix2f::Identity() * 0.09f;
    model::process_noise_covariance(paras) = detail::initQ();
    partition::expansion_2d square(0.f, 2.0, 0.f, 2.0f);
    for (int max_ts = 2; max_ts < 50; ++max_ts) {
        observation_collection obs_col1;
        for (int ts = 0; ts < max_ts; ++ts)
            obs_col1.insert(new_obs(1.f, 1.f, ts));

        boost::shared_ptr<track_collection> tracks(new track_collection());
        clutter_ptr clutter(new clutter_t());
        tracks->insert(track(0, max_ts, obs_col1.begin(), obs_col1.end(), 1.f));
        partition part_tr(tracks, clutter, square);

        boost::shared_ptr<track_collection> tracks2(new track_collection());
        clutter_ptr clutter2(new clutter_t(obs_col1.begin(), obs_col1.end()));
        partition part_cl(tracks2, clutter2, square);

        const boost::shared_ptr<const track>& t_ptr = *part_tr.tracks().begin();
        float log_tr_pdf = model::log_track_given_parameters_density(t_ptr, paras);
        float log_cl_pdf = model::log_clutter_given_parameters_density(part_cl, paras);

        int total = 10000;

        float log_part_pdf_tr = -std::numeric_limits<float>::max();
        partition& part = part_tr;
        for (int i = 0; i < total; ++i) {
            sample_birth_rate_given_partition(part, paras);
            sample_clutter_rate_given_partition(part, paras);
            sample_observation_probability_given_partition(part, paras);
            sample_survival_probability_given_partition(part, paras);
            float log_pdf = model::log_partition_given_parameters_density(part, paras);
            if (log_part_pdf_tr < log_pdf) log_part_pdf_tr = log_pdf;
        }

        float log_part_pdf_cl = -std::numeric_limits<float>::max();
        part = part_cl;
        for (int i = 0; i < total; ++i) {
            sample_birth_rate_given_partition(part, paras);
            sample_clutter_rate_given_partition(part, paras);
            sample_observation_probability_given_partition(part, paras);
            sample_survival_probability_given_partition(part, paras);
            float log_pdf = model::log_partition_given_parameters_density(part, paras);
            if (log_part_pdf_cl < log_pdf) log_part_pdf_cl = log_pdf;
        }
        diag("%d, %f, %f", max_ts, log_part_pdf_tr + log_tr_pdf, log_part_pdf_cl + log_cl_pdf );
    }
    return 0;
}

int main(int argc, char** argv)
{
    general_paras gen_pars;
    gen_pars.total = 100000;
    gen_pars.seed = 0;
    gen_pars.rgp.lambda = 2;
    gen_pars.rgp.p_no = 0.;
    gen_pars.rgp.p_yes = 0.;
    gen_pars.rgp.p_tr = 0.4;
    gen_pars.rgp.min_tracks = 5;
    gen_pars.rgp.max_tracks = 6;
    if (not parse_args(argc, argv, gen_pars))
        return exit_status();
    biggles::detail::seed_prng(gen_pars.seed);

    plan_tests(2);
    clutter_likelihood_test();
    clutter_likelihood_test2();
    /* these have changed
    partition_likelihood_test1();
    partition_likelihood_test2();
    partition_likelihood_test3();
    */
    observation_likelihood1();
    diag("random seed = 0x%x", gen_pars.seed);
    return exit_status();
}
