#include "kalman_filter.hpp"
#include "kalman_filter_cache.hpp"
#include "samplers.hpp"
#include "sampling/simple.hpp"

#include <boost/foreach.hpp>
#include <boost/random.hpp>
//#include <boost/random/exponential_distribution.hpp>
#include <boost/shared_ptr.hpp>
#include <cmath>
#include <ctime>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <limits>
#include <cmath>

extern "C" {
#include "third-party/gen_beta.h"
}

#include "detail/random.hpp"

namespace biggles
{

void sample_birth_rate_given_partition_orig(const partition& part_sample, model::parameters& params) {
    off_t n_total_born(part_sample.tracks().size());
    model::mean_new_tracks_per_frame(params) = sampling::sample_gamma(1.f + n_total_born, 1.f) /
        static_cast<float>(part_sample.duration());
    if (not std::isfinite(model::mean_new_tracks_per_frame(params))) {
        std::cerr << "duration = " << part_sample.duration() << std::endl;
        std::cerr << "number of tracks = " << n_total_born << std::endl;
        std::cerr << "sampled birth rate = " << model::mean_new_tracks_per_frame(params) << std::endl;
        BOOST_ASSERT(std::isfinite(model::mean_new_tracks_per_frame(params)));
    }
    //BOOST_ASSERT(model::mean_new_tracks_per_frame(params)>0.f);
    if (model::mean_new_tracks_per_frame(params) == 0.f)
        model::mean_new_tracks_per_frame(params) = std::numeric_limits<float>::min();
}

void sample_birth_rate_given_partition_alt(const partition& part_sample, model::parameters& params) {
    off_t n_total_born(part_sample.tracks().size());
    const size_t min_surv(1);
    model::mean_new_tracks_per_frame(params) = sampling::sample_gamma(1.f + n_total_born, 1.f) /
        static_cast<float>(part_sample.duration() - min_surv); // no birth in the last frames
    BOOST_ASSERT(std::isfinite(model::mean_new_tracks_per_frame(params)));
    BOOST_ASSERT(model::mean_new_tracks_per_frame(params)>0.f);
}

void sample_birth_rate_given_partition(const partition& part_sample, model::parameters& params) {
    sample_birth_rate_given_partition_orig(part_sample, params);
}

void sample_clutter_rate_given_partition(const partition& part_sample, model::parameters& params) {
    off_t n_total_false(part_sample.clutter().size());
    model::mean_false_observations_per_frame(params) = sampling::sample_gamma(1 + n_total_false, 1.f) /
        static_cast<float>(part_sample.duration());
    if (not std::isfinite(model::mean_false_observations_per_frame(params))) {
        std::cerr << "number of clutter = " << n_total_false << std::endl;
        std::cerr << "partition duration = " << part_sample.duration() << std::endl;
        std::cerr << "sampled clutter rate = " << model::mean_false_observations_per_frame(params) << std::endl;
        throw std::logic_error("clutter rate sampling returned a not finite value");
    }
    if (model::mean_false_observations_per_frame(params) == 0.f)
        model::mean_false_observations_per_frame(params) = std::numeric_limits<float>::min();
}

void sample_observation_probability_given_partition_orig(const partition& part_sample, model::parameters& params) {

    if (part_sample.tracks().empty()) {
        model::generate_observation_probability(params) = sampling::uniform_real(0.f, 1.f);
        return;
    }
    off_t n_total_generated(0);
    off_t sum_of_durations(0);
    // calculate sum_of_durations and n_total_generated
    BOOST_FOREACH(const boost::shared_ptr<const track>& t_ptr, part_sample.tracks())
    {
        sum_of_durations += t_ptr->duration();
        // the track's size is the number of observations within it
        n_total_generated += t_ptr->size();
    }
    //BOOST_ASSERT(n_total_generated > 0); // there would be no tracks otherwise
    // the 1.0001s below should be 1 but the sampling functions for beta distributions hit an infinite loop if the sizes
    // are zero.
    // beta(1,1) is uniform distribution FIXME
    model::generate_observation_probability(params) =
        sampling::sample_beta(1.0001f + n_total_generated, 1.0001f + (sum_of_durations - n_total_generated));

}

void sample_observation_probability_given_partition_alt(const partition& part_sample, model::parameters& params) {

    if (part_sample.tracks().empty()) {
        model::generate_observation_probability(params) = sampling::uniform_real(0.f, 1.f);
        return;
    }
    off_t n_total_generated(0);
    off_t sum_of_durations(0);
    const size_t min_obs(2);
    // calculate sum_of_durations and n_total_generated
    BOOST_FOREACH(const boost::shared_ptr<const track>& t_ptr, part_sample.tracks())
    {
        sum_of_durations += t_ptr->duration() - min_obs;
        n_total_generated += t_ptr->size() - min_obs;
    }
    BOOST_ASSERT(n_total_generated >= 0);
    BOOST_ASSERT(sum_of_durations >= n_total_generated);
    model::generate_observation_probability(params) =
        sampling::sample_beta(1.0001f + n_total_generated, 1.0001f + sum_of_durations - n_total_generated);
}

void sample_observation_probability_given_partition(const partition& part_sample, model::parameters& params) {
    sample_observation_probability_given_partition_orig(part_sample, params);
}


void sample_survival_probability_given_partition_orig(const partition& part_sample, model::parameters& params) {
    if(part_sample.tracks().empty()) {
        model::frame_to_frame_survival_probability(params) = sampling::uniform_real(0.f, 1.f);
        return;
    }

    off_t n_total_died(part_sample.tracks().size()); // everything that has a beginning has an end, Neo. Whatever.
    off_t n_total_survived(0);

    // calculate n_total_survived
    BOOST_FOREACH(const boost::shared_ptr<const track>& t_ptr, part_sample.tracks())
    {
        // a track 'survives' for one fewer ticks than it's duration
        if(t_ptr->duration() >= 1)
            n_total_survived += t_ptr->duration() - 1;
    }
    //BOOST_ASSERT(n_total_survived > 0); // there would be no tracks otherwise
    //BOOST_ASSERT(n_total_died > 0); // there would be no tracks otherwise
    do model::frame_to_frame_survival_probability(params) = sampling::sample_beta(
            1.0001f + n_total_survived, 1.0001f + n_total_died);
    while (model::frame_to_frame_survival_probability(params) == 1.f);
}

void sample_survival_probability_given_partition_alt(const partition& part_sample, model::parameters& params) {
    if(part_sample.tracks().empty()) {
        model::frame_to_frame_survival_probability(params) = sampling::uniform_real(0.f, 1.f);
        return;
    }

    off_t n_total_died(0);
    off_t n_total_survived(0);

    const size_t min_surv(1);

    // calculate n_total_survived
    // a track 'survives' for one fewer ticks than it's duration and one survival is guaranteed
    BOOST_FOREACH(const boost::shared_ptr<const track>& t_ptr, part_sample.tracks()) {
        n_total_survived += t_ptr->duration() - 1 - min_surv ;
        n_total_died += off_t(t_ptr->last_time_stamp() < part_sample.last_time_stamp());
    }
    BOOST_ASSERT(n_total_survived >= 0);
    BOOST_ASSERT(n_total_died >= 0);
    model::frame_to_frame_survival_probability(params) = sampling::sample_beta(1.0001f + n_total_survived, 1.0001f + n_total_died);
}

void sample_survival_probability_given_partition(const partition& part_sample, model::parameters& params) {
    sample_survival_probability_given_partition_orig(part_sample, params);
}

void sample_tracking_control_parameters_given_partition(const partition& partition, model::parameters& params) {
    sample_birth_rate_given_partition(partition, params);
    sample_clutter_rate_given_partition(partition, params);
    sample_observation_probability_given_partition(partition, params);
    sample_survival_probability_given_partition(partition, params);
}

void sample_observation_error_given_partition(const partition& partition, model::parameters& params) {
    // special case: no tracks. sample from a fairly all-encompassing distribution
    if (not is_symmetric(model::observation_error_covariance(params))) {
        std::cout << measure_asymmetry(model::observation_error_covariance(params)) << std::endl;
        throw std::runtime_error("R is not symmetric at the start");
    }

    if(partition.tracks().empty())
    {
        model::observation_error_covariance(params) = sampling::sample_inverse_wishart(Eigen::Matrix2f::Identity() * 2.f, 5);
        //model::observation_error_covariance(params) = Eigen::Matrix2f::Identity() * 0.09f; // FIXME PRIOR calculation
        return;
    }

    if (not is_symmetric(model::observation_error_covariance(params))) {
        std::cout << measure_asymmetry(model::observation_error_covariance(params)) << std::endl;
        throw std::runtime_error("R is not symmetric after track empty.");
    }

    // calculate Wishart Phi parameters for sampling R

    // this is the Wishart prior
    Eigen::Matrix2f Phi = Eigen::Matrix2f::Identity() * 2.f;
    size_t s = 5;

    if (not is_symmetric(Phi)) {
        std::cout << measure_asymmetry(Phi) << std::endl;
        throw std::runtime_error("Phi is not symmetric after init");
    }



    // update parameters for each track
    BOOST_FOREACH(const boost::shared_ptr<const track>& t_ptr, partition.tracks())
    {
        Phi += wishart_phi_parameter_for_track(t_ptr, params);
        s += t_ptr->size(); // number of observations in track
    }

    if (not is_symmetric(Phi)) {
        std::cout << measure_asymmetry(Phi) << std::endl;
        throw std::runtime_error("Phi is not symmetric after update");
    }

    model::observation_error_covariance(params) = sampling::sample_inverse_wishart(Phi, s);
    //model::observation_error_covariance(params) = Eigen::Matrix2f::Identity() * 0.09f; // FIXME PRIOR test
    if (not is_symmetric(model::observation_error_covariance(params))) {
        std::cout << measure_asymmetry(model::observation_error_covariance(params)) << std::endl;
        throw std::runtime_error("R is not symmetric after sampling.");
    }


}

void sample_model_parameters_given_partition(const partition& partition, model::parameters& params) {
    sample_tracking_control_parameters_given_partition(partition, params);
    /*
    model::mean_new_tracks_per_frame(params) = 0.2f; // FIXME PRIOR calculation
    model::mean_false_observations_per_frame(params) = 0.01f; // FIXME PRIOR calculation
    model::generate_observation_probability(params) = .9f; // FIXME PRIOR calculation
    model::frame_to_frame_survival_probability(params) = .9f; // FIXME PRIOR calculation
    */

    // sample constraint radius [removed]

    sample_observation_error_given_partition(partition, params);

}

Eigen::Matrix2f wishart_phi_parameter_for_track(const boost::shared_ptr<const track>& track_p, model::parameters& params)
{
    const track& track(*track_p);
    Eigen::Matrix2f Phi = Eigen::Matrix2f::Zero();
    BOOST_ASSERT(is_symmetric(Phi));

    // create a Kalman filter for the track
    static kalman_filter_cache kf_cache;
    const matrix2f& R(model::observation_error_covariance(params));
    const matrix4f& Q(model::process_noise_covariance(params));
    const kalman_filter& kf(kf_cache.get(track_p, R, Q));

    // ... and a list of smoothed states and covariances
    kalman_filter::states_and_cov_deque smoothed_states_and_covariances;

    // use a _front_ inserter since we generate results in reverse order
    rts_smooth(kf, std::front_inserter(smoothed_states_and_covariances));

    // the observation matrix
    Eigen::Matrix<float, 2, 4> B;
    B << 1, 0, 0, 0,
         0, 0, 1, 0;

    // iterate over state and covariances for the track. keep a track of which time stamp this is
    time_stamp time_stamp(track.first_time_stamp());
    track::const_iterator track_obs_it(track.begin()); // also start iterating over track observations
    BOOST_FOREACH(const kalman_filter::state_covariance_pair& state_and_cov, smoothed_states_and_covariances)
    {
        BOOST_ASSERT(time_stamp < track.last_time_stamp());

        // extract the state and covariance
        const Eigen::Vector4f& state(state_and_cov.first);
        const Eigen::Matrix4f& cov(state_and_cov.second);

        // do we have an associated observation?
        if((track_obs_it != track.end()) && (t(*track_obs_it) == time_stamp))
        {
            // extract observation
            const observation& obs(*track_obs_it);

            // sample state
            Eigen::Vector4f sampled_state(sampling::sample_multivariate_gaussian(state, cov));

            // sampled observation
            Eigen::Vector2f sampled_obs = B * sampled_state;

            // calculate delta
            Eigen::Vector2f obs_vector, delta;
            obs_vector << x(obs), y(obs);
            delta = obs_vector - sampled_obs;

            // update Wishart parameter
            Phi += delta * delta.transpose();
            if (not is_symmetric(Phi)) {
                OK1(cov);
                OK(time_stamp - track.first_time_stamp());
                OK1(Phi);
                OK1(obs_vector);
                OK1(sampled_obs);
                BOOST_ASSERT(is_symmetric(Phi));
            }


            ++track_obs_it;
        }

        ++time_stamp;
    }
    if (not is_symmetric(Phi)) {
        std::cout << "Phi = " << std::endl;
        std::cout << Phi << std::endl;
        throw std::runtime_error("Phi is not symmetric");
    }
    return Phi;
}


} // namespace biggles
