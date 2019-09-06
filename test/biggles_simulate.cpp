#include "biggles/kalman_filter.hpp"
#include "biggles/observation.hpp"
#include "biggles/partition.hpp"
#include "biggles/samplers.hpp"
#include "biggles/simulate.hpp"
#include "biggles/track.hpp"
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <cmath>
#include <ctime>
#include <inttypes.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <utility>
#include <vector>

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

static boost::mt19937 random_generator(std::time(0));

using namespace biggles;

void simulate_one_track(const model::parameters& params)
{
    boost::shared_ptr<track> t_ptr(simulate::generate_track(params, point_2d(10,10), point_2d(), 30));
    boost::shared_ptr<track> t2_ptr(simulate::generate_track(params, point_2d(50,10), point_2d(), t_ptr->last_time_stamp() + 10));

    BOOST_FOREACH(const observation& o, t2_ptr->observations())
    {
        t_ptr->insert(o);
    }

    ok1(t_ptr->size() > 0);

    // save track out to file. This file can be loaded by graph from GNU plotutils
    std::ofstream tx("biggles_simulate_track_tx.txt");
    std::ofstream ty("biggles_simulate_track_ty.txt");

    tx << "#m=-1,S=4\n";
    ty << "#m=-1,S=4\n";

    BOOST_FOREACH(const observation& o, *t_ptr)
    {
        tx << std::setw(16) << t(o) << ' '
           << std::setw(16) << x(o) << " 0" << '\n';
        ty << std::setw(16) << t(o) << ' '
           << std::setw(16) << y(o) << " 0" << '\n';
    }

    // mark end of this dataset
    tx << '\n';
    ty << '\n';

    tx << "#m=2,S=2\n";
    ty << "#m=2,S=2\n";

    Eigen::Matrix2f R(Eigen::Matrix2f::Identity() * 0.1f * 0.1f);
    matrix4f Q(detail::initQ());

    // interpolate missing observations
    kalman_filter kf(*t_ptr, R, Q);
    kalman_filter::states_and_cov_deque smoothed_states_and_covariances;

    // use a _front_ inserter since we generate results in reverse order
    rts_smooth(kf, std::front_inserter(smoothed_states_and_covariances));

    smoothed_states_and_covariances = kf.predictions();

    // the observation matrix
    Eigen::Matrix<float, 2, 4> B;
    B << 1, 0, 0, 0,
         0, 0, 1, 0;

    // save interpolated track out to file
    time_stamp time_stamp(t_ptr->first_time_stamp());
    track::const_iterator track_obs_it(t_ptr->begin());
    const float normalisation(-logf(2.f * M_PI));
    BOOST_FOREACH(const kalman_filter::state_covariance_pair& state_and_covariance,
                  smoothed_states_and_covariances)
    {
        // estimate covariance of pdf
        Eigen::Matrix2f covariance = B * state_and_covariance.second * B.transpose() + R;

        // calculate predicted observation
        Eigen::Vector2f predicted_obs = B * state_and_covariance.first;

        tx << std::setw(16) << time_stamp << ' '
           << std::setw(16) << predicted_obs[0] << ' '
           << std::setw(16) << sqrt(covariance(0,0)) << '\n';

        ty << std::setw(16) << time_stamp << ' '
           << std::setw(16) << predicted_obs[1] << ' '
           << std::setw(16) << sqrt(covariance(1,1)) << '\n';

        // do we have an observation at this time stamp?
        if((track_obs_it != t_ptr->end()) && (t(*track_obs_it) == time_stamp))
        {
            // ... yes

            // find associated observation
            const observation& o(*track_obs_it);
            ++track_obs_it;

            // estimate covariance of pdf
            Eigen::Matrix2f covariance = B * state_and_covariance.second * B.transpose() + R;

            // calculate predicted observation
            Eigen::Vector2f predicted_obs = B * state_and_covariance.first;

            // actual observation
            Eigen::Vector2f obs;
            obs[0] = x(o);
            obs[1] = y(o);

            // this is just a multivariate Gaussian log pdf
            Eigen::Vector2f delta = obs - predicted_obs;
            float a = delta.transpose() * covariance.inverse() * delta;
            float log_pdf = 0.f;
            log_pdf += -0.5f * a;
            log_pdf += -0.5f * logf(covariance.determinant());
            log_pdf += normalisation;

            std::cout << "# " << time_stamp << ' ' << log_pdf << '\n';
        }

        ++time_stamp;
    }
}

void simulate_partition(const model::parameters& params)
{
    partition_ptr_t p;

    // simulate full data
    simulate::generate_partition(p, params, 0, 128 + 128, 128.f, 128.f);

    diag("Generated %zu new tracks.", p->tracks().size());
    diag("Generated %zu clutter observations.", p->clutter().size());

    // crop partition to a specific region
    simulate::crop_partition(p, p, 0.f, 128.f, 0.f, 128.f, 64, 64 + 128);

    const track_collection& simulated_tracks(p->tracks());
    const clutter_t& simulated_clutter(p->clutter());

    diag("Cropped to %zu tracks.", simulated_tracks.size());
    diag("Cropped to %zu clutter observations.", simulated_clutter.size());

    // create a synthetic dataset from the partition
    partition_ptr_t partition_ptr;
    simulate::demote_all_tracks_from_partition(partition_ptr, p);

    // check observations count
    size_t expected_count = p->clutter().size();
    BOOST_FOREACH(const boost::shared_ptr<const track>& p_t, p->tracks())
    {
        expected_count += p_t->size();
    }

    ok1(partition_ptr->clutter().size() == expected_count);
    ok1(partition_ptr->tracks().size() == 0);

}

int main(int argc, char** argv)
{
    plan_tests(1 + 2);

    float sigma_R = 0.1f;

    // initialise model parameters
    model::parameters params;
    model::mean_new_tracks_per_frame(params)           = 1.5f;
    model::mean_false_observations_per_frame(params)   = 2.f;
    model::frame_to_frame_survival_probability(params) = 0.95f;
    model::generate_observation_probability(params)    = 0.8f;
    model::observation_error_covariance(params)       << sigma_R * sigma_R, 0, 0, sigma_R * sigma_R;

    simulate_one_track(params);
    simulate_partition(params);

    return exit_status();
}
