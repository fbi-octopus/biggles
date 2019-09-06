#include "detail/physics.hpp"
#include "detail/fun.hpp"
#include "tools/debug.hpp"

#include "kalman_filter.hpp"
#include "model.hpp"
#include "observation.hpp"
#include "track.hpp"

#include <cmath>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <iostream>
#include <limits>
#include <utility>

#include <Eigen/Dense>

namespace biggles
{

void track::update_boundary_cache()
{
    // if there are no observations in the track, set the bounding box to an extremal value.
    if(observations_.empty())
    {
        boundary_min_location_ = point_2d(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
        boundary_max_location_ = point_2d(std::numeric_limits<float>::min(), std::numeric_limits<float>::min());
        return;
    }

    // initialise the bounding box to the first location
    const observation& first_obs(*begin());
    boundary_min_location_ = point_2d(x(first_obs), y(first_obs));
    boundary_max_location_ = point_2d(x(first_obs), y(first_obs));

    BOOST_FOREACH(const observation& o, std::make_pair(begin(), end()))
    {
        // check observation fits into birth and death time
        if (not (t(o) < first_time_stamp() + duration())) {
            std::cout << "t(o) =               " << t(o) << std::endl;
            std::cout << "first_time_stamp() = " << first_time_stamp() << std::endl;
            std::cout << "duration() =         " << duration() << std::endl;
        }
        BOOST_ASSERT(t(o) < first_time_stamp() + duration());
        BOOST_ASSERT(t(o) >= first_time_stamp());

        // update minimum bound
        x(boundary_min_location_) = std::min(x(boundary_min_location_), x(o));
        y(boundary_min_location_) = std::min(y(boundary_min_location_), y(o));

        // update maximum bound
        x(boundary_max_location_) = std::max(x(boundary_max_location_), x(o));
        y(boundary_max_location_) = std::max(y(boundary_max_location_), y(o));
    }
}

void track::check_common_timestamps()
{
    if(size() == 0)
        return;

    time_stamp last_ts = t(*begin());

    BOOST_FOREACH(const observation &o, std::make_pair(++begin(), end()))
    {
        if(t(o) == last_ts)
            throw std::runtime_error("track::check_common_timestamps: at least one pair of observations share a common time stamp.");
        last_ts = t(o);
    }
}

float track::log_posterior(const model::parameters& parameters) const
{
    const Eigen::Matrix2f& R(model::observation_error_covariance(parameters));

    if (not is_symmetric(R)) {
        std::cout << measure_asymmetry(R) << std::endl;
        throw std::runtime_error("R is not symmetric");
    }

    if((cached_posterior_.first == R) && (cached_posterior_.second < 0.f)) {
        return cached_posterior_.second;
    }


    float log_pdf = 0.f;


    // create a kalman filter for the track
    boost::shared_ptr<const kalman_filter> p_kf(make_kalman_filter(parameters));
    //kalman_filter::states_and_cov_deque smoothed_states_and_covariances;
    //rts_smooth(*p_kf, std::front_inserter(smoothed_states_and_covariances));

    // early out if no states
    if(p_kf->predictions().empty()) {
        return log_pdf;
    }

    // the observation matrix
    Eigen::Matrix<float, 2, 4> B;
    B << 1, 0, 0, 0,
         0, 0, 1, 0;

    // iterate over all state predictions
    time_stamp time_stamp(first_time_stamp());
    track::const_iterator track_obs_it(begin());
    const float normalisation(-logf(2.f * M_PI));
    BOOST_FOREACH(const kalman_filter::state_covariance_pair& state_and_cov, p_kf->predictions())
    {
        BOOST_ASSERT(time_stamp < last_time_stamp());

        // do we have an observation at this time stamp?
        if((track_obs_it != end()) && (t(*track_obs_it) == time_stamp))
        {
            // ... yes

            // find associated observation
            const observation& o(*track_obs_it);
            ++track_obs_it;

            /*
            if (state_and_cov.second.determinant()<=0) {
                std::cerr << "................" << std::endl;
                OK(size());
                OK(duration());
            }
            */

            // estimate covariance of pdf
            Eigen::Matrix2f covariance = B * state_and_cov.second * B.transpose() + R;

            // calculate predicted observation
            Eigen::Vector2f predicted_obs = B * state_and_cov.first;

            // actual observation
            Eigen::Vector2f obs;
            obs[0] = x(o);
            obs[1] = y(o);

            // this is just a multivariate Gaussian log pdf
            Eigen::Vector2f delta = obs - predicted_obs;
            float a = delta.transpose() * covariance.inverse() * delta;
            BOOST_ASSERT(std::isfinite(a));
            log_pdf += -0.5f * a;
            if (not (covariance.determinant()>0)) {
                std::cerr << "====================" << std::endl;
                OK(covariance);
                OK(R);
                OK(state_and_cov.second);
                OK(state_and_cov.second.determinant());
            }
            //BOOST_ASSERT(covariance.determinant()>0);
            log_pdf += -0.5f * logf(covariance.determinant());
            log_pdf += normalisation;
        }

        ++time_stamp;
        //std::cout << boost::format("  ** %.7f") % log_pdf << std::endl;
    }
    BOOST_ASSERT(track_obs_it == end());

    // exponential distribution over maximum offset
    //float lambda = 1.f / model::constraint_radius(parameters);
    //BOOST_ASSERT(std::isfinite(lambda));
    BOOST_ASSERT(std::isfinite(radius()));
    //log_pdf += logf(lambda) - lambda * radius();

    cached_posterior_ = std::make_pair(R, log_pdf);
    //std::cout << boost::format("  ** %.7f") % log_pdf << std::endl;
    return log_pdf;
}

bool verify_track(const track& track)
{
    observation_collection::const_iterator previous(track.observations().begin());

    if(track.observations().end() == previous)
    {
        std::cerr << "track has no observations.\n";
        return false;
    }

    BOOST_FOREACH(const observation& obs,
                  std::make_pair(++(track.observations().begin()), track.observations().end()))
    {
        if(t(obs) <= t(*previous))
        {
            std::cerr << "track has non monotonic time stamps.\n";
            return false;
        }

        time_stamp delta_t = t(obs) - t(*previous);
        BOOST_ASSERT(delta_t > 0);
        /*
        if(delta_t >= detail::last_delta_t_)
        {
            std::cerr << "track has too large a gap in time (" << delta_t << ").\n";
            return false;
        }
        */

        if(not observations_see_each_other(obs, *previous, detail::speed_of_light_))
        {
            float dx(x(*previous) - x(obs));
            float dy(y(*previous) - y(obs));
            float ds_sq = dx*dx + dy*dy;
            std::cerr << "track's displacement is too large (" << sqrt(ds_sq) << ") for time lag " << delta_t << ".\n";
            std::cerr << "Time points: " << t(*previous) << ", " << t(obs) << std::endl;
            time_stamp ts=track.first_time_stamp();
            std::cerr << boost::format("duration: [%d, %d]") % ts % track.last_time_stamp() << std::endl;
            std::cerr << "num obs =" << track.size() << std::endl;
            std::cerr << boost::format("obs range: [%d, %d]")
                % track.observations().first_time_stamp()
                % track.observations().last_time_stamp()
                << std::endl;
            for (; ts < track.last_time_stamp(); ++ts) {
                observation_collection::const_range obs_at_t(track.observations().at_time_stamp(ts));
                if (obs_at_t.first != obs_at_t.second) {
                    while (obs_at_t.first != obs_at_t.second) {
                        observation obs = *obs_at_t.first++;
                        std::cout << boost::format("[%f, %f, %d]") % x(obs) % y(obs) % t(obs) << std::endl;
                    }
                }
                else {
                    std::cout << boost::format("[-, -, %d]") % ts << std::endl;
                }
            }
            return false;
        }

        ++previous;
    }

    return true;
}

float track::radius() const
{
    if(empty())
        return 0.f;

    const point_2d& min(min_location());
    const point_2d& max(max_location());
    float dx = x(min) - x(max);
    float dy = y(min) - y(max);
    return sqrtf(dx*dx + dy*dy);
}

float dist(const track& t1, const track& t2) {
    float result = 0.f;
    if ( x(t1.min_location()) > x(t2.max_location()) ) result += square(x(t1.min_location()) - x(t2.max_location()));
    if ( x(t2.min_location()) > x(t1.max_location()) ) result += square(x(t2.min_location()) - x(t1.max_location()));
    if ( y(t1.min_location()) > y(t2.max_location()) ) result += square(y(t1.min_location()) - y(t2.max_location()));
    if ( y(t2.min_location()) > y(t1.max_location()) ) result += square(y(t2.min_location()) - y(t1.max_location()));
    return sqrtf(result);
}

boost::shared_ptr<const kalman_filter> track::make_kalman_filter(const model::parameters& p) const
{
    return boost::shared_ptr<kalman_filter>(new kalman_filter(*this,
        model::observation_error_covariance(p), model::process_noise_covariance(p)) );
}

}
