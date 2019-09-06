#ifndef WITHIN_BIGGLES_KALMAN_FILTER_HPP__
#error "This file should only be included by kalman_filter.hpp"
#endif // WITHIN_BIGGLES_KALMAN_FILTER_HPP__

#include <boost/assert.hpp>
#include <boost/foreach.hpp>
#include <cmath>
#include <deque>
#include <Eigen/Dense>

#include "detail/physics.hpp"
#include "tools/debug.hpp"

namespace biggles
{

template<typename OutputIterator>
void rts_smooth(const kalman_filter& filter,
                OutputIterator output_reversed_states_and_covariances)
{
    // Rauch-Tung-Striebel smoother:
    //
    // L_t = Sigma_{t|t} A^T P_{t+1|t}^{-1}
    //
    // mu_{t|T} = mu_{t|t} + L_t(mu_{t+1|T} - mu_{t+1|t})
    //
    // P_{t|T} = P_{t|t} + L_t(P_{t+1|T} - P_{t+1|t})L^T_t

    // the state evolution matrix: X_{k+1} = A X_k + <noise>
    Eigen::Matrix4f A(Eigen::Matrix4f::Identity());
    BOOST_ASSERT(filter.dynamic_drag() == 1);
    A(0,1) = A(2,3) = filter.dynamic_drag(); // <- off-diagonal velocity integration

    const kalman_filter::states_and_cov_deque& correction_states_and_covs(filter.corrections());
    const kalman_filter::states_and_cov_deque& prediction_states_and_covs(filter.predictions());

    // final interpolated state is output of forward filter, this is mu_{T|T} and Sigma_{T|T}
    if (not is_symmetric(correction_states_and_covs.back().second) and
            measure_asymmetry(correction_states_and_covs.back().second) > 0.f) {
        OK("foo");
        OK1(correction_states_and_covs.back().second);
        OK1(correction_states_and_covs.back().first);
        BOOST_ASSERT(measure_asymmetry(correction_states_and_covs.back().second)==0.f);
    }
    *output_reversed_states_and_covariances = correction_states_and_covs.back();
    ++output_reversed_states_and_covariances;

    // iterator pointing to mu_{t+1|t}, starting at t = T-1
    kalman_filter::states_and_cov_deque::const_reverse_iterator pred_it(prediction_states_and_covs.rbegin());

    // iterator pointing to mu_{t|t}, starting at t = T-1
    kalman_filter::states_and_cov_deque::const_reverse_iterator corr_it(++correction_states_and_covs.rbegin());

    for(kalman_filter::state_covariance_pair prior_pair(correction_states_and_covs.back());
        corr_it != correction_states_and_covs.rend();
        ++pred_it, ++corr_it, ++output_reversed_states_and_covariances)
    {
        // mu_{t+1|t} and Sigma_{t+1|t}
        const Eigen::Vector4f& pred_state(pred_it->first);
        const Eigen::Matrix4f& pred_cov(pred_it->second);

        // mu_{t|t} and Sigma_{t|t}
        const Eigen::Vector4f& corr_state(corr_it->first);
        const Eigen::Matrix4f& corr_cov(corr_it->second);

        // mu_{t+1|T} and Sigma_{t+1|T}
        const Eigen::Vector4f& last_smoothed_state(prior_pair.first);
        const Eigen::Matrix4f& last_smoothed_cov(prior_pair.second);

        // calculate the RTS smoothed values
        Eigen::Matrix4f L = corr_cov * A.transpose() * pred_cov.inverse();
        Eigen::Vector4f smoothed_state = corr_state + L * (last_smoothed_state - pred_state);
        Eigen::Matrix4f smoothed_cov = corr_cov + L * (last_smoothed_cov - pred_cov) * L.transpose();
        smoothed_cov = enforce_symmetry(smoothed_cov);

        if (not is_symmetric(corr_cov) and measure_asymmetry(corr_cov) > 0.f) {
            OK("rts_smooth");
            OK1(corr_cov);
            OK(measure_asymmetry(corr_cov));
            BOOST_ASSERT(measure_asymmetry(corr_cov)==0.f);
        }

        if (not is_symmetric(smoothed_cov) and measure_asymmetry(smoothed_cov) > 0.f) {
            OK("rts_smooth");
            OK1(smoothed_cov);
            OK(measure_asymmetry(corr_cov));
            OK(measure_asymmetry(smoothed_cov));
            BOOST_ASSERT(measure_asymmetry(smoothed_cov)==0.f);
        }

        // record smoothed state and covariance
        prior_pair = kalman_filter::state_covariance_pair(smoothed_state, smoothed_cov);
        *output_reversed_states_and_covariances = prior_pair;
    }
}

}
