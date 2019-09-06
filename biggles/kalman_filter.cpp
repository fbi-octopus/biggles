#include "kalman_filter.hpp"
#include <Eigen/QR>
namespace biggles
{

matrix4x2f kalman_gain(const matrix2f& R, const matrix4f& Pred, const matrix2x4f& B) {
    matrix4x2f LHS = Pred*B.transpose();
    matrix2f S = B*LHS + R;
    S = enforce_symmetry(S);
    if (not is_symmetric(S)) {
        OK1(S);
        OK1(R);
        OK1(B*LHS);
        BOOST_ASSERT(is_symmetric(S));
    }
    matrix2x4f gain_transposed = S.transpose().colPivHouseholderQr().solve(LHS.transpose());
    return gain_transposed.transpose();
}

//template<typename InputIterator>
void kalman_filter::reinitialise(time_stamp first_time_stamp,
                                 time_stamp last_time_stamp,
                                 track::const_iterator first, track::const_iterator last,
                                 const matrix2f &R, const matrix4f &Q,
                                 float dynamic_drag)
{
    // set the dynamic drag associated with this filter
    dynamic_drag_ = dynamic_drag;

    // remove any existing solution
    prediction_states_and_covs_.clear();
    correction_states_and_covs_.clear();

    BOOST_ASSERT(last_time_stamp >= first_time_stamp);
    BOOST_ASSERT(first != last);

    // shortcut the entire filter if we have nothing to generate
    if(first_time_stamp == last_time_stamp)
        return;

    // the state evolution matrix: X_{k+1} = A X_k + <noise>
    Eigen::Matrix4f A(Eigen::Matrix4f::Identity());
    A(0,1) = A(2,3) = dynamic_drag_; // <- off-diagonal velocity integration


    //// covariance matrix of the state noise
    //Eigen::Matrix4f Q(Eigen::Matrix4f::Zero());
    //Q(0,0) = Q(2,2) = 1e-1f * 1e-1f * detail::speed_of_light_ * detail::speed_of_light_;
    //Q(1,1) = Q(3,3) = 1e-2f * 1e-2f * detail::speed_of_light_ * detail::speed_of_light_;

    // the observation matrix: Y_k = B X_k + <noise>
    Eigen::Matrix<float, 2, 4> B( Eigen::Matrix<float, 2, 4>::Zero());
    B(0,0) = B(1,2) = 1.f;

    // initial state
    Eigen::Vector4f mu;
    mu << 64.f, 0.f, 64.f, 0.f;
    //mu << x(*first), 0.f, y(*first), 0.f;
    //std::cout << mu << std::endl;

    // initial error
    //const float image_size(128.f);
    //const float s2(4.f*image_size*image_size);
    //const float s2(1e-1f * detail::speed_of_light_ * 1e-1f * detail::speed_of_light_);
    //const float s2(powf(2.f, 16.f));
    const float s2(powf(2.f, 16.f));
    Eigen::Matrix4f Sigma(s2 * Eigen::Matrix4f::Identity());
    Sigma(1,1) = Sigma(3,3) = 1e-2f * detail::speed_of_light_ * detail::speed_of_light_;

    // prediction_states_and_covs_ is mu_{t|t-1} and Sigma_{t|t-1}
    // correction_states_and_covs_ is mu_{t|t} and Sigma_{t|t}
    // forward loop
    // BOOST_ASSERT(Q.determinant() > 0.f);
    for(time_stamp time_stamp = first_time_stamp; time_stamp < last_time_stamp; ++time_stamp)
    {
        // prediction
        Eigen::Vector4f predict_mu = A * mu;
        Eigen::Matrix4f predict_Sigma = A * Sigma * A.transpose() + Q;

        // record the predictions mu_{t|t-1} and Sigma_{t|t-1}
        /*
        if (not (predict_Sigma.determinant() > 0)) {
            OK(Sigma.determinant());
            OK(Q.determinant());
        }
        */
        prediction_states_and_covs_.push_back(state_covariance_pair(predict_mu, predict_Sigma));

        // do we have an observation for this time step?
        if((first == last) || (t(*first) != time_stamp))
        {
            // no observation for this time stamp, simply evolve state
            mu = predict_mu;
            Sigma = predict_Sigma;
        }
        else
        {
            // we have an observation, woo! Retrieve it and increment observation iterator.
            const observation& obs(*first);
            ++first;

            // extract observation location
            Eigen::Vector2f obs_vec;
            obs_vec << x(obs), y(obs);

            // update
            matrix4x2f K = kalman_gain(R, predict_Sigma, B);
                //predict_Sigma * B.transpose() * (B * predict_Sigma * B.transpose() + R).inverse();

            mu = predict_mu + K * (obs_vec - B * predict_mu);
            //Eigen::Matrix4f proto_sigma = predict_Sigma - K * B * predict_Sigma;
            Eigen::Matrix4f raw_sigma = predict_Sigma - K * B * predict_Sigma;
            Sigma = enforce_symmetry(raw_sigma);
            if (not is_symmetric(Sigma)) {
                OK1((K * B));
                OK1(K * B * predict_Sigma);
                OK(is_symmetric(K * B * predict_Sigma));
                OK(measure_asymmetry(K * B * predict_Sigma));
            }
        }

        if (not is_symmetric(Sigma) and measure_asymmetry(Sigma) > 0.f) {
            OK("1");
            OK1(Sigma);
            OK(measure_asymmetry(Sigma));
            OK1(predict_Sigma);
            OK1(R);
            BOOST_ASSERT(is_symmetric(Sigma));
        }

        // record the corrections mu_{t|t} and Sigma_{t|t}
        correction_states_and_covs_.push_back(state_covariance_pair(mu, Sigma));
    }

    // check we generated the expected number of states
    BOOST_ASSERT(prediction_states_and_covs_.size() == static_cast<size_t>(last_time_stamp - first_time_stamp));
    BOOST_ASSERT(prediction_states_and_covs_.size() > 0);
    BOOST_ASSERT(correction_states_and_covs_.size() == static_cast<size_t>(last_time_stamp - first_time_stamp));
    BOOST_ASSERT(correction_states_and_covs_.size() > 0);

    // last element of prediction_states_and_covs_ is (mu_{T|T-1}, Sigma_{T|T-1})
    // last element of correction_states_and_covs_ is (mu_{T|T}, Sigma_{T|T})
    // where T is last_time_stamp
}



} // namespace biggles
