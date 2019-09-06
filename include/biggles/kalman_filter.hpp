/// @file kalman_filter.hpp Kalman filter and Rauch–Tung–Striebel smoother for track observations

#ifndef BIGGLES_KALMAN_FILTER_HPP__
#define BIGGLES_KALMAN_FILTER_HPP__

#include <Eigen/Dense>
#include <deque>

#include "observation.hpp"
#include "track.hpp"
#include "detail/physics.hpp"

namespace biggles
{

/// @brief A Kalman filter for predicting missing observations in tracks.
///
/// A Kalman filter is the optimal linear state space evaluator for tracking a hidden state given noisy observations.
/// Specifically we assume some hidden state vector at time \f$ T \f$, \f$ x_t \f$ which evolves with the following
/// model:
///
/// \f[
/// x_t = A x_{t-1} + V_t
/// \f]
///
/// where \f$ A \f$ is some state evolution matrix and \f$ V_t \f$ is a zero-mean Gaussian process with (known)
/// covariance matrix \f$ Q \f$.
///
/// We assume that we can make a noisy observation, \f$ y_t \f$:
///
/// \f[
/// y_t = B x_t + W_t
/// \f]
///
/// where \f$ B \f$ is an observation matrix and \f$ W_t \f$ is a zero-mean Gaussian process with (known) covariance
/// matrix \f$ R \f$.
///
/// A Kalman filter consists of a prediction step followed by an update step. In the prediction step, the current
/// estimate of the hidden state, \f$ \hat{x}_t \f$, is evolved via the state evolution matrix. In the update step, this
/// estimate is refined by fusing any observations made. In this implementation the lack of an observation causes this
/// update step to be skipped and the refined estimate is assumed to be equal to the prediction.
///
/// The prediction step is represented by the following recurrence relations:
///
/// \f[
/// \hat{x}_{t|t-1} = A \hat{x}_{t-1|t-1}, \quad P_{t|t-1} = A P_{t-1|t-1} A^T + Q
/// \f]
///
/// where \f$ P_{t|t-1} \f$ and \f$ P_{t|t} \f$ are, respectively, our prediction of the state estimation error and our
/// refined prediction of the state estimation error.
///
/// The update step is represented by the following recurrence relations:
///
/// \f[
/// \hat{x}_{t|t} = \hat{x}_{t|t-1} + K_t \tilde{z}_t, \quad P_{t|t} = (I - K_t B) P_{t|t-1}
/// \f]
///
/// where
///
/// \f[
/// \tilde{z}_t = y_t - B \hat{x}_{t|t-1}, \quad K_t = P_{t|t-1} B^T + S_t^{-1}, \quad S_t = B P_{t|t-1} B^T + R.
/// \f]
///
/// We initialise the filter by choosing some arbitrary initial state estimate, \f$ \hat{x}_{0|0} \f$, and setting the
/// initial state covariance matrix, \f$ P_{0|0} \f$, to some sufficiently large multiple of \f$ I \f$ so as to specify
/// almost no certainty on the initial estimate. We use the recurrence relations to compute estimates of states and
/// estimation error covariances up until the last time stamp for the track.
///
/// In our system, the state is a position and instantaneous velocity:
///
/// \f[
/// x_t = [ x, x', y, y' ]^T
/// \f]
///
/// States evolve using first order dynamics and the velocity is hidden:
///
/// \f[
/// A = \left[
/// \begin{array}{cccc}
/// 1 & d & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & d \\ 0 & 0 & 0 & 1
/// \end{array}
/// \right], \quad
///
/// B = \left[
/// \begin{array}{cccc}
/// 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0
/// \end{array}
/// \right].
/// \f]
///
/// The matrix \f$ R \f$ and dynamic drag \f$ d \f$ are given as parameters to the constructor. The matrix Q is set to
/// the following by default:
///
/// \f[
/// Q = \left[
/// \begin{array}{cccc}
/// 0.8^2 & 0 & 0 & 0 \\ 0 & 0.2^2 & 0 & 0 \\ 0 & 0 & 0.8^2 & 0 \\ 0 & 0 & 0 & 0.2^2
/// \end{array}
/// \right].
/// \f]
///
/// @sa rts_smooth()
/// @sa http://en.wikipedia.org/wiki/Kalman_filter
/// @sa http://automation.berkeley.edu/resources/KalmanSmoothing.ppt
class kalman_filter
{
public:
    /// @brief The underlying state vector type.
    ///
    /// The state vector represents the instantaneous position and velocity of the molecule as a vector \f$ X \equiv [x,
    /// x', y, y'] \f$.
    typedef Eigen::Vector4f state_vector;

    /// @brief The covariance of the state.
    ///
    /// The Kalman filter maintains an estimate of the instantaneous error in state estimation as a state covariance
    /// matrix \f$ \Sigma \equiv E(XX^T) - E(X)E(X)^T \f$.
    typedef Eigen::Matrix4f covariance_matrix;

    /// @brief A pair holding an interpolated state vector and its associated covariance.
    typedef std::pair<state_vector, covariance_matrix> state_covariance_pair;

    /// @brief The collection type used to hold states and covariances.
    //typedef std::deque<state_covariance_pair, Eigen::aligned_allocator<state_covariance_pair> > states_and_cov_deque;
    typedef std::deque<state_covariance_pair> states_and_cov_deque;

    /// @brief The default observation covariance.
    ///
    /// The default value is
    /// \f[
    /// R = \left[
    /// \begin{array}{cc}
    /// 0.1^2 & 0 \\ 0 & 0.1^2
    /// \end{array}
    /// \right]
    /// \f]
    static const Eigen::Matrix2f default_observation_covariance;

    /// @brief Default constructor.
    kalman_filter() : dynamic_drag_(1.f) { }

    /// @brief Initialise filter from a track's observations.
    ///
    /// This constructor uses the observations from a track to interpolate states for all timestamps within a track.
    ///
    /// @param t The track to initialise from.
    /// @param observation_covariance The observation covariance matrix \f$ R \f$. The default is default_observation_covariance.
    kalman_filter(const track& t, const matrix2f &R, const matrix4f &Q)
    {
        reinitialise(t.first_time_stamp(), t.last_time_stamp(), t.begin(), t.end(), R, Q, t.dynamic_drag());
    }

    /// @brief Initialise filter with a set of observations.
    ///
    /// @tparam InputIterator An InputIterator yielding biggles::observation instances.
    /// @param first_time_stamp The first time stamp of the track.
    /// @param last_time_stamp The time stamp immediately after the last time stamp of the track: \p first_time_stamp +
    /// duration.
    /// @param first The first observation for the track.
    /// @param last Just beyond the last observation for the track.
    /// @param observation_covariance The observation covariance matrix \f$ R \f$.
    /// @param dynamic_drag The dynamic drag factor.
    //template<typename InputIterator>
    kalman_filter(time_stamp first_time_stamp,
                  time_stamp last_time_stamp,
                  track::const_iterator first, track::const_iterator last,
                  const matrix2f &R, const matrix4f &Q,
                  float dynamic_drag)
    {
        reinitialise(first_time_stamp, last_time_stamp, first, last, R, Q, dynamic_drag);
    }

    /// @brief Copy constructor.
    ///
    /// @param kf The biggles::kalman_filter instance to copy.
    kalman_filter(const kalman_filter& kf)
        : prediction_states_and_covs_(kf.prediction_states_and_covs_)
        , correction_states_and_covs_(kf.correction_states_and_covs_)
        , dynamic_drag_(kf.dynamic_drag_)
    { }

    /// @brief Assignment operator.
    ///
    /// @param kf The biggles::kalman_filter instance to copy.
    const kalman_filter& operator = (const kalman_filter& kf)
    {
        prediction_states_and_covs_ = kf.prediction_states_and_covs_;
        correction_states_and_covs_ = kf.correction_states_and_covs_;
        dynamic_drag_ = kf.dynamic_drag_;
        return *this;
    }

    /// @brief Re-initialise filter with a set of observations.
    ///
    /// @tparam InputIterator An InputIterator yielding biggles::observation instances.
    /// @param first_time_stamp The first time stamp of the track.
    /// @param last_time_stamp The time stamp immediately after the last time stamp of the track: \p first_time_stamp +
    /// duration.
    /// @param first The first observation for the track.
    /// @param last Just beyond the last observation for the track.
    /// @param observation_covariance The observation covariance matrix \f$ R \f$.
    /// @param dynamic_drag The dynamic drag factor.
    //template<typename InputIterator>
    void reinitialise(time_stamp first_time_stamp, time_stamp last_time_stamp,
                      track::const_iterator first, track::const_iterator last,
                      const matrix2f &R, const matrix4f &Q,
                      float dynamic_drag);

    /// @brief The predicted states and covariances.
    ///
    /// A collection of predicted states, \f$ x_{t|t-1} \f$, and covariances, \f$ P_{t|t-1} \f$, for all time stamps
    /// covered by the input range.
    const states_and_cov_deque& predictions() const { return prediction_states_and_covs_; }

    /// @brief The corrected states and covariances.
    ///
    /// A collection of corrected states, \f$ x_{t|t} \f$, and covariances, \f$ P_{t|t} \f$, for all time stamps
    /// covered by the input range.
    const states_and_cov_deque& corrections() const { return correction_states_and_covs_; }

    /// @brief The dynamic drag factor associated with this filter.
    float dynamic_drag() const { return dynamic_drag_; }

protected:
    // NOTE: The situation with storing Eigen dense matrices in a std::vector is complex. To avoid this, we use a deque
    // here even though our usage pattern would suggest a vector be more appropriate.
    //
    // See: http://eigen.tuxfamily.org/dox-devel/TopicStlContainers.html

    /// @brief The collection of interpolated states and covariances of their errors.
    states_and_cov_deque prediction_states_and_covs_;

    /// @brief The collection of interpolated states and covariances of their errors.
    states_and_cov_deque correction_states_and_covs_;

    /// @brief The dynamic drag factor to use in state evolution matrices.
    float dynamic_drag_;


};

/// @brief Perform a Rauch–Tung–Striebel backwards smoothing step on the Kalman filter predicted and corrected states.
///
/// @note Since this is a <em>backward</em> process, the results are written to \p
/// output_reversed_states_and_covariances in <em>reverse</em> order; i.e. the smoothed state and covariance for the
/// <em>last</em> time stamp is the first written out.
///
/// Once the forward prediction-update step has been completed for a Kalman filter, we can use a Rauch–Tung–Striebel
/// smoother to refine our earlier state estimates. This is a backwards step which starts from the final estimated state
/// (i.e. the one which has been influenced by all observed observations) and works backwards creating optimal estimates
/// of the hidden state, \f$ \hat{x}_{t|T} \f$, and estimation error covariance, \f$ P_{t|T} \f$. Note that these
/// estimates have been computed given all observations.
///
/// The estimates are computed via the following recurrence relations:
///
/// \f[
/// \hat{x}_{t|T} = \hat{x}_{t|t} + L_t ( \hat{x}_{t+1|T} - \hat{x}_{t+1|t} ), \quad
/// P_{t|T} = P_{t|t} + L_t ( P_{t+1|T} - P_{t+1|t} ) L^T_t
/// \f]
///
/// where \f$ L_t = P_{t|t} A^T P_{t+1|t}^{-1} \f$.
///
/// These smoothed estimates of state and estimation error generated by this function may be used as mean and
/// covariances of a multi-variate Gaussian in order to sample possible state-space configurations for a track. Biggles
/// uses these estimates not only as mean and covariances for evaluation log-likelihoods on track configurations but
/// also for sampling missing data from tracks.
///
/// @sa biggles::kalman_filter
///
/// @tparam OutputIterator Where entries of type kalman_filter::state_covariance_pair are written.
/// @param kalman A Kalman filter to take predicted and corrected states and covariances from.
/// @param output_reversed_states_and_covariances An iterator to write smoothed states and covariances to <em>in reverse
/// order</em>.
template<typename OutputIterator>
void rts_smooth(const kalman_filter& kalman,
                OutputIterator output_reversed_states_and_covariances);

}

#define WITHIN_BIGGLES_KALMAN_FILTER_HPP__
#include "kalman_filter.tcc"
#undef WITHIN_BIGGLES_KALMAN_FILTER_HPP__

#endif // BIGGLES_KALMAN_FILTER_HPP__
