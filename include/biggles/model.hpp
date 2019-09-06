/// @file model.hpp Various log-PDF evaluators for parts of the Bayesian model

#ifndef BIGGLES_MODEL_HPP__
#define BIGGLES_MODEL_HPP__

#include <vector>
#include <boost/tuple/tuple.hpp>
#include <Eigen/Dense>

#include "partition.hpp"
#include "types.hpp"

namespace biggles
{

/// @brief Functions and types defining the Biggles dynamic model
namespace model
{

/// @name Access to the tracking control parameters
/// @{
inline const float& birth_rate(const parameters& p) { return p.birth_rate; }
inline float& birth_rate(parameters& p) { return p.birth_rate; }

inline const float& clutter_rate(const parameters& p) { return p.clutter_rate; }
inline float& clutter_rate(parameters& p) { return p.clutter_rate; }

inline const float& survival_probability(const parameters& p) { return p.survival_probability; }
inline float& survival_probability(parameters& p) { return p.survival_probability; }

inline const float& observation_probability(const parameters& p) { return p.observation_probability; }
inline float& observation_probability(parameters& p) { return p.observation_probability; }

inline const matrix4f& process_noise_covariance(const parameters& p) { return p.process_noise_covariance; }
inline matrix4f& process_noise_covariance(parameters& p) { return p.process_noise_covariance; }
/// @}


/// @brief Access the mean number of tracks per-frame parameter from a biggles::parameters tuple.
///
/// @param p
inline const float& mean_new_tracks_per_frame(const parameters& p) { return p.birth_rate; }

/// @brief Access the mean number of tracks per-frame parameter from a biggles::parameters tuple.
///
/// @param p
inline float& mean_new_tracks_per_frame(parameters& p) { return p.birth_rate; }

/// @brief Access the mean number of spurious observations per-frame parameter from a biggles::parameters tuple.
///
/// @param p
inline const float& mean_false_observations_per_frame(const parameters& p) { return p.clutter_rate; }

/// @brief Access the mean number of spurious observations per-frame parameter from a biggles::parameters tuple.
///
/// @param p
inline float& mean_false_observations_per_frame(parameters& p) { return p.clutter_rate; }

/// @brief Access the frame-to-frame survival probability parameter from a biggles::parameters tuple.
///
/// @param p
inline const float& frame_to_frame_survival_probability(const parameters& p) { return p.survival_probability; }

/// @brief Access the frame-to-frame survival probability parameter from a biggles::parameters tuple.
///
/// @param p
inline float& frame_to_frame_survival_probability(parameters& p) { return p.survival_probability; }

/// @brief Access the probability that a track will generate an observation from a biggles::parameters tuple.
///
/// @param p
inline const float& generate_observation_probability(const parameters& p) { return p.observation_probability; }

/// @brief Access the probability that a track will generate an observation from a biggles::parameters tuple.
///
/// @param p
inline float& generate_observation_probability(parameters& p) { return p.observation_probability; }

/// @brief Access the covariance matrix of the observation error from a biggles::parameters tuple.
///
/// @param p
inline const Eigen::Matrix2f& observation_error_covariance(const parameters& p) { return p.observation_error_covariance; }

/// @brief Access the covariance matrix of the observation error from a biggles::parameters tuple.
///
/// @param p
inline Eigen::Matrix2f& observation_error_covariance(parameters& p) { return p.observation_error_covariance; }

/// @brief Access the constraint radius from a biggles::parameters tuple.
///
/// @param p
inline const float& constraint_radius(const parameters& p) { return p.constraint_radius; }

/// @brief Access the constraint radius from a biggles::parameters tuple.
///
/// @param p
inline float& constraint_radius(parameters& p) { return p.constraint_radius; }

/// @brief Calculate the log-pdf of a partition given a set of model parameters, independent of observed data.
///
/// **The following needs to be reviewed. Binomial coefficients need to be included**
///
/// Calculate the value of \f$ \ell(T | \theta) \f$ from the partition and model parameters given. Since the model
/// parameters are independent of one another, we may factorise the pdf as follows:
///
/// \f[
/// P(T | \theta) = P(T | p_s) P(T | p_d) P(T | \lambda_b) P(T | \lambda_f).
/// \f]
///
/// Note that \f$ R \f$ does not appear here because it depends on data; these terms are all the data-independent parts
/// of the posterior on \f$ T \f$. Each of these terms may defined in terms of the following values:
///
/// - \f$ N_t^s \f$: the number of tracks which survive from \f$ t-1 \f$ to \f$ t \f$;
/// - \f$ N_t^d \f$: the number of tracks which <em>did not</em> survive from \f$ t-1 \f$ to \f$ t \f$;
/// - \f$ N_t^b \f$: the number of tracks which newly appeared at \f$ t \f$;
/// - \f$ N_t^o \f$: the number of observations assigned to a track at \f$ t \f$;
/// - \f$ N_t^f \f$: the number of observations deemed spurious at \f$ t \f$,
///
/// where \f$ N_0^s = N_o^b = 0 \f$ by convention. Each of the individual parameter pdf terms can be written down based
/// on the definition of the corresponding parameters:
///
/// - \f$ P(T | p_s) = \prod_{t=1}^K p_s^{N_t^s} (1-p_s)^{N_t^d}; \f$
/// - \f$ P(T | p_d) = \prod_{t=1}^K p_d^{N_t^o} (1-p_d)^{N_t^s + N_t^b - N_t^o}; \f$
/// - \f$ P(T | \lambda_b) = \prod_{t=1}^K \mathcal{P}(N_t^b ; \lambda_b); \f$
/// - \f$ P(T | \lambda_f) = \prod_{t=1}^K \mathcal{P}(N_t^f ; \lambda_f), \f$
///
/// where \f$ \mathcal{P}(x ; \lambda) \f$ is the Poisson pmf with mean \f$ \lambda \f$ evaluated at \f$ x \f$. The
/// corresponding log terms become:
///
/// - \f$ \ell(T | p_s) = \sum_{t=1}^K N_t^s \log(p_s) + N_t^d \log(1-p_s); \f$
/// - \f$ \ell(T | p_d) = \sum_{t=1}^K N_t^o \log(p_d) + (N_t^s + N_t^b - N_t^o) \log(1-p_d); \f$
/// - \f$ \ell(T | \lambda_b) = \sum_{t=1}^K N_t^b \log(\lambda_b) - \lambda_b - \log(\Gamma(1 + N_t^b)); \f$
/// - \f$ \ell(T | \lambda_f) = \sum_{t=1}^K N_t^f \log(\lambda_f) - \lambda_f - \log(\Gamma(1 + N_t^f)); \f$
///
/// @sa biggles::parameters
///
/// @param part A reference to the partition to consider.
/// @param parameters A reference to the model parameters.
///
/// @return The value of \f$ \ell(T | \theta) \f$.
float log_partition_given_parameters_density(const partition& part, const parameters& parameters);

float log_partition_given_observation_prob_density(const partition& part, const parameters& parameters);
float log_partition_given_survival_prob_density(const partition& part, const parameters& parameters);
float log_partition_given_birth_rate_density(const partition& part, const parameters& parameters);
float log_partition_given_clutter_rate_density(const partition& part, const parameters& parameters);

/// @brief Calculate the log-pdf for the track <em>observations</em> given the parameters.
///
/// This function works by using a biggles::kalman_filter to sample missing states for a track and then to calculate
/// the likelihood of the observations we've seen given the parameters.
///
/// Specifically, suppose we have an observation \f$ y \f$, a predicted state, \f$ \hat{x} \f$ and state estimation
/// error estimate, \f$ \hat{P} \f$. Then the total error in the predicted observation, \f$ \hat{y} = B \hat{x} \f$ is
/// given by \f$ \hat{\Sigma} = B \hat{P} B^T + R \f$. We therefore calculate the likelihood of \f$ y \f$ assuming a
/// Gaussian model:
///
/// \f[
/// P(y | \hat{x}, \hat{P}) = \mathcal{N}(y ; \hat{y}, \hat{\Sigma}).
/// \f]
///
/// We combine all these likelihoods for each observation in the track.
///
/// @sa biggles::kalman_filter
///
/// @param track_p
/// @param parameters
///
/// @return The value of \f$ \ell(d_i|t_i, \theta) \f$.
float log_track_given_parameters_density(const boost::shared_ptr<const track>& track_p,
                                         const parameters& parameters);

/// @brief Compute the log-prior on the model parameters.
///
/// The prior on the model parameters is as follows:
///
/// - \f$ \lambda_b \f$ and \f$ \lambda_s \f$ have improper uninformative priors on them being positive.
/// - \f$ p_s \f$ and \f$ p_d \f$ have uniform priors over [0, 1].
/// - \f$ R \f$ has an inverse Wishart prior with parameters \f$ \Phi = 2I \f$ and \f$ s = 5 \f$. These are the same
/// parameters as in sample_parameters_given_partition().
///
/// @sa sample_parameters_given_partition()
///
/// @param parameters The model parameters whose prior should be calculated.
///
/// @return The value of \f$ \ell(\theta) \f$.
float log_parameters_prior_density(const parameters& parameters);

/// @brief Calculate the log pdf of the clutter observations given the model parameters.
///
/// The clutter observations themselves are independent of the model parameter \f$ R \f$ and so their distribution is
/// very simple:
///
/// \f[
/// P(d_0 | t_0, \theta) = \prod_{t=1}^K \left( \frac{1}{V} \right)^{N_t^f}
/// \f]
///
/// where \f$ V \f$ is the total number of pixels in the image.
///
/// @note Currently it is assumed that \f$ V = 256^2 \f$. This is really a bug but for the moment, this is the only
/// place where the absolute image size matters.
///
/// @param clutter
/// @param parameters
///
/// @return The value of \f$ \ell(d_0 | t_0, \theta) \f$.
float log_clutter_given_parameters_density(const partition& part, const parameters& parameters);

/// @brief Calculate the log pdf of the track observations given the model parameters.
///
/// This uses the Kalman filter. track likelihood depends on the model parameter \f$ R \f$ only.
float log_tracks_given_parameters_density(const partition& part, const model::parameters& parameters);

/** \brief returns log( p(data| partition, parameters) )
 *
 *
 */
float log_likelihood(const partition& part_sample, const model::parameters& para_sample);

/// @brief Evaluate the full log pdf for a particular partition
///
/// In the documentation for biggles::partition_sampler, it was shown that
///
/// \f[
/// \ell(T | \theta, D) = \kappa + \ell(d_0|t_0, \theta) + \sum_{i=1}^K \ell(d_i|t_i, \theta) + \ell(T|\theta) +
/// \ell(\theta).
/// \f]
///
/// This function computes this value by calling other log density calculation functions.
///
/// @sa biggles::partition_sampler
/// @sa log_clutter_given_parameters_density()
/// @sa log_track_given_parameters_density()
/// @sa log_partition_given_parameters_density()
/// @sa log_parameters_prior_density()
///
/// @param part
/// @param parameters
///
/// @return The value of \f$ \ell(T | \theta, D) + \kappa \f$ where \f$ \kappa \f$ is some constant offset.
float log_partition_given_parameters_and_data_density(const partition& part, const parameters& parameters);


class log_factorial {
    typedef float VALUE_TYPE;
    std::deque< VALUE_TYPE > factorial_;
    log_factorial() {
        factorial_.push_back(0.f);
        factorial_.push_back(0.f);
    }
    log_factorial& operator=(log_factorial&);
    log_factorial(const log_factorial&);
public:
    /// \brief the call operator (non-const)
    VALUE_TYPE operator() (size_t n) {
        while (n >= factorial_.size()) {
            factorial_.push_back(factorial_.back() + logf(factorial_.size()));
        }
        return factorial_.at(n);
    }
    /// \brief calculate without creating a reference
    static VALUE_TYPE calc (size_t n) {
        static log_factorial fact;
        return fact(n);
    }
    /// \brief return an instance of the binomial coefficent
    static log_factorial& get() {
        static log_factorial fact;
        return fact;
    }
};

/// \brief calculates the binomial coefficent
///
/// Each value is calculated once and is stored in Pascal's triangle.
/// The calculation includes an overflow control.
/// It is a recursive procedure.
/// This is implemented as a singleton
class binomial_coefficient {
    typedef float VALUE_TYPE;
    /// \brief the row type; the size of each vector is fixed
    typedef std::vector<VALUE_TYPE> row_t;
    std::deque< row_t > triangle_; /// \brief Pascal's triangle
    size_t num_rows_; /// \brief the current number of rows of Pascal's triangle
    /// \brief The calculation routine
    void build_row(size_t n) {
        row_t row(n+1, VALUE_TYPE(1));
        for (size_t k = 1; k < n; ++k) {
            VALUE_TYPE s1 = operator()(n-1, k - 1);
            VALUE_TYPE s2 = operator()(n-1, k);
            if (s1 + s2 < s1) {
                std::stringstream errmsg;
                errmsg << "maximum reached n = " << n << ", k = " << k;
                throw std::overflow_error(errmsg.str());
            }
            row[k] = s1 + s2;
        }
        triangle_.push_back(row);
        num_rows_ = triangle_.size();
    }
    /// \brief The constructor initalises the Pascal's triangle with (0 choose 0)
    binomial_coefficient() {
        triangle_.push_back(row_t(1,VALUE_TYPE(1)));
        num_rows_ = triangle_.size();
    }
    binomial_coefficient& operator=(binomial_coefficient&);
    binomial_coefficient(const binomial_coefficient&);
public:
    /// \brief the call operator (non-const)
    VALUE_TYPE operator() (size_t n, size_t k) {
        if (num_rows_ <= n) build_row(n);
        return triangle_[n][k];
    }
    /// \brief return an instance of the binomial coefficent
    static binomial_coefficient& get() {
        static binomial_coefficient bc;
        return bc;
    }
};

}

}

#endif // BIGGLES_MODEL_HPP__
