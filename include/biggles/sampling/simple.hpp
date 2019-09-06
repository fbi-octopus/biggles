/// @file simple.hpp Sampling for various simple base distributions

#ifndef BIGGLES_SAMPLING_SIMPLE_HPP__
#define BIGGLES_SAMPLING_SIMPLE_HPP__

#include <numeric> // std::partial_sum
#include <vector>
#include <boost/tuple/tuple.hpp>

#include "../detail/random.hpp"
#include "../detail/fun.hpp"
#include "../observation.hpp"
#include "../observation_collection.hpp"

namespace biggles {

/// @brief Utility sampling methods and functions.
namespace sampling
{

/// @brief Sample an integer from a uniform distribution on the interval [\p first, \p last).
///
/// @note This function will never return an integer equal to \p last.
///
/// @param first
/// @param last
int uniform_int(int first, int last);

/// @brief Sample a real uniformly from the interval [\p first, \p last).
///
/// @param first
/// @param last
float uniform_real(float first = 0.f, float last = 1.f);

/// @brief Sample a floating point value from the exponential distribution with rate parameter \p lambda.
///
/// @param lambda
float exponential(float lambda = 1.f);

/// @brief Sample a Poisson variate from a distribution with mean \p lambda.
///
/// @param lambda
int poisson(float lambda);

/// @brief Sample an iterator uniformly from the range [\p first, \p last).
///
/// If \p first = \p last, then return \p last.
///
/// @note This is currently only implemented for \c InputIterator types. Really there should be a specialisation for
/// a \c RandomAccessIterator. This is a FIXME for someone enthused to learn about iterator traits.
///
/// @tparam InputIterator
/// @param first
/// @param last
template<typename InputIterator>
InputIterator from_range(InputIterator first, InputIterator last);

/// @brief Sample an iterator uniformly from the range [\p first, \p last).
///
/// If \p first = \p last, then return \p last. The reference \p log_prob is only written to when \p first != \p last.
///
/// @note This is currently only implemented for \c InputIterator types. Really there should be a specialisation for
/// a \c RandomAccessIterator. This is a FIXME for someone enthused to learn about iterator traits.
///
/// @tparam InputIterator
/// @param first
/// @param last
/// @param[out] log_prob write the log probability of having sampled the return value to this reference
template<typename InputIterator>
InputIterator from_range(InputIterator first, InputIterator last, float& log_prob);

/*
/// @brief Sample an observation from within a light cone.
///
/// Given the definition of a light cone as in biggles::observations_within_light_cone(), sample an observation from it
/// with the following procedure:
///
/// - sample uniformly some time stamp within the light cone
/// - sample uniformly from the observations with that time stamp
///
/// The log probability of drawing the sample is returned in \p output_log_prob and the sample itself is returned in \p
/// output_obs.
///
/// @note The sampling procedure may sometimes fail to sample an observation. In this case the output references are
/// untouched and the function returns \p false.
///
/// @sa observations_within_light_cone()
/// @sa log_prob_of_sampling_from_light_cone()
///
/// @param observations
/// @param start
/// @param speed_of_light
/// @param first_time_stamp
/// @param last_time_stamp
/// @param[out] output_obs
/// @param[out] output_log_prob
///
/// @return \c true iff a sample could be drawn
bool observation_from_light_cone(const observation_collection& observations,
                                 const observation& start,
                                 float speed_of_light,
                                 time_stamp first_time_stamp,
                                 time_stamp last_time_stamp,
                                 observation& output_obs,
                                 float& output_log_prob);
                                 */

/// @brief Return the log probability of having sampled an observation from a light cone.
///
/// This is the distribution from which observation_from_light_cone() actually samples. If \p o could never have
/// been sampled, this returns \c -FLT_MAX as a place holder for negative infinity.
///
/// The parameters are as those in observation_from_light_cone().
///
/// @sa observations_within_light_cone()
/// @sa observation_from_light_cone()
///
/// @param o
/// @param observations
/// @param start
/// @param speed_of_light
/// @param first_time_stamp
/// @param last_time_stamp
float log_prob_of_sampling_from_light_cone(const observation& o,
                                           const observation_collection& observations,
                                           const observation& start,
                                           float speed_of_light,
                                           time_stamp first_time_stamp,
                                           time_stamp last_time_stamp);

/// @brief Sample an observation from an observation collection
///
/// Sample an observation from \p oc using the following procedure:
///
/// - sample a time stamp uniformly from the range covered by \p oc;
/// - sample an observation uniformly from those in \p oc with that timestamp
///
/// It is possible for this procedure to fail to sample an observation in which case this function returns \c false. If
/// the sampling succeeded, the sample is written to \p output_obs and the log probability of having sampled it is
/// written to \p output_log_prob.
///
/// @sa log_prob_of_sampling_observation()
///
/// @param oc Sample an observation from this collection
/// @param[out] output_obs
/// @param[out] output_log_prob
///
/// @return \c true iff an observation was sampled
bool observation_from_collection(const observation_collection& oc,
                                 observation& output_obs,
                                 float& output_log_prob);

/// @brief The distribution from which observation() samples.
///
/// This is the distribution from which observation() actually samples. If \p o could never have been sampled,
/// this returns \c -FLT_MAX as a place holder for negative infinity.
///
/// The parameters are as those in observation().
///
/// @param observation The observation sampled
/// @param oc The collection from which it was sampled
float log_prob_of_sampling_observation(const observation& observation,
                                       const observation_collection& oc);


/// @brief samples from a range of \p items with probabilities proporational to \p weights
///
/// @param items_begin the iterator pointing to the begin of \p items
/// @param weights_begin the iterator pointing to the begin of \p weights
/// @param weights_end the iterator pointing to the end of \p weights
///
/// @return the sampled item
template <class ITEMITER, class WEIGHTITER>
typename ITEMITER::value_type
weighted_choice(ITEMITER items_begin, WEIGHTITER weights_begin, WEIGHTITER weights_end) {
    std::vector<float> partsum(std::distance(weights_begin, weights_end));
    std::partial_sum(weights_begin, weights_end, partsum.begin());
    size_t i = std::distance( partsum.begin(),
        std::lower_bound( partsum.begin(), partsum.end(), uniform_real(0.f, partsum.back()) )
    );
    return *(items_begin+i);
}

/// @brief samples from a range of \p items with probabilities proporational to \p weights
///
/// @param items_begin the iterator pointing to the begin of \p items
/// @param weights_begin the iterator pointing to the begin of \p weights
/// @param weights_end the iterator pointing to the end of \p weights
/// @param log_prob the log-probability to sample the returned item
///
/// @return the sampled item
template <class ITEMITER, class WEIGHTITER>
typename ITEMITER::value_type
weighted_choice(ITEMITER items_begin, WEIGHTITER weights_begin, WEIGHTITER weights_end, float& log_prob) {
    std::vector<float> partsum(std::distance(weights_begin, weights_end));
    float total = std::accumulate(weights_begin, weights_end, 0.0);
    std::partial_sum(weights_begin, weights_end, partsum.begin());
    size_t i = std::distance( partsum.begin(),
        std::lower_bound( partsum.begin(), partsum.end(), uniform_real(0.f, partsum.back()) )
    );
    log_prob = logf(*(weights_begin+i)/total);
    std::advance(items_begin, i);
    return *items_begin;
}

/// @brief samples from a range of \p items with probabilities proporational to weights given by \em fun
///
/// @param items_begin the iterator pointing to the begin of \p items
/// @param items_end the iterator pointing to the end of \p items
/// @param fun the the weight functor
/// @param log_prob the log-probability to sample the returned item
///
/// @return the iterator to the sampled item
template <class ITEMITER, class FUN>
ITEMITER weighted_choice(ITEMITER items_begin, ITEMITER items_end, FUN& fun , float& log_prob) {
    std::vector<float> partsum(std::distance(items_begin, items_end));
    fun_partial_sum(items_begin, items_end, partsum.begin(), fun, 0.f);
    if (partsum.back() == 0) return items_end;
    size_t i = std::distance( partsum.begin(),
        std::lower_bound( partsum.begin(), partsum.end(), uniform_real(0.f, partsum.back()) )
    );
    log_prob = std::log((i == 0 ? *partsum.begin() : partsum[i] - partsum[i-1]) / partsum.back()) ;
    std::advance(items_begin, i);
    return items_begin;
}

template <class WEIGHTITER>
size_t weighted_choice_index(WEIGHTITER weights_begin, WEIGHTITER weights_end, float& log_prob) {
    std::vector<float> partsum(std::distance(weights_begin, weights_end));
    float total = std::accumulate(weights_begin, weights_end, 0.f);
    std::partial_sum(weights_begin, weights_end, partsum.begin());
    size_t i = std::distance( partsum.begin(),
        std::lower_bound( partsum.begin(), partsum.end(), uniform_real(0.f, partsum.back()) )
    );
    log_prob = logf(*(weights_begin+i)/total);
    return i;
}


/// @brief  sample inverse wishart
Eigen::Matrix2f sample_inverse_wishart(const Eigen::Matrix2f& Phi, size_t s);
/// @brief  sample wishart
Eigen::Matrix2f sample_wishart(const Eigen::Matrix2f& Phi, size_t s);

/// @brief A simple functor yielding a uniformly sampled real from the interval [0, 1).
struct uniform_real_functor
{
    typedef float result_type;
    float operator () () const { return sampling::uniform_real(); }
};

float sample_normal();

/// @brief A default uniform real generator.
///
/// This is a convenience functor struct which can be used to get a uniform variate on the interval [0,1]. It is the
/// default source of randomness for metropolis_hastings::sampler<> unless otherwise specified.
///
/// It is implemented using an internal boost Mersenne twister random number generator and a uniform real
/// distribution.
struct boost_uniform
{
    boost_uniform() : dist_(0.,1.), uni_real_(biggles::detail::random_generator, dist_) { }
    boost_uniform(const boost_uniform&) : dist_(0.,1.), uni_real_(biggles::detail::random_generator, dist_) { }
    const boost_uniform& operator = (const boost_uniform&) { return *this; }

    float operator() () { return uni_real_(); }
protected:
    boost::uniform_real<float> dist_;
    boost::variate_generator<boost::mt19937&, boost::uniform_real<float> > uni_real_;
};

/// @brief Sample a Beta variate from a distribution with shape parameters \p alpha and \p beta.
///
/// @param alpha
/// @param beta
float sample_beta(float alpha, float beta);

/** \brief Sample a truncated Beta variate from a distribution with shape parameters \p alpha and \p beta.
 * and limits \p lo_limit and \p hi_limit
 *
 * The truncated beta is defined by:
 * * p(x; alpha, beta) = B(x; alpha, beta)/B([lo_limit, hi_limit]; alpha, beta) if \b x is in [lo_limit, hi_limit]
 * * p(x; alpha, beta) = 0 if \b x is not in  [lo_limit, hi_limit]
 *
 */
float sample_truncated_beta(float alpha, float beta, float lo_limit, float hi_limit);

/// @brief Sample a Gamma variate from a distribution with shape parmater \p k and scale parameter \p theta.
///
/// @param k shape
/// @param theta scale
float sample_gamma(unsigned int k, float theta);

/// @brief Sample a Gamma variate from a distribution with shape parmater \p k and scale parameter \p theta.
///
/// Ensures that the result is > 0.0f
/// @param k shape
/// @param theta scale
float sample_gamma0(unsigned int k, float theta);

template<typename Real, int Dim>
Eigen::Matrix<Real, Dim, 1> sample_multivariate_gaussian(const Eigen::Matrix<Real, Dim, 1>& mean,
                                                         const Eigen::Matrix<Real, Dim, Dim>& covariance)
{
    // this uses http://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution

    typedef Eigen::Matrix<Real, Dim, Dim> Matrix;
    typedef Eigen::Matrix<Real, Dim, 1> Vector;

    Vector sample;

    for(int i=0; i<Dim; ++i)
        sample(i) = sample_normal();

    Matrix cholesky_decomposition = covariance.llt().matrixL();

    return mean + cholesky_decomposition * sample;
}



}

}

#define WITHIN_BIGGLES_SAMPLING_SIMPLE_HPP__
#include "simple.tcc"
#undef WITHIN_BIGGLES_SAMPLING_SIMPLE_HPP__

#endif // BIGGLES_SAMPLING_SIMPLE_HPP__
