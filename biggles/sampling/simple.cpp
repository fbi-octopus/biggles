#include "detail/random.hpp"
#include "observation.hpp"
#include "sampling/simple.hpp"

#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/shared_ptr.hpp>
#include <cmath>
#include <ctime>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <limits>

extern "C" {
#include "../third-party/gen_beta.h"
}

#include "detail/random.hpp"


namespace biggles { namespace sampling
{

using detail::random_generator;

int uniform_int(int first, int last)
{
    BOOST_ASSERT(last > first);
    boost::uniform_int<> dist(first, last-1); // <- Oh, FFS boost::random. Really?
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > uni(detail::random_generator, dist);
    return uni();
}

float uniform_real(float first, float last)
{
    BOOST_ASSERT(last > first);
    boost::uniform_real<float> dist(first, last);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<float> > uni(detail::random_generator, dist);
    return uni();
}

float exponential(float lambda)
{
    return -log(uniform_real()) / lambda;
}

int poisson(float lambda)
{
    // Knuth's algrotihsm
    float L(exp(-lambda)), p(1.f);
    int k(0);
    do {
        ++k;
        p *= uniform_real();
    } while(p > L);
    return k-1;
}

/*
bool observation_from_light_cone(const observation_collection& observations,
                                        const observation& start,
                                        float speed_of_light,
                                        time_stamp first_time_stamp,
                                        time_stamp last_time_stamp,
                                        observation& output_obs,
                                        float& output_log_prob)
{
    BOOST_ASSERT(first_time_stamp <= last_time_stamp);

    // early-out if given a zero-sized time range
    if(first_time_stamp == last_time_stamp)
        return false;

    // find candidate observations
    observation_collection candidates;
    std::set<time_stamp> time_stamp_with_cand;
    observations_within_light_cone(observations,
                                   start,
                                   speed_of_light,
                                   first_time_stamp, last_time_stamp,
                                   time_stamp_with_cand,
                                   std::inserter(candidates, candidates.begin()));

    // did we find any?
    if(candidates.empty())
    {
        // no, return failure
        return false;
    }

    // which time stamp do we choose?
    std::set<time_stamp>::const_iterator ts_iter = time_stamp_with_cand.begin();
    std::advance(ts_iter, uniform_int(0, time_stamp_with_cand.size()));
    float time_stamp_log_prob = -logf(static_cast<float>(time_stamp_with_cand.size()));

    // what candidates are available to us?
    observation_collection::const_range refined_range(candidates.at_time_stamp(*ts_iter));

    // does that leave us with any to choose from?
    BOOST_ASSERT((refined_range.first != refined_range.second));

    // _finally_ let's sample an observation
    float refined_log_prob(1.f);
    output_obs = *from_range(refined_range.first, refined_range.second, refined_log_prob);
    BOOST_ASSERT(refined_log_prob <= 0.f);

    // update the overall log prob
    output_log_prob = time_stamp_log_prob + refined_log_prob;

    return true;
}
*/

float log_prob_of_sampling_from_light_cone(const observation& o,
                                           const observation_collection& observations,
                                           const observation& start,
                                           float speed_of_light,
                                           time_stamp first_time_stamp,
                                           time_stamp last_time_stamp)
{
    BOOST_ASSERT(first_time_stamp <= last_time_stamp);
    BOOST_ASSERT(observations.find(o) != observations.end());

    // early-out if given a zero-sized time range
    if(first_time_stamp == last_time_stamp)
        return -std::numeric_limits<float>::max();

    // early-out if observation is outside of time cone
    if(!observation_is_within_light_cone(o, start, speed_of_light, first_time_stamp, last_time_stamp))
        return -std::numeric_limits<float>::max();

    // the observation *is* somewhere in the light cone

    // find candidate observations
    observation_collection candidates;
    std::set<time_stamp> time_stamp_with_cand;
    observations_within_light_cone(observations,
                                   start,
                                   speed_of_light,
                                   first_time_stamp, last_time_stamp,
                                   time_stamp_with_cand,
                                   std::inserter(candidates, candidates.begin()));
    BOOST_ASSERT(candidates.find(o) != candidates.end());

    // log prob of having sampled the time stamp
    float time_stamp_log_prob = -logf(time_stamp_with_cand.size());

    // the log prob. of sampling the observation
    size_t n_candidates_at_t(candidates.count_at_time_stamp(t(o)));
    BOOST_ASSERT(n_candidates_at_t > 0);
    float refined_log_prob = -logf(static_cast<float>(n_candidates_at_t));

    return time_stamp_log_prob + refined_log_prob;
}

bool observation_from_collection(const observation_collection& oc, observation& output_obs, float& output_log_prob)
{
    throw std::logic_error("this function is no longer in use");
    // early out if we cannot sample
    if(oc.empty())
        return false;

    std::set<time_stamp> ts_with_obs;

    for (observation_collection::const_iterator iter = oc.begin(); iter != oc.end(); ++iter)
        ts_with_obs.insert(t(*iter));

    // choose a time stamp to sample an observation from
    std::set<time_stamp>::iterator b = ts_with_obs.begin();
    std::advance(b, uniform_int(0, ts_with_obs.size()));
    time_stamp ts = *b;
    float log_ts_prob = -logf(ts_with_obs.size());
    /*
    time_stamp ts(uniform_int(oc.first_time_stamp(), oc.last_time_stamp()));
    float log_ts_prob = -logf(oc.last_time_stamp() - oc.first_time_stamp());
    */

    // get the observations at this time stamp
    observation_collection::const_range candidate_range(oc.at_time_stamp(ts));

    // found no observations?
    if(candidate_range.first == candidate_range.second)
        return false;

    // sample an observation
    float log_obs_prob(1.f);
    output_obs = *from_range(candidate_range.first, candidate_range.second, log_obs_prob);
    BOOST_ASSERT(log_obs_prob <= 0.f);

    // update log probability
    output_log_prob = log_ts_prob + log_obs_prob;

    return true;
}

float log_prob_of_sampling_observation(const observation& o,
                                       const observation_collection& oc)
{
    // early out if we cannot sample
    if(oc.empty())
        return -std::numeric_limits<float>::max();

    // early-out if observation is outside of temporal bounds
    if((t(o) < oc.first_time_stamp()) || (t(o) >= oc.last_time_stamp()))
        return -std::numeric_limits<float>::max();

    std::set<time_stamp> ts_with_obs;
    for (observation_collection::const_iterator iter = oc.begin(); iter != oc.end(); ++iter)
        ts_with_obs.insert(t(*iter));

    // log prob of choosing the time stamp
    float log_ts_prob = -logf(ts_with_obs.size());

    // get the observations at the observation's time stamp
    observation_collection::const_range candidate_range(oc.at_time_stamp(t(o)));

    // found no observations?
    if(candidate_range.first == candidate_range.second)
        return -std::numeric_limits<float>::max();

    // how likely was it we sampled this observation?
    size_t n_at_ts(oc.count_at_time_stamp(t(o)));
    BOOST_ASSERT(n_at_ts > 0);
    float log_obs_prob = -logf(static_cast<float>(n_at_ts));
    BOOST_ASSERT(log_obs_prob <= 0.f);

    // return log probability
    return log_ts_prob + log_obs_prob;
}

// sampling functions
static double rand_double();
//static float sample_normal(); // zero-mean, unit variance
inline float sample_exponential(float lambda);



template<typename Real, int Dim>
Eigen::Matrix<Real, Dim, 1> sample_multivariate_gaussian(const Eigen::Matrix<Real, Dim, 1>& mean,
                                                         const Eigen::Matrix<Real, Dim, Dim>& covariance);


// return a random double on the interval [0,1].
static double rand_double()
{
    static boost::uniform_real<double> dist(0.f, 1.f);
    static boost::variate_generator<boost::mt19937&, boost::uniform_real<double> > uni(random_generator, dist);
    return uni();
}

float sample_normal()
{
    static boost::normal_distribution<float> dist(0.f, 1.f);
    static boost::variate_generator<boost::mt19937&, boost::normal_distribution<float> > norm(random_generator, dist);
    return norm();
}

// sample X ~ Exponential(lambda)
inline float sample_exponential(float lambda)
{
    typedef boost::exponential_distribution<float> dist_t;
    dist_t dist(lambda);
    boost::variate_generator<boost::mt19937&, dist_t > sample(random_generator, dist);
    return sample();
}

// sample X ~ Gamma(k, theta)
float sample_gamma(unsigned int k, float theta)
{
    float sum(0.f);
    for(unsigned int i=0; i<k; ++i)
        sum += sample_exponential(1.f / theta);
    return sum;
}

// sample X ~ Gamma(k, theta)
float sample_gamma0(unsigned int k, float theta)
{
    float sum(std::numeric_limits<float>::min());
    for(unsigned int i=0; i<k; ++i)
        sum += sample_exponential(1.f / theta);
    return sum;
}

// sample X ~ Beta(alpha, beta)
float sample_beta(float alpha, float beta)
{
    gen_beta_param gen;
    gen_beta_initialize(&gen, alpha, beta);
    return gen_beta(&gen, rand_double);
}

float sample_truncated_beta(float alpha, float beta, float lo_limit, float hi_limit) {
    gen_beta_param gen;
    gen_beta_initialize(&gen, alpha, beta);
    float sample = -1;
    while (sample < lo_limit or sample > hi_limit)
        sample = gen_beta(&gen, rand_double);
    return sample;
}

Eigen::Matrix2f sample_inverse_wishart(const Eigen::Matrix2f& Phi, size_t s)
{
    return sample_wishart(Phi.inverse(), s).inverse();
}

Eigen::Matrix2f sample_wishart(const Eigen::Matrix2f& Phi, size_t s)
{
    BOOST_ASSERT(s > 1);

    Eigen::Matrix2f sample(Eigen::Matrix2f::Zero());

    // this simply makes use of the definition of the Wishart distribution.
    for(size_t i=0; i<s; ++i)
    {
        Eigen::Vector2f norm(sample_multivariate_gaussian<float,2>(Eigen::Vector2f::Zero(), Phi));
        sample += norm * norm.transpose();
    }

    return sample;
}

//float t_sample_normal() { return sample_normal(); }


} }  // biggles::sampling
