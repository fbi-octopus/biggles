/// @file metropolis_hastings.hpp A generic templated Metropolis-Hastings sampler implementation

#ifndef BIGGLES_SAMPLING_METROPOLIS_HASTINGS_HPP
#define BIGGLES_SAMPLING_METROPOLIS_HASTINGS_HPP

#include <cmath>
#include <boost/tuple/tuple.hpp>

#include "../partition_sampler_items.hpp"
#include "../detail/random.hpp"
#include "../detail/fun.hpp"
#include "../mh_moves/mh_moves.hpp"

namespace biggles { namespace sampling {

/// \brief attributes of the current sampling state
struct mh_state_observer {
    /// @brief The log probability density of the current sample.
    ///
    /// this is either last_sample_log_density (if proposal was rejected)
    /// or last_proposal_log_density (if proposal was accepted)
    float sample_log_density;

    /// @brief last log(alpha) where alpha is uniformly distributed between 0 and 1
    float last_log_alpha;

    /// @brief log density of the last proposed move
    float last_proposal_log_density;

    /// \brief the density of the last sample p(T|theta, y)
    float last_sample_log_density;

    /// @brief the last proposal density ratio Q(T|T*)/Q(T*|T)
    float last_pdr;

    /// @brief The number of samples which have been accepted during the lifetime of this sampler. This will always be
    /// less than or equal to \c n_proposed_.
    size_t n_accepted;

    /// @brief The number of samples which have been proposed during the lifetime of this sampler.
    size_t n_proposed;

    /// @brief was the last proposal accepted?
    bool accepted;

    mh_state_observer() :
        sample_log_density(-std::numeric_limits<float>::max()),
        last_log_alpha(-std::numeric_limits<float>::max()),
        last_proposal_log_density(-std::numeric_limits<float>::max()),
        last_sample_log_density(-std::numeric_limits<float>::max()),
        last_pdr(-std::numeric_limits<float>::max()),
        n_accepted(1),
        n_proposed(1),
        accepted(false) {}

    mh_state_observer(const mh_state_observer& o) :
        sample_log_density(o.sample_log_density),
        last_log_alpha(o.last_log_alpha),
        last_proposal_log_density(o.last_proposal_log_density),
        last_sample_log_density(o.last_sample_log_density),
        last_pdr(o.last_pdr),
        n_accepted(o.n_accepted),
        n_proposed(o.n_proposed),
        accepted(o.accepted) {}

    const mh_state_observer& operator = (const mh_state_observer& o) {
        if (&o == this) return *this;
        sample_log_density = o.sample_log_density;
        last_log_alpha = o.last_log_alpha;
        last_proposal_log_density = o.last_proposal_log_density;
        last_sample_log_density = o.last_sample_log_density;
        last_pdr = o.last_pdr;
        n_accepted = o.n_accepted;
        n_proposed = o.n_proposed;
        accepted = o.accepted;
        return *this;
    }

    /// \brief the total accptance rate
    float acceptance_rate() const { return static_cast<float>(n_accepted) / static_cast<float>(n_proposed); }

};

/// @brief A Metropolis-Hastings sampler.
///
/// Metropolis-Hastings is an algorithm which can sample from an arbitrary probability density function, \f$ P(x) \f$,
/// given only that the PDF can be evaluated. In fact it is even more general: the function may return some constant
/// multiple of the true PDF meaning that M-H is exceptionally useful when the PDF is known up to some normalisation
/// factor. The PDF \f$ P(x) \f$ is termed the <em>target</em> distribution.
///
/// In addition to being able to evaluate the target distribution, one must be able to draw samples from a proposal
/// distribution, \f$ Q(x' | x) \f$. In this case \f$ x' \f$ is a new sample and \f$ x \f$ is the sample we currently
/// are at.
///
/// To use the sampler you must provide at a minimum a functor which can evaluate the logarithm of the target PDF at a
/// given sample and a proposal functor. The proposal functor takes a current sample, \f$ x \f$, and returns a new
/// sample, \f$ x' \f$, conditioned on it. In addition, it returns
///
/// \f[
/// \log\left(\frac{Q(x|x')}{Q(x'|x)}\right) = \log(Q(x|x')) - \log(Q(x'|x)).
/// \f]
///
/// This logarithm may of course be zero to indicate a symmetric proposal distribution.
///
/// In addition you may optionally provide the C++ type of the sample and a functor which can sample uniformly from the
/// interval [0,1]. The boost library is used to provide a Mersenne twister-based uniform sampler if you do not specify
/// an alternative.
///
/// As an example of use, consider sampling the Rosenbrock function:
///
/// @code
/// typedef boost::tuples::tuple<float, float> point;
///
/// inline float rosenbrock(const point& p)
/// {
///     float x, y;
///     boost::tuples::tie(x,y) = p;
///     return (1.f-x)*(1.f-x) + 100.f*(y-x*x)*(y-x*x);
/// }
/// @endcode
///
/// Firstly, we need to convert this function which has a minimum at (1,1) to a log-PDF with a maximum at (1,1). Since
/// we do not need to worry about normalisation, this is quite easy:
///
/// @code
/// // a cooked-up PDF which has a maximum where the Rosenbrock function has a minimum.
/// inline float rosenbrock_log_pdf(const point& p)
/// {
///     return -rosenbrock(p);
/// }
/// @endcode
///
/// We also need a proposal function. The proposal function should accept a sample, in this case a \c point, and return
/// a tuple containing a new sample and the log proposal density ration outlined above. We shall use a simple Guassian
/// proposal function which moves the sample according to a Gaussian distribution centred on the current sample with a
/// standard deviation of 0.2. Since this proposal distribution is symmetric, we can return 0 as the log proposal ratio:
///
/// @code
/// typedef boost::tuples::tuple<point, float> proposal_result_t;
///
/// inline proposal_result_t propose_gaussian(const point& p)
/// {
///     static boost::mt19937 rng;
///     static boost::normal_distribution<float> norm(0.f, 0.2f);
///     static boost::variate_generator<boost::mt19937&, boost::normal_distribution<float> > variate(rng, norm);
///
///     point new_p(p.get<0>() + variate(), p.get<1>() + variate());
///
///     return proposal_result_t(new_p, 0.f);
/// }
/// @endcode
///
/// Since these functions are plain C++ functions, we shall wrap them in the functor classes available in the standard
/// library's \c functor header. We firstly need to define a type for them:
///
/// @code
/// typedef std::pointer_to_unary_function<const point&, float> pdf_t;
/// typedef std::pointer_to_unary_function<const point&, proposal_result_t> proposal_func_t;
/// @endcode
///
/// Using the standard library \c std::ptr_fun function, we can now declare our sampler.
///
/// @code
/// // MH sampler for Rosenbrock function. An example of using a function pointer directly.
/// sampler<pdf_t, proposal_func_t> s(std::ptr_fun(rosenbrock_log_pdf),
///                                   std::ptr_fun(propose_gaussian));
/// @endcode
///
/// In this case, the sample type has been inferred directly from the types of the proposal and target functions. You
/// may also specify it directly if necessary.
///
/// Calling draw() on \c s will now draw point locations from the Rosenbrock PDF we defined earlier. In real-world use,
/// you should keep an eye on the value returned from acceptance_rate() to check it is around 0.25.
///
/// @sa http://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm
///
class metropolis_hastings_sampler {
public:
    /// @brief Initialise the sampler.
    ///
    /// @param target A functor representing the target distribution. This is copied internally.
    /// @param propose A functor for proposing new samples. This is copied internally.
    /// @param initial_sample An initial starting sample.
    /// @param uniform A functor for generating uniform reals on the interval [0,1].
    metropolis_hastings_sampler(const partition_distribution& target, const partition_proposal& propose,
            const partition_sampler_sample& initial_sample )
        : target_(target)
        , propose_(propose)
        , sample_(initial_sample)
        , last_proposal_(partition_proposal::result_type(initial_sample, -std::numeric_limits<float>::max()))
        , temperature_(1.0f)
        , acceptance_seq_(5000)
        , internal_records_ptr_(new capability_recorder(initial_sample.partition_sample_ptr->first_time_stamp(),
                            initial_sample.partition_sample_ptr->last_time_stamp(),
                            initial_sample.partition_sample_ptr->tracks()))
    {
        //state_.sample_log_density = target(initial_sample);
    }

    /// @brief Draw the next sample from the sampler.
    ///
    /// @return A reference to the new sample. This can be retrieved before the next call to draw() via last_sample().
    const partition_sampler_sample& draw();

    /// @brief Obtain a reference to the last sample drawn from the sampler.
    const partition_sampler_sample& last_sample() const { return sample_; }

    /// @brief Obtain the log probability density of the last sample drawn from the sampler. log(P(x| theta, D)
    float current_sample_log_density() const { return state_.sample_log_density; }

    /// @brief Obtain the current acceptance rate of the sampler.
    ///
    /// The acceptance rate is defined as the ratio of accepted proposals to total number of proposals. The usual rule
    /// of thumb for Metropolis-Hastings is that this value, often termed \f$ \alpha \f$, should be around 0.25.
    float acceptance_rate() const { return state_.acceptance_rate(); }

    /// @brief Obtain the recent acceptance rate of the sampler.
    ///
    /// This acceptance rate only takes the last samples into account instead of the whole sampling history
    float recent_acceptance_rate() const {
        return std::accumulate(acceptance_seq_.begin(), acceptance_seq_.end(), 0.f)/float(acceptance_seq_.size());
    }

    /// @brief The last proposed sample
    partition_sampler_sample proposed_sample() const { return last_proposal_.partition_sample; }

    /*
    /// @brief The last proposed move
    mh_moves::move_type proposed_move() const { return proposed_sample().executed_move; }
    */

    /// @brief last log alpha
    float last_log_alpha() const { return state_.last_log_alpha; }

    /// @brief log density of the last proposed partition log(P(x'|theta, D))
    float last_log_density() const { return state_.last_proposal_log_density; }

    /// @brief last proposal density ratio Q(T|T*)/Q(T*|T)
    float last_pdr() const {return state_.last_pdr; }

    bool last_proposal_accepted() const { return state_.accepted; }

    mh_state_observer current_state() const { return state_; }
    internal_containers_observer current_internals() const { return intern_; }

    void set_temperature(float temperature) {
        if (temperature > 1.f or temperature <= 0.f) {
            throw std::runtime_error("temperature must be in (0, 1].");
        }
        temperature_ = temperature;
    }

protected:
    /// @brief The target distribution functor.
    partition_distribution target_;

    /// @brief The proposal function and PDF. See class documentation.
    partition_proposal propose_;

    /// @brief The current sample.
    partition_sampler_sample sample_;

    /// @brief The last (possibly not accepted) proposal
    partition_proposal::result_type last_proposal_;

    mh_state_observer state_;
    internal_containers_observer intern_;

    float temperature_;

    fqueue<int> acceptance_seq_;

    capability_recorder_ptr internal_records_ptr_;

};


} }

#endif // BIGGLES_SAMPLING_METROPOLIS_HASTINGS_HPP
