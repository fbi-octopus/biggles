/// @file gibbs.hpp A generic templated Gibbs sampler implementation

#ifndef SAMPLING_GIBBS_HPP__
#define SAMPLING_GIBBS_HPP__

#include <boost/tuple/tuple.hpp>
#include <functional>
#include <stdexcept>

#include "../partition_sampler.hpp"



namespace biggles {

typedef boost::tuple<partition_sampler_sample, model::parameters> tracker_sample;

namespace sampling {

/// @brief An implementation of a Gibbs sampler.
///
/// A Gibbs sampler can sample from some conditional joint posterior distribution, \f$ P(x_1, x_2, ..., x_N | Y) \f$
/// given only the ability to sample from individual conditional probability distributions \f$ P(x_i | x_1, ...,
/// x_{i-1}, x_{i+1}, ..., x_N, Y) \f$.
///
/// Each sample from the sampler is a tuple \f$(x_1, x_2, ..., x_N)\f$ of samples from the joint distribution. This
/// tuple has type \p SampleTuple which must be an appropriate form of \c boost::tuple<>. The sampler maintains a copy
/// of the last sample tuple drawn which may be accessed via the last_sample() member function.
///
/// The sampler draws samples from the individual conditional distributions via the \p ConditionalSampler type. The
/// sampler maintains an internal instance of this type which is used to draw samples. If \c cs is an instance of \p
/// ConditionalSampler and \c st is an instance of \c SampleTuple, then the following expression should draw a sample
/// for \f$x_{i+1}\f$ where \f$0 \le i < N\f$:
///
/// @code
/// cs.template draw<i>(st);
/// @endcode
///
/// The sampler may use values within \c st not corresponding to index \c i as conditional values for sampling. As a
/// concrete example, the following \c struct would be a minimal declaration matching the \c ConditionalSampler concept.
///
/// @code
/// #include <boost/tuple.hpp>
///
/// // ...
///
/// template<typename SampleTuple>
/// struct ConditionalSampler
/// {
///     // Draw a sample for x_{idx+1} from P(x_{idx+1} | x_1, ..., x_{idx}, x_{idx+2}, ..., x_N).
///     template<int idx>
///     boost::tuples::element<idx, SampleTuple>::type draw(const SampleTuple&);
/// };
/// @endcode
///
/// @sa two_var_cond_sampler An example of a \p ConditionalSampler which assumes that the Gibbs sampler is sampling from
/// the joint distribution of just two random variables.
///
/// @tparam SampleTuple A boost::tuple which specifies the types of \f$ x_1 \f$, \f$ x_2 \f$, etc.
/// @tparam ConditionalSampler A sampler which can sample from each conditional distribution. See
/// two_var_cond_sampler for an example.
template<typename ConditionalSampler, typename SampleTuple = typename ConditionalSampler::sample_tuple_type>
class gibbs_sampler
{
public:
    /// @brief Specify initial sample tuple and conditional sampler instance.
    ///
    /// @param cond_sampler A copy of this conditional sampler is taken to initialise the internal instance.
    /// @param initial_sample The initial sample tuple to use for the sampler.
    gibbs_sampler(const ConditionalSampler& cond_sampler = ConditionalSampler(),
            const SampleTuple& initial_sample = SampleTuple())
        : sample_(initial_sample)
        , cond_sampler_(cond_sampler)
    { }

    /// @brief Draw a new sample from the conditional joint distribution.
    ///
    /// @return A reference to the sample drawn by the sampler which can be queried before the next invocation of draw()
    /// via the last_sample() member function.
    const SampleTuple& draw()
    {
        typedef gibbs_sampler<ConditionalSampler, SampleTuple> self_t;

        // start the sample drawing process.
        sample_drawer<self_t, boost::tuples::length<SampleTuple>::value, 0>()(*this);
        return last_sample();
    }

    /// @brief Obtain a reference to the last sample drawn by the sampler.
    const SampleTuple& last_sample() const { return sample_; }

protected:
    /// @brief The current sample tuple for the sampler.
    SampleTuple sample_;

    /// @brief The sample evolution functor.
    ConditionalSampler cond_sampler_;

private:
    /// These implement a template metaprogramming recursion to draw one sample in turn from the conditional PDF.
    /// @{

    template<typename Parent, int maxidx, int idx>
    struct sample_drawer
    {
        inline void operator() (Parent& p) const
        {
            boost::tuples::get<idx>(p.sample_) = p.cond_sampler_.template draw<idx>(p.sample_);
            sample_drawer<Parent, maxidx, idx+1>()(p);
        };
    };

    template<typename Parent, int maxidx>
    struct sample_drawer<Parent, maxidx, maxidx>
    {
        inline void operator() (Parent&) const { }
    };

    /// @}
};

class new_gibbs_sampler {
    //new_gibbs_sampler();
    new_gibbs_sampler(const new_gibbs_sampler&);
    void operator=(const new_gibbs_sampler&);
public:
    new_gibbs_sampler(const partition_ptr_t& initial_partition_ptr = partition_ptr_t(new partition),
        const model::parameters& initial_parameters = model::parameters());
    tracker_sample draw();
    tracker_sample last_sample() const;
    const partition_sampler_sample& last_partition_sample() const;
    const model::parameters& last_parameter_sample() const;
    /** set the observation error to the passed values and stop sampling it
     *
     */
    void fix_observation_error(const float R00, const float R01, const float R11);
    /** \brief draw a sample of the oservation error every *lag* samples
     *
     * lag==0 does the same as lag==1
     */
    void sample_observation_error(const size_t lag = 1);
protected:
    /// \brief read access to the Metropolis-Hastings sampler
    const partition_sampler& underlying_partition_sampler() const { return partition_sampler_; }
    /// \brief write access to the Metropolis-Hastings sampler
    partition_sampler& underlying_partition_sampler() { return partition_sampler_; }
private:
    size_t observation_error_lag_; ///< \brief the gap between two observation error samples
    size_t gibbs_sample_count_; ///< \brief the internal counter for the observation error samples
    partition_sampler_sample partition_sample_; ///< \brief the current partition sample incl. proposed/executed moves
    model::parameters parameter_sample_; ///< \brief the current parameter sample
    partition_sampler partition_sampler_; ///< \brief the Metropolis-Hastings sampler
};



namespace gibbs {

/// @brief An base class for a two-variable sampler.
///
/// This class serves as an implementation of the \c ConditionalSampler concept which can sample from distributions of
/// the form \f$ P(x_1, x_2 | Y) \f$. It does this by taking two \c UnaryFunction functors which can draw a sample of
/// \f$ x_1 \f$ given \f$ x_2 \f$ and vice-versa. Any static data, \f$Y\f$, should be hidden within the functor state.
///
/// The sampler is templated on two types implementing the \c UnaryFunction concept, \p UnaryFunction1 and \p
/// UnaryFunction2. The return type of \p UnaryFunction1 should correspond to \f$ x_1 \f$ and the return type of \p
/// UnaryFunction2 should correspond to \f$ x_2 \f$. The argument types should correspond to \f$ x_2 \f$ and \f$ x_1 \f$
/// respectively.
///
/// The instances of these functors provided should sample values of \f$ x_1 \f$ and \f$ x_2 \f$ appropriately
/// conditioned on their arguments.
///
/// @sa mem_fun_two_var_cond_sampler A convenience derived class capable of using its member functions as sampling
/// functions.
///
/// @tparam UnaryFunction1 A std::unary_function to draw samples of \f$x_1\f$.
/// @tparam UnaryFunction2 A std::unary_function to draw samples of \f$x_2\f$.
template<typename UnaryFunction1, typename UnaryFunction2>
class two_var_cond_sampler
{
public:
    /// @brief The C++ type corresponding to \f$ x_1 \f$.
    typedef typename UnaryFunction1::result_type first_sample_type;

    /// @brief The C++ type corresponding to \f$ x_2 \f$.
    typedef typename UnaryFunction2::result_type second_sample_type;

    /// @brief The C++ type corresponding to the tuple \f$ (x_1, x_2) \f$.
    typedef boost::tuple<first_sample_type, second_sample_type> sample_tuple_type;

    /// @brief Initialise the sampling functors.
    ///
    /// @param first_sampler An instance of \p UnaryFunction1 which is copied internally.
    /// @param second_sampler An instance of \p UnaryFunction2 which is copied internally.
    two_var_cond_sampler(
        const UnaryFunction1 first_sampler = UnaryFunction1(),
        const UnaryFunction2 second_sampler = UnaryFunction2())
        :   first_sampler_(first_sampler)
        ,   second_sampler_(second_sampler)
    { }

    /// @brief Draw a sample from one conditional distribution.
    ///
    /// Draw a sample of \f$ x_{i+1} \f$ from the individual conditional distributions.
    ///
    /// @tparam i Which random variable to sample.
    /// @param sample_tuple A tuple giving the conditional values.
    ///
    /// @return A new sample of \f$ x_{i+1} \f$.
    template<int i>
    inline typename boost::tuples::element<i, sample_tuple_type>::type draw(const sample_tuple_type& sample_tuple)
    {
        typedef two_var_cond_sampler<UnaryFunction1, UnaryFunction2> self_t;
        return sampler_func<self_t, i>()(*this, boost::tuples::get<1-i>(sample_tuple));
    }

protected:
    /// @brief A functor for sampling \f$x_1\f$ from \f$P(x_1 | x_2, Y)\f$.
    UnaryFunction1 first_sampler_;

    /// @brief A functor for sampling \f$x_2\f$ from \f$P(x_2 | x_1, Y)\f$.
    UnaryFunction2 second_sampler_;

private:
    /// These implement a templated functor for sampling from either the first or second sampler.
    /// @{

    template<typename Parent, int>
    struct sampler_func { };

    template<typename Parent>
    struct sampler_func<Parent, 0>
    {
        inline first_sample_type operator() (Parent& p, const second_sample_type& a) const
        {
            return p.first_sampler_(a);
        }
    };

    template<typename Parent>
    struct sampler_func<Parent, 1>
    {
        inline second_sample_type operator() (Parent& p, const first_sample_type& a) const
        {
            return p.second_sampler_(a);
        }
    };

    /// @}
};

/// @brief An abstract base class for a class which implements two-variable conditional sampling in member functions.
///
/// Derived classes should implement the following the two member functions \c sample_variable_1() and \c
/// sample_variable_2(). For example:
///
/// @code
/// template<typename Variable1, typename Variable2>
/// struct derived_two_var_cond_sampler
///     : public mem_fun_two_var_cond_sampler<derived_two_var_cond_sampler, Variable1, Variable2>
/// {
///     Variable1 sample_variable_1(Variable2);
///
///     Variable2 sample_variable_2(Variable1);
/// };
/// @endcode
///
/// @tparam Derived The type of the derived class.
/// @tparam Variable1 The type corresponding to \f$ x_1 \f$.
/// @tparam Variable2 The type corresponding to \f$ x_2 \f$.
template<typename Derived, typename Variable1, typename Variable2>
class mem_fun_two_var_cond_sampler
    : public two_var_cond_sampler<
      std::binder1st<std::mem_fun1_t<Variable1, Derived, Variable2> >,
      std::binder1st<std::mem_fun1_t<Variable2, Derived, Variable1> > >
{
public:
    /// @brief A \c UnaryFunctor type which returns a sample of \f$ x_1 \f$ given a value of \f$ x_2 \f$.
    typedef std::binder1st<std::mem_fun1_t<Variable1, Derived, Variable2> > sampler_1_func;

    /// @brief A \c UnaryFunctor type which returns a sample of \f$ x_2 \f$ given a value of \f$ x_1 \f$.
    typedef std::binder1st<std::mem_fun1_t<Variable2, Derived, Variable1> > sampler_2_func;

    mem_fun_two_var_cond_sampler()
        : two_var_cond_sampler<sampler_1_func, sampler_2_func>(sampler_1(), sampler_2())
    { }

private:
    inline sampler_1_func sampler_1() const
    {
        return std::bind1st(std::mem_fun(&Derived::sample_variable_1), this);
    }

    inline sampler_2_func sampler_2() const
    {
        return std::bind1st(std::mem_fun(&Derived::sample_variable_2), this);
    }
};

} } } // biggles::sampling::gibbs

#endif // SAMPLING_GIBBS_HPP__
