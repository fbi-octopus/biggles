#ifndef BIGGLES_PARTITION_SAMPLER_HPP__
#define BIGGLES_PARTITION_SAMPLER_HPP__

#include <boost/tuple/tuple.hpp>

#include "model.hpp"
#include "partition.hpp"
#include "mh_moves/mh_moves.hpp"
#include "sampling/metropolis_hastings.hpp"
#include "sampling/simple.hpp"

namespace biggles {

/// @brief A class which can sample a set of tracks and clutter given a set of model parameters.
///
/// Internally this class is implemented via the Metropolis-Hastings algorithm. See the documentation for sample() for
/// more information.
///
/// This samples from the posterior distribution \f$ P(T | \theta, D) \f$ where \f$ T \f$ is a set of tracks and
/// clutter which encompass all observations, \f$ \theta \f$ is a set of model parameters and \f$ D \f$ are the
/// observations.
///
/// We can re-arrange this posterior using the law of conditional probability which states \f$ P(T | \theta, D) =
/// P(T, \theta | D) / P(\theta | D) \f$ combined with Bayes' theorem:
///
/// \f[
/// P(T, \theta | D)
///     = \frac{P(D|T, \theta) P(T, \theta)}{P(D)}
///     = \frac{P(D|T, \theta) P(T | \theta) P(\theta)}{P(D)}
/// \f]
///
/// and hence
///
/// \f[
/// P(T | \theta, D)
///     = \frac{P(D|T, \theta) P(T | \theta) P(\theta)}{P(\theta | D)P(D)}
///     = \frac{P(D|T, \theta) P(T | \theta) P(\theta)}{P(\theta, D)}.
/// \f]
///
/// Internally we use the Metropolis-Hastings algorithm which requires only a value proportional to this density.
/// Removing terms independent of \f$ T \f$ we obtain:
///
/// \f[
/// P(T | \theta, D) \propto P(D|T, \theta) P(T | \theta) P(\theta).
/// \f]
///
/// The \f$ P(\theta) \f$ term is simply the prior on the model parameters. The \f$ P(T | \theta) \f$ term is a
/// little more subtle; it is the likelihood of the tracks independent of the observations within them. This is
/// a function of the model parameters and birth and death times for the tracks only. Since the data is partitioned
/// between tracks and the clutter with no overlap, we may factorise the data likelihood term as follows:
///
/// \f[
/// P(D | T, \theta) = P(d_0 | \theta) \prod_{i = 1}^K P(d_i | t_i, \theta)
/// \f]
///
/// where \f$ d_0 \f$ is used to represent the clutter observations and \f$ d_1, ..., d_K \f$ are the \f$ K \f$
/// sets of observations corresponding to the \f$ K \f$ tracks in the partition. \f$ P(d_0 | \theta) \f$ is the
/// likelihood of having seen the clutter observations we have given the parameters. \f$ P(d_i | t_i, \theta) \f$ is
/// the likelihood of track \f$ i \f$ having generated the data we saw.
///
/// Taking logarithms, we obtain the final log density function used within the Metropolis-Hastings sampler:
///
/// \f[
/// \ell(T | \theta, D) = \kappa + \ell(d_0 | t_0, \theta) + \sum_{i = 1}^K \ell(d_i | t_i, \theta) + \ell(T | \theta) + \ell(\theta)
/// \f]
///
/// where \f$ \ell(\cdot) \f$ is used to denote the log-PDF and \f$ \kappa \f$ is some arbitrary normalising offset.
/// Without loss of generality, we set it to zero.
///
/// Each log-PDF term in the expansion above is calculated by one of the log-PDF functions in biggles/model.hpp.
class partition_sampler : public sampling::metropolis_hastings_sampler
{
private:
    typedef sampling::metropolis_hastings_sampler base_type;
public:
    /// @brief Construct a sampler from an initial partition and set of model parameters.
    ///
    /// @param initial_partition
    /// @param initial_parameters
    partition_sampler(const partition_ptr_t& initial_partition_ptr,
                      const model::parameters& initial_parameters = model::parameters())
        : base_type(partition_distribution(initial_parameters),
                    partition_proposal(),
                    partition_sampler_sample(initial_partition_ptr, mh_moves::NONE, mh_moves::NONE))
    {
        // Merge partition sampler and MH sampler?
    }

    /// @brief The model parameters associated with this sampler.
    const model::parameters& parameters() const { return target_.params; }

    /// @brief Return a reference to the partition of the last sample drawn so far.
    const partition_ptr_t last_partition() const { return last_sample().partition_sample_ptr; }

    /// @brief Return a reference to the move type of the last sample drawn so far.
    const mh_moves::move_type& last_move() const { return last_sample().executed_move; }

    /// @brief Modify the parameters associated with this sampler.
    ///
    /// @param params A new set of parameters for the sampler.
    void set_parameters(const model::parameters& params) { target_.params = params; }

    /// @brief set the partition_is_valid value
    void set_partition_is_valid(bool v) { propose_.partition_is_valid = v; }
};


} // namespace biggles

#endif //BIGGLES_PARTITION_SAMPLER_HPP__
