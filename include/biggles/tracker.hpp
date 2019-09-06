/// @file tracker.hpp Tracking features in input partitions

#ifndef BIGGLES_TRACKER_HPP___
#define BIGGLES_TRACKER_HPP___
#include <boost/tuple/tuple.hpp>
#include <boost/assert.hpp>

#include "mh_moves/mh_moves.hpp"
#include "model.hpp"
#include "partition.hpp"
#include "samplers.hpp"
#include "partition_sampler.hpp"
#include "sampling/gibbs.hpp"
#include "sampling/metropolis_hastings.hpp"

#include "tracker_detail.hpp"

namespace biggles
{

/** \brief a record for the move statistics
 *
 * this stores how often a proposed move was accepted, rejected or produced an identity move
 * the birth move has index 0, the number of *proper* moves is IDENTITY - BIRTH
 */
struct move_stats_t {
    std::deque<uint64_t> rejected;
    std::deque<uint64_t> accepted;
    std::deque<uint64_t> identity;
    move_stats_t () {
        size_t len = static_cast<size_t>(mh_moves::IDENTITY - mh_moves::BIRTH);
        rejected.assign(len, 0);
        accepted.assign(len, 0);
        identity.assign(len, 0);
    }
};

/// @brief An implementation of the Biggles tracking algorithm.
///
/// The Biggles tracking algorithm is at it's highest level a sampler which maintains an optimum over the samples drawn.
/// More specifically, the Biggles model specifies that there are three main components to the mathematical model used
/// to represent the state of the world:
///
/// - the input observations from the feature detector, \f$ D \f$;
/// - the association between these observations and one or more tracks, \f$ T \f$; and
/// - parameters for the underlying dynamic model of molecule evolution, \f$ \theta \f$.
///
/// The track association also includes a special 'track' known as the clutter which holds all observations we assert
/// are due to spurious detection rather than a tracked molecule. This is the only track which may hold more than one
/// observation per time stamp.
///
/// The value of \f$ T \f$ is represented in this code by biggles::partition which stores a collection of
/// biggles::observation instances for the clutter and a collection of biggles::track instances for each track. The
/// collection of tracks is implemented by biggles::track_collection which implicitly shares tracks between collections.
///
/// The tracking algorithm proceeds by drawing samples of \f$ T \f$ and \f$ \theta \f$ from the joint distribution \f$
/// P(T, \theta | D) \f$. By definition of the mode, we expect more of these samples to be close to the optimum than far
/// from it. We therefore maintain an 'optimal' sample 2-tuple, \f$ ( \hat{T}, \hat{\theta} ) \f$, which serve as our
/// current estimate for the 'true' tracks. After drawing each sample from the joint distribution, we compare the value
/// of the distribution at that point to that of our 'best' sample. If the newly drawn sample has a higher posterior
/// than the current best one, replace the best sample with the one just drawn.
///
/// Details of precisely how these samples are drawn and the mathematical model we use can be found elsewhere in the
/// documentation as indicated by the 'see also' section below.
///
/// @sa biggles::model
/// @sa biggles::partition
/// @sa biggles::sampling
/// @sa biggles::track
/// @sa sampling::gibbs_sampler
///
class tracker : public sampling::new_gibbs_sampler {
public:
    /// @brief Default constructor.
    tracker();

    tracker(const tracker& t);
    /// @brief Initialise the tracker with an initial set of model parameters and partition.
    ///
    /// @param initial_parameters
    /// @param initial_partition
    tracker(const model::parameters& initial_parameters,
            const partition_ptr_t& initial_partition);

    /// @brief Copy constructor.
    ///
    /// @param t

    /// @brief Advance the tracking algorithm one step.
    ///
    /// This draws a single sample for \f$ (T, \theta) \f$ and performs the necessary bookkeeping to maintain the 'best'
    /// sample.
    void advance();

    /// @brief Equivalent to calling advance().
    void operator () () { advance(); }

    /// @name Access to the current tracker state
    /// @{

    /// @brief The number of iterations performed for the algorithm.
    size_t n_iterations() const { return n_iterations_; }

    /// @brief The last value of \f$ T \f$ drawn by the sampler.
    const partition_ptr_t last_partition() const { return last_partition_sample().partition_sample_ptr; }

    /// @brief The last move type accepted by the partition sampler.
    const mh_moves::move_type last_move_type() const {
        return last_proposal_accepted() ? last_partition_sample().executed_move : mh_moves::NONE;
    }

    /// @brief The last value of \f$ \theta \f$ drawn by the sampler.
    const model::parameters& last_parameters() const { return last_sample().get<1>(); }

    /// @brief The last value of \f$ P(T | \theta) \f$ for the last sample drawn.
    float last_log_pdf() const
        { return underlying_partition_sampler().current_sample_log_density(); }

    /// \brief the current state of the Metropolis Hastings sampler
    sampling::mh_state_observer mh_state() const {
        return underlying_partition_sampler().current_state();
    }

    internal_containers_observer mh_internals() const {
        return underlying_partition_sampler().current_internals();
    }

    /// @brief The current best estimate of the true partition.
    const partition_ptr_t best_partition() const { return best_partition_; }

    /// @brief The current best estimate of the true model parameters.
    const model::parameters& best_parameters() const { return best_parameters_; }

    /// @brief The current value of \f$ P(\hat{T} | \hat{\theta}) \f$.
    float best_log_pdf() const { return best_log_pdf_; }

    float acceptance_rate() const
        { return underlying_partition_sampler().acceptance_rate(); }

    void set_partition_is_valid(bool v)
        { underlying_partition_sampler().set_partition_is_valid(v); }

    bool last_proposal_accepted() const
        { return underlying_partition_sampler().last_proposal_accepted(); }

    /// \brief the partition that was proposed in the latest sample (accepted or not)
    const partition_ptr_t last_proposal_partition() const {
        return underlying_partition_sampler().proposed_sample().partition_sample_ptr;
    }
    /// \brief the move that was proposed in the latest sample (accepted or not)
    mh_moves::move_type last_proposal_move() const {
        return underlying_partition_sampler().proposed_sample().proposed_move;
    }
    partition_sampler_sample last_proposed_sample() const {
        return underlying_partition_sampler().proposed_sample();
    }

    /// \brief the proposal density of the last proposal
    float last_proposal_density() const {
        return underlying_partition_sampler().last_log_density();
    }

    /// \brief the PDR of the last MH comparision
    float last_pdr() const {
        return underlying_partition_sampler().last_pdr();
    }

    /// \brief the recent acceptance rate
    float recent_acceptance_rate() const {
        return underlying_partition_sampler().recent_acceptance_rate();
    }

    /// \brief the log alpha of the last MH comparision
    float last_log_alpha() const {
        return underlying_partition_sampler().last_log_alpha();
    }

    /// @brief Return a reference to the histogram of accepted move types.
    const std::deque<uint64_t>& move_histogram() const { return move_hist_; }

    /// @brief Return a reference to the data of move type statistics (a/r/i).
    const move_stats_t& move_statistics() const { return move_stats_; }

    void set_temperature(float temperature) { underlying_partition_sampler().set_temperature(temperature); }

    /// @}

    /// @brief Assignment operator.
    ///
    /// @param t The tracker to copy state from.
    ///
    /// @return A reference to \c *this.
    const tracker& operator = (const tracker& t);

protected:
    /// @brief The number of iterations which have currently been performed.
    size_t n_iterations_;

    /// @brief Our current estimate of \f$ \hat{T} \f$.
    partition_ptr_t best_partition_;

    /// @brief Our current estimate of \f$ \hat{\theta} \f$.
    model::parameters best_parameters_;

    /// @brief The value of \f$ P(T|\theta) \f$ for our current best estimate.
    float best_log_pdf_;

    /// @brief Histogram of accepted move types.
    std::deque<uint64_t> move_hist_;

    /// @brief data for move statistics ([a]ccepted/[r]ejected/[i]dentity)
    move_stats_t move_stats_;

    /// @brief Record an accepted move type in the move type histogram.
    ///
    /// @param mt
    void record_move_type_(mh_moves::move_type mt);

    /// @brief Record the the statistics ([a]ccepted/[r]ejected/[i]dentity) for each move.
    ///
    /// @param sample
    void record_move_stats_(const partition_sampler_sample& sample);
};

}

#endif // BIGGLES_TRACKER_HPP___
