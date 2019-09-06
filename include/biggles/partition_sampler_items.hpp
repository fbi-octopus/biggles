#ifndef BIGGLES_PARTITION_SAMPLER_ITEMS_HPP__
#define BIGGLES_PARTITION_SAMPLER_ITEMS_HPP__

#include <boost/tuple/tuple.hpp>

#include "model.hpp"
#include "partition.hpp"
#include "mh_moves/mh_moves.hpp"


namespace biggles {
/// @brief The sample type for the partition sampler.
///
/// Something as a hack, the samples drawn by the partition sampler are the partition itself <em>and</em> the proposal
/// move which generated it.
struct partition_sampler_sample {
    partition_ptr_t partition_sample_ptr;
    mh_moves::move_type proposed_move; /// \brief the move that has been attempted
    mh_moves::move_type executed_move; /// \brief could be the proposed move, identity or maybe none
    partition_sampler_sample () :
         partition_sample_ptr(new partition()), proposed_move(mh_moves::NONE), executed_move(mh_moves::NONE) {}
    partition_sampler_sample (
        const partition_ptr_t& part, const mh_moves::move_type prop_move, const mh_moves::move_type exec_move)
        : partition_sample_ptr(part), proposed_move(prop_move), executed_move(exec_move) {};
    explicit partition_sampler_sample (const partition_ptr_t& part)
        : partition_sample_ptr(part), proposed_move(mh_moves::NONE), executed_move(mh_moves::NONE) {};
};

struct partition_sampler_result_t {
    partition_sampler_sample partition_sample;
    float proposal_mass_ratio;
    partition_sampler_result_t(const partition_sampler_sample& sample, const float& pmr) :
        partition_sample(sample), proposal_mass_ratio(pmr) {}
};

/// @brief A proposal function for the Biggles sampler.
struct partition_proposal : public std::unary_function<const partition&, partition_sampler_result_t >
{
    /// @brief we assume that the partition is valid
    partition_proposal() : partition_is_valid(true) {}

    /// @brief Propose a new sample given an input partition.
    ///
    /// @param p
    partition_sampler_result_t operator () (const partition_sampler_sample& p) const;

    /// @brief is the partition OK for making valid proposals
    ///
    /// if it is not posibble to make a valid proposal for partition an infinite loop will occur
    bool partition_is_valid;
};

/// @brief Evaluate the partition posterior ditribution.
///
struct partition_distribution : public std::unary_function<const partition_sampler_sample&, float>
{
    partition_distribution(const model::parameters& p = model::parameters()) : params(p) { }

    float operator () (const partition_sampler_sample& sample) const {
        return model::log_partition_given_parameters_and_data_density(*sample.partition_sample_ptr, params);
    }

    model::parameters params;
};


}

#endif
