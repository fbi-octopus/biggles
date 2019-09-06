#include "sampling/metropolis_hastings.hpp"
namespace biggles { namespace sampling {

const partition_sampler_sample& metropolis_hastings_sampler::draw()
{
    typedef partition_proposal::result_type proposal_type;


    // calc P(omega|theta, data)
    float old_log_density = target_(sample_);

    // draw a proposal given our current sample
    sample_.partition_sample_ptr->set_capability_recorder(internal_records_ptr_);
    proposal_type proposal(propose_(sample_));
    last_proposal_ = proposal;

    // evaluate the target PDF at the proposed sample
    // PP(omega*|theta, data)
    float new_log_density = target_(proposal.partition_sample);
    if (not std::isfinite(new_log_density)) {
        std::cerr << new_log_density << std::endl;
        throw  std::runtime_error("new_log_density not finite");
    }

    if (not std::isfinite(old_log_density)) {
        std::cerr << old_log_density << std::endl;
        throw  std::runtime_error("old_log_density not finite");
    }

    // calculate threshold
    //float log_threshold = new_log_density - sample_log_density_ + proposal.template get<1>();
    float log_threshold = new_log_density - old_log_density + proposal.proposal_mass_ratio;

    // draw some uniform sample
    //float alpha = uniform_();
    float alpha = sampling::uniform_real();

    state_.last_log_alpha = logf(alpha);
    state_.last_proposal_log_density = new_log_density;
    state_.last_sample_log_density = old_log_density;
    state_.last_pdr = proposal.proposal_mass_ratio;

    // do we accept?
    if(logf(alpha) < log_threshold)
    {
        // ... yes, update state
        // storing the proposed sample
        sample_ = proposal.partition_sample;
        state_.sample_log_density = new_log_density;

        // increment accepted count
        ++state_.n_accepted;
        state_.accepted = true;
        acceptance_seq_.push(1);
        if (sample_.executed_move != mh_moves::IDENTITY) {
            internal_records_ptr_->commit_changes(sample_.partition_sample_ptr->tracks());
        } else {
            BOOST_ASSERT(new_log_density = old_log_density);
        }

    }
    else
    {
        state_.sample_log_density = old_log_density;
        state_.accepted = false;
        acceptance_seq_.push(0);
    }

    // FIXME: refacture that. put it to the proper place
    //intern_.update(*sample_.partition_sample_ptr);

    // increment proposal count
    ++state_.n_proposed;
    intern_.update(*internal_records_ptr_);

    //std::cout << mh_moves::move_name((proposal.template get<0>()).proposed_move) << " @@ ";

    return last_sample();
}

} }
