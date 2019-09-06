#include "partition_sampler.hpp"

namespace biggles {

partition_sampler_result_t partition_proposal::operator () (
    const partition_sampler_sample& sample) const
{
    partition_ptr_t new_partition_ptr;
    float log_prob_ratio(0.f);

    mh_moves::move_type proposed_move_type(mh_moves::select_move_type(*sample.partition_sample_ptr));
    bool succeeded(mh_moves::propose(proposed_move_type, sample.partition_sample_ptr, new_partition_ptr, log_prob_ratio));
    if(!succeeded)
    {
        //log_prob_ratio = -std::numeric_limits<float>::max();
        log_prob_ratio = 0.f;
        new_partition_ptr = sample.partition_sample_ptr;
    }

    if(succeeded && !verify_partition(*new_partition_ptr))
    {
        std::cerr << "verification failed after move type " << mh_moves::move_sign(proposed_move_type) << std::endl;
        throw std::runtime_error("partition_proposal::operator(): Partition verification failed");
    }

    partition_sampler_sample new_sample(new_partition_ptr, proposed_move_type,
        succeeded ? proposed_move_type : mh_moves::IDENTITY);

    return partition_sampler_result_t(new_sample, log_prob_ratio);

}

} // namespace biggles
