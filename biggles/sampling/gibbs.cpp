#include "biggles/sampling/gibbs.hpp"
#include "biggles/detail/constants.hpp"
#include "biggles/samplers.hpp"

namespace biggles{
namespace sampling {
	const size_t zillion(std::numeric_limits<size_t>::max());
    new_gibbs_sampler::new_gibbs_sampler(const partition_ptr_t& initial_partition_ptr,
                                         const model::parameters& initial_parameters)
        : observation_error_lag_(1), gibbs_sample_count_(0)
        , partition_sample_(initial_partition_ptr), parameter_sample_(initial_parameters)
        , partition_sampler_(initial_partition_ptr, initial_parameters)
    {
        BOOST_ASSERT(model::process_noise_covariance(initial_parameters) != matrix4f::Zero());
    }
    tracker_sample new_gibbs_sampler::draw() {
        try {
            partition_sample_ = partition_sampler_.draw();
        } catch (const std::exception& e) {
            std::cerr << "partition_sampler" << std::endl;
            throw;
        }
        try {
            sample_tracking_control_parameters_given_partition(*partition_sample_.partition_sample_ptr, parameter_sample_);
        } catch (const std::exception& e) {
            std::cerr << "parameter sampler" << std::endl;
            throw;
        }
        try {
            if (gibbs_sample_count_ >= observation_error_lag_) {
                sample_observation_error_given_partition(*partition_sample_.partition_sample_ptr, parameter_sample_);
                gibbs_sample_count_ = 0;
            }
        } catch (const std::exception& e) {
            std::cerr << "R sampler" << std::endl;
            throw;
        }
        ++gibbs_sample_count_;
        partition_sampler_.set_parameters(parameter_sample_);
        return last_sample();
    }
    /// \brief returns the complete last sample
    tracker_sample new_gibbs_sampler::last_sample() const {
        return tracker_sample(partition_sample_, parameter_sample_);
    }
    /// \brief returns the last partition sample
    const partition_sampler_sample& new_gibbs_sampler::last_partition_sample() const {
        return partition_sample_;
    }
    /// \brief returns the last parameter sample
    const model::parameters& new_gibbs_sampler::last_parameter_sample() const {
        return parameter_sample_;
    }
    void new_gibbs_sampler::fix_observation_error(const float R00, const float R01, const float R11) {
        observation_error_lag_ = zillion;
        gibbs_sample_count_ = 0;
        model::observation_error_covariance(parameter_sample_) << R00, R01, R01, R11;
    }
    void new_gibbs_sampler::sample_observation_error(const size_t lag) {
        observation_error_lag_ = lag;
        gibbs_sample_count_ = 0;
    }

} } // namespace biggles::sampling
