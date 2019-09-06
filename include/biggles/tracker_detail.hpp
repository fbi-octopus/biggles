/// @file tracker_detail.hpp tracker details

#ifndef BIGGLES_TRACKER_DETAIL_HPP___
#define BIGGLES_TRACKER_DETAIL_HPP___

namespace biggles {

    namespace detail
    {


        /// @brief A conditional sampler for the Biggles Gibbs sampler.
        class tracker_cond_sampler
            : public sampling::gibbs::mem_fun_two_var_cond_sampler<tracker_cond_sampler, partition_sampler_sample, model::parameters>
        {
        public:
            //tracker_cond_sampler() { }

            tracker_cond_sampler(const model::parameters& initial_parameters,
                                 const partition_ptr_t& initial_partition_ptr)
                : partition_sampler_(initial_partition_ptr, initial_parameters)
                , sampled_params_(initial_parameters)
            { }

            tracker_cond_sampler(const tracker_cond_sampler& s)
                : partition_sampler_(s.partition_sampler_)
                , sampled_params_(s.sampled_params_)
            { }

            const tracker_cond_sampler& operator = (const tracker_cond_sampler& s)
            {
                partition_sampler_ = s.partition_sampler_;
                sampled_params_ = s.sampled_params_;
                return *this;
            }

            partition_sampler_sample sample_variable_1(model::parameters params)
            {
                partition_sampler_.set_parameters(params);
                partition_sampler_.draw();
                return partition_sampler_.last_sample();
            }

            model::parameters sample_variable_2(partition_sampler_sample sample)
            {
                sample_model_parameters_given_partition(*sample.partition_sample_ptr, sampled_params_);
                return sampled_params_;
            }

            const partition_sampler& underlying_partition_sampler() const { return partition_sampler_; }
            partition_sampler& underlying_partition_sampler() { return partition_sampler_; }

        protected:
            partition_sampler partition_sampler_;

            model::parameters sampled_params_;
        };
    }
}

#endif
