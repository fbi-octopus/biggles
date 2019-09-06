#include "biggles/model.hpp"
#include "biggles/samplers.hpp"
#include "biggles/partition_sampler.hpp"
#include "biggles/simulate.hpp"
#include <boost/foreach.hpp>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>

#include "biggles/tools/sundries.hpp"

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

int main(int argc, char** argv)
{
    plan_tests(2);

    float sigma_R = 0.1f;

    // initialise model parameters
    model::parameters params;
    model::mean_new_tracks_per_frame(params)           = 3.f;
    model::mean_false_observations_per_frame(params)   = 1.f;
    model::frame_to_frame_survival_probability(params) = 0.95f;
    model::generate_observation_probability(params)    = 0.9f;
    model::observation_error_covariance(params)       << sigma_R * sigma_R, 0, 0, sigma_R * sigma_R;
    model::process_noise_covariance(params) = detail::initQ();

    float image_size(128.f);
    time_stamp time_extent(128);

    // simulate full data
    partition_ptr_t partition_ptr;
    simulate::generate_partition(partition_ptr, params, 0, (time_extent>>1) + time_extent, image_size * 1.5f, image_size * 1.5f);
    diag("Partition volume = %.2f", partition_ptr->volume());

    diag("Generated %zu new tracks.", partition_ptr->tracks().size());
    diag("Generated %zu clutter observations.", partition_ptr->clutter().size());

    // crop partition to a specific region
    simulate::crop_partition(partition_ptr, partition_ptr,
                             image_size * 0.25f, image_size * 1.25f,
                             image_size * 0.25f, image_size * 1.25f,
                             time_extent >> 2, (time_extent >> 2) + time_extent);

    diag("Cropped to %zu tracks.", partition_ptr->tracks().size());
    diag("Cropped to %zu clutter observations.", partition_ptr->clutter().size());
    diag("Ground truth log pdf: %f", model::log_partition_given_parameters_and_data_density(*partition_ptr, params));
    diag("Partition volume = %.2f", partition_ptr->volume());

    // create a synthetic dataset from the partition
    simulate::demote_all_tracks_from_partition(partition_ptr, partition_ptr);

    // check we have some data and no tracks
    ok1(partition_ptr->clutter().size() > 0);
    ok1(partition_ptr->tracks().size() == 0);

    diag("Final data set consists of %zu observations.", partition_ptr->clutter().size());

    partition_sampler sample_partition(partition_ptr, params);
    partition_ptr_t best_partition_ptr(sample_partition.last_partition());
    float best_log_pdf(sample_partition.current_sample_log_density());

    int n_it = 4000;
    if(argc > 1)
    {
        n_it = atoi(argv[1]);
    }

    std::ofstream stats("biggles_partition_sampler_stats.txt");
    for(int i=0; i<n_it; ++i)
    {
        sample_partition.draw();

        if(sample_partition.current_sample_log_density() > best_log_pdf)
        {
            best_partition_ptr = sample_partition.last_partition();
            best_log_pdf = sample_partition.current_sample_log_density();
        }

        if(0 == (i & 0xff))
        {
            diag("i: %i, alpha: %f, track count: %zu, log PDF: %f, current: %f",
                 i,
                 sample_partition.acceptance_rate(),
                 best_partition_ptr->tracks().size(),
                 best_log_pdf,
                 sample_partition.current_sample_log_density());
            diag("sample partition volume = %.2f", sample_partition.last_partition()->volume());
        }

        if(0 == (i & 0xfff))
        {
            write_partition("partition_sampler", *best_partition_ptr, "post_sampling");
        }

        stats << std::setw(16) << sample_partition.acceptance_rate()
              << std::setw(16) << best_log_pdf
              << std::setw(16) << best_partition_ptr->tracks().size()
              << '\n';
    }

    // write out post-sampling dataset
    write_partition("partition_sampler", *best_partition_ptr, "post_sampling");

    return exit_status();
}
