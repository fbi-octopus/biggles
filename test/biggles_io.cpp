#include "biggles/io.hpp"
#include "biggles/observation.hpp"
#include "biggles/partition.hpp"
#include "biggles/simulate.hpp"
#include "biggles/track.hpp"
#include "biggles/detail/physics.hpp"
#include <boost/foreach.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>
#include <iostream>
#include <sstream>
#include <stdlib.h>

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

boost::shared_ptr<const partition> simulate_partition(const model::parameters& params)
{

    // simulate full data
    partition_ptr_t p(new partition());
    simulate::generate_partition(p, params, 0, 128 + 128, 128.f, 128.f);

    diag("Generated %zu new tracks.", p->tracks().size());
    diag("Generated %zu clutter observations.", p->clutter().size());

    // crop partition to a specific region
    simulate::crop_partition(p, p, 0.f, 128.f, 0.f, 128.f, 64, 64 + 128);

    diag("Cropped to %zu tracks.", p->tracks().size());
    diag("Cropped to %zu clutter observations.", p->clutter().size());

    return p;
}

int main(int argc, char** argv)
{
    plan_tests(2 + 2 + 1 + 3);

    // initialise model parameters
    model::parameters params;
    model::mean_new_tracks_per_frame(params)           = 1.5f;
    model::mean_false_observations_per_frame(params)   = 2.f;
    model::frame_to_frame_survival_probability(params) = 0.95f;
    model::generate_observation_probability(params)    = 0.8f;
    float sigma_R = 0.1f;
    model::observation_error_covariance(params)       << sigma_R * sigma_R, 0, 0, sigma_R * sigma_R;
    model::process_noise_covariance(params) = detail::initQ();

    {
        // serialise it to a string
        std::stringstream write_stream;
        io::write_model_parameters_to_json_stream(write_stream, params);
        std::string serialised_model_parameters(write_stream.str());

        // check length of string is sensible
        diag("Serialised representation is %zu characters.", serialised_model_parameters.size());
        ok1(serialised_model_parameters.size() > 10);

        // de-serialise
        std::stringstream read_stream(serialised_model_parameters);
        boost::shared_ptr<const model::parameters> recovered_model_parameters_ptr(
            io::read_model_parameters_from_json_stream(read_stream));

        // check equality (dangerous!)
        ok1(*recovered_model_parameters_ptr == params);
    }

    // create a simulated partition
    boost::shared_ptr<const partition> source_partition_ptr(simulate_partition(params));
    ok1(source_partition_ptr->clutter().size() > 0);
    ok1(source_partition_ptr->tracks().size() > 0);

    {
        // serialise it to a string
        std::stringstream write_stream;
        io::write_partition_to_json_stream(write_stream, *source_partition_ptr);
        std::string serialised_partition(write_stream.str());

        // check length of string is sensible
        diag("Serialised representation is %zu characters.", serialised_partition.size());
        ok1(serialised_partition.size() > 10);

        // de-serialise
        std::stringstream read_stream(serialised_partition);
        boost::shared_ptr<const partition> recovered_partition_ptr(
            io::read_partition_from_json_stream(read_stream));

        // check that this is a different partition in memory but has the same shape
        ok1(source_partition_ptr != recovered_partition_ptr);
        ok1(recovered_partition_ptr->clutter().size() == source_partition_ptr->clutter().size());
        ok1(recovered_partition_ptr->tracks().size() == source_partition_ptr->tracks().size());

        // FIXME: be a little more thorough in the checking
    }

    return exit_status();
}
