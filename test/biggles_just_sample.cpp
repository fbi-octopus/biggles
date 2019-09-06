#include "biggles/io.hpp"
#include "biggles/model.hpp"
#include "biggles/partition.hpp"
#include "biggles/tracker.hpp"
#include "biggles/simulate.hpp"
#include "biggles/detail/physics.hpp"
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <fstream>
#include <string>
#include <cmath>

#include "biggles/tools/sundries.hpp"

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

struct part_para {
    float b, c, o, s;
    part_para (float bb, float cc, float oo, float ss) : b(bb), c(cc), o(oo), s(ss) {}
    part_para (const model::parameters& params) :
        b(params.get<0>()), c(params.get<1>()), o(params.get<3>()), s(params.get<2>()) {}
};

bool is_in_01_incl(float f) { return 0.0f <= f and f <= 1.0f; }
bool is_positive(float f) { return f > 0.0f; }

typedef std::pair<float, float> data_pair;

float quartiles(data_iter b, data_iter e, float &lo, float &hi) {
    std::deque<float> data_copy(b, e);
    std::sort(data_copy.begin(), data_copy.end());
    data_iter median_ptr = data_copy.begin();
    int step = data_copy.size()/4;
    std::advance(median_ptr, step);
    lo = *median_ptr;
    std::advance(median_ptr, step);
    float result = *median_ptr;
    std::advance(median_ptr, step);
    hi = *median_ptr;
    return result;
}

bool dump_data(data_iter b, data_iter e, const std::string& fname) {
    try {
        std::ofstream dumpsite(fname.c_str());
        while (b != e) dumpsite << *b++ << std::endl;
    }
    catch (const std::exception &e) {
        return false;
    }
    return true;
}

data_pair generate_partition(const model::parameters& params, partition& dataset_partition, bool shutup) {

    int image_size = sampling::uniform_int(4, 32);
    int time_extent = sampling::uniform_int(5, 200);
    simulate::generate_partition(dataset_partition, params, 0, time_extent, image_size, image_size);

    if (not shutup) {
        diag("Generated %zu new tracks.", dataset_partition.tracks().size());
        diag("Generated %zu clutter observations.", dataset_partition.clutter().size());
    }

    return std::make_pair(image_size, time_extent);
}

bool generate_parameters2(model::parameters& params) {
    model::mean_new_tracks_per_frame(params)           = sampling::uniform_real(5, 10);
    model::mean_false_observations_per_frame(params)   = sampling::uniform_real(5, 10);
    model::frame_to_frame_survival_probability(params) = sqrtf(sampling::uniform_real());
    model::generate_observation_probability(params)    = sqrtf(sampling::uniform_real());
    float sigma_R = sampling::uniform_real();
    model::observation_error_covariance(params)       << sigma_R * sigma_R, 0, 0, sigma_R * sigma_R;
    model::constraint_radius(params) = sampling::uniform_real(0, 10);
    return true;
}

const observation_collection::const_range penultimate_observations(const observation_collection &oc) {
    if (oc.first_time_stamp() == oc.last_time_stamp()) {
        return oc.at_time_stamp(oc.first_time_stamp());
    }
    if (oc.first_time_stamp() == t(oc.back())) {
        return oc.at_time_stamp(oc.first_time_stamp());
    }
    const observation_collection::const_range& end_obs = oc.at_time_stamp(t(oc.back()));
    observation_collection::const_iterator obs_iter = end_obs.first;
    assert(obs_iter != oc.begin());
    --obs_iter;
    return oc.at_time_stamp(t(*obs_iter));
}

void count_obs(const track_collection &tracks, std::map<observation, int> &obs_count) {
    for (track_collection::const_iterator tp = tracks.begin(); tp != tracks.end(); ++tp) {
        const observation_collection &oc = (*tp)->observations();
        for (observation_collection::const_iterator op = oc.begin(); op != oc.end(); ++op) {
            obs_count[*op]++;
        }
    }
}

int main(int argc, char** argv)
{
    plan_no_plan();

    //bool dump = (argc > 1) && (0 == strcmp(argv[1], "-d"));

    biggles::detail::seed_prng(0xfacedead);

    // initialise model parameters
    model::parameters params;


    int n_it(1);
    if(argc > 1 and 0 != strcmp(argv[1], "-d"))
        n_it = atoi(argv[1]);

    diag("number of iterations = %d", n_it);

    //const bool shutup = true;
    const bool tellme = false;
    bool issues = false;

    for (int i=0; i< n_it; ++i) {
        diag(" *** iteration %d *** ", i);
        model::parameters params;
        partition original_partition;
        generate_parameters2(params);
        std::pair<int, int> dims;
        do {
            dims = generate_partition(params, original_partition, tellme);
            simulate::demote_all_tracks_from_partition(original_partition, original_partition);
        } while (not observation_collection_does_allow_tracks(original_partition.clutter(), detail::speed_of_light_));
        tracker tracking(params, original_partition);
        diag("image size=%d, time extent=%d", dims.first, dims.second);
        diag("b=%f, c=%f, o=%f, s=%f",
            model::mean_new_tracks_per_frame(params),
            model::mean_false_observations_per_frame(params),
            model::frame_to_frame_survival_probability(params),
            model::generate_observation_probability(params));
        int bcount(0);
        int dcount(0);
        int scount(0);
        int mcount(0);
        int accepted(0);
        std::map<observation, int> obs_count;
        const int total = 200;
        for(int j=0; j<total; ++j)
        {
            try {
                tracking();
                if (tracking.last_move_type() == mh_moves::BIRTH) {
                    bcount++;
                }
                if (tracking.last_move_type() == mh_moves::DEATH) {
                    dcount++;
                }
                if (tracking.last_move_type() == mh_moves::SPLIT) {
                    scount++;
                }
                if (tracking.last_move_type() == mh_moves::MERGE) {
                    mcount++;
                }
                accepted += int(tracking.last_proposal_accepted());
                count_obs(tracking.last_partition().tracks(), obs_count);

            }
            catch (const std::exception &e) {
                std::cout << e.what() << std::endl;
                std::cout << "iteration " << j << std::endl;
                std::ofstream write_stream("/tmp/biggles_just_sample_parameters.json");
                io::write_model_parameters_to_json_stream(write_stream, params);
                std::ofstream part_stream("/tmp/biggles_just_sample_partition.json");
                io::write_partition_to_json_stream(part_stream, original_partition);
                std::ofstream cur_part_stream("/tmp/biggles_just_sample_last_partition.json");
                std::ofstream cur_para_stream("/tmp/biggles_just_sample_last_parameter.json");
                model::parameters last_params = tracking.last_parameters();
                diag("b=%f, c=%f, o=%f, s=%f",
                    model::mean_new_tracks_per_frame(last_params),
                    model::mean_false_observations_per_frame(last_params),
                    model::frame_to_frame_survival_probability(last_params),
                    model::generate_observation_probability(last_params));
                io::write_partition_to_json_stream(cur_part_stream, tracking.last_partition());
                io::write_model_parameters_to_json_stream(cur_para_stream, tracking.last_parameters());
                issues = true;
                break;
            }

        }
        const std::deque< uint64_t > &move_hist(tracking.move_histogram());
        ok(int(tracking.last_partition().tracks().size()) == bcount + scount - dcount - mcount,
            "tracks.size() = birth + splits - deaths - merges");
        ok1((total - accepted) == int(move_hist[mh_moves::NONE]));
        float acceptance_rate = tracking.acceptance_rate();
        ok1(float(accepted+1)/float(total+1) == acceptance_rate);
        std::multimap<time_stamp, float> time_to_prob;
        for (std::map<observation, int>::const_iterator it = obs_count.begin(); it != obs_count.end(); ++it) {
            time_to_prob.insert(std::make_pair(t(it->first), float(it->second)/float(total)));
        }
        //for (std::multimap<time_stamp, float>::const_iterator it = time_to_prob.begin(); it != time_to_prob.end(); ++it) {
        //    diag("obs at %ld: %.4f", it->first, it->second);
        //}
        if (issues) break;
    }

    ok(not issues, "no issues");


    //return exit_status();
    return 0;
}
