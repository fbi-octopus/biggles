#include "biggles/model.hpp"
#include "biggles/partition.hpp"
#include "biggles/tracker.hpp"
#include "biggles/simulate.hpp"
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/tuple/tuple_io.hpp>
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

std::vector<std::string> move_names(mh_moves::MOVE_COUNT);

struct part_para {
    float b, c, o, s;
    part_para (float bb, float cc, float oo, float ss) : b(bb), c(cc), o(oo), s(ss) {}
};

int test1(int n_it, model::parameters &params)
{
    float image_size(32.f);
    time_stamp time_extent(128);

    // simulate full data
    partition_ptr_t partition_ptr;
    simulate::generate_partition(partition_ptr, params, 0, (time_extent>>1) + time_extent, image_size * 1.5f, image_size * 1.5f);

    diag("Generated %zu new tracks.", partition_ptr->tracks().size());
    diag("Generated %zu clutter observations.", partition_ptr->clutter().size());

    // crop partition to a specific region
    simulate::crop_partition(partition_ptr, partition_ptr,
                             image_size * 0.25f, image_size * 1.25f,
                             image_size * 0.25f, image_size * 1.25f,
                             time_extent >> 2, (time_extent >> 2) + time_extent);

    partition_ptr_t ground_truth(new partition(*partition_ptr));

    diag("Cropped to %zu tracks.", partition_ptr->tracks().size());
    diag("Cropped to %zu clutter observations.", partition_ptr->clutter().size());
    diag("Ground truth log pdf: %f", model::log_partition_given_parameters_and_data_density(*partition_ptr, params));

    BOOST_FOREACH(const boost::shared_ptr<const track>& t_ptr, partition_ptr->tracks())
    {
        diag("P(t|theta): %f", model::log_track_given_parameters_density(t_ptr, params));
    }

    diag("P(clutter|theta): %f", model::log_clutter_given_parameters_density(*partition_ptr, params));
    diag("P(T|theta): %f", model::log_partition_given_parameters_density(*partition_ptr, params));
    diag("Priors: %f", model::log_parameters_prior_density(params));


    // create a synthetic dataset from the partition
    simulate::demote_all_tracks_from_partition(partition_ptr, partition_ptr);

    // check we have some data and no tracks
    ok1(partition_ptr->clutter().size() > 0);
    ok1(partition_ptr->tracks().size() == 0);

    diag("P(T|theta): %f", model::log_partition_given_parameters_density(*partition_ptr, params));

    diag("Final data set consists of %zu observations.", partition_ptr->clutter().size());

    // create tracker
    tracker tracking(params, partition_ptr);
    //tracker tracking(params, ground_truth);
    for(int i=0; i<n_it; ++i)
    {
        tracking();

        if(0 == (i&0xff))
        {
            diag("i: %i, track count: %zu, log PDF: %f, last PDF: %f",
                 i,
                 tracking.best_partition()->tracks().size(),
                 tracking.best_log_pdf(),
                 tracking.last_log_pdf());

            model::parameters theta(tracking.best_parameters());
            diag("parameters: lambda_b = %f, lambda_f = %f, p_s = %f, p_d = %f",
                model::birth_rate(theta),
                model::clutter_rate(theta),
                model::observation_probability(theta),
                model::survival_probability(theta));
        }

        if(0 == (i&0x1ff))
            write_partition("tracker", *tracking.best_partition(), "estimate");
    }

    write_partition("tracker", *tracking.best_partition(), "estimate");
    return 0;
}

int test2(int n_it, model::parameters &params)
{
    observation_collection track1_obs;
    track1_obs.insert(new_obs(0, 0, 0));
    track1_obs.insert(new_obs(0, 0, 2));
    track1_obs.insert(new_obs(0, 0, 4));
    track1_obs.insert(new_obs(0, 0, 6));
    track1_obs.insert(new_obs(0, 0, 8));

    observation_collection track2_obs;
    track2_obs.insert(new_obs(0, 0, 1));
    track2_obs.insert(new_obs(0, 0, 3));
    track2_obs.insert(new_obs(0, 0, 5));
    track2_obs.insert(new_obs(0, 0, 7));
    track2_obs.insert(new_obs(0, 0, 9));

    observation_collection clutter_obs;
    clutter_obs.insert(new_obs(9, 9, 1));

    boost::shared_ptr<track_collection> tracks2(new track_collection());
    clutter_ptr clutter(new clutter_t(clutter_obs.begin(), clutter_obs.end()));

    track track1(0, 10, track1_obs.begin(), track1_obs.end(), 1.0);
    track track2(0, 10, track2_obs.begin(), track2_obs.end(), 1.0);

    tracks2->insert(track1);
    tracks2->insert(track2);

    partition_ptr_t partition_ptr(new partition(tracks2, clutter));

    model::mean_new_tracks_per_frame(params)           = 0.2f;
    model::mean_false_observations_per_frame(params)   = 0.1f;
    model::frame_to_frame_survival_probability(params) = 0.95f;
    model::generate_observation_probability(params)    = 0.5f;

    BOOST_FOREACH(const boost::shared_ptr<const track>& t_ptr, partition_ptr->tracks())
    {
        diag("P(t|theta): %f", model::log_track_given_parameters_density(t_ptr, params));
    }

    diag("P(clutter|theta): %f", model::log_clutter_given_parameters_density(*partition_ptr, params));
    diag("P(T|theta): %f", model::log_partition_given_parameters_density(*partition_ptr, params));
    diag("Priors: %f", model::log_parameters_prior_density(params));

    diag("Original partition log pdf: %f", model::log_partition_given_parameters_and_data_density(*partition_ptr, params));

    diag("** starting tracking **");

    tracker trackit(params, partition_ptr);

    std::deque<mh_moves::move_type> move_seq;
    std::deque<mh_moves::move_type> prop_seq;
    std::deque<int> num_tr_seq;
    std::deque<int> longest_dur_seq;
    std::deque<int> longest_obs_seq;
    std::deque<int> sample_id_seq;
    std::deque<int> clutter_size_seq;
    std::deque<float> log_pdf_seq;
    std::deque<part_para> para_seq;

    for(int i=0; i<n_it; ++i)
    {
        try {
            trackit();
        }
        catch (const std::exception& e) {
            for (size_t k = 0; k < move_seq.size(); ++k ) {
                std::cerr
                    << "  "
                    << mh_moves::move_sign(prop_seq[k])
                    << "  "
                    << mh_moves::move_sign(move_seq[k])
                    << std::endl
                ;
            }
            std::cerr << "iteration i=" << i << std::endl;
            throw;
        }
        prop_seq.push_back(trackit.last_proposal_move());
        if (trackit.last_proposal_accepted()) {
            //diag("%s", mh_moves::move_name(trackit.last_move_type()).c_str());
            sample_id_seq.push_back(i+1);
            move_seq.push_back(trackit.last_move_type());
            track_collection tracks = trackit.last_partition()->tracks();
            clutter_size_seq.push_back(trackit.last_partition()->clutter().size());
            num_tr_seq.push_back(tracks.size());
            track_collection::const_iterator iter = tracks.begin();
            int maxdur = 0;
            int maxobs = 0;
            for (; iter != tracks.end(); ++iter) {
                if ((*iter)->duration() > maxdur) {
                    maxdur = (*iter)->duration();
                    maxobs = (*iter)->size();
                }
            }
            longest_dur_seq.push_back(maxdur);
            longest_obs_seq.push_back(maxobs);
            log_pdf_seq.push_back(trackit.last_log_pdf());
            model::parameters params(trackit.last_parameters());
            para_seq.push_back(part_para(model::birth_rate(params), model::clutter_rate(params),
                model::observation_probability(params), model::survival_probability(params)));
        } else {
            move_seq.push_back(mh_moves::NONE);
        }

    }
    partition_ptr_t last_partition = trackit.last_partition();
    track_collection::const_iterator iter = last_partition->tracks().begin();
    for (; iter != last_partition->tracks().end(); ++iter) {
        diag("duration %zu, observations %zu", (*iter)->duration(), (*iter)->size());
    }
    const std::deque< uint64_t >& move_histogram = trackit.move_histogram();
    for (size_t i = 0; i < move_histogram.size(); ++i)
        diag("move: %s, count: %zu", move_names[i].c_str(), move_histogram[i]);

    /*
    for (int i = 0; i < move_seq.size(); ++i) {
        std::stringstream sstr;
        sstr << boost::format("%4d #%1d/%1d (%2d, %2d) log PDF %.2f ; b=%.2f, c=%.2f, o=%.2f, s=%.2f")
            % sample_id_seq[i]
            % num_tr_seq[i]
            % clutter_size_seq[i]
            % longest_dur_seq[i]
            % longest_obs_seq[i]
            % log_pdf_seq[i]
            % para_seq[i].b % para_seq[i].c % para_seq[i].o % para_seq[i].s
        ;
        diag("%s", sstr.str().c_str());
    }
    */
    return 0;
}

int main(int argc, char** argv)
{
    plan_tests(1);

    const float sigma_R = 0.1f;

    //biggles::detail::seed_prng(0xfacedead);

    // initialise model parameters
    model::parameters params;
    model::mean_new_tracks_per_frame(params)           = 1.f;
    model::mean_false_observations_per_frame(params)   = 1.f;
    model::frame_to_frame_survival_probability(params) = 0.95f;
    model::generate_observation_probability(params)    = 0.9f;
    model::observation_error_covariance(params)       << sigma_R * sigma_R, 0, 0, sigma_R * sigma_R;
    model::constraint_radius(params) = 100.0f;
    model::process_noise_covariance(params) = detail::initQ();


    int n_it(4000);
    if(argc > 1)
        n_it = atoi(argv[1]);

    //test1(n_it, params);
    test2(n_it, params);
    ok1(true);

    return exit_status();
}
