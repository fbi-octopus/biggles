#include "biggles/detail/random.hpp"
#include "biggles/mh_moves/mh_moves.hpp"
#include "biggles/simulate.hpp"
#include "biggles/tools/sundries.hpp"
#include <boost/foreach.hpp>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

int main(int argc, char** argv)
{
    // seed the PRNG to ensure consistent results across runs.
    biggles::detail::seed_prng(0xfacedead);

    plan_tests(4 + 28);

    float sigma_R = 0.1f;

    // initialise model parameters
    model::parameters params;
    model::mean_new_tracks_per_frame(params)           = 1.5f;
    model::mean_false_observations_per_frame(params)   = 2.f;
    model::frame_to_frame_survival_probability(params) = 0.95f;
    model::generate_observation_probability(params)    = 0.8f;
    model::observation_error_covariance(params)       << sigma_R * sigma_R, 0, 0, sigma_R * sigma_R;

    // simulate full data
    partition_ptr_t partition_ptr;
    partition_ptr_t target_ptr;
    simulate::generate_partition(partition_ptr, params, 0, 128 + 128, 128.f, 128.f);

    diag("Generated %zu new tracks.", partition_ptr->tracks().size());
    diag("Generated %zu clutter observations.", partition_ptr->clutter().size());

    // crop partition to a specific region
    simulate::crop_partition(partition_ptr, partition_ptr, 0.f, 128.f, 0.f, 128.f, 64, 64 + 128);

    diag("Cropped to %zu tracks.", partition_ptr->tracks().size());
    diag("Cropped to %zu clutter observations.", partition_ptr->clutter().size());

    // create a synthetic dataset from the partition
    simulate::demote_all_tracks_from_partition(partition_ptr, partition_ptr);

    // check we have some data and no tracks
    ok1(partition_ptr->clutter().size() > 0);
    ok1(partition_ptr->tracks().size() == 0);

    diag("Final data set consists of %zu observations.", partition_ptr->clutter().size());
    size_t input_obs_count(partition_ptr->clutter().size());

    // try some birth moves
    int n_succeeded(0), n_tried(0);
    for(int birth_idx = 0; birth_idx < 300; ++birth_idx)
    {
        float log_pdf(1337.f);
        ++n_tried;
        if(mh_moves::birth(partition_ptr, partition_ptr, log_pdf))
        {
            ++n_succeeded;
            diag("birth move succeeded with log( P(backward) / P(forward) ): %f", log_pdf);
        }
    }
    diag("%i birth move(s) succeeded out of %i attempt(s) (~%i%%)",
         n_succeeded, n_tried, (100 * n_succeeded) / n_tried);
    ok(n_succeeded >= 2, "at least two birth moves succeeded");


    // check that total number of observations hasn't changed
    size_t new_obs_count(partition_ptr->clutter().size());
    BOOST_FOREACH(const boost::shared_ptr<const track>& t_p, partition_ptr->tracks())
    {
        new_obs_count += t_p->size();
    }
    ok(new_obs_count == input_obs_count, "birth moves did not change overall observation count");

    // try some split moves
    //size_t n_to_try = std::min(2 * partition_ptr->tracks().size(), size_t(200));
    size_t n_to_try = 50;
    ok1(n_to_try > 0);

    size_t pre_split_track_count(partition_ptr->tracks().size());
    ok1(pre_split_track_count > 0);

    n_succeeded = 0;
    n_tried = 0;
    for(size_t split_idx = 0; split_idx < n_to_try; ++split_idx)
    {
        float log_pdf(1337.f);
        ++n_tried;
        if(mh_moves::split(partition_ptr, partition_ptr, log_pdf))
        {
            ++n_succeeded;
            diag("split move succeeded with log( P(backward) / P(forward) ): %f", log_pdf);
        }
    }
    diag("%i split move(s) succeeded out of %i attempt(s) (~%i%%)",
         n_succeeded, n_tried, (100 * n_succeeded) / n_tried);
    ok(n_succeeded >= 2, "at least two split moves succeeded");


    ok(partition_ptr->tracks().size() == pre_split_track_count + n_succeeded,
       "each split move added one track");

    // try some death moves
    //n_to_try = std::min(partition_ptr->tracks().size() >> 1, size_t(200));
    ok1(n_to_try > 0);

    size_t pre_death_track_count(partition_ptr->tracks().size());
    ok1(pre_death_track_count > 0);

    n_succeeded = 0;
    n_tried = 0;
    for(size_t death_idx = 0; death_idx < n_to_try; ++death_idx)
    {
        float log_pdf(1337.f);
        ++n_tried;
        if(mh_moves::death(partition_ptr, partition_ptr, log_pdf))
        {
            ++n_succeeded;
            diag("death move succeeded with log( P(backward) / P(forward) ): %f", log_pdf);
        }
    }
    diag("%i death move(s) succeeded out of %i attempt(s) (~%i%%)",
         n_succeeded, n_tried, (100 * n_succeeded) / n_tried);
    ok(n_succeeded >= 2, "at least two death moves succeeded");


    ok(partition_ptr->tracks().size() == pre_death_track_count - n_succeeded,
       "each death move removed one track");

    // check that total number of observations hasn't changed
    new_obs_count = partition_ptr->clutter().size();
    BOOST_FOREACH(const boost::shared_ptr<const track>& t_p, partition_ptr->tracks())
    {
        new_obs_count += t_p->size();
    }
    ok(new_obs_count == input_obs_count, "split moves did not change overall observation count");

    // try some merge moves
    //n_to_try = std::min(partition_ptr->tracks().size() >> 1, size_t(200));
    ok1(n_to_try > 0);

    size_t pre_merge_track_count(partition_ptr->tracks().size());
    ok1(pre_merge_track_count > 0);

    n_succeeded = 0;
    n_tried = 0;
    for(size_t merge_idx = 0; merge_idx < n_to_try; ++merge_idx)
    {
        float log_pdf(1337.f);
        ++n_tried;
        if(mh_moves::merge(partition_ptr, partition_ptr, log_pdf))
        {
            ++n_succeeded;
            diag("merge move succeeded with log( P(backward) / P(forward) ): %f", log_pdf);
        }
    }
    diag("%i merge move(s) succeeded out of %i attempt(s) (~%i%%)",
         n_succeeded, n_tried, (100 * n_succeeded) / n_tried);
    ok(n_succeeded >= 2, "at least two merge moves succeeded");


    ok(partition_ptr->tracks().size() == pre_merge_track_count - n_succeeded,
       "each merge move removed one track");

    // check that total number of observations hasn't changed
    new_obs_count = partition_ptr->clutter().size();
    BOOST_FOREACH(const boost::shared_ptr<const track>& t_p, partition_ptr->tracks())
    {
        new_obs_count += t_p->size();
    }
    ok(new_obs_count == input_obs_count, "merge moves did not change overall observation count");

    // try some extend moves
    //n_to_try = 100;

    size_t pre_extend_track_count = partition_ptr->tracks().size();
    size_t pre_extend_clutter_count = partition_ptr->clutter().size();
    ok1(pre_extend_track_count > 0);

    n_succeeded = 0;
    n_tried = 0;
    for(size_t extend_idx = 0; extend_idx < n_to_try; ++extend_idx)
    {
        float log_pdf(1337.f);
        ++n_tried;
        if(mh_moves::extend(partition_ptr, partition_ptr, log_pdf))
        {
            ++n_succeeded;
            diag("extend move succeeded with log( P(backward) / P(forward) ): %f", log_pdf);
        }
    }
    diag("%i extend move(s) succeeded out of %i attempt(s) (~%i%%)",
         n_succeeded, n_tried, (100 * n_succeeded) / n_tried);
    ok(n_succeeded >= 2, "at least two extend moves succeeded");


    ok(partition_ptr->tracks().size() == pre_extend_track_count, "extend moves retained track count");

    // check that total number of observations hasn't changed
    new_obs_count = partition_ptr->clutter().size();
    BOOST_FOREACH(const boost::shared_ptr<const track>& t_p, partition_ptr->tracks())
    {
        new_obs_count += t_p->size();
    }
    ok(new_obs_count == input_obs_count, "extend moves did not change overall observation count");

    ok(partition_ptr->clutter().size() <= pre_extend_clutter_count - n_succeeded,
       "each extend move updated the clutter by at least one");

    // try some reduce moves
    //n_to_try = 100;

    size_t pre_reduce_track_count = partition_ptr->tracks().size();
    size_t pre_reduce_clutter_count = partition_ptr->clutter().size();
    ok1(pre_reduce_track_count > 0);

    n_succeeded = 0;
    n_tried = 0;
    partition_ptr_t proposed_part_ptr;
    for(size_t reduce_idx = 0; reduce_idx < n_to_try; ++reduce_idx)
    {
        float log_pdf(1337.f);
        ++n_tried;
        if(mh_moves::reduce(partition_ptr, proposed_part_ptr, log_pdf))
        {
            ++n_succeeded;
            diag("reduce move succeeded with log( P(backward) / P(forward) ): %f", log_pdf);
            partition_ptr = proposed_part_ptr;
        }
    }
    diag("%i reduce move(s) succeeded out of %i attempt(s) (~%i%%)",
         n_succeeded, n_tried, (100 * n_succeeded) / n_tried);
    ok(n_succeeded >= 2, "at least two reduce moves succeeded");


    ok(partition_ptr->tracks().size() == pre_reduce_track_count, "reduce moves retained track count");

    // check that total number of observations hasn't changed
    new_obs_count = partition_ptr->clutter().size();
    BOOST_FOREACH(const boost::shared_ptr<const track>& t_p, partition_ptr->tracks())
    {
        new_obs_count += t_p->size();
    }
    ok(new_obs_count == input_obs_count, "reduce moves did not change overall observation count");

    ok(partition_ptr->clutter().size() >= pre_reduce_clutter_count + n_succeeded,
       "each reduce move increased the clutter by at least one");

    // try some update moves
    //n_to_try = 100;

    size_t pre_update_track_count = partition_ptr->tracks().size();
    ok1(pre_update_track_count > 0);

    n_succeeded = 0;
    n_tried = 0;
    for(size_t update_move_idx = 0; update_move_idx < n_to_try; ++update_move_idx)
    {
        float log_pdf(1337.f);
        ++n_tried;
        if(mh_moves::update(partition_ptr, partition_ptr, log_pdf))
        {
            ++n_succeeded;
        }
    }
    diag("%i update move(s) succeeded out of %i attempt(s) (~%i%%)",
         n_succeeded, n_tried, (100 * n_succeeded) / n_tried);
    ok(n_succeeded >= 2, "at least two update moves succeeded");


    ok(partition_ptr->tracks().size() == pre_update_track_count, "update moves retained track count");

    // check that total number of observations hasn't changed
    new_obs_count = partition_ptr->clutter().size();
    BOOST_FOREACH(const boost::shared_ptr<const track>& t_p, partition_ptr->tracks())
    {
        new_obs_count += t_p->size();
    }
    ok(new_obs_count == input_obs_count, "update moves did not change overall observation count");

    diag("Using propose() to perform one move of each type.");
    float ratio(0.f);
    //mh_moves::propose(mh_moves::NONE, partition_ptr, partition_ptr, ratio);

    for (int i = 1; i != int(mh_moves::IDENTITY); ++i) {
        mh_moves::move_type mv = mh_moves::move_type(i);
        diag("Move type = %s", mh_moves::move_name(mv).c_str());
        mh_moves::propose(mv, partition_ptr, proposed_part_ptr, ratio);
    }
    /*
    mh_moves::propose(mh_moves::BIRTH, partition_ptr, partition_ptr, ratio);
    mh_moves::propose(mh_moves::DEATH, partition_ptr, partition_ptr, ratio);
    mh_moves::propose(mh_moves::EXTEND, partition_ptr, partition_ptr, ratio);
    mh_moves::propose(mh_moves::REDUCE, partition_ptr, partition_ptr, ratio);
    mh_moves::propose(mh_moves::SPLIT, partition_ptr, partition_ptr, ratio);
    mh_moves::propose(mh_moves::MERGE, partition_ptr, partition_ptr, ratio);
    mh_moves::propose(mh_moves::UPDATE, partition_ptr, partition_ptr, ratio);
    */

    return exit_status();
}
