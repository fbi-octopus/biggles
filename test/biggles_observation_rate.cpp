#include "biggles/observation.hpp"
#include "biggles/detail/random.hpp"
#include "biggles/mh_moves/mh_moves.hpp"
#include "biggles/mh_moves/utility.hpp"
#include "biggles/tracker.hpp"
#include "biggles/simulate.hpp"
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <string>
#include "biggles/tools/sundries.hpp"

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

class select_move_from {
    std::deque<mh_moves::move_type> m_data;
    int m_size;
public:
    select_move_from &operator<<(const mh_moves::move_type &type) {
        m_data.push_back(type);
        m_size = int(m_data.size());
        return *this;
    }
    mh_moves::move_type operator()() const {
        return m_data[biggles::sampling::uniform_int(0, m_size)];
    }
};

mh_moves::move_type simple_metropolis_hastings(partition_ptr_t &partition_ptr, const model::parameters& para,
        const select_move_from &selector)
{
    const capability_recorder_ptr &cap_rec_ptr(partition_ptr->get_capability_recorder());
    mh_moves::move_type this_move(selector());
    biggles_move proposal(mh_moves::move_name(this_move));
    float pdr;
    partition_ptr_t proposed_sample_ptr;
    bool success = proposal(partition_ptr, proposed_sample_ptr, pdr);
    if (not success) {
        return mh_moves::IDENTITY;
    }

    float proposed_density = model::log_partition_given_parameters_density(*proposed_sample_ptr, para);
    float current_density = model::log_partition_given_parameters_density(*partition_ptr, para);
    float mh_ratio = proposed_density - current_density + pdr;
    bool accept = logf(sampling::uniform_real()) < mh_ratio;
    if (accept) {
        partition_ptr = proposed_sample_ptr;
        cap_rec_ptr->commit_changes(partition_ptr->tracks());
    }
    else {
        this_move = mh_moves::NONE;
    }
    return this_move;
}

std::ostream &get_ostream(const std::string &fname, std::ofstream &fout) {
    if (fname.size() > 0) {
        fout.open(fname.c_str());
        if (fout.good()) {
            return fout;
        }
    }
    return std::cout;
}


void log_part_for_para_sampler(partition &partition_ptr, const general_paras &opts) {
    model::parameters para_sample;
    std::string fname(boost::filesystem::temp_directory_path().string());
    std::ofstream fh((fname + "/cpp_biggles.log_part_for_para.txt").c_str());
    for (int i = 0; i < opts.total; ++i) {
        sample_model_parameters_given_partition(partition_ptr, para_sample);
        fh
            << model::log_partition_given_parameters_density(partition_ptr, para_sample)
            << std::endl
        ;
    }
}

void parameter_sampler(partition &partition_ptr, const general_paras &opts) {
    model::parameters para_sample;
    std::ofstream fh;
    std::ostream &out = get_ostream(opts.rgp.output, fh);
    for (int i = 0; i < opts.total; ++i) {
        sample_model_parameters_given_partition(partition_ptr, para_sample);
        out << model::mean_new_tracks_per_frame(para_sample) << ","
            << model::mean_false_observations_per_frame(para_sample) << ","
            << model::generate_observation_probability(para_sample) << ","
            << model::frame_to_frame_survival_probability(para_sample)
            << std::endl
        ;
    }
}

bool mh_only(partition_ptr_t &partition_ptr, model::parameters& para_sample, const general_paras &opts) {
    int half_time = opts.total/2;
    mh_moves::move_type last_move;
    select_move_from selector;
    std::map<mh_moves::move_type, size_t> move_hist;
    std::map<mh_moves::move_type, size_t>::iterator it;
    selector << mh_moves::BIRTH << mh_moves::DEATH << mh_moves::UPDATE
        << mh_moves::MERGE << mh_moves::SPLIT << mh_moves::EXTEND << mh_moves::REDUCE;
    diag("parameters: %s", parameters_to_string(para_sample).c_str());

    for (int i = 0; i < half_time; ++i) {
        last_move = simple_metropolis_hastings(partition_ptr, para_sample, selector);
        ++move_hist[last_move];
    }
    for (it = move_hist.begin(); it not_eq move_hist.end(); ++it) {
        diag("%8s => %6.2f%%", mh_moves::move_name(it->first).c_str(), 100.f*(it->second)/half_time);
    }
    return true;
}

bool simple_gibbs(partition_ptr_t partition_ptr, model::parameters& para_sample, const general_paras &opts) {
    int half_time = opts.total/2;
    mh_moves::move_type last_move;
    select_move_from selector;
    std::map<mh_moves::move_type, size_t> move_hist;
    std::map<mh_moves::move_type, size_t>::iterator it;
    //selector << mh_moves::BIRTH << mh_moves::DEATH << mh_moves::UPDATE;
    selector << mh_moves::BIRTH << mh_moves::DEATH << mh_moves::UPDATE
        << mh_moves::MERGE << mh_moves::SPLIT << mh_moves::EXTEND << mh_moves::REDUCE;
    diag("parameters: %s", parameters_to_string(para_sample).c_str());
    capability_recorder_ptr cap_rec_ptr(new_cap_recorder(*partition_ptr));

    for (int i = 0; i < half_time; ++i) {
        partition_ptr->set_capability_recorder(cap_rec_ptr);
        last_move = simple_metropolis_hastings(partition_ptr, para_sample, selector);
        ++move_hist[last_move];
        sample_model_parameters_given_partition(*partition_ptr, para_sample);

    }
    for (it = move_hist.begin(); it not_eq move_hist.end(); ++it) {
        diag("%8s => %6.2f%%", mh_moves::move_name(it->first).c_str(), 100.f*(it->second)/half_time);
    }
    move_hist.clear();
    bool write_parameters = opts.rgp.output.size() > 0;
    std::ofstream fh;
    std::ostream &out = get_ostream(opts.rgp.output, fh);
    int generic_ctr = 0;
    for (int i = 0; i < half_time; ++i) {
        partition_ptr->set_capability_recorder(cap_rec_ptr);
        last_move = simple_metropolis_hastings(partition_ptr, para_sample, selector);
        sample_model_parameters_given_partition(*partition_ptr, para_sample);
        ++move_hist[last_move];
        if (partition_ptr->tracks().size() == 0)
            generic_ctr++;
        if (write_parameters) out << model::mean_new_tracks_per_frame(para_sample) << ","
            << model::mean_false_observations_per_frame(para_sample) << ","
            << model::generate_observation_probability(para_sample) << ","
            << model::frame_to_frame_survival_probability(para_sample)
            << std::endl
        ;
    }
    for (it = move_hist.begin(); it not_eq move_hist.end(); ++it) {
        diag("%8s => %6.2f%%", mh_moves::move_name(it->first).c_str(), 100.f*(it->second)/half_time);
    }
    diag("trivial proportion %6.2f%%", 100.0f*float(generic_ctr)/float(half_time));
    return true;
}

int main(int argc, char** argv)
{
    general_paras gen_pars;
    gen_pars.total = 200000;
    gen_pars.seed = 0;
    gen_pars.opts = "gibbs";
    gen_pars.rgp.lambda = 2.5;
    gen_pars.rgp.p_no = 0.2;
    gen_pars.rgp.p_yes = 0.2;
    gen_pars.rgp.p_tr = 0.5;
    gen_pars.rgp.min_tracks = 1;
    gen_pars.rgp.max_tracks = 2;
    if (not parse_args(argc, argv, gen_pars))
        return exit_status();
    biggles::detail::seed_prng(gen_pars.seed);

    plan_no_plan();

    model::parameters para_sample;
    generate_parameters(para_sample);
    diag("paramters: %s", parameters_to_string(para_sample).c_str());
    partition_ptr_t partition_ptr(random_debug_partition(gen_pars.rgp));
    diag("opts = %s, total = %d", gen_pars.opts.c_str(), gen_pars.total);
    if (gen_pars.opts == "part")
        log_part_for_para_sampler(*partition_ptr, gen_pars);
    else if (gen_pars.opts == "part")
        parameter_sampler(*partition_ptr, gen_pars);
    else if (gen_pars.opts == "gibbs")
        simple_gibbs(partition_ptr, para_sample, gen_pars);
    else if (gen_pars.opts == "moves") {
        model::mean_new_tracks_per_frame(para_sample)           = 0.1f;
        model::mean_false_observations_per_frame(para_sample)   = 0.25f;
        model::frame_to_frame_survival_probability(para_sample) = 0.9f;
        model::generate_observation_probability(para_sample)    = 0.9f;
        mh_only(partition_ptr, para_sample, gen_pars);
    }
    else
        std::cout << "--opts is required"
            << std::endl
            << "Possible values: \"para\", \"part\", \"gibbs\", \"moves\""
            << std::endl
            ;
    ok1(true);

    diag("random seed = 0x%x", gen_pars.seed);
    return exit_status();
}
