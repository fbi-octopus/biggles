#include <limits>

#include "tracker.hpp"
#include "tools/debug.hpp"

namespace biggles
{

tracker::tracker() : n_iterations_(0) { std::cerr << "WARNING this is unexpected." << std::endl; }

tracker::tracker(const model::parameters& initial_parameters,
                 const partition_ptr_t& initial_partition_ptr)
    : sampling::new_gibbs_sampler(initial_partition_ptr, initial_parameters)
    , n_iterations_(0)
    , best_partition_(initial_partition_ptr)
    , best_parameters_(initial_parameters)
    , best_log_pdf_(-std::numeric_limits<float>::max())
    , move_hist_(mh_moves::MOVE_COUNT, 0)
{ }


tracker::tracker(const tracker& t)
    : n_iterations_(t.n_iterations_)
    , best_partition_(t.best_partition_)
    , best_parameters_(t.best_parameters_)
    , best_log_pdf_(t.best_log_pdf_)
    , move_hist_(mh_moves::MOVE_COUNT, 0)
{ std::cerr << "WARNING Thou shall not copy." << std::endl; }


const tracker& tracker::operator = (const tracker& t)
{
    n_iterations_ = t.n_iterations_;
    best_partition_ = t.best_partition_;
    best_parameters_ = t.best_parameters_;
    best_log_pdf_ = t.best_log_pdf_;
    move_hist_ = t.move_hist_;
    return *this;
}

void tracker::advance()
{
    draw();

    record_move_type_(last_move_type());
    record_move_stats_(last_proposed_sample());

    if(last_log_pdf() > best_log_pdf())
    {
        best_log_pdf_ = last_log_pdf();
        best_partition_ = last_partition();
        best_parameters_ = last_parameters();
    }

    ++n_iterations_;
}

void tracker::record_move_type_(mh_moves::move_type mt) {
    BOOST_ASSERT(mt >= 0);
    BOOST_ASSERT(mt < mh_moves::MOVE_COUNT);
    if (last_proposal_accepted())
        ++move_hist_[mt];
    else
        ++move_hist_[mh_moves::NONE];
}

void tracker::record_move_stats_(const partition_sampler_sample& sample) {
    int index = static_cast<int>(sample.proposed_move - mh_moves::BIRTH);
    if (not last_proposal_accepted()) {
        move_stats_.rejected[index]++;
    } else if (sample.executed_move == mh_moves::IDENTITY) {
        move_stats_.identity[index]++;
    } else {
        move_stats_.accepted[index]++;
    }
}

}
