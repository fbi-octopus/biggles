#include "model.hpp"
#include "partition.hpp"
#include "server/engine.hpp"
#include "tracker.hpp"
#include "tools/debug.hpp"

#include <boost/assert.hpp>
#include <boost/thread/locks.hpp>
#include <functional>
#include <stdexcept>
#include <iostream>

namespace biggles
{

namespace server
{

engine::engine()
    : should_stop_(true)
    , tracking_state_ptr_(new tracking_state())
{
    //std::cout << "random seed = " << std::hex << std::uppercase << detail::get_random_seed() << std::endl;
    tracking_state_ptr_->current_partition = partition_ptr_t(new partition);
    tracking_state_ptr_->best_partition = tracking_state_ptr_->current_partition;
    tracking_state_ptr_->last_proposed_partition = tracking_state_ptr_->current_partition;
}

engine::~engine()
{
    stop();
}

engine::engine(const engine& e)
    : should_stop_(true)
    , initial_parameters_(e.initial_parameters_)
    , initial_partition_(e.initial_partition_)
    , tracking_state_ptr_(e.tracking_state_ptr_)
{ }

void engine::set_initial_conditions(
    const biggles::model::parameters& params,
    const biggles::partition_ptr_t& partition_ptr)
{
    // stop any tracking which is currently in progress
    stop();

    // update the initial conditions
    initial_parameters_ = params;
    initial_partition_ = partition_ptr;
}

void engine::set_initial_conditions_part(const biggles::model::parameters& params,const biggles::partition& part) {
    // stop any tracking which is currently in progress
    stop();

    // update the initial conditions
    initial_parameters_ = params;
    initial_partition_ = partition_ptr_t(new partition(part));
}

void engine::start()
{
    // starting multiple times is a NOP
    if(tracking_thread_ptr_)
        return;

    // reset tracking state
    tracking_state_ptr_ = boost::shared_ptr<tracking_state>(new tracking_state());

    // kick off a new thread
    should_stop_ = false;
    tracking_thread_ptr_ = boost::shared_ptr<boost::thread>(
        new boost::thread(std::mem_fun(&engine::track_loop), this));
}

void engine::stop()
{
    // stopping multiple times is a NOP
    if(!tracking_thread_ptr_)
        return;

    // stop thread
    should_stop_ = true;
    tracking_thread_ptr_->join();
    tracking_thread_ptr_.reset();
}

boost::shared_ptr<tracking_state> engine::tracking_state_ptr() const
{
    boost::lock_guard<boost::mutex> guard(tracking_state_mutex_);
    return tracking_state_ptr_;
}

void engine::record_state(const biggles::tracker& track) {
        // record tracker state
        boost::shared_ptr<tracking_state> state(new tracking_state());

        state->sample_count = track.n_iterations();

        state->current_model_parameters = track.last_parameters();
        state->current_partition = track.last_partition();
        state->current_log_pdf = track.last_log_pdf();
        state->current_move_type = track.last_move_type();
        state->last_proposal_density = track.last_proposal_density();
        state->last_sample_density = track.mh_state().last_sample_log_density;

        state->best_model_parameters = track.best_parameters();
        state->best_partition = track.best_partition();
        state->best_log_pdf = track.best_log_pdf();

        state->acceptance_rate = track.acceptance_rate();
        state->recent_acceptance_rate = track.recent_acceptance_rate();
        state->accepted = track.mh_state().accepted;

        state->last_proposed_partition = track.last_proposal_partition();

        state->last_proposed_move = track.last_proposal_move();


        state->move_histogram = track.move_histogram();
        state->moves_accepted = track.move_statistics().accepted;
        state->moves_rejected = track.move_statistics().rejected;
        state->moves_identity = track.move_statistics().identity;

        state->internals = track.mh_internals().get();
        /* */
        // update state record
        {
            boost::lock_guard<boost::mutex> guard(tracking_state_mutex_);
            tracking_state_ptr_ = state;
        }
}

void engine::track_loop()
{
    using namespace biggles;

    // create the tracker
    biggles::tracker track(initial_parameters_, initial_partition_);

    // check partition
    if (track.last_partition()->tracks().size() == 0 and track.last_partition()->clutter().size() < 2)
        track.set_partition_is_valid(false);
    else
        track.set_partition_is_valid(true);

    // keep tracking until told to stop
    while(!should_stop_)
    {
        track.advance();
        record_state(track);
    }
}

//tracking_state::tracking_state() {}
    //: current_partition(new partition), best_partition(current_partition), last_proposed_partition(current_partition) {}

} // biggles::server

} // biggles
