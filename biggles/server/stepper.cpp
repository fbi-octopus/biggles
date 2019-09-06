#include "model.hpp"
#include "partition.hpp"
#include "server/stepper.hpp"
#include "tools/debug.hpp"

#include <boost/assert.hpp>
#include <functional>
#include <stdexcept>

namespace biggles { namespace server {

stepper::stepper(const biggles::model::parameters& params, const biggles::partition_ptr_t& partition_ptr)
    : initial_parameters_(params)
    , initial_partition_(partition_ptr)
    , tracking_state_ptr_(new tracking_state())
    , track_(params, partition_ptr)
{        record_state();  }

stepper::stepper(const biggles::model::parameters& params, const biggles::partition& part)
    : initial_parameters_(params)
    , initial_partition_(partition_ptr_t(new partition(part)))
    , tracking_state_ptr_(new tracking_state())
    , track_(params, initial_partition_)
{         record_state(); }

stepper::~stepper () {}

boost::shared_ptr<tracking_state> stepper::tracking_state_ptr() const {
    return tracking_state_ptr_;
}

void stepper::record_state() {
        //boost::shared_ptr<tracking_state> state(new tracking_state());

        tracking_state_ptr_->sample_count = track_.n_iterations();
        tracking_state_ptr_->current_model_parameters = track_.last_parameters();
        tracking_state_ptr_->current_partition = track_.last_partition();
        tracking_state_ptr_->current_log_pdf = track_.last_log_pdf();
        tracking_state_ptr_->current_move_type = track_.last_move_type();
        tracking_state_ptr_->last_proposal_density = track_.last_proposal_density();
        tracking_state_ptr_->last_sample_density = track_.mh_state().last_sample_log_density;
        tracking_state_ptr_->best_model_parameters = track_.best_parameters();
        tracking_state_ptr_->best_partition = track_.best_partition();
        tracking_state_ptr_->best_log_pdf = track_.best_log_pdf();
        tracking_state_ptr_->acceptance_rate = track_.acceptance_rate();
        tracking_state_ptr_->recent_acceptance_rate = track_.recent_acceptance_rate();
        tracking_state_ptr_->accepted = track_.mh_state().accepted;
        tracking_state_ptr_->last_proposed_partition = track_.last_proposal_partition();
        tracking_state_ptr_->last_proposed_move = track_.last_proposal_move();
        tracking_state_ptr_->last_pdr = track_.last_pdr();
        tracking_state_ptr_->move_histogram = track_.move_histogram();
        tracking_state_ptr_->moves_accepted = track_.move_statistics().accepted;
        tracking_state_ptr_->moves_rejected = track_.move_statistics().rejected;
        tracking_state_ptr_->moves_identity = track_.move_statistics().identity;
        tracking_state_ptr_->internals = track_.mh_internals().get();

}

void stepper::step(const size_t num_steps) {
    /* this must be done somewhere else
    // check partition
    if (track.last_partition().tracks().size() == 0 and track.last_partition().clutter().size() < 2)
        track.set_partition_is_valid(false);
    else
        track.set_partition_is_valid(true);
     */

    // keep tracking until told to stop
    for (size_t i = 0; i < num_steps; ++i ) {
        track_.advance();
        record_state();
    }

}

void stepper::fix_observation_error(const float R00, const float R01, const float R11) {
    track_.fix_observation_error(R00, R01, R11);
}

void stepper::sample_observation_error(const size_t lag) {
    track_.sample_observation_error(lag);
}




}} // biggles::server
