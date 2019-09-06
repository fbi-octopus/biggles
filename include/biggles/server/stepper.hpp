#ifndef BIGGLES_STEPPER_HPP__
#define BIGGLES_STEPPER_HPP__

#include "../mh_moves/mh_moves.hpp"
#include "../model.hpp"
#include "../partition.hpp"

#include "../tracker.hpp"
#include "engine.hpp"

#include <boost/shared_ptr.hpp>

namespace biggles { namespace server {
    /// \brief alternative to the Engine class
    ///
    /// allows to do the sampling of p(omega, theta | data) step by step
    /// or a given number of steps at a time
    class stepper {
    public:
        stepper() { throw std::logic_error("undefined behaviour"); }
        stepper(const biggles::model::parameters& params, const biggles::partition_ptr_t& partition_ptr_t);
        stepper(const biggles::model::parameters& params, const biggles::partition& part);

        void record_state();

        ~stepper();

        /* maybe later
        /// @brief Copy constructor
        ///
        /// @param s
        stepper(const stepper& s);
        */

        /// @brief Obtain the initial model parameters used by the tracker.
        const biggles::model::parameters& initial_model_parameters() const { return initial_parameters_; }

        /// @brief Obtain the initial partition used by the tracker.
        const biggles::partition_ptr_t& initial_partition() const { return initial_partition_; }

        /// @brief Return a shared pointer to the current tracking state.
        ///
        /// The contents are only useful if the tracking is currently running. If the tracking is not running, this is the
        /// tracking state at the time the last tracker stopped or an empty initial value if there was no prior tracking.
        boost::shared_ptr<server::tracking_state> tracking_state_ptr() const;

        /* maybe not. this could be done in the constructor
        /// @brief Start the tracker thread running.
        ///
        /// @note Calls to start() when the thread is already running are NOPs.
        void start();
        */

        /* maybe not. This could be done in the destructor
        /// @brief Stop the tracker thread running.
        ///
        /// This member function will block until the thread has stopped.
        ///
        /// @note Calls to stop() when the thread is not running are NOPs.
        void stop();
        */

        /// \brief fix the observation error
        void fix_observation_error(const float R00, const float R01, const float R11);

        /// \brief sample the observation error every *lag* samples
        void sample_observation_error(const size_t lag);

        /// @brief The tracking loop implementation.
        ///
        /// This does \c num_steps number of sampling steps.
        void step(const size_t num_steps = 1);
    protected:
        /// @brief Initial model parameters for the tracker.
        biggles::model::parameters initial_parameters_;

        /// @brief An initial partition for the tracker.
        biggles::partition_ptr_t initial_partition_;

        /// @brief The current state of the tracks.
        boost::shared_ptr<tracking_state> tracking_state_ptr_;

        biggles::tracker track_;

    };

}} // biggles::server


#endif // BIGGLES_STEPPER_HPP__
