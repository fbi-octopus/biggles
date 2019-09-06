#ifndef BIGGLES_SERVER_STATE_HPP__
#define BIGGLES_SERVER_STATE_HPP__

#include "../mh_moves/mh_moves.hpp"
#include "../model.hpp"
#include "../partition.hpp"
#include "../tracker.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/thread/thread.hpp>

namespace biggles
{

/// @brief Support for writing a Bigglers tracking server.
///
/// This namespace contains a tracking engine suitable for embedding into a server. The engine itself runs the tracking
/// in a separate thread and allows querying the current state of the tracker.
namespace server
{

/// @brief A snapshot of the state of a Biggles tracking process.
///
/// @sa server::engine
struct tracking_state
{
    /// @brief The number of samples which have been drawn.
    size_t sample_count;

    /// @brief The model parameters corresponding to the last drawn sample.
    biggles::model::parameters current_model_parameters;

    /// @brief The partition corresponding to the last drawn sample.
    biggles::partition_ptr_t current_partition;
    /// \brief This is for the python interface
    biggles::partition current_partition_fun() const {return *current_partition;}

    /// @brief The proposal move type by which the current partition was proposed.
    biggles::mh_moves::move_type current_move_type;

    /// @brief The computed log PDF of the last drawn sample.
    /// is either last_proposal_density (if proposal was accepted)
    /// or last_sample_density (if proposal was rejected)
    float current_log_pdf;

    /// @brief the target density of the last proposal
    float last_proposal_density;

    /// @brief the target density of the last sample in comparison
    float last_sample_density;

    /// @brief The model parameters associated with the sample with lowest PDF value.
    biggles::model::parameters best_model_parameters;

    /// @brief The partition associated with the sample with lowest PDF value.
    biggles::partition_ptr_t best_partition;
    /// \brief This is for the python interface
    biggles::partition best_partition_fun() const {return *best_partition;}

    /// @brief The lowest PDF value of any sample drawn so far.
    float best_log_pdf;

    /// @brief The acceptance rate for all samples of the Metropolis Hastings part of the tracker.
    float acceptance_rate;

    /// @brief The acceptance rate for the last samples of the Metropolis Hastings part of the tracker.
    float recent_acceptance_rate;

    bool accepted; // was the last move accepted

    /// \brief the last proposal partition; 1st component of last_proposal
    biggles::partition_ptr_t last_proposed_partition;
    /// \brief This is for the python interface
    biggles::partition last_proposed_partition_fun() const {return *last_proposed_partition;}

    /// \brief last PDR
    float last_pdr;

    /// \brief the last proposal move; 2st component of last_proposal
    biggles::mh_moves::move_type last_proposed_move;

    /// @brief A histogram of accepted move types.
    std::deque<uint64_t> move_histogram;

    /// @brief A record for the move statisicis (a/r/i).
    std::deque<uint64_t> moves_accepted;
    std::deque<uint64_t> moves_rejected;
    std::deque<uint64_t> moves_identity;

    std::deque<uint64_t> internals;

    //tracking_state();
};

/// @brief A threaded engine for a Biggles server.
///
/// This class encapsulates the logic required to implement Biggles as a server keeping the tracking process within a
/// separate thread. It can be used to expose the Biggles tracker via almost any RPC mechanism with an appropriate
/// adaptor.
class engine
{
public:
    engine();

    ~engine();

    /// @brief Copy constructor
    ///
    /// @param e
    engine(const engine& e);

    /// @brief Set the tracker's initial conditions.
    ///
    /// @note This will implicitly cause any current tracking to be stopped via stop().
    ///
    /// @param params
    /// @param partition
    void set_initial_conditions(const biggles::model::parameters& params, const biggles::partition_ptr_t& partition_ptr);
    void set_initial_conditions_part(const biggles::model::parameters& params,const biggles::partition& part);

    /// @brief Obtain the initial model parameters used by the tracker.
    const biggles::model::parameters& initial_model_parameters() const { return initial_parameters_; }

    /// @brief Obtain the initial partition used by the tracker.
    const biggles::partition_ptr_t& initial_partition() const { return initial_partition_; }

    /// @brief Return a shared pointer to the current tracking state.
    ///
    /// The contents are only useful if the tracking is currently running. If the tracking is not running, this is the
    /// tracking state at the time the last tracker stopped or an empty initial value if there was no prior tracking.
    boost::shared_ptr<server::tracking_state> tracking_state_ptr() const;

    /// @brief Start the tracker thread running.
    ///
    /// @note Calls to start() when the thread is already running are NOPs.
    void start();

    /// @brief Stop the tracker thread running.
    ///
    /// This member function will block until the thread has stopped.
    ///
    /// @note Calls to stop() when the thread is not running are NOPs.
    void stop();

    void record_state(const biggles::tracker& track);

protected:
    /// @brief Set to \c true iff the tracking thread should stop.
    bool should_stop_;

    /// @brief Initial model parameters for the tracker.
    biggles::model::parameters initial_parameters_;

    /// @brief An initial partition for the tracker.
    biggles::partition_ptr_t initial_partition_;

    /// @brief The current tracking thread. Can be NULL if no tracking is currently running.
    boost::shared_ptr<boost::thread> tracking_thread_ptr_;

    /// @brief The current state of the tracks.
    boost::shared_ptr<tracking_state> tracking_state_ptr_;

    /// @brief A mutex protecting tracking_state_ptr_
    ///
    /// This is needed because also shared_ptr allows simultaneous reads *or* simultaneous writes, a simultaneous write
    /// and read is undefined.
    mutable boost::mutex tracking_state_mutex_;

    /// @brief The tracking loop implementation.
    ///
    /// This does not return until \c should_stop_ is \c true.
    void track_loop();
};

}

}

#endif // BIGGLES_SERVER_STATE_HPP__
