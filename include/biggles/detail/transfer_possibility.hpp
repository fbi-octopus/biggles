#ifndef TRANSFER_POSSIBILITY_HPP__
#define TRANSFER_POSSIBILITY_HPP__

#include <deque>

#include "../track.hpp"
#include "../sampling/simple.hpp"
#include "../detail/physics.hpp"

namespace biggles { namespace detail {

class transfer_possibility_test {
    const shared_const_track_ptr& donor_;
    const shared_const_track_ptr& acceptor_;
    time_stamp t_beg, t_end;
    bool passed_;
    observation_collection::const_iterator obs_beg, obs_end;
    std::deque<time_stamp> possible_times_;

    bool no_transfer_time();
    void get_possible_times();

public:

    transfer_possibility_test(const shared_const_track_ptr& donor, const shared_const_track_ptr& acceptor)
        : donor_(donor), acceptor_(acceptor)
    {
        if  (donor->size() < 3)
            passed_ = false;
        else if (no_transfer_time())
            passed_ = false;
        else {
            get_possible_times();
            passed_ = possible_times_.size() > 0;
        }
    }

    bool passed() const {
        return passed_;
    }

    time_stamp sample_transfer_time() const {
        return possible_times_.at(sampling::uniform_int(0, possible_times_.size()));
    }

    float log_prob_transfer_time() const {
        return -logf(possible_times_.size());
    }
};

typedef boost::shared_ptr<const transfer_possibility_test> transfer_possibility_ptr;

}} //biggles::detail

#endif //TRANSFER_POSSIBILITY_HPP__
