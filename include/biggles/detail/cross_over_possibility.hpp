#ifndef CROSS_OVER_POSSIBILITY_HPP__
#define CROSS_OVER_POSSIBILITY_HPP__

#include <deque>

#include "../track.hpp"
#include "../sampling/simple.hpp"
#include "../detail/physics.hpp"

namespace biggles { namespace detail {

/// \brief A test if a cross-over move for 2 given tracks is possible.
class cross_over_possibility_test {
    const shared_const_track_ptr& bro;
    const shared_const_track_ptr& sis;
    time_stamp t_beg;
    time_stamp t_end;
    bool passed_;
    std::deque<time_stamp> possible_times;
    bool no_cross_over_time();
    void get_possible_times();
public:
    cross_over_possibility_test(const shared_const_track_ptr& brother, const shared_const_track_ptr& sister)
        : bro(brother), sis(sister)
    {
        if (no_cross_over_time())
            passed_ = false;
        else
            get_possible_times();
        passed_ = possible_times.size() > 0;
    }
    bool passed() const {
        return passed_;
    }

    time_stamp sample_cross_over_time() const {
        return possible_times.at(sampling::uniform_int(0, possible_times.size()));
    }
};

typedef boost::shared_ptr<const cross_over_possibility_test> cross_over_possibility_ptr;

}} //biggles::detail

#endif // CROSS_OVER_POSSIBILITY_HPP__
