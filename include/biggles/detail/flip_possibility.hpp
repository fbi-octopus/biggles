#ifndef CHANGE_POSSIBILITY_HPP__
#define CHANGE_POSSIBILITY_HPP__

#include <deque>

#include "../track.hpp"
#include "../sampling/simple.hpp"
#include "../detail/physics.hpp"

namespace biggles { namespace detail {

enum track_side {NO_SIDE, FRONT_END, BACK_END};

class flip_possibility_test {
    const shared_const_track_ptr& bro;
    const shared_const_track_ptr& sis;
    bool front_ok; ///< \brief is a flip possible at the back end?
    bool back_ok; ///< \brief is a flip possible at the back end?
    bool passed_;
    void get_possible_times();
public:
    flip_possibility_test(const shared_const_track_ptr& brother, const shared_const_track_ptr& sister)
        : bro(brother), sis(sister)
    {
        get_possible_times();
        passed_ = front_ok or back_ok;
    }

    bool passed() const { return passed_; }

    track_side sample_flip_side() const {
        if (front_ok and back_ok) {
            if (sampling::uniform_int(0, 2) == 1)
                return FRONT_END;
            return BACK_END;
        }
        if (front_ok)
            return FRONT_END;
        if (back_ok)
            return BACK_END;
        return NO_SIDE;
    }
    int num_flips() const {
        return int(front_ok) + int(back_ok);
    }
};

typedef boost::shared_ptr<const flip_possibility_test> flip_possibility_ptr;

}} //biggles::detail

#endif // CHANGE_POSSIBILITY_HPP__
