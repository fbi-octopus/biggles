#include "detail/transfer_possibility.hpp"

namespace biggles { namespace detail {

bool transfer_possibility_test::no_transfer_time() {
    time_stamp time_obs_end( t(donor_->observations().back())+1 ); // the time stamp *after* the last observation
    int d_before_a = int(time_obs_end) - int(acceptor_->first_time_stamp());
    int d_after_a = int(int(acceptor_->last_time_stamp() - t(donor_->observations().front())));
    if (d_before_a < 1 or d_after_a < 1)
        return true;
    t_beg = std::max(t(donor_->observations().front()), acceptor_->first_time_stamp());
    t_end = std::min(time_obs_end, acceptor_->last_time_stamp());
    obs_beg = donor_->observations().lower_bound_for_time_stamp(t_beg);
    obs_end = donor_->observations().upper_bound_for_time_stamp(t_end-1);
    return false;
}

void transfer_possibility_test::get_possible_times() {
    const observation_collection& acceptor_obs = acceptor_->observations();
    observation_collection::const_iterator pred = acceptor_obs.begin();
    observation_collection::const_iterator scout = acceptor_obs.begin();
    scout++;
    for (observation_collection::const_iterator donor_it = obs_beg; donor_it != obs_end; ++donor_it) {
        time_stamp time = t(*donor_it);
        while(scout != acceptor_obs.end() and t(*scout) < time) {
            pred++;
            scout++;
        }
        if (t(*pred) >= time) { // the donor obs is earlier than all acceptor obs
            if (observations_see_each_other(*pred, *donor_it, detail::speed_of_light_))
                possible_times_.push_back(time);
            continue;
        }
        if (scout == acceptor_obs.end()) { // the donor obs is later than all acceptor obs
            if (observations_see_each_other(*pred, *donor_it, detail::speed_of_light_))
                possible_times_.push_back(time);
            continue;
        }
        if (t(*scout) == time) // there is a acceptor obs at the same time
            continue;
        if (observations_see_each_other(*pred, *donor_it, detail::speed_of_light_) and
            observations_see_each_other(*donor_it, *scout, detail::speed_of_light_))
                possible_times_.push_back(time);
    }
}

}} //biggles::detail
