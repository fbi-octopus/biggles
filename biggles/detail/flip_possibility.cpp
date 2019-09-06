#include "detail/flip_possibility.hpp"

namespace biggles { namespace detail {

bool check_front(shared_const_track_ptr bro, shared_const_track_ptr sis) {
    time_stamp t_bro_front = bro->first_time_stamp();
    time_stamp t_sis_front = sis->first_time_stamp();
    if (t_sis_front == t_bro_front)
        return false;
    shared_const_track_ptr acceptor_cand(bro);
    shared_const_track_ptr donor_cand(sis);
    if (t_sis_front > t_bro_front) {
        acceptor_cand = sis;
        donor_cand = bro;
    }
    observation accept_first = acceptor_cand->observations().front();
    time_stamp ts = acceptor_cand->first_time_stamp();
    observation_collection::const_iterator lower_bound = donor_cand->observations().lower_bound_for_time_stamp(ts);
    //size_t num_obs = std::distance(donor_cand->observations().begin(), lower_bound);
    if (lower_bound == donor_cand->observations().begin() or lower_bound == donor_cand->observations().end()) // there are no observations in the donated branch
        return false; // no empty exchange
    if (std::distance(lower_bound, donor_cand->observations().end()) < 2) // not enough left
        return false;
    lower_bound--;
    return observations_see_each_other(*lower_bound, accept_first, detail::speed_of_light_);
}

bool check_back(shared_const_track_ptr bro, shared_const_track_ptr sis) {
    time_stamp t_bro_back = bro->last_time_stamp();
    time_stamp t_sis_back = sis->last_time_stamp();
    if (t_sis_back == t_bro_back)
        return false;
    shared_const_track_ptr acceptor_cand(bro);
    shared_const_track_ptr donor_cand(sis);
    if (t_sis_back < t_bro_back) {
        acceptor_cand = sis;
        donor_cand = bro;
    }
    observation accept_back = acceptor_cand->observations().back();
    time_stamp ts = acceptor_cand->last_time_stamp();
    observation_collection::const_iterator lower_bound = donor_cand->observations().lower_bound_for_time_stamp(ts);
    if (lower_bound == donor_cand->observations().end() or lower_bound == donor_cand->observations().begin())
        return false; // no empty exchange
    if (std::distance(donor_cand->observations().begin(), lower_bound)<2)
        return false;
    return observations_see_each_other(accept_back, *lower_bound, detail::speed_of_light_);
}

void flip_possibility_test::get_possible_times() {
    front_ok = check_front(bro, sis);
    back_ok = check_back(bro, sis);
}

}}  // biggles::detail
