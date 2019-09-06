#include "detail/cross_over_possibility.hpp"

namespace biggles { namespace detail {

bool cross_over_possibility_test::no_cross_over_time() {
    time_stamp sis_first_obs = t(sis->observations().front());
    time_stamp bro_first_obs = t(bro->observations().front());
    time_stamp sis_last_obs = t(sis->observations().back());
    time_stamp bro_last_obs = t(bro->observations().back());
    if (bro_first_obs >= sis_last_obs or sis_first_obs >= bro_last_obs)
            return true;
    t_beg = std::max(sis_first_obs, bro_first_obs)+1;
    t_end = std::min(sis_last_obs, bro_last_obs);
    return false;
}

void cross_over_possibility_test::get_possible_times() {
    // the observations just before the begin; no iterator will be *end*
    observation_collection::const_iterator bro_pred_obs = bro->observations().begin();
    observation_collection::const_iterator bro_scout = bro_pred_obs;
    observation_collection::const_iterator sis_pred_obs = sis->observations().begin();
    observation_collection::const_iterator sis_scout = sis_pred_obs;
    sis_scout++;
    bro_scout++;
    // the observations just after and incl. the end; no iterator will be *end*
    for (time_stamp time = t_beg; time <= t_end; ++time) {
        while (t(*bro_scout) < time) {
            bro_pred_obs++;
            bro_scout++;
        }
        while (t(*sis_scout) < time) {
            sis_pred_obs++;
            sis_scout++;
        }
        BOOST_ASSERT(sis_scout != sis->observations().end());
        BOOST_ASSERT(bro_scout != bro->observations().end());
        // now the xxx_pred_obs points to an obs that is earlier than *time*
        // and the xxx_scout points to an obs that at the same time or later than *time*
        if (observations_see_each_other(*sis_pred_obs, *bro_scout, detail::speed_of_light_) and
            observations_see_each_other(*bro_pred_obs, *sis_scout, detail::speed_of_light_)) {
                possible_times.push_back(time);
        }
    }
}

}}  // biggles::detail
