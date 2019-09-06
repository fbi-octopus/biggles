#include <boost/assert.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/format.hpp>

#include "mh_moves/utility.hpp"
//#include "sampling/simple.hpp"
#include "tools/debug.hpp"
#include "track.hpp"

namespace biggles
{

namespace mh_moves
{


namespace cross_over_fun {

    /// \brief this just does a pure time check; no distances are checked
    ///
    /// a cross over cuts each track in 2 sub-tracks and crosses them over
    /// (make sure that the crossed over tracks have enough observations)
    bool no_cross_over_time(shared_const_track_ptr_pair& tracks_to_cross, cross_over_data_t& data) {
        const track& bro = *(tracks_to_cross.first);
        const track& sis = *(tracks_to_cross.second);
        time_stamp sis_first_obs = t(sis.observations().front());
        time_stamp bro_first_obs = t(bro.observations().front());
        time_stamp sis_last_obs = t(sis.observations().back());
        time_stamp bro_last_obs = t(bro.observations().back());
        if (bro_first_obs >= sis_last_obs or sis_first_obs >= bro_last_obs)
                return true;
        data.t_beg = std::max(sis_first_obs, bro_first_obs)+1;
        data.t_end = std::min(sis_last_obs, bro_last_obs);
        return false;
    }

    void get_possible_times(shared_const_track_ptr_pair& tracks_to_cross, std::deque<time_stamp>& possible_times,
        const cross_over_data_t& data)
    {
        shared_const_track_ptr& bro = tracks_to_cross.first;
        shared_const_track_ptr& sis = tracks_to_cross.second;
        // the observations just before the begin; no iterator will be *end*
        observation_collection::const_iterator bro_pred_obs = bro->observations().begin();
        observation_collection::const_iterator bro_scout = bro_pred_obs;
        observation_collection::const_iterator sis_pred_obs = sis->observations().begin();
        observation_collection::const_iterator sis_scout = sis_pred_obs;
        sis_scout++;
        bro_scout++;
        // the observations just after and incl. the end; no iterator will be *end*
        for (time_stamp time = data.t_beg; time <= data.t_end; ++time) {
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
            if (observations_could_be_neighbours(*sis_pred_obs, *bro_scout) and
                observations_could_be_neighbours(*bro_pred_obs, *sis_scout)) {
                    possible_times.push_back(time);
            }
        }
    }

}

track_collection::const_iterator sample_track_with_minimum_size(
    const track_collection& tc,
    size_t minimum_size,
    float& output_log_prob)
{
    size_t n_in_set(0);
    track_collection::const_iterator current_sample(tc.end());
    track_collection::const_iterator first(tc.begin());
    track_collection::const_iterator last(tc.end());

    for(; first != last; ++first)
    {
        if((*first)->size() < minimum_size)
            continue;
        if(0 == sampling::uniform_int(0, n_in_set+1))
        {
            current_sample = first;
        }
        ++n_in_set;
    }

    output_log_prob = (n_in_set == 0)
        ? (-std::numeric_limits<float>::max())
        : (-logf(static_cast<float>(n_in_set)));

    return current_sample;
}

shared_const_track_ptr_pair split_track_at_time_stamp(
    const shared_const_track_ptr& track_to_split,
    time_stamp split_ts, float& log_splitting_prob)
{
    // check that the split point will not leave a zero-duration track
    BOOST_ASSERT(split_ts > track_to_split->first_time_stamp());
    BOOST_ASSERT(split_ts < track_to_split->last_time_stamp());
    BOOST_ASSERT(!track_to_split->empty());

    // find the splitting point: the first observation with a time stamp >= split_ts
    observation_collection::const_iterator first(track_to_split->observations().begin());
    observation_collection::const_iterator last(track_to_split->observations().end());
    observation_collection::const_iterator split(track_to_split->observations().lower_bound_for_time_stamp(split_ts));

    BOOST_ASSERT(t(*first) >= track_to_split->first_time_stamp());
    BOOST_ASSERT(t(*first) < track_to_split->last_time_stamp());

    // HACK: track a's final time stamp is set after the final observation within it
    time_stamp track_a_last(track_to_split->first_time_stamp() + 1);
    if(split != first)
    {
        observation_collection::const_iterator prior_to_split(split);
        --prior_to_split;
        track_a_last = t(*prior_to_split) + 1;
    }

    /*
    time_stamp track_b_first(track_to_split->last_time_stamp() - 1);
    if(split != last)
    {
        BOOST_ASSERT(t(*split) >= track_to_split->first_time_stamp());
        BOOST_ASSERT(t(*split) < track_to_split->last_time_stamp());
        track_b_first = t(*split);
    }
    */

    // get rid of some noobs at the *end* of the *first* track
    log_splitting_prob = 0.f;
    if (track_a_last != split_ts) {
        log_splitting_prob += -logf(split_ts + 1 - track_a_last);
        track_a_last = sampling::uniform_int(track_a_last, split_ts + 1);
    }
    // create the new tracks
    shared_const_track_ptr track_a(new track(
            track_to_split->first_time_stamp(), track_a_last,
            first, split, track_to_split->dynamic_drag()));
    shared_const_track_ptr track_b(new track(
            split_ts, track_to_split->last_time_stamp(),
            split, last, track_to_split->dynamic_drag()));

    return std::make_pair(track_a, track_b);
}

std::deque<float> get_split_weights(observation_collection::const_iterator fst,
    const observation_collection::const_iterator& lst)
{
    observation_collection::const_iterator lead = fst;
    std::deque<float> result;
    const float inc = 0.1f;
    while (lead != lst) {
        ++lead;
        time_stamp dt = t(*lead) - t(*fst);
        float dobs = dist(*fst, *lead) + inc;
        for (time_stamp t = 1; t <= dt; ++t) result.push_back(dobs);
        ++fst;
    }
    return result;
}

std::deque<float> get_split_weights_unif(observation_collection::const_iterator fst,
    observation_collection::const_iterator& lst)
{
    return std::deque<float>(size_t(t(*lst)-t(*fst)), 1.f);
}

float get_time_stamp_to_split_track_prob(const shared_const_track_ptr& track_to_split,
        const time_stamp& ts_to_split_at)
{
    // FIXME: check that time stamp is in range
    BOOST_ASSERT(ts_to_split_at >= track_to_split->first_time_stamp() and
                 ts_to_split_at  < track_to_split->last_time_stamp());
    BOOST_ASSERT(track_to_split->size() > 3);
    // find the second and second to last observations
    observation_collection::const_iterator second_obs(track_to_split->observations().begin());
    std::advance(second_obs, 1);

    observation_collection::const_iterator second_to_last_obs(track_to_split->observations().end());
    std::advance(second_to_last_obs, -2);

    // how much slack do we have for the splitting point
    time_stamp slack(t(*second_to_last_obs) - t(*second_obs));
    BOOST_ASSERT(slack > 0);

    std::deque<float> weights = get_split_weights(second_obs, second_to_last_obs);
    //OK(slack);
    //OK(weights.size());
    BOOST_ASSERT(slack == time_stamp(weights.size()));
    const float total = std::accumulate(weights.begin(), weights.end(), 0.f);
    return logf(weights[ts_to_split_at - t(*second_obs) - 1]/total);
}

bool sample_time_stamp_to_split_track(
        const shared_const_track_ptr& track_to_split,
        time_stamp& ts_to_split_at,
        float& log_prob_ts)
{
    // check it is big enough to split
    if(track_to_split->size() < 4)
        return false;

    // find the second and second to last observations
    observation_collection::const_iterator second_obs(track_to_split->observations().begin());
    std::advance(second_obs, 1);

    observation_collection::const_iterator second_to_last_obs(track_to_split->observations().end());
    std::advance(second_to_last_obs, -2);

    // how much slack do we have for the splitting point
    time_stamp slack(t(*second_to_last_obs) - t(*second_obs));
    BOOST_ASSERT(slack > 0);

    std::deque<float> weights = get_split_weights(second_obs, second_to_last_obs);
    //OK(slack);
    //OK(weights.size());
    BOOST_ASSERT(slack == time_stamp(weights.size()));

    // choose a splitting point
    ts_to_split_at = t(*second_obs);
    /*
    time_stamp time_incr = sampling::uniform_int(0, slack);
    log_prob_ts = -logf(static_cast<float>(slack));
    */
    time_stamp time_incr = sampling::weighted_choice_index(weights.begin(), weights.end(), log_prob_ts);
    ts_to_split_at += time_incr + 1;
    return true;
}

bool observations_could_be_neighbours(const observation& o1, const observation& o2)
{
    //if (std::abs(t(o2) - t(o1)) >= detail::last_delta_t_)
    //    return false;
    return observations_see_each_other(o1, o2, detail::speed_of_light_);
}

bool tracks_could_be_merged(const track& t1, const track& t2) {
    if(t1.empty() || t2.empty())
        return false;

    if (t2.first_time_stamp() < t1.last_time_stamp() and t1.first_time_stamp() < t2.last_time_stamp())
        return false;

    if(t(t1.observations().back()) < t(t2.observations().front()))
        return observations_could_be_neighbours(t1.observations().back(), t2.observations().front());

    if(t(t2.observations().back()) < t(t1.observations().front()))
        return observations_could_be_neighbours(t2.observations().back(), t1.observations().front());

    return false;
}

float tracks_could_be_merged2(const track& t1, const track& t2) {
    if(t1.empty() || t2.empty())
        return -1.f;

    if (t2.first_time_stamp() < t1.last_time_stamp())
        return -1.f;

    return observations_see_each_other2(t1.observations().back(), t2.observations().front(), detail::speed_of_light_);
}

shared_const_track_ptr_pair choose_two_tracks(const track_collection& tracks) {
    size_t n_donor = sampling::uniform_int(0, tracks.size());
    size_t n_acceptor = sampling::uniform_int(0, tracks.size());
    while (n_acceptor == n_donor) n_acceptor = sampling::uniform_int(0, tracks.size());

    track_collection::const_iterator it_don(tracks.begin());
    std::advance(it_don, n_donor);
    track_collection::const_iterator it_acc = tracks.begin();
    std::advance(it_acc, n_acceptor);
    return shared_const_track_ptr_pair(*it_don, *it_acc);
}

bool track_could_transfer_obs_to_other(const track& t1,  const track& t2) {
    throw std::logic_error("track_could_transfer_obs_to_other: not implemented.");
}

bool tracks_could_cross(const track& t1,  const track& t2) {
    throw std::logic_error("tracks_could_cross: not implemented.");
}

void print_track(const track& tr) {
    for (time_stamp ts(tr.first_time_stamp()); ts < tr.last_time_stamp(); ++ts) {
        observation_collection::const_range obs_at_t(tr.observations().at_time_stamp(ts));
        if (obs_at_t.first != obs_at_t.second) {
            while (obs_at_t.first != obs_at_t.second) {
                observation obs = *obs_at_t.first++;
                std::cout << boost::format("[%f, %f, %d]") % x(obs) % y(obs) % t(obs) << std::endl;
            }
        }
        else {
            std::cout << boost::format("[-, -, %d]") % ts << std::endl;
        }
    }
}



}

}
