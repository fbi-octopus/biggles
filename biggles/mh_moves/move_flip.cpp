#include <boost/assert.hpp>
#include <cmath>
#include <limits>
#include <utility>

#include "mh_moves/mh_moves.hpp"
#include "sampling/simple.hpp"
#include "mh_moves/utility.hpp"
#include "detail/flip_possibility.hpp"
#include "tools/debug.hpp"

namespace biggles { namespace mh_moves {

shared_const_track_ptr_pair flip_at_front(shared_const_track_ptr_pair& tracks_to_flip) {
    shared_const_track_ptr donor(tracks_to_flip.first);
    shared_const_track_ptr acceptor(tracks_to_flip.second);

    if (acceptor->first_time_stamp() < donor->first_time_stamp()) {
        donor = tracks_to_flip.second;
        acceptor = tracks_to_flip.first;
    }

    time_stamp flip_time = acceptor->first_time_stamp();

    observation_collection::const_iterator lower_bound = donor->observations().lower_bound_for_time_stamp(flip_time);
    /* this should not be necessary
    size_t num_obs = std::distance(donor->observations().begin(), lower_bound);
    if (num_obs == 0) // there are no observations in the donated branch
        return false; // no empty exchange
        */
    BOOST_ASSERT(std::distance(donor->observations().begin(), lower_bound)>0);

    shared_track_ptr donated(new track(flip_time, donor->last_time_stamp(),
                                    lower_bound, donor->observations().end() ));
    shared_track_ptr accepted(new track(donor->first_time_stamp(), acceptor->last_time_stamp(),
                                    donor->observations().begin(), lower_bound));

    observation_collection::const_iterator acceptor_first = acceptor->observations().begin();
    for ( ; acceptor_first != acceptor->observations().end(); ++acceptor_first)
        accepted->insert(*acceptor_first);

    return shared_const_track_ptr_pair(donated, accepted);
}

shared_const_track_ptr_pair flip_at_back(shared_const_track_ptr_pair& tracks_to_flip) {
    shared_const_track_ptr donor(tracks_to_flip.first);
    shared_const_track_ptr acceptor(tracks_to_flip.second);

    if (acceptor->last_time_stamp() > donor->last_time_stamp()) {
        donor = tracks_to_flip.second;
        acceptor = tracks_to_flip.first;
    }
    BOOST_ASSERT(acceptor->last_time_stamp() < donor->last_time_stamp());

    time_stamp flip_time = acceptor->last_time_stamp();

    observation_collection::const_iterator lower_bound = donor->observations().lower_bound_for_time_stamp(flip_time);
    /* this should not be necessary
    size_t num_obs = std::distance(donor->observations().begin(), lower_bound);
    if (num_obs == 0) // there are no observations in the donated branch
        return false; // no empty exchange
        */
    BOOST_ASSERT(std::distance(lower_bound, donor->observations().end())>0);

    shared_track_ptr donated(new track(donor->first_time_stamp(), flip_time,
                                    donor->observations().begin(), lower_bound ));
    shared_track_ptr accepted(new track(acceptor->first_time_stamp(), donor->last_time_stamp(),
                                    acceptor->observations().begin(), acceptor->observations().end()));

    for ( ; lower_bound != donor->observations().end(); lower_bound++)
        accepted->insert(*lower_bound);

    return shared_const_track_ptr_pair(donated, accepted);
}

bool flip(const partition_ptr_t& start_partition_ptr, partition_ptr_t& end_partition_ptr,
    float& proposal_mass_ratio)
{
    const track_collection& tracks(start_partition_ptr->tracks());

    //proposal_mass_ratio = 0.f;

    capability_recorder_ptr cap_rec_ptr = start_partition_ptr->get_capability_recorder();
    BOOST_ASSERT(cap_rec_ptr.get() not_eq 0); // there actually is something assigned

    const flippables_t& flippables = cap_rec_ptr->flip_pairs();
    if (flippables.empty())
        return false;

    const float log_prob_sample_pair = -std::log(float(flippables.size()));

    const flippables_t::item_type& flip_item = *sampling::from_range(flippables.begin(), flippables.end());
    const flippables_t::value_type possibility_test_ptr(flip_item.value);
    if (not possibility_test_ptr->passed())
        return false;
    shared_const_track_ptr_pair tracks_to_flip(std::make_pair(flip_item.first, flip_item.second));

    if (not tracks.contains(tracks_to_flip.first)) { throw std::runtime_error("first track not found"); }
    if (not tracks.contains(tracks_to_flip.second)) { throw std::runtime_error("second track not found"); }
    if (tracks_to_flip.first == tracks_to_flip.second) { throw std::runtime_error("tracks identical"); }

    // flip over tracks
    detail::track_side flip_side = possibility_test_ptr->sample_flip_side();
    BOOST_ASSERT(flip_side != detail::NO_SIDE);
    shared_const_track_ptr_pair flipped_tracks = (flip_side == detail::FRONT_END) ? flip_at_front(tracks_to_flip) :
                                                        flip_at_back(tracks_to_flip);

    //float log_prob_flip_side = -std::log(possibility_test_ptr->num_flips());
    const float num_flips = possibility_test_ptr->num_flips();
    //    std::cout << num_flips << std::endl;

    BOOST_ASSERT(flipped_tracks.first->size()>1);
    BOOST_ASSERT(flipped_tracks.second->size()>1);

    cap_rec_ptr->add_erase_track(tracks_to_flip.first);
    cap_rec_ptr->add_erase_track(tracks_to_flip.second);
    cap_rec_ptr->add_insert_track(flipped_tracks.first);
    cap_rec_ptr->add_insert_track(flipped_tracks.second);
    cap_rec_ptr->set_editing_finished();

    // create new track partition
    shared_track_collection_ptr new_tracks(new track_collection(tracks));

    if (not new_tracks->contains(tracks_to_flip.first)) { throw std::runtime_error("NEW first track not found"); }
    if (not new_tracks->contains(tracks_to_flip.second)) { throw std::runtime_error("NEW second track not found"); }

    new_tracks->remove(tracks_to_flip.first);
    new_tracks->remove(tracks_to_flip.second);
    new_tracks->insert(flipped_tracks.first);
    new_tracks->insert(flipped_tracks.second);

    const_clutter_ptr new_clutter(new clutter_t(start_partition_ptr->clutter()));
    end_partition_ptr = partition_ptr_t(
        new partition(start_partition_ptr->pool(), new_tracks, new_clutter, start_partition_ptr->expansion()));

    detail::flip_possibility_test new_test(flipped_tracks.first, flipped_tracks.second);
    float new_num_flips = new_test.num_flips();
    //if (new_num_flips != num_flips)
    //    std::cout << "number of flips different" << std::endl;


    capability_recorder_ptr new_cap_rec_ptr = new_cap_recorder(*end_partition_ptr);

    const flippables_t& new_flippables = new_cap_rec_ptr->flip_pairs();
    float new_log_prob_sample_pair = -std::log(float(new_flippables.size()));

    float log_back_prob = new_log_prob_sample_pair -std::log(new_num_flips);
    float log_forw_prob = log_prob_sample_pair -std::log(num_flips);
    proposal_mass_ratio = log_back_prob - log_forw_prob;


    return true;
}

} } // biggles::mh_moves
