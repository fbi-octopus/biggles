#include <boost/assert.hpp>
#include <cmath>
#include <limits>
#include <utility>

#include "mh_moves/mh_moves.hpp"
#include "sampling/simple.hpp"
#include "mh_moves/utility.hpp"
#include "detail/flip_possibility.hpp"

namespace biggles { namespace mh_moves {

shared_const_track_ptr_pair cross_over_at(shared_const_track_ptr_pair& tracks_to_cross, time_stamp time) {
    const track& bro(*(tracks_to_cross.first));
    const track& sis(*(tracks_to_cross.second));

    time_stamp bro_first = bro.first_time_stamp();
    time_stamp bro_last = bro.last_time_stamp();
    time_stamp sis_first = sis.first_time_stamp();
    time_stamp sis_last = sis.last_time_stamp();

    observation_collection::const_iterator bro_front_first = bro.observations().begin();
    observation_collection::const_iterator bro_front_last = bro.observations().lower_bound_for_time_stamp(time);
    observation_collection::const_iterator bro_back_first = bro_front_last;
    observation_collection::const_iterator bro_back_last = bro.observations().end();

    observation_collection::const_iterator sis_front_first = sis.observations().begin();
    observation_collection::const_iterator sis_front_last = sis.observations().lower_bound_for_time_stamp(time);
    observation_collection::const_iterator sis_back_first = sis_front_last;
    observation_collection::const_iterator sis_back_last = sis.observations().end();

    shared_track_ptr brosis(new track(bro_first, sis_last, bro_front_first, bro_front_last));
    for ( ; sis_back_first != sis_back_last; sis_back_first++)
        brosis->insert(*sis_back_first);
    shared_track_ptr sisbro(new track(sis_first, bro_last, sis_front_first, sis_front_last));
    for ( ; bro_back_first != bro_back_last; bro_back_first++)
        sisbro->insert(*bro_back_first);
    return shared_const_track_ptr_pair(brosis, sisbro);
}

bool cross_over(const partition_ptr_t& start_partition_ptr, partition_ptr_t& end_partition_ptr,
    float& proposal_density_ratio)
{
    const track_collection& tracks(start_partition_ptr->tracks());

    proposal_density_ratio = 0.f; // hmm ... this needs to be reconsidered; does a cross-over create/remove cross-over options?

    capability_recorder_ptr cap_rec_ptr = start_partition_ptr->get_capability_recorder();
    BOOST_ASSERT(cap_rec_ptr.get() not_eq 0); // there actually is something assigned

    const crossables_t& crossables = cap_rec_ptr->cross_over_pairs();
    if (crossables.empty())
        return 0;

    int num_crossables_before = crossables.size();

    const crossables_t::item_type& cross_over_item = *sampling::from_range(crossables.begin(), crossables.end());
    const crossables_t::value_type possibility_test_ptr(cross_over_item.value);
    if (not possibility_test_ptr->passed())
        return false;
    shared_const_track_ptr_pair tracks_to_cross(std::make_pair(cross_over_item.first, cross_over_item.second));
    time_stamp cross_over_time = possibility_test_ptr->sample_cross_over_time();

    if (not tracks.contains(tracks_to_cross.first)) { throw std::runtime_error("first track not found"); }
    if (not tracks.contains(tracks_to_cross.second)) { throw std::runtime_error("second track not found"); }
    if (tracks_to_cross.first == tracks_to_cross.second) { throw std::runtime_error("tracks identical"); }

    // cross over tracks
    shared_const_track_ptr_pair crossed_tracks = cross_over_at(tracks_to_cross, cross_over_time);

    BOOST_ASSERT(crossed_tracks.first->size()>1);
    BOOST_ASSERT(crossed_tracks.second->size()>1);

    cap_rec_ptr->add_erase_track(tracks_to_cross.first);
    cap_rec_ptr->add_erase_track(tracks_to_cross.second);
    cap_rec_ptr->add_insert_track(crossed_tracks.first);
    cap_rec_ptr->add_insert_track(crossed_tracks.second);
    cap_rec_ptr->set_editing_finished();

    // create new track partition
    shared_track_collection_ptr new_tracks(new track_collection(tracks));

    if (not new_tracks->contains(tracks_to_cross.first)) { throw std::runtime_error("NEW first track not found"); }
    if (not new_tracks->contains(tracks_to_cross.second)) { throw std::runtime_error("NEW second track not found"); }

    new_tracks->remove(tracks_to_cross.first);
    new_tracks->remove(tracks_to_cross.second);
    new_tracks->insert(crossed_tracks.first);
    new_tracks->insert(crossed_tracks.second);

    int num_crossables_after = cap_rec_ptr->num_crossables_if_commited(*new_tracks);
    proposal_density_ratio = -log(num_crossables_after) + log(num_crossables_before);

    const_clutter_ptr new_clutter(new clutter_t(start_partition_ptr->clutter()));
    end_partition_ptr = partition_ptr_t(
        new partition(start_partition_ptr->pool(), new_tracks, new_clutter, start_partition_ptr->expansion()));
    return true;
}

}

}
