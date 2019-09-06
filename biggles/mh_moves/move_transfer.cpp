#include <boost/assert.hpp>
#include <cmath>
#include <limits>
#include <utility>

#include "mh_moves/mh_moves.hpp"
#include "sampling/simple.hpp"
#include "mh_moves/utility.hpp"
#include "detail/transfer_possibility.hpp"
#include "tools/debug.hpp"

namespace biggles { namespace mh_moves {

bool transfer(const partition_ptr_t& start_partition_ptr, partition_ptr_t& end_partition_ptr, float& proposal_density_ratio) {
    capability_recorder_ptr cap_rec_ptr = start_partition_ptr->get_capability_recorder();
    BOOST_ASSERT(cap_rec_ptr.get() not_eq 0); // there actually is something assigned

    // copy the start partition so that start_partition can equal end_partition and we can modify end_partition
    const track_collection& tracks(start_partition_ptr->tracks());

    const transferables_t& transferables(cap_rec_ptr->transfer_pairs());
    if (transferables.empty())
        return false;
    const transferables_t::item_type& transfer_item = *sampling::from_range(transferables.begin(), transferables.end());
    const transferables_t::value_type possibility_test_ptr(transfer_item.value);
    if (not possibility_test_ptr->passed())
        return false;
    shared_const_track_ptr_pair donor_acceptor(std::make_pair(transfer_item.first, transfer_item.second));
    float log_prob_time_selection = possibility_test_ptr->log_prob_transfer_time();
    time_stamp transfer_time = possibility_test_ptr->sample_transfer_time();
    float log_pair_selection = - std::log(float(transferables.size()));

    // transfer observation
    shared_const_track_ptr donor = donor_acceptor.first;
    shared_const_track_ptr acceptor = donor_acceptor.second;
    observation transfer_obs = *(donor->observations().at_time_stamp(transfer_time).first);
    shared_observation_collection_ptr new_donor_obs(new observation_collection(donor->observations()));
    new_donor_obs->erase(transfer_obs);
    shared_observation_collection_ptr new_acceptor_obs(new observation_collection(acceptor->observations()));
    new_acceptor_obs->insert(transfer_obs);

    shared_track_collection_ptr new_tracks(new track_collection(tracks));
    new_tracks->remove(donor);
    new_tracks->remove(acceptor);

    shared_track_ptr new_donor = shared_track_ptr(new track(
        donor->first_time_stamp(), donor->last_time_stamp(),
        new_donor_obs->begin(), new_donor_obs->end(), 1.f));
    shared_track_ptr new_acceptor = shared_track_ptr(new track(
        acceptor->first_time_stamp(), acceptor->last_time_stamp(),
        new_acceptor_obs->begin(), new_acceptor_obs->end(), 1.f));

    BOOST_ASSERT(new_donor->size()>1);
    BOOST_ASSERT(new_acceptor->size()>1);

    new_tracks->insert(new_donor);
    new_tracks->insert(new_acceptor);

    cap_rec_ptr->add_erase_track(donor);
    cap_rec_ptr->add_erase_track(acceptor);
    cap_rec_ptr->add_insert_track(new_donor);
    cap_rec_ptr->add_insert_track(new_acceptor);
    cap_rec_ptr->set_editing_finished();

    const_clutter_ptr new_clutter(new clutter_t(start_partition_ptr->clutter()));
    end_partition_ptr = partition_ptr_t(
        new partition(start_partition_ptr->pool(), new_tracks, new_clutter, start_partition_ptr->expansion())
    );

    float forward_log_prob = log_prob_time_selection + log_pair_selection;

    // backward probability
    detail::transfer_possibility_test backward_test(new_acceptor, new_donor);
    BOOST_ASSERT(backward_test.passed());
    // FIXME calc the new transfer pair size
    long transfer_pairs_size = 0;
    {
        transfer_pairs_size = transferables.size();
        detail::ContainsTrack<transferables_t::item_type> contains_donor(donor);
        detail::ContainsTrack<transferables_t::item_type> contains_acceptor(acceptor);
        transfer_pairs_size -= std::count_if(transferables.begin(), transferables.end(), contains_donor);
        transfer_pairs_size -= std::count_if(transferables.begin(), transferables.end(), contains_acceptor);
        transfer_pairs_size += 1; // the donor-acceptor item has been subtracted twice
        if (detail::transfer_possibility_test(acceptor, donor).passed())
            transfer_pairs_size += 1;
        for (track_collection::const_iterator tit = tracks.begin(); tit != tracks.end(); ++tit) {
            if (*tit == donor or *tit == acceptor)
                continue;
            detail::transfer_possibility_test new_donor_test(new_donor, *tit);
            if (new_donor_test.passed())
                transfer_pairs_size++;
            detail::transfer_possibility_test new_donor_test2(*tit, new_donor);
            if (new_donor_test2.passed())
                transfer_pairs_size++;
            detail::transfer_possibility_test new_acceptor_test(new_acceptor, *tit);
            if (new_acceptor_test.passed())
                transfer_pairs_size++;
            detail::transfer_possibility_test new_acceptor_test2(*tit, new_acceptor);
            if (new_acceptor_test2.passed())
                transfer_pairs_size++;
        }
        if (detail::transfer_possibility_test(new_donor, new_acceptor).passed())
            transfer_pairs_size += 1; // new_donor-new_acceptor is transferable
        detail::transfer_possibility_test acc_don_test(new_acceptor, new_donor);
        if (acc_don_test.passed())
                transfer_pairs_size++;
    }
    //OK(transfer_pairs_size);
    float backward_log_prob = backward_test.log_prob_transfer_time() - logf(transfer_pairs_size);

    proposal_density_ratio = backward_log_prob - forward_log_prob;
    return true;

}

}}
