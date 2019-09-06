#include <boost/assert.hpp>

#include "mh_moves/mh_moves.hpp"
#include "detail/fun.hpp"

namespace biggles { namespace mh_moves {

std::string track_to_str(const track & track0, time_stamp end_ts) {
    time_stamp ts=track0.first_time_stamp();
    std::stringstream track_ss;
    track_ss << std::string(ts, ' ');
    for (; ts < track0.last_time_stamp(); ++ts) {
        observation_collection::const_range obs_at_t(track0.observations().at_time_stamp(ts));
        if (obs_at_t.first != obs_at_t.second)
            track_ss << x(*obs_at_t.first);
        else
           track_ss << "_";
    }
    while (ts++ < end_ts) {
        track_ss << " ";
    }
    return track_ss.str();
}



bool propose(move_type move,
             const partition_ptr_t& start_partition,
             partition_ptr_t& end_partition,
             float& proposal_density_ratio)
{
    proposal_density_ratio = -std::numeric_limits<float>::max();
    try {
        switch(move) {
            case mh_moves::BIRTH:
                return mh_moves::birth(start_partition, end_partition, proposal_density_ratio);
                break;
            case mh_moves::DEATH:
                return mh_moves::death(start_partition, end_partition, proposal_density_ratio);
                break;
            case mh_moves::EXTEND:
                return mh_moves::extend(start_partition, end_partition, proposal_density_ratio);
                break;
            case mh_moves::REDUCE:
                return mh_moves::reduce(start_partition, end_partition, proposal_density_ratio);
                break;
            case mh_moves::UPDATE:
                return mh_moves::update(start_partition, end_partition, proposal_density_ratio);
                break;
            case mh_moves::SPLIT:
                return mh_moves::split(start_partition, end_partition, proposal_density_ratio);
                break;
            case mh_moves::MERGE:
                return mh_moves::merge(start_partition, end_partition, proposal_density_ratio);
                break;
            case mh_moves::TRANSFER:
                return mh_moves::transfer(start_partition, end_partition, proposal_density_ratio);
                break;
            case mh_moves::CROSS_OVER:
                return mh_moves::cross_over(start_partition, end_partition, proposal_density_ratio);
                break;
            case mh_moves::FLIP:
                return mh_moves::flip(start_partition, end_partition, proposal_density_ratio);
                break;
            case mh_moves::NONE:
                // NOP
                throw std::logic_error("you shouldn't try to make an NONE move");
                proposal_density_ratio = 0.f;
                return true;
                break;
            case mh_moves::IDENTITY:
                throw std::logic_error("you shouldn't try to make an identity move");
                proposal_density_ratio = 0.f;
                return true;
                break;
            default:
                BOOST_ASSERT(false);
                break;
        }
    } catch (const std::exception& e) {
        std::cerr << "culprit: " << move_name(move) << std::endl;
        throw;
    }

    // should be unreachable
    BOOST_ASSERT(false);
    return false;
}

move_fun get_move(move_type type) {
    switch (type) {
        case BIRTH      : return birth;
        case DEATH      : return death;
        case SPLIT      : return split;
        case MERGE      : return merge;
        case EXTEND     : return extend;
        case REDUCE     : return reduce;
        case UPDATE     : return update;
        case TRANSFER   : return transfer;
        case CROSS_OVER : return cross_over;
        case FLIP : return flip;
        default : throw std::logic_error("Not a executable move type");
    }
}

move_type select_move_type(const partition& part) {
    if (part.tracks().size() == 0)
        return BIRTH;
    const capability_recorder_ptr& cr_ptr = part.get_capability_recorder();
    BOOST_ASSERT(cr_ptr.get() != 0);
    weighted_move_sampler selector;
    selector.add(DEATH).add(REDUCE).add(SPLIT).add(UPDATE).add(MERGE);
    selector.add(EXTEND);
    const float num_clutter(part.clutter().size());
    const float num_obs(part.observation_count());
    if (cr_ptr->transfer_pairs().size() > 0)
        selector.add(TRANSFER, std::min(1.f, 4.f/3.f * (num_obs-num_clutter)/num_obs));
    if (not cr_ptr->cross_over_pairs().empty())
        selector.add(CROSS_OVER, std::min(1.f, 4.f/3.f * (num_obs-num_clutter)/num_obs));
    // there is still somewhere a bug.
    //if (not cr_ptr->flip_pairs().empty())
    //    selector.add(FLIP, std::min(1.f, 4.f/3.f * (num_obs-num_clutter)/num_obs));
    if (part.clutter().size() > 1)
        selector.add(BIRTH);
    return selector();
}

} } // biggles::mh_moves
