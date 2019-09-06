#include "detail/track_pair_data.hpp"

namespace biggles { namespace detail {
    /*
    void track_pair_data_t::erase_items_with_track(const shared_const_track_ptr& track_p) {
        for(;;) {
            // This erases a item pair with track_p as first or second element in the key
            iterator i(track_pair_data_.begin());
            for(; (i != track_pair_data_.end()) and (i->first.first != track_p) and (i->first.second != track_p); ++i);
            if(i != track_pair_data_.end()) {
                track_pair_data_.erase(i);
            } else {
                break;
            }
        }
    }
    */
}} // biggles::detail
