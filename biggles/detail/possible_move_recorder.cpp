#include "detail/cross_over_possibility.hpp"
#include "detail/possible_move_recorder.hpp"
#include "mh_moves/utility.hpp"
#include "tools/debug.hpp"
namespace biggles {

void capability_recorder::remove_track_from_merge_pairs(const shared_const_track_ptr& t)
{
    //mergeable_pairs_.erase_pairs_with_element(t);
    //mergeable_pairs_.erase_items_with_track(t);
}

void capability_recorder::insert_track_into_merge_pairs(
    const shared_const_track_ptr& new_p, const track_collection &tracks)
{
    if(new_p->empty())
        return;

    BOOST_FOREACH(const shared_const_track_ptr& t_p, tracks) {
        if(t_p == new_p)
            continue;
        insert_track_pair_into_merge_pairs(t_p, new_p);
    }
}



void capability_recorder::insert_track_pair_into_transferables(
        const shared_const_track_ptr& t1, const shared_const_track_ptr& t2)
{
    detail::transfer_possibility_ptr new_track_as_donor(new detail::transfer_possibility_test(t1, t2));

    if(new_track_as_donor->passed()) {
        transfer_pairs_.insert(t1, t2, new_track_as_donor);
    }

    detail::transfer_possibility_ptr new_track_as_acceptor(new detail::transfer_possibility_test(t2, t1));

    if (new_track_as_acceptor->passed()) {
        transfer_pairs_.insert(t2, t1, new_track_as_acceptor);
    }
}

void capability_recorder::insert_track_into_transferables(
    const shared_const_track_ptr& new_p, const track_collection &tracks)
{
    if(new_p->empty())
        return;

    BOOST_FOREACH(const shared_const_track_ptr& t_p, tracks) {
        if(t_p == new_p)
            continue;
        insert_track_pair_into_transferables(new_p, t_p);
    }
}

bool capability_recorder::verify_against_track(const track_collection& tracks) const {
    throw std::logic_error("this is not implemented yet");
    return true;
}


/*
void capability_recorder::initial_setup(const track_collection& tc) {
    typedef track_collection::const_iterator tc_iter;
    if (tc.empty())
        return;
    // now there is at least one element in tc
    tc_iter one(tc.begin());
    bool done = false;
    while (not done) {
        insert_track_into_extensibles(*one);
        tc_iter two(one);
        two++;
        for (; two != tc.end(); two++) {
            if (not (*one < *two))
                continue;
            insert_track_pair_into_crossables(*one, *two);
            insert_track_pair_into_merge_pairs(*one, *two);
            insert_track_pair_into_transferables(*one, *two);
        }
        one++;
        if (one == tc.end())
            done = true;
    }
};
*/

void capability_recorder::insert_track_pair_into_crossables(
        const shared_const_track_ptr& t1, const shared_const_track_ptr& t2)
{
    detail::cross_over_possibility_ptr possibility_test(new detail::cross_over_possibility_test(t1, t2));

    if(possibility_test->passed())
        cross_over_pairs_.insert(t1, t2, possibility_test);
}

void capability_recorder::insert_track_pair_into_flippables(
        const shared_const_track_ptr& t1, const shared_const_track_ptr& t2)
{
    detail::flip_possibility_ptr possibility_test(new detail::flip_possibility_test(t1, t2));

    if(possibility_test->passed())
        flip_pairs_.insert(t1, t2, possibility_test);
}

void capability_recorder::insert_track_pair_into_merge_pairs(
        const shared_const_track_ptr& t1, const shared_const_track_ptr& t2)
{
    //if(mh_moves::tracks_could_be_merged(*t1, *t2))
    //    mergeable_pairs_.insert(t1, t2, 0);
}

void capability_recorder::insert_track_into_crossables(const shared_const_track_ptr& new_p,
    const track_collection &tracks)
{
    if(new_p->empty())
        return;

    BOOST_FOREACH(const shared_const_track_ptr& t_p, tracks) {
        if(t_p == new_p)
            continue;
        insert_track_pair_into_crossables(t_p, new_p);
    }
}

size_t capability_recorder::count_inserts_into_crossables(const shared_const_track_ptr& new_p,
    const track_collection &tracks)
{
    size_t result = 0;

    if(new_p->empty()) return result;

    BOOST_FOREACH(const shared_const_track_ptr& t_p, tracks) {
        if(t_p == new_p)
            continue;
        detail::cross_over_possibility_ptr possibility_test(new detail::cross_over_possibility_test(t_p, new_p));
        if(possibility_test->passed())  result++;
    }
    return result;
}

void capability_recorder::insert_track_into_flippables(const shared_const_track_ptr& new_p,
    const track_collection &tracks)
{
    if(new_p->empty())
        return;

    BOOST_FOREACH(const shared_const_track_ptr& t_p, tracks) {
        if(t_p == new_p)
            continue;
        insert_track_pair_into_flippables(t_p, new_p);
    }
}

void capability_recorder::remove_track_from_transferables(const shared_const_track_ptr& track_p) {
    transfer_pairs_.erase_items_with_track(track_p);
}

void capability_recorder::remove_track_from_crossables(const shared_const_track_ptr& track_p) {
    cross_over_pairs_.erase_items_with_track(track_p);
}

void capability_recorder::remove_track_from_flippables(const shared_const_track_ptr& track_p) {
    flip_pairs_.erase_items_with_track(track_p);
}

void capability_recorder::initial_setup(const track_collection& tracks) {
    BOOST_FOREACH(const shared_const_track_ptr& t_p, tracks) {
        insert_into_records(t_p, tracks);
    }
};


} // namespace biggles
