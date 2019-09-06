#include "biggles/kalman_filter.hpp"
#include "biggles/observation.hpp"
#include "biggles/partition.hpp"
#include "biggles/samplers.hpp"
#include "biggles/simulate.hpp"
#include "biggles/track.hpp"
#include "biggles/server/engine.hpp"
#include "biggles/server/stepper.hpp"
#include "biggles/detail/track_pair_data.hpp"
#include "biggles/tools/debug.hpp"
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <cmath>
#include <ctime>
#include <inttypes.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <utility>
#include <vector>

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}


const float PI = 3.141592653589793238462643383279502884197169399375f;

static boost::mt19937 random_generator(std::time(0));

using namespace biggles;

typedef std::pair<float, float> float_pair;

float_pair get_dx_dy(const float& rmax) {
    float phi = sampling::uniform_real(0.f, 2.f*PI);
    float r = sampling::uniform_real(0.f, rmax);
    return std::make_pair(r * std::cos(phi), r * std::sin(phi));
}

bool operator>>(char needle, const std::string& haystack) { return haystack.find(needle) != std::string::npos; }

shared_track_ptr make_track(time_stamp begin_ts, time_stamp end_ts, const float x0, const float y0,
    const float& rmax, const std::string& obs_pos)
{
    BOOST_ASSERT(0<=begin_ts);
    BOOST_ASSERT(begin_ts<end_ts);
    BOOST_ASSERT(time_stamp(obs_pos.size()) + begin_ts == end_ts);
    std::deque<observation> obs;
    const std::string gap(". -");
    float x = x0;
    float y = y0;
    float r = rmax;
    if (not (obs_pos.at(0) >> gap))
        obs.push_back(new_obs(x, y, 0+ begin_ts));
    for (time_stamp ts = 1; ts < end_ts - begin_ts; ++ts) {
        if (not (obs_pos.at(ts) >> gap)) {
            float_pair dxdy = get_dx_dy(r);
            r = rmax;
            x += dxdy.first;
            y += dxdy.second;
            obs.push_back(new_obs(x, y, ts + begin_ts));
        } else {
            r+= rmax;
        }
    }
    //OK(std::distance(obs.begin(), obs.end()));
    return shared_track_ptr(new track(begin_ts, end_ts, obs.begin(), obs.end()));
}

void test1(const model::parameters& para) {
    //boost::shared_ptr<track_collection> tracks(new track_collection());
    shared_track_ptr t1 = make_track(0, 10, 0.f, 0.f, 0.5f,   "012.456.89");
    shared_track_ptr t2 = make_track(5, 15, 0.5f, 0.5f, 0.5f,      "012.456.89");
    shared_track_ptr t3 = make_track(0, 5, 0.f, 0.5f, 0.5f,   ".12.3");
    shared_track_ptr t4 = make_track(10, 15, 0.5f, 0.f, 0.5f,           "01..3");
    //boost::shared_ptr<observation_collection> clutter;
    typedef detail::track_pair_item_t<int> item_t;
    typedef detail::sym_track_pair_item_t<int> sym_item_t;
    detail::new_track_pair_data_t < item_t > list1;
    detail::new_track_pair_data_t < sym_item_t > list2;
    list1.insert(t2, t4, 0);
    list1.insert(t1, t4, 0);
    list1.insert(t3, t4, 0);
    list1.insert(t2, t1, 0);
    list1.insert(t2, t3, 0);
    list1.insert(t4, t1, 0);
    list1.insert(t2, t4, 0);
    diag("size = %zu", list1.size());
    ok1(list1.size() == 6);
    list2.insert(t2, t4, 0);
    list2.insert(t1, t4, 0);
    list2.insert(t3, t4, 0);
    list2.insert(t2, t1, 0);
    list2.insert(t2, t3, 0);
    list2.insert(t4, t1, 0);
    list2.insert(t2, t4, 0);
    diag("size = %zu", list2.size());
    ok1(list2.size() == 5);
    ok1(is_sorted(list1.begin(), list1.end()));
    ok1(is_sorted(list2.begin(), list2.end()));
    list1.erase_items_with_track(t3);
    list2.erase_items_with_track(t3);
    ok1(list1.size() == 4);
    ok1(list2.size() == 3);
    ok1(is_sorted(list1.begin(), list1.end()));
    ok1(is_sorted(list2.begin(), list2.end()));
    shared_track_ptr t5 = make_track(10, 15, 0.5f, 0.6f, 0.5f,           "0.234");
    list1.insert(t5, t2, 0);
    list1.insert(t5, t3, 0);
    list2.insert(t5, t2, 0);
    list2.insert(t5, t3, 0);
    ok1(is_sorted(list1.begin(), list1.end()));
    ok1(is_sorted(list2.begin(), list2.end()));
}

void test2() {
    shared_track_ptr t1 = make_track(0, 10, 0.f, 0.f, 0.5f,   "0124.56.89");
    shared_track_ptr t2 = make_track(5, 15, 0.5f, 0.5f, 0.5f,      "012.456.89");
    shared_track_ptr t3 = make_track(0, 5, 0.f, 0.5f, 0.5f,   "0.2.3");
    shared_track_ptr t4 = make_track(10, 15, 0.5f, 0.f, 0.5f,           "01..3");
    shared_track_ptr t5 = make_track(10, 15, 0.5f, 0.6f, 0.5f,          "0.234");
    boost::shared_ptr<track_collection> tracks(new track_collection());
    boost::shared_ptr<observation_collection> clutter;
    tracks->insert(t1);
    tracks->insert(t2);
    diag("size t1 = %zu", t1->size());
    diag("size t2 = %zu", t2->size());
    tracks->insert(t3);
    tracks->insert(t4);
    //tracks->insert(t5);
    capability_recorder cap_rec(0, 15, *tracks);
    diag("cross-over %zu", cap_rec._size_cross_over_pairs_());
    diag("mergeables %zu", cap_rec._size_mergeable_pairs_());
    diag("transferable %zu", cap_rec._size_transfer_pairs_());
    diag("extendibles %zu", cap_rec._size_extendible_tracks_());
    const transferables_t &trans_pairs = cap_rec.transfer_pairs();
    for (transferables_t::const_iterator it = trans_pairs.begin(); it != trans_pairs.end(); ++it) {
        diag("size 1st %zu, size 2nd %zu", it->first->size(), it->second->size());
    }
    cap_rec.add_erase_track(t3);
    cap_rec.set_editing_finished();
    cap_rec.commit_changes(*tracks);
    diag("remove t3");
    diag("cross-over %zu", cap_rec._size_cross_over_pairs_());
    diag("mergeables %zu", cap_rec._size_mergeable_pairs_());
    diag("transferable %zu", cap_rec._size_transfer_pairs_());
    diag("extendibles %zu", cap_rec._size_extendible_tracks_());
    ok1(is_sorted(cap_rec.transfer_pairs().begin(), cap_rec.transfer_pairs().end()));
    ok1(is_sorted(cap_rec.mergeable_pairs().begin(), cap_rec.mergeable_pairs().end()));
    ok1(is_sorted(cap_rec.cross_over_pairs().begin(), cap_rec.cross_over_pairs().end()));
    ok1(is_sorted(cap_rec.extendible_tracks().begin(), cap_rec.extendible_tracks().end()));
}

void test3() {
    shared_track_ptr t1 = make_track(0, 10, 0.f, 0.f, 0.5f,   "0124.56.89");
    shared_track_ptr t2 = make_track(5, 15, 0.5f, 0.5f, 0.5f,      "012.456.89");
    shared_track_ptr t3 = make_track(0, 5, 0.f, 0.5f, 0.5f,   "0.2.3");
    shared_track_ptr t4 = make_track(10, 15, 0.5f, 0.f, 0.5f,           "01..3");
    shared_track_ptr t5 = make_track(10, 15, 0.5f, 0.6f, 0.5f,          "0.234");
    boost::shared_ptr<track_collection> tracks(new track_collection());
    boost::shared_ptr<observation_collection> clutter;
    tracks->insert(t1);
    tracks->insert(t2);
    diag("size t1 = %zu", t1->size());
    diag("size t2 = %zu", t2->size());
    tracks->insert(t3);
    tracks->insert(t4);
    tracks->insert(t5);
    capability_recorder cap_rec(0, 15, *tracks);
    int num_crossable_before = cap_rec._size_cross_over_pairs_();
    cap_rec.add_erase_track(t3);
    cap_rec.set_editing_finished();
    int num_predicted = cap_rec.num_crossables_if_commited(*tracks);
    cap_rec.commit_changes(*tracks);
    int num_crossable_after = cap_rec._size_cross_over_pairs_();
    diag("num_crossable_before %d", num_crossable_before);
    diag("num_predicted %d", num_predicted);
    diag("num_crossable_after %d", num_crossable_after);
    //crossables_t crossables = cap_rec.cross_over_pairs();
    ok1(num_predicted == num_crossable_after);
    num_crossable_before = cap_rec._size_cross_over_pairs_();
    cap_rec.add_erase_track(t5);
    cap_rec.set_editing_finished();
    num_predicted = cap_rec.num_crossables_if_commited(*tracks);
    cap_rec.commit_changes(*tracks);
    num_crossable_after = cap_rec._size_cross_over_pairs_();
    diag("num_crossable_before %d", num_crossable_before);
    diag("num_predicted %d", num_predicted);
    diag("num_crossable_after %d", num_crossable_after);
    //crossables_t crossables = cap_rec.cross_over_pairs();
    ok1(num_predicted == num_crossable_after);
}


int main(int argc, char** argv)
{
    plan_no_plan();

    float sigma_R = 0.1f;

    // initialise model parameters
    model::parameters params;
    model::mean_new_tracks_per_frame(params)           = 0.1f;
    model::mean_false_observations_per_frame(params)   = 0.05f;
    model::frame_to_frame_survival_probability(params) = 0.95f;
    model::generate_observation_probability(params)    = 0.8f;
    model::observation_error_covariance(params)       << sigma_R * sigma_R, 0, 0, sigma_R * sigma_R;

    test1(params);
    test2();
    test3();
}
