#include <boost/format.hpp>
#include "biggles/tools/max_creation.hpp"
#include "biggles/model.hpp"
#include <iostream>
#include <sstream>
//#include <stdlib.h>

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

void test_proc1() {
    clutter_ptr clutter_p(new clutter_t());
    const time_stamp max_time = 20;
    size_t total_obs = 0;
    for (time_stamp ts = 0; ts < max_time; ++ts) {
        int num_obs = sampling::uniform_int(1, 10);
        total_obs += size_t(num_obs);
        for (int i = 0; i < num_obs; ++i) {
            clutter_p->insert(new_obs(sampling::uniform_real(0.f, 10.f), sampling::uniform_real(0.f, 10.f), ts));
        }
    }
    partition::track_collection_ptr tracks_p(new track_collection());
    partition_ptr_t start_partition_ptr(new partition(tracks_p, clutter_p));
    partition_ptr_t end_partition_ptr;
    max_creation(start_partition_ptr, end_partition_ptr);
    diag("num obs = %zu", total_obs);
    ok1(end_partition_ptr->observation_count() == total_obs);
    diag("%zu clutter, %zu tracks", end_partition_ptr->clutter().size(), end_partition_ptr->tracks().size());
    BOOST_FOREACH(const boost::shared_ptr<const track>& t, end_partition_ptr->tracks()) {
        diag("duration = %zu", t->duration());
    }

    ok1(end_partition_ptr->volume() == start_partition_ptr->volume());

    model::parameters paras;
    model::clutter_rate(paras) = 2.0f;
    model::birth_rate(paras) = 0.2f;
    model::observation_probability(paras) = 0.8f;
    model::survival_probability(paras) = 0.98;
    model::observation_error_covariance(paras) = Eigen::Matrix2f::Identity() * 0.04f;
    model::process_noise_covariance(paras) = detail::initQ();

    model::log_partition_given_parameters_and_data_density(*end_partition_ptr, paras);
}

int main(int argc, char** argv)
{
    //plan_tests(15);
    plan_no_plan();

    test_proc1();

    return exit_status();
}
