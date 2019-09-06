#include "biggles/observation.hpp"
#include "biggles/samplers.hpp"
#include "biggles/observation_collection.hpp"
#include "biggles/track.hpp"
#include "biggles/track_collection.hpp"
#include "biggles/partition.hpp"
#include "biggles/detail/random.hpp"
#include "biggles/sampling/simple.hpp"
#include <boost/foreach.hpp>
#include <boost/math/distributions.hpp>
#include <boost/timer.hpp>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <string>

#include "biggles/tools/sundries.hpp"

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

void print_obs(observation obs) {
    std::cout << "x=" << x(obs) << ", y= " << y(obs) << ", t= " << t(obs);
}

typedef std::deque<float> data_container;
typedef std::pair<float, float> data_pair;

float sample_beta_alt(float alpha, float beta) {
    boost::math::beta_distribution<float> dist(alpha, beta);
    return quantile(dist, sampling::uniform_real());
}

int test_normal(int total) {
    diag("Normal (%d)",total);
    data_container samples(total);
    for (int i = 0; i < total; ++i) {
        samples[i] = sampling::sample_normal();
    }
    data_pair ms = mean_stdev(samples.begin(), samples.end());
    diag(" mean %.2f, stdev %.2f", ms.first, ms.second);
    return 0;
}

int test_uniform() {
    diag("Uniform 1 ");
    float minimum = 1.f;
    float maximum = 0.f;
    float current = 0.f;
    size_t ctr = 0;
    while (minimum > 0.f and maximum < 1.f) {
        current = sampling::uniform_real();
        if (current < minimum) {
            minimum = current;
            diag("  (%g, %g)", minimum, 1.f - maximum);
        }
        if (current > maximum) {
            maximum = current;
            diag("  (%g, %g)", minimum, 1.f - maximum);
        }
        ++ctr;
    }
    diag(" max == 1 ? %d",  maximum == 1);
    diag(" min == 0 ? %d",  minimum == 0);
    diag(" nloops = %zu", ctr);
    return 0;
}

int test_uniform_inv() {
    diag("Uniform inverse ");
    float minimum = 1.f;
    float maximum = 0.f;
    float current = 0.f;
    int ctr = 0;
    while (minimum > 0.f and maximum < 1.f) {
        current = 1.f - sampling::uniform_real();
        if (current < minimum) {
            minimum = current;
        }
        if (current > maximum) {
            maximum = current;
        }
        ++ctr;
    }
    diag(" max == 1 ? %d",  maximum == 1);
    diag(" min == 0 ? %d",  minimum == 0);
    diag(" nloops = %d", ctr);
    return 0;
}

int test_factorial() {
    model::log_factorial &fact = model::log_factorial::get();
    diag("  6! = %f", fact(6));
    diag(" 12! = %f", fact(12));
    diag("  8! = %f", fact(8));
    diag("  3! = %f", fact(3));
    diag("  0! = %f", fact(0));
    diag("  1! = %f", fact(1));
    diag("113! = %f", fact(113));
    return 0;
}

int test_trunc_beta(int total) {
    diag("Truncated Beta (%d)",total);
    data_container samples(total);
    for (int i = 0; i < total; ++i) {
        samples[i] = sampling::sample_truncated_beta(1.f+0.0001f,1.f + 40.0001f, 0.8, 1.0);
    }
    data_pair ms = mean_stdev(samples.begin(), samples.end());
    diag(" mean %.6f, stdev %.6f, min %.6f, max %.6f", ms.first, ms.second,
        *std::min_element(samples.begin(), samples.end()),
        *std::max_element(samples.begin(), samples.end()));
    return 0;
}

int test_beta(int total) {
    diag("Beta (%d)",total);
    data_container samples(total);
    for (int i = 0; i < total; ++i) {
        samples[i] = sampling::sample_beta(1.f+10.f,1.f + 10.f);
    }
    data_pair ms = mean_stdev(samples.begin(), samples.end());
    diag(" mean %.6f, stdev %.6f", ms.first, ms.second);
    return 0;
}

int test_beta2(int total) {
    diag("Beta2 (%d)",total);
    data_container samples(total);
    for (int i = 0; i < total; ++i) {
        samples[i] = sampling::sample_beta(1.f+10.f,1.f + 10.f);
    }
    data_pair ms = mean_stdev(samples.begin(), samples.end());
    diag(" mean %.6f, stdev %.6f", ms.first, ms.second);
    return 0;
}

int test_beta3(int total) {
    diag("Beta3 (%d)",total);
    data_container samples(total);
    for (int i = 0; i < total; ++i) {
        samples[i] = sampling::sample_beta(1.0001f+2785.f,1.0001f);
    }
    data_pair ms = mean_stdev(samples.begin(), samples.end());
    diag(" mean %.6f, stdev %.6f", ms.first, ms.second);
    return 0;
}

int test_gamma(int total) {
    diag("Gamma (%d)",total);
    data_container samples(total);
    for (int i = 0; i < total; ++i) {
        samples[i] = sampling::sample_gamma(1.f+71.f,1.f);
    }
    data_pair ms = mean_stdev(samples.begin(), samples.end());
    ok1(size_t(std::count_if(samples.begin(), samples.end(), std::isfinite<float>)) == samples.size());
    diag(" mean %.2f, stdev %.2f", ms.first, ms.second);
    return 0;
}

int test_gamma0(int total) {
    diag("Gamma 0 (%d)",total);
    data_container samples(total);
    for (int i = 0; i < total; ++i) {
        samples[i] = sampling::sample_gamma(1.f,1.f);
        if (not std::isfinite(samples[i])) {
            diag("value is not finite");
        }
    }
    data_pair ms = mean_stdev(samples.begin(), samples.end());
    ok1(size_t(std::count_if(samples.begin(), samples.end(), std::isfinite<float>)) == samples.size());
    diag(" mean %.5f, stdev %.5f", ms.first, ms.second);
    return 0;
}

int test_gamma3(int total) {
    diag("Gamma 3 (%d)",total);
    data_container samples(total);
    for (int i = 0; i < total; ++i) {
        samples[i] = sampling::sample_gamma(1.f,1.f);
        if (not std::isfinite(samples[i])) {
            diag("value is not finite");
        }
    }
    data_pair ms = mean_stdev(samples.begin(), samples.end());
    ok1(size_t(std::count_if(samples.begin(), samples.end(), std::isfinite<float>)) == samples.size());
    diag(" mean %.5f, stdev %.5f", ms.first, ms.second);
    return 0;
}

int test_gamma1() {
    diag("Gamma 1");
    float minimum = 1.f;
    float current = 0.f;
    size_t ctr = 0;
    while (minimum > 0.f) {
        current = sampling::sample_gamma(1.f,1.f);
        if (current < minimum) {
            minimum = current;
            //diag("  min = %g", minimum);
        }
        ++ctr;
    }
    diag(" nloops = %zu", ctr);
    return 0;
}

int test_gamma2() {
    diag("Gamma 2");
    float minimum = 1.f;
    float current = 0.f;
    int ctr = 0;
    while (minimum > 0.f) {
        if (ctr == std::numeric_limits<int>::max())
            break;
        current = sampling::sample_gamma0(1.f,1.f);
        if (current < minimum) {
            minimum = current;
            diag("min = %g", minimum);
        }
        ++ctr;
    }
    diag(" nloops = %d", ctr);
    return 0;
}

int test_clutter_rate_sampling(int total) {
    observation_collection obs_coll;
    for (time_stamp ts = 0; ts < 50; ++ts)
        obs_coll.insert(new_obs(0.f, 0.f, ts));
    boost::shared_ptr<track_collection> tracks(new track_collection());
    clutter_ptr clutter(new clutter_t(obs_coll.begin(), obs_coll.end()));
    const partition part(tracks, clutter);
    diag("partition duration = %zu", part.duration());
    diag("clutter size = %zu", part.clutter().size());
    model::parameters para;
    for (int i = 0; i < total; ++i) {
        try {
            sample_clutter_rate_given_partition(part, para);
        }
        catch (const std::exception &e) {
            diag("loop %d", i);
        }
    }
    return 0;
}

int main(int argc, char** argv)
{
    unsigned int seed=std::time(0);
    // seed the PRNG to ensure consistent results across runs.
    // biggles::detail::seed_prng(0xfacedead);
    //seed = 0x55fa8e3d;
    biggles::detail::seed_prng(seed);

    bool verbose = (argc > 1) && (0 == strcmp(argv[1], "-v"));

    plan_tests(13+3-1);

#ifdef NDEBUG
    diag("NDEBUG is defined");
#else
    diag("NDEBUG is **not** defined");
#endif

    observation_collection cl_obs;

    cl_obs.insert(new_obs(-3, 0, 0));
    cl_obs.insert(new_obs(-2, 0, 1));
    //cl_obs.insert(new_obs(-1, 0, 2));
    cl_obs.insert(new_obs( 0, 0, 3));
    cl_obs.insert(new_obs( 1, 0, 4));
    cl_obs.insert(new_obs( 2, 0, 5));
    cl_obs.insert(new_obs( 3, 0, 6));

    observation_collection collection2;
    collection2.insert(new_obs( 0, 0, 1));
    collection2.insert(new_obs( 0, 0.1, 3));
    collection2.insert(new_obs( 0, 0.2, 3));
    collection2.insert(new_obs( 0, 0.3, 3));
    collection2.insert(new_obs( .1, 0, 3));
    collection2.insert(new_obs( .2, 0, 3));

    //const float sqrt2 = sqrtf(2.0f);

    const int total = 10000;
    const float lambda = 2.5f;

    data_container poisson_samples(total);

    for (int i = 0; i < total; ++i) {
        poisson_samples[i] = sampling::poisson(lambda);
    }


    data_pair pois_ms = mean_stdev(poisson_samples.begin(), poisson_samples.end());
    data_pair pois_minmax = std::make_pair(*std::min_element(poisson_samples.begin(), poisson_samples.end()),
        *std::max_element(poisson_samples.begin(), poisson_samples.end()));
    if (verbose) {
        std::cout
            << "Poisson distribution: " << std::endl
            << "mean = " << pois_ms.first
            << ", "
            << "stdev = " << pois_ms.second
            << ", "
            << "stdev*stdev = " << (pois_ms.second * pois_ms.second)
            << ", "
            << "min = " << pois_minmax.first
            << ", "
            << "max = " << pois_minmax.second
            << std::endl;
        ;
    }

    ok(*std::min_element(poisson_samples.begin(), poisson_samples.end()) == 0, "min Poisson == 0");
    ok(fabsf(lambda - pois_ms.first) < 0.1, "Poisson: sample mean near lambda");

    data_container exponential_samples(total);

    for (int i = 0; i < total; ++i) {
        exponential_samples[i] = sampling::exponential(lambda);
    }


    data_pair exp_ms = mean_stdev(exponential_samples.begin(), exponential_samples.end());
    data_pair exp_minmax = std::make_pair(*std::min_element(exponential_samples.begin(), exponential_samples.end()),
        *std::max_element(exponential_samples.begin(), exponential_samples.end()));
    if (verbose) {
        std::cout
            << "Exponential distribution: " << std::endl
            << "mean = " << exp_ms.first
            << ", "
            << "stdev = " << exp_ms.second
            << ", "
            << "stdev*stdev = " << (exp_ms.second * exp_ms.second)
            << ", "
            << "min = " << exp_minmax.first
            << ", "
            << "max = " << exp_minmax.second
            << std::endl;
        ;
    }

    ok(*std::min_element(exponential_samples.begin(), exponential_samples.end()) >= 0.0f, "min Exponential >= 0");
    ok(*std::min_element(exponential_samples.begin(), exponential_samples.end()) < 0.01f, "min Exponential near 0");
    ok(fabsf(1.0f/lambda - exp_ms.first) < 0.1, "Exponential: sample mean near lambda");


    float half_pois = ceilf(pois_minmax.second/2.0);
    bool all_ok = true;
    for (int num = 0 ; num<= int(half_pois) ; ++num) {
        bool done = false;
        for (int i = 0; i < total; ++i) {
            if (num == *sampling::from_range(poisson_samples.begin(), poisson_samples.end())) {
                done = true;
                if (verbose) {
                    std::cout
                        << "from_range "
                        << num << " -- " << i
                        << std::endl
                    ;
                }
                break;
            }
        }
        all_ok = all_ok and done;
    }
    ok(all_ok, "from_range sampling");

    const int repeats = 100;
    int count_ok = 0;
    const float from_range_prob = -logf(float(total));
    for (int i=0; i < repeats; ++i) {
        float prob = 0.f;
        sampling::from_range(poisson_samples.begin(), poisson_samples.end(), prob);
        count_ok += int(prob == from_range_prob);
    }
    ok(count_ok == repeats, "from_range_prob");

    std::vector<int> unif_ints(repeats);

    const int i_upper = 5;
    for (int i=0; i < repeats; ++i) {
        unif_ints[i] = sampling::uniform_int(0, i_upper);
    }
    bool unif_sampling_ok = true;
    for (int i=0; i < i_upper; ++i) {
        unif_sampling_ok = unif_sampling_ok and (std::count(unif_ints.begin(), unif_ints.end(), i) > 0);
    }
    unif_sampling_ok = unif_sampling_ok and (std::count(unif_ints.begin(), unif_ints.end(), i_upper) == 0);
    ok(unif_sampling_ok, "uniform integer sampling");

    float count_0 = 0;
    float count_1 = 0;
    //const int total = 100000;
    for (int i = 0; i < total; ++i) {
        int isam = sampling::uniform_int(0, 2);
        count_1 += isam == 1 ? 1.f : 0.f;
        count_0 += isam == 0 ? 1.f : 0.f;
    }
    count_0 /= float(total);
    count_1 /= float(total);
    ok1(count_0 > 0.48 and count_1 > 0.48 and count_0 + count_1 == 1.f);

    const float f_upper = 5.0f;
    unif_sampling_ok = true;
    data_container unif_floats(repeats);
    for (int i=0; i < repeats; ++i) {
        unif_floats[i] = sampling::uniform_real(0, f_upper);
    }
    unif_sampling_ok = unif_sampling_ok and *std::min_element(unif_floats.begin(), unif_floats.end()) < 1.0f;
    unif_sampling_ok = unif_sampling_ok and *std::min_element(unif_floats.begin(), unif_floats.end()) >= 0.0f;
    unif_sampling_ok = unif_sampling_ok and *std::max_element(unif_floats.begin(), unif_floats.end()) < f_upper;
    ok(unif_sampling_ok, "uniform float sampling");

    observation centre(new_obs(0, 0, 3));
    observation sample(new_obs(0, 0, 0));
    //float prob;

    const float eps = std::numeric_limits<float>::epsilon();

    const float speed_of_light(1.0f + eps);

    std::set<time_stamp> time_stamp_with_obs;
    observation_collection candidates;
    observations_within_light_cone(cl_obs, centre, speed_of_light,
                               1, 6,
                               time_stamp_with_obs,
                               std::inserter(candidates, candidates.begin()));
    if (verbose)
        std::cout << "candidates size = " << candidates.size() << std::endl;
    ok(candidates.size() == 4, "observations_within_light_cone");

    //int num_sampling_success = 0;
    /*
    observation start(new_obs(0, 0, 0));
    for (int i = 0 ; i < total; ++i) {
        sample = new_obs(-1, -1, -1);
        prob = 0.0f;
        bool success = sampling::observation_from_light_cone(collection2, start, speed_of_light, 1, 4, sample, prob);
        if (not (t(sample) == 1 and prob == logf(0.5f)) and not (t(sample) == 3 and prob == logf(0.1f))) {
            success = false;
        }
        num_sampling_success += success ? 1 : 0;
    }
    ok(num_sampling_success == total, "observation_from_light_cone successful");
    */

    /*
    num_sampling_success = 0;
    for (int i = 0 ; i < total; ++i) {
        sample = new_obs(-1, -1, -1);
        prob = 0.0f;
        bool success = sampling::observation_from_collection(collection2, sample, prob);
        if (not (t(sample) == 1 and prob == logf(0.5f)) and not (t(sample) == 3 and prob == logf(0.1f))) {
            //std::cout << t(sample) << "; p=" << expf(prob) << std::endl;
            success = false;
        }
        num_sampling_success += success ? 1 : 0;
    }
    ok(num_sampling_success == total, "observation_from_collection");
    */
    ok(true, "observation_from_collection not longer in use");

    test_normal(total);

    /*
    boost::timer timer;
    int nloops = 10000000;
    test_beta(nloops);
    diag("time 1st beta test: %f", timer.elapsed());
    timer.restart();
    test_beta2(nloops);
    diag("time 2nd beta test: %f", timer.elapsed());
    */
    test_beta3(1000);
    test_gamma(1000);
    test_gamma0(1000);
    test_gamma3(1000);
    test_gamma1();
    //test_trunc_beta(10000);
    test_factorial();
    //test_uniform();
    //test_uniform_inv();
    //test_gamma2();
    //test_clutter_rate_sampling(1000000);
    diag("random seed = 0x%x", seed);
    return exit_status();
}
