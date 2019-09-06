#include <functional>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string.h>

#include <boost/random.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

#include "biggles/detail/random.hpp"
#include "biggles/sampling/metropolis_hastings.hpp"

extern "C" {
#include <ccan/tap/tap.h>
}

typedef boost::tuples::tuple<float, float> point;

inline float rosenbrock(const point& p)
{
    float x, y;
    boost::tuples::tie(x,y) = p;
    return (1.f-x)*(1.f-x) + 100.f*(y-x*x)*(y-x*x);
}

// a cooked-up PDF which has a maximum where the Rosenbrock function has a minimum.
inline float rosenbrock_log_pdf(const point& p)
{
    return -rosenbrock(p);
}

typedef boost::tuples::tuple<point, float> proposal_result_t;

inline proposal_result_t propose_gaussian(const point& p)
{
    static boost::normal_distribution<float> norm(0.f, 0.2f);
    static boost::variate_generator<boost::mt19937&, boost::normal_distribution<float> > variate(
        biggles::detail::random_generator, norm);

    point new_p(p.get<0>() + variate(), p.get<1>() + variate());

    return proposal_result_t(new_p, 0.f);
}

using namespace biggles::sampling::metropolis_hastings;

int main(int argc, char** argv)
{
    biggles::detail::seed_prng(0xdeadbeef);

    bool verbose = (argc > 1) && (0 == strcmp(argv[1], "-v"));

    plan_tests(4);

    typedef std::pointer_to_unary_function<const point&, float> pdf_t;
    typedef std::pointer_to_unary_function<const point&, proposal_result_t> proposal_func_t;

    // MH sampler for Rosenbrock function. An example of using a function pointer directly.
    sampler<pdf_t, proposal_func_t> s(std::ptr_fun(rosenbrock_log_pdf), std::ptr_fun(propose_gaussian));

    point optimum(s.last_sample());
    float optimum_log_density(s.current_sample_log_density());

    for(int i=0; i<4096; ++i)
    {
        s.draw();

        if(s.current_sample_log_density() > optimum_log_density)
        {
            optimum_log_density = s.current_sample_log_density();
            optimum = s.last_sample();
        }

        if(verbose)
        {
            std::stringstream ss;
            ss << "sample: " << s.last_sample();
            diag("%s", ss.str().c_str());
        }
    }

    {
        std::stringstream ss;
        ss << "optimum: " << optimum << ", log-pdf: " << optimum_log_density
        << ", alpha: " << s.acceptance_rate() << ", rosenbrock: " << rosenbrock(optimum);
        diag("%s", ss.str().c_str());
    }

    diag("acceptance rate = %.2f%%", s.acceptance_rate()*100.f);
    //ok(fabs(s.acceptance_rate() - 0.25f) < 0.2f, "acceptance rate within 0.2 of 0.25");
    ok(fabs(optimum.get<0>() - 1.f) < 0.1f, "x-co-ordinate within 10%% of optimum");
    ok(fabs(optimum.get<1>() - 1.f) < 0.1f, "y-co-ordinate within 10%% of optimum");

    // the same thing using the optimise function
    point optimum2 = optimise(std::ptr_fun(rosenbrock_log_pdf), std::ptr_fun(propose_gaussian));

    {
        std::stringstream ss;
        ss << "optimum via optimise(): " << optimum2 << ", log-pdf: " << rosenbrock_log_pdf(optimum2)
        << ", rosenbrock: " << rosenbrock(optimum2);
        diag("%s", ss.str().c_str());
    }

    ok(fabs(optimum2.get<0>() - 1.f) < 0.1f, "x-co-ordinate within 10%% of optimum");
    ok(fabs(optimum2.get<1>() - 1.f) < 0.1f, "y-co-ordinate within 10%% of optimum");

    return exit_status();
}
