#include <iostream>
#include <functional>
#include <sstream>
#include <stdlib.h>

#include <boost/random.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

#include "biggles/sampling/gibbs.hpp"

extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles::sampling;

class uniform_sampler : public gibbs::mem_fun_two_var_cond_sampler<uniform_sampler, int, float>
{
public:
    uniform_sampler()
        : uni_int_(1,6)
        , uni_real_(10.f,20.f)
    { }

    int sample_variable_1(float)
    {
        return uni_int_(rng_);
    }

    float sample_variable_2(int)
    {
        return uni_real_(rng_);
    }

protected:
    boost::mt19937 rng_;
    boost::uniform_int<> uni_int_;
    boost::uniform_real<> uni_real_;
};

int main(int argc, char** argv)
{
    plan_tests(4*20);

    // uniform independent sampler
    uniform_sampler es;
    gibbs_sampler<uniform_sampler> s(es);
    for(int i=0; i<20; ++i)
    {
        std::stringstream ss;
        ss << "drawn sample: " << s.draw();
        diag("%s", ss.str().c_str());

        ok(boost::get<0>(s.last_sample()) <= 6, "int sample is not too big");
        ok(boost::get<0>(s.last_sample()) >= 1, "int sample is not too small");

        ok(boost::get<1>(s.last_sample()) <= 20.f, "float sample is not too big");
        ok(boost::get<1>(s.last_sample()) >= 10.f, "float sample is not too small");
    }

    return exit_status();
}
