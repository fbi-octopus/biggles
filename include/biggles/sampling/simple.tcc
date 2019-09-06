#ifndef WITHIN_BIGGLES_SAMPLING_SIMPLE_HPP__
#error "this file must only be included from samplers.hpp"
#endif // WITHIN_BIGGLES_SAMPLING_SIMPLE_HPP__

#include <cmath>
#include <limits>

namespace biggles { namespace sampling
{

template<typename InputIterator>
InputIterator from_range(InputIterator first, InputIterator last)
{
    float ignored(0.f);
    return from_range(first, last, ignored);
}

template<typename InputIterator>
InputIterator from_range(InputIterator first, InputIterator last, float& log_prob) {
    int len(std::distance(first, last));
    log_prob = -logf(static_cast<float>(len));
    InputIterator result = first;
    std::advance(result, uniform_int(0, len));
    return result;
}

/*
{
    InputIterator current_sample(first);

    // this implements the following algorithm: suppose we already have a sample, s, from a range of length N. If we
    // sample from a range of length N+1 which is the original range plus a new element, the probability that the
    // new element is in this new sample is 1/(N+1) and the probability we retain the original sample is N/(N+1).

    size_t n_in_previous_range = 0;
    for(; first != last; ++first, ++n_in_previous_range)
    {
        // accept this new sample with prob. 1 / (n_in_previous_range + 1)
        if(0 == uniform_int(0,n_in_previous_range+1))
            current_sample = first;
    }

    if(n_in_previous_range > 0)
    {
        log_prob = -logf(static_cast<float>(n_in_previous_range));
    }
    else
    {
        log_prob = -std::numeric_limits<float>::max();
    }

    return current_sample;
}
*/

} }
