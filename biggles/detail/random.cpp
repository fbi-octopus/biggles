#include <ctime>

#include "detail/random.hpp"

namespace biggles { namespace detail {

const unsigned int random_seed(std::time(0));
//const unsigned int random_seed(0x56cd94d3);

boost::mt19937 random_generator(random_seed);

unsigned int get_random_seed() {
    return random_seed;
}

void seed_prng(unsigned int seed)
{
    random_generator.seed(seed);
}

} } // biggles::detail
