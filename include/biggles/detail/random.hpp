#include <boost/random.hpp>

namespace biggles { namespace detail {

/// @brief The underlying boost random number generator used by Biggles' sampling functions.
extern boost::mt19937 random_generator;

unsigned int get_random_seed();

/// @brief Reset the seed of the PRNG used by Biggles.
///
/// On startup the PRNG function used by Biggles is seeded with a value based on the system clock. This function allows
/// the seed to be set explicitly.
///
/// @param seed The new seed for the PRNG.
void seed_prng(unsigned int seed);

} } // biggles::detail
