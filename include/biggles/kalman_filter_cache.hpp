/// @file kalman_filter_cache.hpp LRU cache of Kalman filter parameters for tracks

#ifndef BIGGLES_KALMAN_FILTER_CACHE_HPP__
#define BIGGLES_KALMAN_FILTER_CACHE_HPP__

#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>
#include <boost/weak_ptr.hpp>
#include <map>

#include "kalman_filter.hpp"

namespace biggles
{

/// @brief A periodically purged cache of Kalman filters associated with a track.
///
/// This class maintains an internal cache mapping weak track pointers to instanced of biggles::kalman_filter and the
/// values of \f$ R \f$ which those Kalman filters were initialised with. Periodically (currently every 1024 insertions
/// into the cache) entries associated with tracks no-longer in use are purged.
///
class kalman_filter_cache : boost::noncopyable
{
public:
    kalman_filter_cache() : purge_timer_(0) { }

    /// Return a biggles::kalman_filter instance for the specified track and covariance matrix.
    ///
    /// @param track_ptr
    /// @param R11
    /// @param R12
    /// @param R22
    const kalman_filter& get(boost::shared_ptr<const track> track_ptr, const matrix2f &R, const matrix4f &Q);

protected:
    /// @brief The type of an entry in the cache.
    typedef boost::tuples::tuple<kalman_filter, float, float, float> cache_entry;

    /// @brief The type of the cache itself.
    typedef std::map<boost::weak_ptr<const track>, cache_entry> cache_map;

    /// @brief The cache.
    cache_map cache_;

    /// @brief Incremented on every call to make_entry().
    size_t purge_timer_;

    /// @brief Force an entry to be updated.
    ///
    /// @param track_ptr
    /// @param R11
    /// @param R12
    /// @param R22
    const kalman_filter& make_entry(boost::weak_ptr<const track> track_ptr, const matrix2f &R, const matrix4f &Q);
};

}

#endif // BIGGLES_KALMAN_FILTER_CACHE_HPP__
