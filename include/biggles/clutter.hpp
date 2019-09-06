#ifndef BIGGLES_CLUTTER_HPP__
#define BIGGLES_CLUTTER_HPP__

#include "detail/constants.hpp"
#include "observation.hpp"
#include "point_2d.hpp"
#include <boost/shared_ptr.hpp>
#include <algorithm>
#include <map>
#include <numeric>
#include <deque>


namespace biggles {

class clutter_t {
public:
    typedef std::deque<observation> observation_container;
    typedef std::map<time_stamp, observation_container>::value_type value_type;
private:
    std::map<time_stamp, observation_container> clutter_map_;
    observation_container emptiness_;
    struct adder {
        size_t operator()(const size_t& a, const value_type& b) {
            return a + b.second.size();
        }
    };
    struct upper_bound_ {
        point_2d operator()(const point_2d& a, const value_type& b) {
            float xhi = -detail::myriad;
            float yhi = -detail::myriad;
            for (observation_container::const_iterator it = b.second.begin(); it != b.second.end(); ++it) {
                xhi = std::max(xhi, x(*it));
                yhi = std::max(yhi, y(*it));
            }
            return point_2d(std::max(x(a), xhi), std::max(y(a), yhi));
        }
    };
    struct lower_bound_ {
        point_2d operator()(const point_2d& a, const value_type& b) {
            float xlo = detail::myriad;
            float ylo = detail::myriad;
            for (observation_container::const_iterator it = b.second.begin(); it != b.second.end(); ++it) {
                xlo = std::min(xlo, x(*it));
                ylo = std::min(ylo, y(*it));
            }
            return point_2d(std::min(x(a), xlo), std::min(y(a), ylo));
        }
    };
public:
    typedef std::map<time_stamp, observation_container>::iterator iterator;
    typedef std::map<time_stamp, observation_container>::const_iterator const_iterator;
    clutter_t() {}
    clutter_t(const clutter_t& c) : clutter_map_(c.clutter_map_) {}
    template<class ITER> clutter_t(ITER b, ITER e) { insert(b, e); }
    /// \brief inserts a single observation clutter
    void insert(const observation& o) {
        iterator it = clutter_map_.find(t(o));
        if (it == end()) {
            std::pair<iterator, bool> iter_bool = clutter_map_.insert(std::make_pair(t(o), observation_container()));
            it = iter_bool.first;
        }
        it->second.push_back(o);
    }
    /// \brief inserts observations from range into clutter
    template<class ITER> void insert(ITER b, ITER e) { for( ; b != e; ++b) insert(*b); }
    /// \brief remove an observation from clutter
    void erase(const observation& o) {
        iterator it = clutter_map_.find(t(o));
        if (it == end()) return;
        observation_container& obs_cont = it->second;
        observation_container::iterator to_remove = std::find(obs_cont.begin(), obs_cont.end(), o);
        if (to_remove != obs_cont.end()) obs_cont.erase(to_remove);
    }
    /// \brief the time stamp of the earliest observation
    time_stamp first_time_stamp() const {
        time_stamp ts = 0;
        for (const_iterator it = clutter_map_.begin(); it != clutter_map_.end(); ++it) {
            if (it->second.size() > 0) {
                ts = it->first;
                break;
            }
        }
        return ts;
    }
    /// \brief the time stamp following the latest observation
    time_stamp last_time_stamp() const {
        time_stamp ts = 0;
        std::map<time_stamp, observation_container>::const_reverse_iterator rit;
        for (rit = clutter_map_.rbegin(); rit != clutter_map_.rend(); ++rit) {
            if (rit->second.size() > 0) {
                ts = rit->first + 1;
                break;
            }
        }
        return ts;
    }
    /// \brief the difference between last and first time stamp
    time_stamp duration() const { return last_time_stamp() - first_time_stamp(); }
    /// \brief the number of observations in clutter
    size_t size() const {
        return std::accumulate(begin(), end(), size_t(0), adder());
    }
    /// \brief is the clutter empty?
    bool empty() const { return size() == size_t(0); }
    /// \brief number of observations at time stamp
    size_t count_at_time_stamp(const time_stamp& ts) const {
        const_iterator it = clutter_map_.find(ts);
        if (it == clutter_map_.end())
            return size_t(0);
        return it->second.size();
    }
    /// \brief iterator to the first <time_stamp, observation_container> pair
    inline const_iterator begin() const { return clutter_map_.begin(); }
    /// \brief end iterator of the <time_stamp, observation_container> pairs
    inline const_iterator   end() const { return clutter_map_.end(); }
    /// \brief the global (left, bottom) location boundary of the clutter
    point_2d lower_bound() const {
        return std::accumulate(clutter_map_.begin(), clutter_map_.end(), point_2d(detail::myriad, detail::myriad),
            lower_bound_());
    }
    /// \brief the global (right, top) location boundary of the clutter
    point_2d upper_bound() const {
        return std::accumulate(clutter_map_.begin(), clutter_map_.end(), point_2d(-detail::myriad, -detail::myriad),
            upper_bound_());
    }
    /// \brief locate near proxy
    std::deque<observation> locate_near(time_stamp ts, const observation& obs, const float radius) const {
        const_iterator it = clutter_map_.find(ts);
        std::deque<observation> result;
        if (it == clutter_map_.end())
            return result;
        const std::deque<observation>& obs_cont = it->second;
        for (std::deque<observation>::const_iterator iter = obs_cont.begin(); iter != obs_cont.end(); ++iter) {
            if (dist(obs, *iter)<=radius) result.push_back(*iter);
        }
        return result;
    }
    /// \brief observations at time point
    const observation_container& observations(time_stamp ts) const {
        const_iterator it = clutter_map_.find(ts);
        if (it == clutter_map_.end())
            return emptiness_;
        return it->second;
    }
    /// \brief all observations
    observation_container observations() const {
        observation_container result;
        for (const_iterator it = begin(); it not_eq end(); ++it) {
            result.insert(result.end(), it->second.begin(), it->second.end());
        }
        return result;
    }
    /// \brief number of observations at time point
    size_t size(time_stamp ts) const { return count_at_time_stamp(ts); }
    /// \brief
    bool contains(const observation& o) const {
        const_iterator it = clutter_map_.find(t(o));
        if (it == clutter_map_.end())
            return false;
        return std::find(it->second.begin(), it->second.end(), o) != it->second.end();
    }

};

template <class OBS_CONTAINER>
void ingest_observations(OBS_CONTAINER& container, const clutter_t& clutter) {
    for (clutter_t::const_iterator it = clutter.begin(); it not_eq clutter.end(); ++it) {
        container.insert(container.end(), it->second.begin(), it->second.end());
    }
}


} // biggles

#endif // BIGGLES_CLUTTER_HPP__
