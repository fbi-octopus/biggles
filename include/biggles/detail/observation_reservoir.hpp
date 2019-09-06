#ifndef BIGGLES_OBSERVATION_RESERVOIR_HPP__
#define BIGGLES_OBSERVATION_RESERVOIR_HPP__

#include <deque>
#include <map>
#include "../observation.hpp"

namespace biggles {

/** \brief Holds observations
 *
 * This class holds observations in a specific order and allows to identify them via an index.
 * It is intended that the order is preserved
 *
 */
class ObservationReservoir {
    typedef std::deque<observation> container_t;
    typedef std::map<observation, size_t> index_t;
    container_t container_;
    index_t index_;
public:
    typedef container_t::const_iterator const_iterator;
    typedef container_t::const_iterator iterator; ///< \brief Cheating boost::python
    ObservationReservoir() {}
    template <class OBS_ITER>
    ObservationReservoir(OBS_ITER obegin, OBS_ITER oend) { while (obegin != oend) push(*obegin++); }
    //ObservationReservoir(const ObservationReservoir& oc) : container_(oc.container_), index_(oc.index_) {}
    size_t size() const {
        return container_.size();
    }
    void push(const observation& o) {
        std::pair<index_t::iterator, bool> result = index_.insert(std::make_pair(o, size()));
        if (result.second) { container_.push_back(o); }
    }
    void push(const float x, const float y, const time_stamp t) {
        observation o = new_obs(x, y, t);
        push(o);
    }
    /// \brief the index of an observation
    size_t idx(const observation& o) const { return index_.at(o); }
    /// \brief writes the indices of the observations given by \e obegin and \e oend to \e ibegin
    template <class OBS_ITER, class IDX_ITER>
    IDX_ITER idx(OBS_ITER obegin, OBS_ITER oend, IDX_ITER ibegin) const {
        while (obegin != oend) *ibegin++ = index_.at(*obegin++);
        return ibegin;
    }
    const observation& operator[](const size_t i) const { return container_[i]; }
    template <class OBS_ITER, class IDX_ITER>
    OBS_ITER obs(IDX_ITER ibegin, IDX_ITER iend, OBS_ITER obegin) const {
        while (ibegin != iend) *obegin++ = container_[*ibegin++];
        return obegin;
    }
    const_iterator begin() const { return container_.begin(); }
    const_iterator end()   const { return container_.end(); }

}; // ObservationReservoir

typedef boost::shared_ptr<const ObservationReservoir> observation_reservoir_ptr;

} // namespace biggles
#endif // BIGGLES_OBSERVATION_RESERVOIR_HPP__
