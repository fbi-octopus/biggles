#ifndef BIGGLES_MH_MOVE_WEIGHTS_HPP___
#define BIGGLES_MH_MOVE_WEIGHTS_HPP___

#include "../track_collection.hpp"
#include "../detail/fun.hpp"
#include "../detail/constants.hpp"
#include "../sampling/simple.hpp"

namespace biggles { namespace mh_moves {

typedef std::unary_function<shared_const_track_ptr, float> track_weight_function;
typedef std::unary_function<observation, float> obs_weight_function;
typedef std::deque<float> weights_t;

/// \brief the weight to select the *first* track in the merge move
template<int max_weight, int max_ratio_numer, int max_ratio_denom>
class front_weight_t : public track_weight_function {
    time_stamp min_weight_ts;
    time_stamp last_max_weight_ts;
    float m_;
    float n_;
public:
    front_weight_t(const time_stamp first_ts, const time_stamp last_ts) :
        min_weight_ts(last_ts - 2),
        last_max_weight_ts(time_stamp(float(last_ts - first_ts - 2) * float(max_ratio_numer)/float(max_ratio_denom))
            + first_ts + 2),
        m_(float(max_weight - 1.f) / float(last_max_weight_ts - min_weight_ts)),
        n_(1.f - min_weight_ts * m_)
    {}

    float operator ()(const shared_const_track_ptr& t_ptr) const {
        time_stamp last_ts = t_ptr->last_time_stamp();
        if (last_ts <= last_max_weight_ts)
            return float(max_weight);
        if (last_ts > min_weight_ts)
            return 0.f;
        return m_* last_ts + n_;
    }
};

/// \brief the weight to select the *second* track in the merge move
template<int offset>
class back_weight_t : public track_weight_function {
    shared_const_track_ptr front_piece_;
public:
    back_weight_t(const shared_const_track_ptr& front_piece) : front_piece_(front_piece) {}
    float operator ()(const shared_const_track_ptr& t_ptr) const {
        float track_dist = tracks_could_be_merged2(*front_piece_, *t_ptr);
        if (track_dist<0.f)
            return 0.f;
        return 1.f/square(track_dist + static_cast<float>(offset));
    }
};

/// \brief the weight to select the track for the extend move
class extend_weight_t : public track_weight_function {
    float max_duration_;
public:
    explicit extend_weight_t(size_t max_duration) : max_duration_(float(max_duration)) {}
    float operator ()(const shared_const_track_ptr& t_ptr) const {
        return (max_duration_ - float(t_ptr->duration()) > 0)
            ? 3.f/(max_duration_ - 2.f)*(float(max_duration_) - float(t_ptr->duration()))
            : 0.f;
    }
};

template<int sigma_numer, int sigma_denom>
class gaussian_weight_t : public obs_weight_function { // this is still not right
    float x_;
    float y_;
    float t_;
    float sigma_;
public:
    explicit gaussian_weight_t(const observation& mean) : x_(x(mean)), y_(y(mean)), t_(t(mean)),
        sigma_(float(sigma_numer)/float(sigma_denom)) {}
    float operator ()(const observation& obs) const {
        const float sig = std::abs(t_ - float(t(obs))) * sigma_;
        return std::exp(-(square(x_ - x(obs)) + square(y_ - y(obs)))/(2.f*sig*sig)) * detail::gauss_factor/sig;
    }
    float sigmas(size_t count, time_stamp ts) const {
        const float sig = std::abs(t_ - float(ts))*sigma_;
        return std::exp(-float(square(count))*0.5f) * detail::gauss_factor/sig;
    }
};

/// \brief the weight to select the track for the reduce move
class reduce_weight_t : public track_weight_function {
public:
    float operator ()(const shared_const_track_ptr& t_ptr) const {
        return t_ptr->duration() > 2 ? 1.f : 0.f;
    }
};

/// \brief the weight to select the track for the split move
class split_weight_t : public track_weight_function {
public:
    float operator ()(const shared_const_track_ptr& t_ptr) const {
        //return t_ptr->size() > 3 ? 1.f : 0.f;
        return t_ptr->size() > 3 ? std::exp(float(t_ptr->duration())/float(t_ptr->size())-1.f) : 0.f;
    }
};

template <class FUN>
float get_track_log_prob(const FUN& fun, const track_collection& tracks, const shared_const_track_ptr& track_ptr) {
    float target_weight = fun(track_ptr);
    return logf(target_weight/fun_sum(tracks.begin(), tracks.end(), fun, 0.f));
}

template <class FUN>
bool select_track(const FUN& fun, const track_collection& tracks, shared_const_track_ptr& track_ptr, float& log_prob) {
    track_collection::const_iterator it(sampling::weighted_choice(tracks.begin(), tracks.end(), fun, log_prob));
    if (it == tracks.end())
        return false;
    track_ptr = *it;
    return true;
}

typedef front_weight_t<10, 1, 3> front_weight_merging;
typedef back_weight_t<0> back_weight_merging;
typedef gaussian_weight_t<1, 1> std_normal_weight;

} } // \brief biggles::mh_moves

#endif // BIGGLES_MH_MOVE_WEIGHTS_HPP___
