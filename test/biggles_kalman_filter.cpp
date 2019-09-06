#include "biggles/observation.hpp"
#include "biggles/detail/random.hpp"
#include "biggles/mh_moves/mh_moves.hpp"
#include "biggles/mh_moves/utility.hpp"
#include "biggles/track.hpp"
#include "biggles/tracker.hpp"
#include "biggles/kalman_filter.hpp"
#include "biggles/simulate.hpp"
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <string>
#include "biggles/tools/sundries.hpp"

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

typedef boost::shared_ptr< const kalman_filter > shared_const_kalman_filter_ptr;
//typedef std::pair<float, float> data_pair;

std::ostream& operator<<(std::ostream& os, const kalman_filter::state_vector& state) {
    os << boost::format("[% .4f, % .4f, % .4f, % .4f]") % state[0] % state[1] % state[2] % state[3];
    return os;
}

float max_deviation_from_line( const observation_collection &oc) {
    time_stamp start = oc.first_time_stamp();
    time_stamp end = oc.last_time_stamp();
    float x_off = x(*oc.begin());
    float y_off = y(*oc.begin());
    float delta_x = (x(oc.back()) - x_off)/float(end-start-1);
    float delta_y = (y(oc.back()) - y_off)/float(end-start-1);
    float result = 0.f;
    for (time_stamp ts = start; ts != end; ++ts) {
        observation_collection::const_range obs_at_t = oc.at_time_stamp(ts);
        float line_x = x_off + (ts - start) * delta_x;
        float line_y = y_off + (ts - start) * delta_y;
        for (; obs_at_t.first != obs_at_t.second; ++(obs_at_t.first)) {
            float diff = hypotf(line_x - x(*obs_at_t.first), line_y - y(*obs_at_t.first));
            if (diff > result) result = diff;
        }
    }
    return result;
}

shared_const_track_ptr generate_rotating_track(int duration) {
    const float delta_angle = 2.0f*3.1415f/float(duration);
    observation_collection oc;
    for (time_stamp t = 0; t < duration; ++t) {
        float angle = t*delta_angle;
        oc.insert(new_obs(std::cos(angle), std::sin(angle), t));
    }
    return shared_const_track_ptr(new track(0, duration, oc.begin(), oc.end(), 1.0f));
}

shared_const_track_ptr generate_linear_track(int duration, float dx, float dy) {
    observation_collection oc;
    for (time_stamp t = 0; t < duration; ++t) oc.insert(new_obs( t*dx, t*dy, t));
    return shared_const_track_ptr(new track(0, duration, oc.begin(), oc.end(), 1.0f));
}

shared_const_track_ptr generate_linear_track(
    time_stamp begin_ts, time_stamp end_ts, float x0, float y0, float dx, float dy, float loc_err)
{
    observation_collection oc;
    for (time_stamp ts = begin_ts; ts < end_ts; ++ts) oc.insert(new_obs(
        x0 + (ts - begin_ts)*dx + sampling::sample_normal() * loc_err,
        y0 + (ts - begin_ts)*dy + sampling::sample_normal() * loc_err,
        ts));
    return shared_const_track_ptr(new track(begin_ts, end_ts, oc.begin(), oc.end(), 1.0f));
}



void print_kalman_filter_states(const kalman_filter::states_and_cov_deque& sc) {
    kalman_filter::states_and_cov_deque::const_iterator sc_iter;
    int i = 0;
    for (sc_iter = sc.begin(); sc_iter != sc.end(); ++sc_iter) {
        std::stringstream sstr;
        sstr
            << boost::format("%3d: ") % i++
            << sc_iter->first
            << ", "
            << boost::format("v = % .4f") % hypotf(sc_iter->first[1], sc_iter->first[3])
            << ", "
            << boost::format("log(||C||) = % .4f") % logf(sc_iter->second.determinant());
        diag("%s", sstr.str().c_str());
    }
}

std::string print_parameters(const model::parameters& params) {
    std::stringstream para_ss;
    para_ss << boost::format("b = %.2f, c = %.2f, o = %.2f%%, s = %.2f%%, cr = %.2f")
            % model::mean_new_tracks_per_frame(params)
            % model::mean_false_observations_per_frame(params)
            % (model::generate_observation_probability(params) * 100.f)
            % (model::frame_to_frame_survival_probability(params) * 100.f)
            % model::constraint_radius(params)
            << ", "
            << boost::format("R = ((%.2f, %.2f), (%.2f, %.2f))")
            % model::observation_error_covariance(params)(0, 0)
            % model::observation_error_covariance(params)(0, 1)
            % model::observation_error_covariance(params)(1, 0)
            % model::observation_error_covariance(params)(1, 1)
    ;
    return para_ss.str();
}

/*
data_pair generate_partition(const model::parameters& params, partition& dataset_partition, bool shutup) {

    int image_size = sampling::uniform_int(14, 17);
    int time_extent = sampling::uniform_int(20, 30);
    simulate::generate_partition(dataset_partition, params, 0, time_extent, image_size, image_size);

    if (not shutup) {
        diag("Generated %zu new tracks.", dataset_partition.tracks().size());
        diag("Generated %zu clutter observations.", dataset_partition.clutter().size());
    }

    return std::make_pair(image_size, time_extent);
}
*/

std::string partition_properties(const partition &const_partition) {
    const track_collection &tracks(const_partition.tracks());
    size_t num_tracks = tracks.size();
    track_collection::const_iterator tr_iter;
    float avg_track_duration = 0.0f;
    float avg_track_size = 0.0f;
    for (tr_iter = tracks.begin(); tr_iter != tracks.end(); ++tr_iter) {
        avg_track_duration += float((*tr_iter)->duration());
        avg_track_size += float((*tr_iter)->size());
    }
    avg_track_duration /= float(num_tracks);
    avg_track_size /= float(num_tracks);
    float o_obs = avg_track_size/avg_track_duration;
    float b_obs = float(num_tracks)/float(const_partition.duration());
    float c_obs = float(const_partition.clutter().size())/float(const_partition.duration());
    float s_obs = 1.f - 1.f/avg_track_duration;
    std::stringstream pp_ss;
    pp_ss << boost::format("num.tracks = %d, avg.length = %.1f, avg.obs.count = %.1f, "
        "observed: b = %.2f, c = %.2f, o = %.2f, s = %.2f")
        % num_tracks % avg_track_duration % avg_track_size % b_obs % c_obs % o_obs % s_obs
    ;
    return pp_ss.str();
}

int plot_tracks(const track_collection &tracks, const model::parameters &const_parameters) {
    track_collection::const_iterator tr_iter;
    for (tr_iter = tracks.begin(); tr_iter != tracks.end(); ++tr_iter) {
        const track & track0(**tr_iter);
        shared_const_kalman_filter_ptr kf = track0.make_kalman_filter(const_parameters);
        float log_posterior = track0.log_posterior(const_parameters);
        const kalman_filter::states_and_cov_deque & sncs(kf->corrections());
        const observation_collection &oc(track0.observations());
        std::deque <float> diffs;
        time_stamp ts=track0.first_time_stamp();
        std::stringstream diff_ss;
        std::stringstream cov_ss;
        std::stringstream track_ss;
        diff_ss <<  "| ";
        track_ss << "| ";
        track_ss << std::string(ts, ' ');
        diff_ss << std::string(ts, ' ');
        cov_ss   << "log(cov) = ";
        for (size_t i=0; i < sncs.size(); ++i, ++ts) {
            kalman_filter::state_vector state = sncs[i].first;
            kalman_filter::covariance_matrix cov = sncs[i].second;
            observation_collection::const_range obs_at_t(oc.at_time_stamp(ts));
            cov_ss << boost::format(" %.1f") % log10f(cov.determinant());
            if (obs_at_t.first != obs_at_t.second) {
                observation obs = *obs_at_t.first;
                float diff = hypotf(x(obs)-state[0], y(obs)-state[2]);
                diff_ss << boost::format("%d") % int(diff);
                diffs.push_back(diff);
                track_ss << "O";
            }
            else {
               diff_ss << "_";
               track_ss << "_";
            }
        }
        diag("log(PDF) = %8.2f %s", log_posterior, track_ss.str().c_str());
        diag("                    %s", diff_ss.str().c_str());
        diag("%s", cov_ss.str().c_str());
        float max_diff = *std::max_element(diffs.begin(), diffs.end());
        float radius = track0.radius();
        float max_dev = max_deviation_from_line(oc);
        if (max_diff >= max_dev or max_dev >= radius) {
            //diag("max diff = %.2f, dev from line = %.2f, radius = %.2f", max_diff, max_dev, radius);
            //diag("%s", diff_ss.str().c_str());
            //diag("length = %ld, observations = %.2f%%", track0.duration(), 100.f*float(track0.size())/float(track0.duration()));
        }
    }
    return 0;
}

int plot_states_of_tracks(const track_collection &tracks, const model::parameters &const_parameters) {
    track_collection::const_iterator tr_iter;
    for (tr_iter = tracks.begin(); tr_iter != tracks.end(); ++tr_iter) {
        const track& track0(**tr_iter);
        shared_const_kalman_filter_ptr kf = track0.make_kalman_filter(const_parameters);
        diag("track [%ld, %ld]", track0.first_time_stamp(), track0.last_time_stamp());
        print_kalman_filter_states(kf->corrections());
    }
    return 0;
}

int test2() {
    observation_collection t1_obs;
    //observation_collection cl_obs;
    t1_obs.insert(new_obs(0, 0, 0));
    t1_obs.insert(new_obs(0, 0, 1));
    t1_obs.insert(new_obs(0, 0, 2));
    t1_obs.insert(new_obs(0, 0, 3));
    track track1(0, 5, t1_obs.begin(), t1_obs.end(), 1.0);

    /*
    cl_obs.insert(new_obs(0, 0, 4));
    cl_obs.insert(new_obs(1, 1, 5));
    cl_obs.insert(new_obs(1, 1, 6));
    cl_obs.insert(new_obs(1, 1, 7));
    cl_obs.insert(new_obs(1, 1, 8));
    cl_obs.insert(new_obs(1, 1, 9));
    boost::shared_ptr<track_collection> tracks(new track_collection());
    tracks->insert(track1);
    boost::shared_ptr<const observation_collection> clutter(new observation_collection(cl_obs.begin(), cl_obs.end()));
    const partition original_partition(tracks, clutter);
    */

    model::parameters params;
    model::mean_new_tracks_per_frame(params)           = 2.f;
    model::mean_false_observations_per_frame(params)   = 1.f;
    model::frame_to_frame_survival_probability(params) = 0.95f;
    model::generate_observation_probability(params)    = 0.9f;
    model::constraint_radius(params) = 100.0f;
    model::observation_error_covariance(params) << 0.001f, 0.0f, 0.0f, 0.001f;

    print_kalman_filter_states( track1.make_kalman_filter(params)->corrections() );
    diag("forward:");
    print_kalman_filter_states( track1.make_kalman_filter(params)->predictions() );
    diag("backward");
    print_kalman_filter_states( generate_rotating_track(20)->make_kalman_filter(params)->corrections() );
    diag("forward");
    print_kalman_filter_states( generate_rotating_track(20)->make_kalman_filter(params)->predictions() );
    diag("backward linear:");
    print_kalman_filter_states( generate_linear_track(40, 1.f, 0.f)->make_kalman_filter(params)->corrections() );
    diag("forward linear:");
    print_kalman_filter_states( generate_linear_track(40, 1.f, 0.f)->make_kalman_filter(params)->predictions() );
    return 0;
}



int test3() {
    const float loc_err(.01f);

    std::deque<float> sing_pdf;
    for (int i=0; i<1; ++i) {
        boost::shared_ptr<track_collection> tracks(new track_collection());
        tracks->insert(*generate_linear_track(0, 50, 0.f, 0.f, 1.f, 0.f, loc_err));
        clutter_ptr clutter(new clutter_t());
        const partition original_partition(tracks, clutter);
        model::parameters params;
        float log_pdf = find_good_parameters(original_partition, params, 1000);
        model::observation_error_covariance(params) << loc_err * loc_err, 0, 0, loc_err * loc_err;
        diag("%s, log(PDF)=%.2f", print_parameters(params).c_str(), log_pdf);
        plot_tracks(*tracks, params);
        plot_states_of_tracks(*tracks, params);
        sing_pdf.push_back(log_pdf);
    }
    std::deque<float> quad_pdf;
    for (int i=0; i<1; ++i) {
        boost::shared_ptr<track_collection> tracks(new track_collection());
        int step = 5;
        for (int j=0; j < 50; j+= step)
            tracks->insert(*generate_linear_track(j, j+step, float(j), 0.f, 1.f, 0.f, loc_err));
        clutter_ptr clutter(new clutter_t());
        const partition original_partition(tracks, clutter);
        model::parameters params;
        float log_pdf = find_good_parameters(original_partition, params, 1000);
        model::observation_error_covariance(params) << loc_err * loc_err, 0, 0, loc_err * loc_err;
        diag("%s, log(PDF)=%.2f", print_parameters(params).c_str(), log_pdf);
        plot_tracks(*tracks, params);
        plot_states_of_tracks(*tracks, params);
        quad_pdf.push_back(log_pdf);
    }
    /*
    std::pair<float, float> sing_stat = mean_stdev(sing_pdf.begin(), sing_pdf.end());
    std::pair<float, float> quad_stat = mean_stdev(quad_pdf.begin(), quad_pdf.end());
    diag("sing %.2f +/- %.2f", sing_stat.first, sing_stat.second);
    diag("quad %.2f +/- %.2f", quad_stat.first, quad_stat.second);
    */
    return 0;
}

struct extract_state : public std::unary_function<kalman_filter::state_covariance_pair, state_t> {
    state_t operator()(const kalman_filter::state_covariance_pair &sc_pair) const {
        return sc_pair.first;
    }
};

struct vdist : public std::binary_function<state_t, state_t, float> {
    float operator()(const state_t &s1, const state_t &s2) { return sqrtf((s1-s2).transpose()*(s1-s2)); }
};

struct mahalanobis : public std::binary_function<kalman_filter::state_covariance_pair, state_t, float> {
    float operator()(const kalman_filter::state_covariance_pair &sc, const state_t &s) {
        return sqrtf( (sc.first-s).transpose() * sc.second.inverse() * (sc.first-s) );
    }
};

Eigen::Vector2f loc(const observation& o) { return Eigen::Vector2f(x(o), y(o)); }

struct maha_obs : public std::binary_function<observation, state_t, float> {
    matrix2f invR;
    Eigen::Matrix<float, 2, 4> B;
    maha_obs (const matrix2f &R) : invR(R.inverse()) {
        B << 1, 0, 0, 0,
             0, 0, 1, 0;
    }
    float operator()(const observation& o, const state_t& s) {
        Eigen::Vector2f diff(loc(o)-B*s);
        return sqrtf(diff.transpose()*invR*diff);
    }
};

void test6() {
    diag("Test 6:");
    model::parameters paras;
    model::observation_probability(paras) = 1.0f;
    const float mean_track_length = 100.f;
    model::survival_probability(paras) = 1.f - 1.f/mean_track_length;
    model::birth_rate(paras) = 0.5f;
    model::clutter_rate(paras) = 0.5f;
    model::constraint_radius(paras) = 1.0f;
    float loc_error = 0.3f;
    model::observation_error_covariance(paras) << loc_error * loc_error, 0, 0, loc_error * loc_error;
    model::process_noise_covariance(paras) = detail::initQ();
    float position_sigma = sqrt(model::process_noise_covariance(paras)(0,0));
    float velocity_sigma = sqrt(model::process_noise_covariance(paras)(1,1));

    const float thresh = sqrt(-2.f*logf(1.0-0.68));
    float total_len = 0.f;
    float dist_pred = 0.f;
    float dist_corr = 0.f;
    float dist_smth = 0.f;
    float lt1_pred = 0.f;
    float lt1_corr = 0.f;
    float lt1_smth = 0.f;
    float obs_stat = 0.f;
    float obs_pred = 0.f;
    float obs_corr = 0.f;
    float obs_smth = 0.f;
    maha_obs maha_dist(model::observation_error_covariance(paras));
    for (int loop = 0; loop < 1000; ++loop) {
        std::deque<state_t> states;
        shared_const_track_ptr syntrack = simulate::generate_track(states, paras, point_2d(), point_2d(), 0,
                                            position_sigma, velocity_sigma);

        boost::shared_ptr<kalman_filter> kf_p(new kalman_filter(*syntrack,
            model::observation_error_covariance(paras), model::process_noise_covariance(paras)));
        kalman_filter::states_and_cov_deque predictions = kf_p->predictions();
        kalman_filter::states_and_cov_deque corrections = kf_p->corrections();
        kalman_filter::states_and_cov_deque smoothed;
        rts_smooth(*kf_p, std::front_inserter(smoothed));
        std::deque<state_t> pred_states(predictions.size());
        std::deque<state_t> corr_states(corrections.size());
        std::deque<state_t> smth_states(smoothed.size());
        std::transform(predictions.begin(), predictions.end(), pred_states.begin(), extract_state());
        std::transform(corrections.begin(), corrections.end(), corr_states.begin(), extract_state());
        std::transform(smoothed.begin(), smoothed.end(), smth_states.begin(), extract_state());
        std::deque<float> dist_to_gt(states.size());
        std::transform(states.begin(), states.end(), pred_states.begin(), dist_to_gt.begin(), vdist());
        total_len += states.size();
        dist_pred += std::accumulate(dist_to_gt.begin(), dist_to_gt.end(), 0.f);
        std::transform(states.begin(), states.end(), corr_states.begin(), dist_to_gt.begin(), vdist());
        dist_corr += std::accumulate(dist_to_gt.begin(), dist_to_gt.end(), 0.f);
        std::transform(states.begin(), states.end(), smth_states.begin(), dist_to_gt.begin(), vdist());
        dist_smth += std::accumulate(dist_to_gt.begin(), dist_to_gt.end(), 0.f);
        std::transform(predictions.begin(), predictions.end(), states.begin(), dist_to_gt.begin(), mahalanobis());
        lt1_pred += std::accumulate(dist_to_gt.begin(), dist_to_gt.end(), 0.f);
        std::transform(corrections.begin(), corrections.end(), states.begin(), dist_to_gt.begin(), mahalanobis());
        lt1_corr += std::accumulate(dist_to_gt.begin(), dist_to_gt.end(), 0.f);
        std::transform(smoothed.begin(), smoothed.end(), states.begin(), dist_to_gt.begin(), mahalanobis());
        lt1_smth += std::accumulate(dist_to_gt.begin(), dist_to_gt.end(), 0.f);

        std::transform(syntrack->begin(), syntrack->end(), states.begin(), dist_to_gt.begin(), maha_dist);
        obs_stat += count_if(dist_to_gt.begin(), dist_to_gt.end(), std::bind2nd(std::less<float>(), thresh));
        std::transform(syntrack->begin(), syntrack->end(), pred_states.begin(), dist_to_gt.begin(), maha_dist);
        obs_pred += count_if(dist_to_gt.begin(), dist_to_gt.end(), std::bind2nd(std::less<float>(), thresh));
        std::transform(syntrack->begin(), syntrack->end(), corr_states.begin(), dist_to_gt.begin(), maha_dist);
        obs_corr += count_if(dist_to_gt.begin(), dist_to_gt.end(), std::bind2nd(std::less<float>(), thresh));
        std::transform(syntrack->begin(), syntrack->end(), smth_states.begin(), dist_to_gt.begin(), maha_dist);
        obs_smth += count_if(dist_to_gt.begin(), dist_to_gt.end(), std::bind2nd(std::less<float>(), thresh));

    }
    diag("average Euclidian GT-prediction %.2f", dist_pred/total_len);
    diag("average Euclidian GT-correction %.2f", dist_corr/total_len);
    diag("average Euclidian GT-smoothed %.2f", dist_smth/total_len);
    diag("--");
    diag("average Mahalanobis GT-prediction = %.2f", lt1_pred/total_len);
    diag("average Mahalanobis GT-correction = %.2f", lt1_corr/total_len);
    diag("average Mahalanobis GT-smooted    = %.2f", lt1_smth/total_len);
    diag("--");
    diag("Mahalanobis obs-GT         less than thresh = %.2f%%", 100.f*obs_stat/total_len);
    diag("Mahalanobis obs-prediction less than thresh = %.2f%%", 100.f*obs_pred/total_len);
    diag("Mahalanobis obs-correction less than thresh = %.2f%%", 100.f*obs_corr/total_len);
    diag("Mahalanobis obs-smoothed   less than thresh = %.2f%%", 100.f*obs_smth/total_len);
}

void test5() {
    diag("Test 5:");
    model::parameters paras;
    model::observation_probability(paras) = 1.0f;
    const float mean_track_length = 20.f;
    model::survival_probability(paras) = 1.f - 1.f/mean_track_length;
    model::birth_rate(paras) = 0.5f;
    model::clutter_rate(paras) = 0.5f;
    model::constraint_radius(paras) = 1.0f;
    float loc_error = 0.3f;
    model::observation_error_covariance(paras) << loc_error * loc_error, 0, 0, loc_error * loc_error;
    model::process_noise_covariance(paras) = detail::initQ();
    float position_sigma = sqrt(model::process_noise_covariance(paras)(0,0));
    float velocity_sigma = sqrt(model::process_noise_covariance(paras)(1,1));
    std::deque<state_t> states;
    shared_const_track_ptr syntrack = simulate::generate_track(states, paras, point_2d(), point_2d(), 0,
                                        position_sigma, velocity_sigma);

    boost::shared_ptr<kalman_filter> kf_p(new kalman_filter(*syntrack,
        model::observation_error_covariance(paras), model::process_noise_covariance(paras)));
    kalman_filter::states_and_cov_deque predictions = kf_p->predictions();
    kalman_filter::states_and_cov_deque corrections = kf_p->corrections();
    kalman_filter::states_and_cov_deque smoothed;
    rts_smooth(*kf_p, std::front_inserter(smoothed));
    std::deque<state_t> pred_states(predictions.size());
    std::deque<state_t> corr_states(corrections.size());
    std::deque<state_t> smth_states(smoothed.size());
    std::transform(predictions.begin(), predictions.end(), pred_states.begin(), extract_state());
    std::transform(corrections.begin(), corrections.end(), corr_states.begin(), extract_state());
    std::transform(smoothed.begin(), smoothed.end(), smth_states.begin(), extract_state());
    std::deque<float> dist_to_gt(states.size());
    diag("num of GT states = %zu", states.size());
    diag("num of GT observations = %zu", syntrack->size());
    diag("num of predicted states = %zu", pred_states.size());
    diag("num of corrected states = %zu", corr_states.size());
    diag("num of smoothed states = %zu", smth_states.size());
    ok(states.size() == pred_states.size(), "predicted number of states == ground truth states");
    ok(states.size() == corr_states.size(), "corrected number of states == ground truth states");
    ok(states.size() == smth_states.size(), " smoothed number of states == ground truth states");
    diag("--");
    BOOST_ASSERT(states.size() == pred_states.size());
    std::transform(states.begin(), states.end(), pred_states.begin(), dist_to_gt.begin(), vdist());
    float len(states.size());
    float sum = std::accumulate(dist_to_gt.begin(), dist_to_gt.end(), 0.f);
    diag("residual GT-prediction %.2f", sum/len);
    std::transform(states.begin(), states.end(), corr_states.begin(), dist_to_gt.begin(), vdist());
    sum = std::accumulate(dist_to_gt.begin(), dist_to_gt.end(), 0.f);
    diag("residual GT-correction %.2f", sum/len);
    std::transform(states.begin(), states.end(), smth_states.begin(), dist_to_gt.begin(), vdist());
    sum = std::accumulate(dist_to_gt.begin(), dist_to_gt.end(), 0.f);
    diag("residual GT-smoothed   %.2f", sum/len);
    float lt1;
    float thresh = sqrt(-2.f*logf(1.0-0.68));
    std::transform(predictions.begin(), predictions.end(), states.begin(), dist_to_gt.begin(), mahalanobis());
    lt1 = count_if(dist_to_gt.begin(), dist_to_gt.end(), std::bind2nd(std::less<float>(), thresh));
    sum = std::accumulate(dist_to_gt.begin(), dist_to_gt.end(), 0.f);
    diag("Mahalanobis GT-prediction %.2f, less than %.2f = %.2f%%", sum/len, thresh, 100.f*lt1/len);
    std::transform(corrections.begin(), corrections.end(), states.begin(), dist_to_gt.begin(), mahalanobis());
    lt1 = count_if(dist_to_gt.begin(), dist_to_gt.end(), std::bind2nd(std::less<float>(), thresh));
    sum = std::accumulate(dist_to_gt.begin(), dist_to_gt.end(), 0.f);
    diag("Mahalanobis GT-correction %.2f, less than %.2f = %.2f%%", sum/len, thresh, 100.f*lt1/len);
    std::transform(smoothed.begin(), smoothed.end(), states.begin(), dist_to_gt.begin(), mahalanobis());
    lt1 = count_if(dist_to_gt.begin(), dist_to_gt.end(), std::bind2nd(std::less<float>(), thresh));
    sum = std::accumulate(dist_to_gt.begin(), dist_to_gt.end(), 0.f);
    diag("Mahalanobis GT-smoothed   %.2f, less than %.2f = %.2f%%", sum/len, thresh, 100.f*lt1/len);
    diag("--");
    //clutter_ptr clutter = clutter_ptr(new clutter_t());
    //mh_moves::shared_track_collection_ptr tracks = mh_moves::shared_track_collection_ptr(new track_collection());
    //tracks->insert(syntrack);
    //partition synpart(tracks, clutter);
    maha_obs maha_dist(model::observation_error_covariance(paras));
    std::transform(syntrack->begin(), syntrack->end(), states.begin(), dist_to_gt.begin(), maha_dist);
    lt1 = count_if(dist_to_gt.begin(), dist_to_gt.end(), std::bind2nd(std::less<float>(), thresh));
    sum = std::accumulate(dist_to_gt.begin(), dist_to_gt.end(), 0.f);
    diag("Mahalanobis obs-GT %.2f, less than %.2f = %.2f%%", sum/len, thresh, 100.f*lt1/len);
    std::transform(syntrack->begin(), syntrack->end(), pred_states.begin(), dist_to_gt.begin(), maha_dist);
    lt1 = count_if(dist_to_gt.begin(), dist_to_gt.end(), std::bind2nd(std::less<float>(), thresh));
    sum = std::accumulate(dist_to_gt.begin(), dist_to_gt.end(), 0.f);
    diag("Mahalanobis obs-predicted %.2f, less than %.2f = %.2f%%", sum/len, thresh, 100.f*lt1/len);
    std::transform(syntrack->begin(), syntrack->end(), corr_states.begin(), dist_to_gt.begin(), maha_dist);
    lt1 = count_if(dist_to_gt.begin(), dist_to_gt.end(), std::bind2nd(std::less<float>(), thresh));
    sum = std::accumulate(dist_to_gt.begin(), dist_to_gt.end(), 0.f);
    diag("Mahalanobis obs-corrected %.2f, less than %.2f = %.2f%%", sum/len, thresh, 100.f*lt1/len);
    std::transform(syntrack->begin(), syntrack->end(), smth_states.begin(), dist_to_gt.begin(), maha_dist);
    lt1 = count_if(dist_to_gt.begin(), dist_to_gt.end(), std::bind2nd(std::less<float>(), thresh));
    sum = std::accumulate(dist_to_gt.begin(), dist_to_gt.end(), 0.f);
    diag("Mahalanobis obs-smoothed %.2f, less than %.2f = %.2f%%", sum/len, thresh, 100.f*lt1/len);

    //
    diag("--");
    diag("track length %zu, number of obs %zu", syntrack->duration(), syntrack->size());
    float pdf = log_track_given_parameters_density(syntrack, paras);
    diag("Q[0, 0] = %g, track likelihood %f", model::process_noise_covariance(paras)(0, 0), pdf);
    diag("--");
    float base = 1.8257418583505538f;
    for (int i = -6; i < 6; ++i) {
        model::process_noise_covariance(paras)(0, 0) = pow(base, i) * pow(base, i);
        model::process_noise_covariance(paras)(2, 2) = pow(base, i) * pow(base, i);
        pdf = log_track_given_parameters_density(syntrack, paras);
        diag("Q[0, 0] = %g, track likelihood %f", model::process_noise_covariance(paras)(0, 0), pdf);
    }
}

void test4() {
    general_paras gen_pars;
    gen_pars.total = 200000;
    gen_pars.seed = 0;
    gen_pars.rgp.lambda = 2.5;
    gen_pars.rgp.p_no = 0.2;
    gen_pars.rgp.p_yes = 0.2;
    gen_pars.rgp.p_tr = 0.5;
    gen_pars.rgp.min_tracks = 1;
    gen_pars.rgp.max_tracks = 2;
    partition_ptr_t part = random_debug_partition(gen_pars.rgp);
    model::parameters paras;
    generate_parameters(paras);
    /*
    std::cerr << "Q=" << std::endl;
    std::cerr << model::process_noise_covariance(paras) << std::endl;
    std::cerr << "initQ=" << std::endl;
    std::cerr << detail::initQ() << std::endl;
    */
    //find_good_parameters(*part, paras, 1000);
    diag("track likelihood %f", model::log_tracks_given_parameters_density(*part, paras));
    diag("track likelihood %f", model::log_tracks_given_parameters_density(*part, paras));
    partition_ptr_t final_part;
    float pdr;
    bool success = false;
    int attempts = 0;
    capability_recorder_ptr cap_rec_ptr(new_cap_recorder(*part));
    part->set_capability_recorder(cap_rec_ptr);
    while (not success and attempts++ < 50) {
        success = mh_moves::update(part, final_part, pdr);
    }
    diag("success: %d, attempts: %d", success, attempts);
    if (success)
        diag("track likelihood %f", model::log_tracks_given_parameters_density(*final_part, paras));
}

int main(int argc, char** argv) {
    // seed the PRNG to ensure consistent results across runs.
    general_paras gen_pars;
    gen_pars.total = 200000;
    gen_pars.seed = 0;
    gen_pars.rgp.lambda = 2.5;
    gen_pars.rgp.p_no = 0.2;
    gen_pars.rgp.p_yes = 0.2;
    gen_pars.rgp.p_tr = 0.5;
    gen_pars.rgp.min_tracks = 1;
    gen_pars.rgp.max_tracks = 2;
    if (not parse_args(argc, argv, gen_pars))
        return exit_status();
    biggles::detail::seed_prng(gen_pars.seed);
    diag("random seed = 0x%x", gen_pars.seed);
    plan_no_plan();

    //test2();
    //test3();
    test4();
    test5();
    test6();
    //diag("dev-from-line-test %.4f", max_deviation_from_line(generate_linear_track(10, 1.f, 0.f)->observations()));
    ok1(true);
    return exit_status();
}
