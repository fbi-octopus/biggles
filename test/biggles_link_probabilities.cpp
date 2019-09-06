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

// include this last to stop pre-processor macros breaking things
extern "C" {
#include <ccan/tap/tap.h>
}

using namespace biggles;

typedef boost::shared_ptr< const kalman_filter > shared_const_kalman_filter_ptr;

std::ostream& operator<<(std::ostream& os, const kalman_filter::state_vector& state) {
    os << boost::format("[% .4f, % .4f, % .4f, % .4f]") % state[0] % state[1] % state[2] % state[3];
    return os;
}

typedef std::deque<float>::const_iterator data_iter;
std::pair<float, float> mean_stdev(data_iter b, data_iter e) {
    float sum = std::accumulate(b, e, 0.0f);
    float mean = sum / std::distance(b, e);

    float sq_sum = std::inner_product(b, e, b, 0.0f);
    float stdev = std::sqrt(sq_sum / std::distance(b, e) - mean * mean);
    return std::make_pair(mean, stdev);
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


boost::shared_ptr<const observation_collection> generate_observations(
    time_stamp begin_ts, time_stamp end_ts, float x0, float y0, float stdev, float lambda)
{
    boost::shared_ptr<observation_collection> oc_ptr(new observation_collection);
    for (time_stamp ts = begin_ts; ts < end_ts; ++ts)
    {
        int num_obs = sampling::poisson(lambda);
        for (int i = 0; i < num_obs; ++i) {
            oc_ptr->insert(
                new_obs(x0 + sampling::sample_normal() * stdev, y0 + sampling::sample_normal() * stdev, ts));
        }
    }
    return oc_ptr;
}


void print_kalman_filter_states(const kalman_filter::states_and_cov_deque& sc) {
    kalman_filter::states_and_cov_deque::const_iterator sc_iter;
    int i = 0;
    for (sc_iter = sc.begin(); sc_iter != sc.end(); ++sc_iter)
        std::cout
            << boost::format("%3d: ") % i++
            << sc_iter->first
            << ", "
            << boost::format("v = % .4f") % hypotf(sc_iter->first[1], sc_iter->first[3])
            << ", "
            << boost::format("log(||C||) = % .4f") % logf(sc_iter->second.determinant())
            << std::endl;
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

float find_good_parameters(partition_ptr_t& const_partition, model::parameters& good_parameters, const int nloops) {
    float best_pdf = -std::numeric_limits<int>::max();
    model::parameters target_parameters;
    model::process_noise_covariance(target_parameters) = detail::initQ();
    for (int i = 0; i < nloops; ++i) {
        sample_model_parameters_given_partition(*const_partition, target_parameters);
        float log_pdf = model::log_partition_given_parameters_and_data_density(*const_partition, target_parameters);
        if (log_pdf > best_pdf) {
            best_pdf = log_pdf;
            good_parameters = target_parameters;
        }
    }
    return best_pdf;
}

bool generate_parameters(model::parameters& params) {
    model::mean_new_tracks_per_frame(params)           = 1.2f;
    model::mean_false_observations_per_frame(params)   = 0.001f;
    model::frame_to_frame_survival_probability(params) = 0.9;
    model::generate_observation_probability(params)    = 0.9f;
    float sigma_R = 0.3f;
    model::observation_error_covariance(params)       << sigma_R * sigma_R, 0, 0, sigma_R * sigma_R;
    model::process_noise_covariance(params) = detail::initQ();
    model::constraint_radius(params) = sampling::uniform_real(0, 10);
    return true;
}

std::string partition_properties(partition_ptr_t &const_partition) {
    const track_collection &tracks(const_partition->tracks());
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
    float b_obs = float(num_tracks)/float(const_partition->duration());
    float c_obs = float(const_partition->clutter().size())/float(const_partition->duration());
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


struct olink {
    observation to, from;
    olink(const observation &t, const observation &f) : to(t), from(f) {}
    bool operator<(const olink &o) const { return to < o.to or (to == o.to and from < o.from); }
};

int test1() {
    const int total = 1;
    const time_stamp last_ts(80);
    boost::shared_ptr<const observation_collection> clutter_obs(generate_observations(0, last_ts, 0.f, 0.f, 1.5f, 3.f));
    clutter_ptr clutter(new clutter_t(clutter_obs->begin(), clutter_obs->end()));
    boost::shared_ptr<const track_collection> tracks(new track_collection);
    partition_ptr_t original_part_ptr(new partition(tracks, clutter));
    for (int i=0; i< total; ++i) {
        diag(" *** iteration %d *** ", i);
        model::parameters params;
        model::process_noise_covariance(params) = detail::initQ();
        std::map<olink, float> link_count;
        generate_parameters(params);
        const std::string bsp("\b\b\b\b\b\b");
        const int nloops(5000);
        diag("number of loops: = %d", nloops);
        tracker tracking(params, original_part_ptr);
        for(int loop=0; loop<nloops; ++loop) {
            //if (loop % 500 == 0) {
            //    std::clog << bsp << boost::format("%3d%%") % (100*loop/nloops)  << std::flush;
            //}
            tracking();
            BOOST_FOREACH(const shared_const_track_ptr& t_p, tracking.last_partition()->tracks()) {
                const track tr0(*t_p);
                BOOST_ASSERT(tr0.size()>1);
                observation_collection::const_iterator from_it = tr0.observations().begin();
                observation_collection::const_iterator to_it(from_it);
                std::advance(to_it, 1);
                for (; to_it != tr0.observations().end(); ++to_it, ++from_it) link_count[olink(*to_it, *from_it)]++;
            }
        }
        //std::clog << bsp;
        //
        std::map<olink, float>::iterator it = link_count.begin();
        std::map<int, float> dist_count;
        std::map<int, float> from_count;
        std::map<int, float> to_count;
        std::map<int, float> dist_occur;
        std::map<int, float> from_occur;
        std::map<int, float> to_occur;
        for (; it != link_count.end(); ++it) {
              to_count[t(it->first.to)] += it->second;
            dist_count[t(it->first.to) - t(it->first.from)] += it->second;
            from_count[t(it->first.from)] += it->second;
              to_occur[t(it->first.to)] ++;
            dist_occur[t(it->first.to) - t(it->first.from)] ++;
            from_occur[t(it->first.from)] ++;
        }
        std::map<int, float>::iterator fp;
        diag("distance");
        for (int i = 0; i < int(last_ts); ++i) {
            float d_occur = dist_occur[i] == 0.f ? 1.f : float(nloops) * dist_occur[i];
            float f_occur = from_occur[i] == 0.f ? 1.f : float(nloops) * from_occur[i];
            float t_occur =   to_occur[i] == 0.f ? 1.f : float(nloops) *   to_occur[i];
            std::stringstream sstr;
            sstr
                << boost::format("%2d => %5.2f%%, %5.2f%%, %5.2f%%")
                % i
                % (dist_count[i]*100.f/ d_occur)
                % (from_count[i]*100.f/ f_occur)
                % (  to_count[i]*100.f/ t_occur)
                ;
            diag("%s", sstr.str().c_str());
        }
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
    model::process_noise_covariance(params) = detail::initQ();


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
    const float loc_err(.3f);

    std::deque<float> sing_pdf;
    for (int i=0; i<1; ++i) {
        boost::shared_ptr<track_collection> tracks(new track_collection());
        tracks->insert(*generate_linear_track(0, 50, 0.f, 0.f, 1.f, 0.f, loc_err));
        clutter_ptr clutter(new clutter_t());
        partition_ptr_t original_part_ptr(new partition(tracks, clutter));
        model::parameters params;
        model::process_noise_covariance(params) = detail::initQ();
        float log_pdf = find_good_parameters(original_part_ptr, params, 1000);
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
        partition_ptr_t original_part_ptr(new partition(tracks, clutter));
        model::parameters params;
        model::process_noise_covariance(params) = detail::initQ();
        float log_pdf = find_good_parameters(original_part_ptr, params, 1000);
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

int main(int argc, char** argv) {
    // seed the PRNG to ensure consistent results across runs.
    //biggles::detail::seed_prng(0xfacedead);

    plan_no_plan();

    test1();

    ok1(true);
    return exit_status();
}
