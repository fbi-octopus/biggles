#include <boost/assert.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
//#include <boost/math/special_functions/binomial.hpp>
#include <cmath>
#include <Eigen/Dense>
#include <limits>
#include <vector>

#include <fstream>

#include "kalman_filter.hpp"
#include "kalman_filter_cache.hpp"
#include "model.hpp"
#include "partition.hpp"
#include "io.hpp"

//using boost::math::binomial_coefficient;

namespace biggles
{

namespace model
{

namespace detail
{

// log multivariate gamma function with p == 2
inline float lgamma_2(float x)
{
    // see http://en.wikipedia.org/wiki/Multivariate_gamma_function
    return 0.5f*logf(M_PI) + lgammaf(x) + lgammaf(x-0.5f);
}

// log of the inverse Wishart PDF.
inline float log_inverse_wishart_pdf(const Eigen::Matrix2f& B, const Eigen::Matrix2f& Phi, int m)
{
    // see http://en.wikipedia.org/wiki/Inverse-Wishart_distribution

    float log_det_B = logf(B.determinant());
    float log_det_Phi = logf(Phi.determinant());
    float trace_Phi_B_inv = (Phi * B.inverse()).trace();

    return 0.5f*m*log_det_Phi - 0.5f*(m+2+1)*log_det_B - 0.5f*trace_Phi_B_inv - 0.5f*(2*m)*logf(2.f) - lgamma_2(0.5f*m);
}

}

/// \brief calculates p(partition | clutter rate)
///
/// This function seems to be uncontroversial
float log_partition_given_clutter_rate_density(const partition& part_sample, const model::parameters& para_sample) {
    const time_stamp part_first_ts = part_sample.first_time_stamp();
    const time_stamp part_last_ts = part_sample.last_time_stamp();
    const time_stamp part_duration = part_last_ts - part_first_ts;
    /*
    std::vector<size_t> clutter_counts(part_duration);
    BOOST_FOREACH(const observation& obs, part_sample.clutter()) {
        ++clutter_counts[t(obs) - part_first_ts];
    }
    */
    float lambda_f = clutter_rate(para_sample);
    /*
    if (not std::isfinite(lambda_f)) {
        std::stringstream sstr;
        sstr << "clutter rate is not finite: " << lambda_f;
        throw std::runtime_error(sstr.str()) ;
    }
    */
    BOOST_ASSERT(std::isfinite(lambda_f));
    BOOST_ASSERT(lambda_f >= 0);
    float log_part_clutter = part_sample.clutter().size() * logf(lambda_f) - part_duration * lambda_f;
    for (clutter_t::const_iterator it = part_sample.clutter().begin(); it not_eq part_sample.clutter().end(); ++it) {
        log_part_clutter += -lgammaf(1.f + it->second.size());
    }
    return log_part_clutter;
}

float log_partition_given_survival_prob_density(const partition& part_sample, const model::parameters& para_sample) {
    binomial_coefficient& binom_coeff = binomial_coefficient::get();
    const time_stamp part_first_ts = part_sample.first_time_stamp();
    const time_stamp part_last_ts = part_sample.last_time_stamp();
    const time_stamp part_duration = part_last_ts - part_first_ts;
    const int min_surv = 0; // number of guaranteed survivals
    std::vector<size_t> survival_events(part_duration);
    std::vector<size_t> death_events(part_duration);
    BOOST_FOREACH(const boost::shared_ptr<const track>& t_ptr, part_sample.tracks()) {
        BOOST_ASSERT(t_ptr->duration() > min_surv);
        time_stamp track_first_ts = t_ptr->first_time_stamp();
        time_stamp track_last_ts = t_ptr->last_time_stamp();
        for (time_stamp ts = track_first_ts - part_first_ts + 1 + min_surv; ts < track_last_ts - part_first_ts; ++ts)
            survival_events[ts]++;
        if (track_last_ts < part_last_ts)
            death_events[track_last_ts - part_first_ts]++;
    }
    float log_part_surv = 0.f;
    float log_p_s = logf(survival_probability(para_sample));
    float log_1m_p_s = logf(1.f - survival_probability(para_sample));
    for ( time_stamp i = 0; i < part_duration; ++i) {
        log_part_surv += logf(binom_coeff(survival_events[i] + death_events[i], survival_events[i])) +
            survival_events[i] * log_p_s + death_events[i] * log_1m_p_s;
    }
    return log_part_surv;
}

float log_partition_given_observation_prob_density(const partition& part_sample, const model::parameters& para_sample) {
    binomial_coefficient& binom_coeff = binomial_coefficient::get();
    const time_stamp part_first_ts = part_sample.first_time_stamp();
    const time_stamp part_last_ts = part_sample.last_time_stamp();
    const time_stamp part_duration = part_last_ts - part_first_ts;
    const size_t min_obs = 0; //number of guaranteed observations
    std::vector<size_t> track_events(part_duration);
    std::vector<size_t> observation_events(part_duration);
    BOOST_FOREACH(const boost::shared_ptr<const track>& t_ptr, part_sample.tracks()) {
        time_stamp track_first_ts = t_ptr->first_time_stamp();
        time_stamp track_last_ts = t_ptr->last_time_stamp();
        for (time_stamp ts = track_first_ts - part_first_ts; ts < track_last_ts - part_first_ts; ++ts)
            track_events[ts]++;
        BOOST_ASSERT(t_ptr->observations().size() >= min_obs);
        observation_collection::const_iterator obs_iter = t_ptr->observations().begin();
        for (size_t i = 0; i < min_obs; ++i) {
            track_events[t(*obs_iter) - part_first_ts]--;
            obs_iter++;
        }
        for (; obs_iter != t_ptr->observations().end(); ++obs_iter)
            observation_events[t(*obs_iter) - part_first_ts]++;
    }
    float log_part_obs = 0.f;
    float log_p_o = logf(observation_probability(para_sample));
    float log_1m_p_o = logf(1.f - observation_probability(para_sample));
    for (time_stamp i = 0; i < part_duration; ++i) {
        log_part_obs += logf(binom_coeff(track_events[i], observation_events[i])) +
            observation_events[i] * log_p_o + (track_events[i] - observation_events[i]) * log_1m_p_o;
    }
    return log_part_obs;
}

float log_partition_given_birth_rate_density(const partition& part_sample, const model::parameters& para_sample) {
    const time_stamp part_first_ts = part_sample.first_time_stamp();
    const time_stamp part_last_ts = part_sample.last_time_stamp();
    const time_stamp part_duration = part_last_ts - part_first_ts;
    const int min_surv = 0; // number of guaranteed survivals
    std::vector<size_t> birth_events(part_duration);
    BOOST_FOREACH(const boost::shared_ptr<const track>& t_ptr, part_sample.tracks()) {
        BOOST_ASSERT(t_ptr->duration() > min_surv);
        time_stamp track_first_ts = t_ptr->first_time_stamp();
        birth_events[track_first_ts - part_first_ts]++;
    }
    float lambda_b = birth_rate(para_sample);
    size_t n_total_born(part_sample.tracks().size());
    float log_part_birth = n_total_born * logf(lambda_b) - (part_duration - min_surv) * lambda_b;
    // TODO: the value for lgammaf may be stored and recalled rather than recalculated; could be faster.
    for (time_stamp i = 0; i < part_duration; ++i)
        log_part_birth += -lgammaf(1.f + birth_events[i]);
    return log_part_birth;
}

float log_partition_given_parameters_density_orig(const partition& partition, const model::parameters& parameters)
{
    float log_pdf = 0.f;

    off_t n_frames = partition.duration();

    // ensemble statistics over all time stamps
    off_t n_total_born(partition.tracks().size());
    off_t n_total_false(partition.clutter().size());
    off_t n_total_died(partition.tracks().size()); // everything that has a beginning has an end, Neo.
    off_t n_total_survived(0);
    off_t n_total_generated(0);

    // calculate n_total_survived and n_total_generated. record birth times for tracks
    std::vector<size_t> birth_counts(partition.last_time_stamp() - partition.first_time_stamp());
    BOOST_FOREACH(const boost::shared_ptr<const track>& t_ptr, partition.tracks())
    {
        BOOST_ASSERT(t_ptr->duration() > 1);
        // a track 'survives' for one fewer ticks than it's duration
        n_total_survived += t_ptr->duration() - 1;


        // the track's size is the number of observations within it
        n_total_generated += t_ptr->size();

        // record birth time
        BOOST_ASSERT(t_ptr->first_time_stamp() >= partition.first_time_stamp());
        BOOST_ASSERT(t_ptr->first_time_stamp() - partition.first_time_stamp() < static_cast<time_stamp>(birth_counts.size()));
        ++birth_counts[t_ptr->first_time_stamp() - partition.first_time_stamp()];
    }

    // record clutter counts
    std::vector<size_t> clutter_counts(partition.last_time_stamp() - partition.first_time_stamp());
    for (clutter_t::const_iterator it = partition.clutter().begin(); it not_eq partition.clutter().end(); ++it) {
        clutter_counts[it->first - partition.first_time_stamp()] = it->second.size();
    }

    // p_s term
    float p_s = frame_to_frame_survival_probability(parameters);
    BOOST_ASSERT(std::isfinite(p_s) and p_s >= 0);
    log_pdf += n_total_survived * logf(p_s) + n_total_died * logf(1.f - p_s);

    // p_d term
    float p_d = generate_observation_probability(parameters);
    BOOST_ASSERT(std::isfinite(p_d) and p_d >= 0);
    log_pdf += n_total_generated * logf(p_d) + (n_total_survived + n_total_born - n_total_generated) * logf(1.f - p_d);

    // start of lambda_b term
    float lambda_b = mean_new_tracks_per_frame(parameters);
    BOOST_ASSERT(std::isfinite(lambda_b) and lambda_b >= 0);
    log_pdf += n_total_born * logf(lambda_b) - n_frames * lambda_b; // ... still to calculate log gamma terms

    // start of lambda_f term
    float lambda_f = mean_false_observations_per_frame(parameters);
    if (not std::isfinite(lambda_f)) {
        std::stringstream sstr;
        sstr << "clutter rate is not finite: " << lambda_f;
        throw std::runtime_error(sstr.str()) ; // TODO debug cleaning up
    }
    BOOST_ASSERT(std::isfinite(lambda_f));
    BOOST_ASSERT(lambda_f >= 0);
    log_pdf += n_total_false * logf(lambda_f) - n_frames * lambda_f; // ... still to calculate log gamma terms

    // finish lambda_b and lambda_f terms
    for(time_stamp ts = partition.first_time_stamp(); ts < partition.last_time_stamp(); ++ts)
    {
        // finish lambda_b term
        log_pdf += -lgammaf(1.f + birth_counts[ts - partition.first_time_stamp()]);

        // finish lambda_f term
        log_pdf += -lgammaf(1.f + clutter_counts[ts - partition.first_time_stamp()]);
    }

    return log_pdf; // the rest is just a repeat of P(T|clutter_rate)

}

float log_partition_given_parameters_density(
        const partition& partition_sample, const model::parameters& parameter_sample)
{
    return log_partition_given_parameters_density_orig(partition_sample, parameter_sample);
}

float log_clutter_given_parameters_density(const partition& part, const model::parameters& parameters) {
    return -logf(part.volume())*part.clutter().size();
}

float log_track_given_parameters_density(const boost::shared_ptr<const track>& track_p, const model::parameters& parameters)
{
    return track_p->log_posterior(parameters);
}

float log_parameters_prior_density(const model::parameters& parameters)
{
    // finite regions of support:
    if(mean_new_tracks_per_frame(parameters) <= 0.f)
        return -std::numeric_limits<float>::max();
    if(mean_false_observations_per_frame(parameters) <= 0.f)
        return -std::numeric_limits<float>::max();
    if(frame_to_frame_survival_probability(parameters) < 0.f)
        return -std::numeric_limits<float>::max();
    if(frame_to_frame_survival_probability(parameters) > 1.f)
        return -std::numeric_limits<float>::max();
    if(generate_observation_probability(parameters) < 0.f)
        return -std::numeric_limits<float>::max();
    if(generate_observation_probability(parameters) > 1.f)
        return -std::numeric_limits<float>::max();


    /*
    if (mean_false_observations_per_frame(parameters) != 0.1f) { // FIXME PRIOR calculation
        return -std::numeric_limits<float>::max();
    }
    if (mean_new_tracks_per_frame(parameters) != .5f) { // FIXME PRIOR calculation
        return -std::numeric_limits<float>::max();
    }
    if (generate_observation_probability(parameters) != .9f) { // FIXME PRIOR calculation
        return -std::numeric_limits<float>::max();
    }
    if (frame_to_frame_survival_probability(parameters) != .9f) { // FIXME PRIOR calculation
        return -std::numeric_limits<float>::max();
    }
    */
    // prior on R:
    // FIXME PRIOR calculation
    /*
    if ( observation_error_covariance(parameters) == Eigen::Matrix2f::Identity()* 0.09f ) {
        return 0.f;
    } else {
        return -std::numeric_limits<float>::max();
    }
    */
    return detail::log_inverse_wishart_pdf(observation_error_covariance(parameters), 2.f * Eigen::Matrix2f::Identity(), 5);
}

float log_likelihood(const partition& part_sample, const model::parameters& para_sample) {
    float log_pdf(0.f);

    // track PDFs
    BOOST_FOREACH(const boost::shared_ptr<const track>& t_ptr, part_sample.tracks())
    {
        log_pdf += log_track_given_parameters_density(t_ptr, para_sample);
    }

    // clutter PDF
    log_pdf += log_clutter_given_parameters_density(part_sample, para_sample);
    return log_pdf;
}

float log_tracks_given_parameters_density(const partition& part, const model::parameters& parameters) {
    // track PDFs
    float log_pdf(0.f);
    BOOST_FOREACH(const boost::shared_ptr<const track>& t_ptr, part.tracks())
    {
        float track_log_pdf = log_track_given_parameters_density(t_ptr, parameters);
        if (not std::isfinite(track_log_pdf)) {
            std::cerr << "track duration " << t_ptr->duration() << std::endl;
            std::cerr << "track size " << t_ptr->size() << std::endl;
        }
        log_pdf += track_log_pdf;
    }
    return log_pdf;
}

void dump_part_para(const partition& part, const model::parameters& para) {
    size_t ctr = 0;
    std::string fname = boost::str(boost::format("no_finite_log_pdf.%04d.json") % ctr);
    while (boost::filesystem::exists(fname)) {
        ctr++;
        if (ctr > 9999) {
            std::cerr << "nothing dumped :(" << std::endl;
            return;
        }
        fname = boost::str(boost::format("no_finite_log_pdf.%04d.json") % ctr);
    }
    std::ofstream fh(fname.c_str());
    fh << "{ \"parameters\": ";
    io::write_model_parameters_to_json_stream(fh, para);
    fh << ", \"partition\": " << std::endl;
    io::write_partition_to_json_stream(fh, part);
    fh << "}" << std::endl;
}

float log_partition_given_parameters_and_data_density(const partition& part, const model::parameters& parameters)
{
    float log_pdf(0.f);
    float temp(0.f);

    // clutter PDF
    temp = log_clutter_given_parameters_density(part, parameters);
    if (not std::isfinite(temp)) {
        std::cerr << " *** clutter log density not finite" << std::endl;
        std::cerr << "     " << part.volume() << std::endl;
        dump_part_para(part, parameters);
    }
    log_pdf += temp;

    temp = log_tracks_given_parameters_density(part, parameters);
    if (not std::isfinite(temp)) {
        std::cerr << " *** track log density not finite" << std::endl;
        std::cerr << temp << std::endl;
        std::cerr << model::observation_error_covariance(parameters) << std::endl;
        dump_part_para(part, parameters);
    }
    log_pdf += temp;


    //log_pdf = 0.f; // FIXME PRIOR test get rid of this line again

    // data-independent track PDF
    temp = log_partition_given_parameters_density(part, parameters);
    if (not std::isfinite(temp)) {
        std::cerr << " *** partition log density not finite" << std::endl;
        dump_part_para(part, parameters);
    }
    log_pdf += temp;

    // prior
    // the prior is only need to record the log pdf and not for the sampling itself.
    // log_pdf += log_parameters_prior_density(parameters);

    return log_pdf;
}

}

}
