#include "detail/random.hpp"
#include "model.hpp"
#include "observation.hpp"
#include "simulate.hpp"
#include "track.hpp"

#include <boost/random.hpp>
#include <cmath>
#include <ctime>
#include <deque>
#include <Eigen/Dense>

using namespace Eigen;

namespace biggles { namespace simulate {

boost::shared_ptr<track> generate_track(
    const model::parameters& parameters,
    const point_2d& initial_position,
    const point_2d& initial_velocity,
    time_stamp first_time_stamp,
    float position_sigma,
    float velocity_sigma)
{
    std::deque<state_t> states;
    return generate_track(states, parameters, initial_position, initial_velocity, first_time_stamp,
        position_sigma, velocity_sigma);
}

boost::shared_ptr<track> generate_track(
    std::deque<state_t> &states,
    const model::parameters& parameters,
    const point_2d& initial_position,
    const point_2d& initial_velocity,
    time_stamp first_time_stamp,
    float position_sigma,
    float velocity_sigma)
{
    boost::uniform_real<float> uni_real_dist(0,1);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<float> >
        uni_real(detail::random_generator, uni_real_dist);

    boost::normal_distribution<float> norm_dist(0,1);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<float> >
        norm(detail::random_generator, norm_dist);

    // the observation vector type: [x, y]
    typedef Vector2f obs_t;

    // the state evolution matrix: X_{k+1} = A X_k + <noise>
    Matrix4f A;
    //float d(uni_real() < 0.5 ? 1 : 0); // the dynamic drag
    const float d(1.0f);
    A << 1, d, 0, 0,
         0, 1, 0, 0,
         0, 0, 1, d,
         0, 0, 0, 1;

    // the observation matrix: Y_k = B X_k + <noise>
    Matrix<float, 2, 4> B;
    B << 1, 0, 0, 0,
         0, 0, 1, 0;

    // initial state
    state_t X;
    X << x(initial_position), x(initial_velocity), y(initial_position), y(initial_velocity);

    // collection we'll store generated observations in
    std::deque<observation> observations;

    // extract the parameters we're interested in
    Matrix2f R(model::observation_error_covariance(parameters));
    float p_s(model::frame_to_frame_survival_probability(parameters));
    float p_d(model::generate_observation_probability(parameters));

    float sqrt_R_00(sqrt(R(0,0))), sqrt_R_11(sqrt(R(1,1)));

    // make sure we generate at least one observation. Keep time stamp in global scope because we use it later as the
    // duration.
    time_stamp time_stamp(0);
    for(time_stamp = first_time_stamp; observations.empty() || (uni_real() <= p_s); ++time_stamp)
    {
        states.push_back(X);
        // generate observation
        obs_t Y = B * X;

        // noise observation - FIXME: does not use cross-covariance
        Y[0] += sqrt_R_00 * norm();
        Y[1] += sqrt_R_11 * norm();

        // record observation if an observation is made
        if(uni_real() <= p_d)
            observations.push_back(new_obs(Y[0], Y[1], time_stamp));

        // evolve state
        X = A * X;

        // noise state
        X[0] += position_sigma * norm();
        X[1] += velocity_sigma * norm();
        X[2] += position_sigma * norm();
        X[3] += velocity_sigma * norm();
    }

    return boost::shared_ptr<track>(
        new track(first_time_stamp, time_stamp,
                  observations.begin(),
                  observations.end(), d));
}

void generate_partition(partition_ptr_t& output_part_ptr,
                        const model::parameters& params,
                        time_stamp start_time_stamp,
                        time_stamp end_time_stamp,
                        float image_width,
                        float image_height,
                        float position_sigma,
                        float velocity_sigma,
                        float initial_velocity_sigma)
{
    // the various random variate generators.
    boost::uniform_real<float> uni_real_dist(0,1);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<float> >
        uni_real(detail::random_generator, uni_real_dist);

    boost::poisson_distribution<int, float> birth_count_dist(model::mean_new_tracks_per_frame(params));
    boost::variate_generator<boost::mt19937&, boost::poisson_distribution<int, float> >
        birth_count(detail::random_generator, birth_count_dist);

    boost::poisson_distribution<int, float> clutter_count_dist(model::mean_false_observations_per_frame(params));
    boost::variate_generator<boost::mt19937&, boost::poisson_distribution<int, float> >
        clutter_count(detail::random_generator, clutter_count_dist);

    boost::normal_distribution<float> norm_dist(0,1);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<float> >
        norm(detail::random_generator, norm_dist);

    boost::shared_ptr<track_collection> tracks(new track_collection());
    clutter_ptr clutter(new clutter_t());

    // the simulation work proper
    size_t expected_clutter_size(0);
    for(time_stamp time_stamp = start_time_stamp; time_stamp < end_time_stamp; ++time_stamp)
    {
        size_t n_new_tracks(birth_count());
        for(size_t track_idx = 0; track_idx < n_new_tracks; ++track_idx)
        {
            boost::shared_ptr<track> new_track = generate_track(
                params,
                point_2d(image_width * uni_real(), image_height * uni_real()),
                point_2d(initial_velocity_sigma * norm(), initial_velocity_sigma * norm()),
                time_stamp, position_sigma, velocity_sigma);

            if(new_track->size() >= 2)
                tracks->insert(new_track);
        }

        size_t n_clutter_obs(clutter_count());
        for(size_t clutter_idx = 0; clutter_idx < n_clutter_obs; ++clutter_idx)
        {
            clutter->insert( new_obs(image_width * uni_real(), image_height * uni_real(), time_stamp));
        }
        BOOST_ASSERT(clutter->count_at_time_stamp(time_stamp) == n_clutter_obs);

        expected_clutter_size += n_clutter_obs;
        BOOST_ASSERT(clutter->size() == expected_clutter_size);
    }

    // set output
    output_part_ptr = partition_ptr_t(new partition(tracks, clutter, partition::expansion_2d(0.f, image_width, 0.f, image_height)));
}

void crop_partition(partition_ptr_t& output_part_ptr,
                    const partition_ptr_t& input_part_ptr,
                    float first_x, float last_x,
                    float first_y, float last_y,
                    time_stamp first_time_stamp,
                    time_stamp last_time_stamp)
{
    clutter_ptr new_clutter(new clutter_t());
    boost::shared_ptr<track_collection> new_tracks(new track_collection());

    // crop clutter observations
    for (clutter_t::const_iterator it = input_part_ptr->clutter().begin();
                    it not_eq input_part_ptr->clutter().end(); ++it) {
        BOOST_FOREACH(const observation& o, it->second) {
            if((x(o) < first_x) || (x(o) >= last_x))
                continue;
            if((y(o) < first_y) || (y(o) >= last_y))
                continue;
            if((t(o) < first_time_stamp) || (t(o) >= last_time_stamp))
                continue;
            new_clutter->insert(o);
        }
    }

    // create new cropped tracks and insert them into the track collection
    BOOST_FOREACH(const boost::shared_ptr<const track>& p_track, input_part_ptr->tracks())
    {
        // create a new set of cropped observations
        std::deque<observation> cropped_observations;
        BOOST_FOREACH(const observation& o, *p_track)
        {
            if((x(o) < first_x) || (x(o) >= last_x))
                continue;
            if((y(o) < first_y) || (y(o) >= last_y))
                continue;
            if((t(o) < first_time_stamp) || (t(o) >= last_time_stamp))
                continue;
            cropped_observations.push_back(o);
        }

        // if no observations remain, remove this track
        if(0 == cropped_observations.size())
            continue;

        // create a new track based on old track.
        boost::shared_ptr<track> new_track(
                new track(std::max(p_track->first_time_stamp(), first_time_stamp),
                          std::min(p_track->last_time_stamp(), last_time_stamp),
                          cropped_observations.begin(),
                          cropped_observations.end(), p_track->dynamic_drag()));

        if(new_track->size() >= 2)
            new_tracks->insert(new_track);
    }

    // assign new tracks and clutter to output partition
    output_part_ptr = partition_ptr_t(
        new partition(new_tracks, new_clutter, partition::expansion_2d(first_x, last_x, first_y, last_y)));
}

void demote_all_tracks_from_partition(partition_ptr_t& output_part_ptr,
                                      const partition_ptr_t& input_part_ptr)
{
    // copy the clutter ready to demote observations
    clutter_ptr clutter(new clutter_t(input_part_ptr->clutter()));

    // a new, empty set of tracks
    boost::shared_ptr<track_collection> tracks(new track_collection());

    // for each input track...
    BOOST_FOREACH(const boost::shared_ptr<const track>& p_track, input_part_ptr->tracks())
    {
        // copy it's observations to the clutter
        clutter->insert(p_track->begin(), p_track->end());
    }

    // return new partition
    output_part_ptr = partition_ptr_t(new partition(tracks, clutter, input_part_ptr->expansion()));
}

} }
