#include "biggles/observation.hpp"
#include "biggles/observation_collection.hpp"
#include "biggles/clutter.hpp"
#include "biggles/partition.hpp"
#include "biggles/track.hpp"
#include "biggles/track_collection.hpp"
#include "biggles/io.hpp"

#include <iostream>

#include <boost/foreach.hpp>

#include "utility/jansson_wrapper.hpp"

using namespace jansson;

namespace biggles { namespace io {

namespace detail
{


    /// @brief Write and array of observations to the specified JSON array.
    ///
    /// @param observations
    /// @param output
    value write_observations_to_array(const clutter_t& clutter)
    {
        value output(json_array());
        for (clutter_t::const_iterator it = clutter.begin(); it != clutter.end(); ++it) {
            //std::deque<observation> obs = it->second->observations();
            const std::deque<observation>& obs = it->second;
            for (std::deque<observation>::const_iterator ot = obs.begin(); ot not_eq obs.end(); ++ot) {
                value observation_array(json_pack("[ffi]",
                    (double) x(*ot), (double) y(*ot), (int) t(*ot)));
                json_array_append(output.get(), observation_array.get());
            }
        }
        return output;
    }

    /// @brief Write and array of observations to the specified JSON array.
    ///
    /// @param observations
    /// @param output
    value write_observations_to_array(const observation_collection& observations)
    {
        value output(json_array());
        BOOST_FOREACH(const observation& obs, observations)
        {
            value observation_array(json_pack("[ffi]",
                (double) x(obs), (double) y(obs), (int) t(obs)));
            json_array_append(output.get(), observation_array.get());
        }
        return output;
    }

    /// @brief Write a track object into a JSON object.
    ///
    /// @param t
    value write_track_to_object(const track& t)
    {
        return value(json_pack("{s[ii]sOsf}",
            "time_span", (int) t.first_time_stamp(), (int) t.last_time_stamp(),
            "observations", write_observations_to_array(t.observations()).get(),
            "dynamic_drag", (double) t.dynamic_drag()
            ));
    }

    boost::shared_ptr<const observation_collection>
    observation_collection_from_json_array(const value& observations)
    {
        if(!observations || !json_is_array(observations.get())) {
            throw std::invalid_argument("Expected a JSON array");
        }

        boost::shared_ptr<observation_collection> obs_col(new observation_collection());

        for(size_t idx=0; idx < json_array_size(observations.get()); ++idx)
        {
            double x, y;
            int t;

            value obs_value(json_incref(json_array_get(observations.get(), idx)));
            unpack(obs_value, "[FFi]", &x , &y, &t);
            obs_col->insert(new_obs(x, y, t));
        }

        return obs_col;
    }

    clutter_ptr clutter_from_json_array(const value& observations) {
        if(!observations || !json_is_array(observations.get())) {
            throw std::invalid_argument("Expected a JSON array");
        }

        clutter_ptr clutter(new clutter_t());

        for(size_t idx=0; idx < json_array_size(observations.get()); ++idx)
        {
            double x, y;
            int t;

            value obs_value(json_incref(json_array_get(observations.get(), idx)));
            unpack(obs_value, "[FFi]", &x , &y, &t);
            clutter->insert(new_obs(x, y, t));
        }

        return clutter;
    }

    boost::shared_ptr<const track>
    track_from_json_object(const value& track_obj)
    {
        if(!track_obj || !json_is_object(track_obj.get())) {
            throw std::invalid_argument("Expected a JSON object");
        }

        boost::shared_ptr<const observation_collection> obs(new observation_collection());
        int first_t(0), last_t(0);
        double dynamic_drag = 1.f;

        for(void* iter = json_object_iter(track_obj.get()); iter != NULL; iter = json_object_iter_next(track_obj.get(), iter))
        {
            std::string key(json_object_iter_key(iter));
            value val(json_incref(json_object_iter_value(iter)));

            if(key == "observations")
            {
                obs = observation_collection_from_json_array(val);
            }
            else if(key == "time_span")
            {
                unpack(val, "[ii]", &first_t, &last_t);
                BOOST_ASSERT(last_t >= first_t);
            }
            else if(key == "dynamic_drag")
            {
                unpack(val, "F", &dynamic_drag);
            }
        }

        return boost::shared_ptr<const track>(new track(first_t, last_t, obs->begin(), obs->end(), dynamic_drag));
    }
}

value write_partition_to_json_object(const biggles::partition& p)
{
    using namespace detail;

    value clutter(write_observations_to_array(p.clutter()));

    value tracks(json_array());
    BOOST_FOREACH(const shared_const_track_ptr& t, p.tracks())
    {
        value track_obj(write_track_to_object(*t));
        json_array_append(tracks.get(), track_obj.get());
    }

    return value(json_pack("{sOsO}",
            "clutter", clutter.get(),
            "tracks", tracks.get()));
}

void write_partition_to_json_stream(std::ostream& os, const partition& p)
{
    write(os, write_partition_to_json_object(p));
}

boost::shared_ptr<biggles::partition> read_partition_from_json_object(const value& obj)
{
    using namespace detail;

    if(!obj || !json_is_object(obj.get())) {
        throw std::invalid_argument("Expected a JSON object");
    }

    const_clutter_ptr clutter(new clutter_t());
    boost::shared_ptr<track_collection> tracks(new track_collection());

    for(void* iter = json_object_iter(obj.get()); iter != NULL; iter = json_object_iter_next(obj.get(), iter))
    {
        std::string key(json_object_iter_key(iter));
        value val(json_incref(json_object_iter_value(iter)));

        if(key == "clutter")
        {
            clutter = clutter_from_json_array(val);
        }
        else if(key == "tracks")
        {
            if(!val || !json_is_array(val.get())) {
                throw std::runtime_error("Expected track array");
            }

            for(size_t idx=0; idx < json_array_size(val.get()); ++idx)
            {
                value track_obj(json_incref(json_array_get(val.get(), idx)));
                tracks->insert(track_from_json_object(track_obj));
            }
        }
    }

    return boost::shared_ptr<biggles::partition>(new partition(tracks, clutter));
}

boost::shared_ptr<biggles::partition> read_partition_from_json_stream(std::istream& is)
{
    return read_partition_from_json_object(read(is));
}

value write_model_parameters_to_json_object(const biggles::model::parameters& params)
{
    using namespace biggles::model;

    value obs_cov(json_array());
    for(int row_idx=0; row_idx<2; ++row_idx)
    {
        json_array_append_new(obs_cov.get(), json_pack("[ff]",
            static_cast<double>(observation_error_covariance(params)(row_idx,0)),
            static_cast<double>(observation_error_covariance(params)(row_idx,1))));
    }

    value proc_cov(json_array());
    for(int row_idx=0; row_idx<4; ++row_idx)
    {
        json_array_append_new(proc_cov.get(), json_pack("[ffff]",
            static_cast<double>(process_noise_covariance(params)(row_idx,0)),
            static_cast<double>(process_noise_covariance(params)(row_idx,1)),
            static_cast<double>(process_noise_covariance(params)(row_idx,2)),
            static_cast<double>(process_noise_covariance(params)(row_idx,3))));
    }


    return value(json_pack("{sfsfsfsfsfsOsO}",
        "births_per_frame", (double) mean_new_tracks_per_frame(params),
        "clutter_per_frame", (double) mean_false_observations_per_frame(params),
        "survival_probability", (double) frame_to_frame_survival_probability(params),
        "observation_probability", (double) generate_observation_probability(params),
        "constraint_radius", (double) constraint_radius(params),
        "observation_covariance", obs_cov.get(),
        "process_noise_covariance", proc_cov.get()
        ));
}

void write_model_parameters_to_json_stream(std::ostream& os,
                const biggles::model::parameters& model_parameters)
{
    write(os, write_model_parameters_to_json_object(model_parameters));
}

boost::shared_ptr<biggles::model::parameters>
    read_model_parameters_from_json_object(const value& obj)
{
    using namespace biggles::model;
    if(!obj || !json_is_object(obj.get())) {
        throw std::invalid_argument("expected object");
    }

    boost::shared_ptr<model::parameters> model_ptr(new model::parameters());

    for(void* iter = json_object_iter(obj.get()); iter != NULL; iter = json_object_iter_next(obj.get(), iter))
    {
        std::string key(json_object_iter_key(iter));
        value val(json_incref(json_object_iter_value(iter)));
        double num;

        if(key == "births_per_frame")
        {
            unpack(val, "F", &num);
            mean_new_tracks_per_frame(*model_ptr) = num;
        }
        else if(key == "clutter_per_frame")
        {
            unpack(val, "F", &num);
            mean_false_observations_per_frame(*model_ptr) = num;
        }
        else if(key == "survival_probability")
        {
            unpack(val, "F", &num);
            frame_to_frame_survival_probability(*model_ptr) = num;
        }
        else if(key == "observation_probability")
        {
            unpack(val, "F", &num);
            generate_observation_probability(*model_ptr) = num;
        }
        else if(key == "observation_covariance")
        {
            double c11, c12, c21, c22;
            unpack(val, "[[FF][FF]]", &c11, &c12, &c21, &c22);
            observation_error_covariance(*model_ptr)(0, 0) = c11;
            observation_error_covariance(*model_ptr)(0, 1) = c12;
            observation_error_covariance(*model_ptr)(1, 0) = c21;
            observation_error_covariance(*model_ptr)(1, 1) = c22;
        }
        else if(key == "process_noise_covariance")
        {
            double c11, c12, c13, c14, c21, c22, c23, c24;
            double c31, c32, c33, c34, c41, c42, c43, c44;
            unpack(val, "[[FFFF][FFFF][FFFF][FFFF]]",
                &c11, &c12, &c13, &c14,
                &c21, &c22, &c23, &c24,
                &c31, &c32, &c33, &c34,
                &c41, &c42, &c43, &c44
                );
            process_noise_covariance(*model_ptr)(0, 0) = c11;
            process_noise_covariance(*model_ptr)(0, 1) = c12;
            process_noise_covariance(*model_ptr)(0, 2) = c13;
            process_noise_covariance(*model_ptr)(0, 3) = c14;
            process_noise_covariance(*model_ptr)(1, 0) = c21;
            process_noise_covariance(*model_ptr)(1, 1) = c22;
            process_noise_covariance(*model_ptr)(1, 2) = c23;
            process_noise_covariance(*model_ptr)(1, 3) = c24;
            process_noise_covariance(*model_ptr)(2, 0) = c31;
            process_noise_covariance(*model_ptr)(2, 1) = c32;
            process_noise_covariance(*model_ptr)(2, 2) = c33;
            process_noise_covariance(*model_ptr)(2, 3) = c34;
            process_noise_covariance(*model_ptr)(3, 0) = c41;
            process_noise_covariance(*model_ptr)(3, 1) = c42;
            process_noise_covariance(*model_ptr)(3, 2) = c43;
            process_noise_covariance(*model_ptr)(3, 3) = c44;
        }
        else if(key == "constraint_radius")
        {
            unpack(val, "F", &num);
            constraint_radius(*model_ptr) = num;
        }
    }

    return model_ptr;
}

boost::shared_ptr<biggles::model::parameters>
    read_model_parameters_from_json_stream(std::istream& is)
{
    return read_model_parameters_from_json_object(read(is));
}

} }
