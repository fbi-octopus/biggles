#include "biggles/python_api/py_biggles.hpp"

using namespace biggles::python;

BOOST_PYTHON_MODULE(_biggles)
{
    // only implement 'const' members
    typedef std::deque<uint64_t> deque_unit64;
    py::class_<deque_unit64>("_UInt64Deque")
        .def("__len__", &deque_unit64::size)
        .def("clear", &deque_unit64::clear)
        .def("__getitem__", &const_std_random_access<deque_unit64>::get,
             py::return_value_policy<py::copy_const_reference>())
        .def("__iter__", py::iterator<deque_unit64>())
    ;

    typedef std::vector<observation> obs_vector;
    py::class_<obs_vector>("VectorOfObservations")
        .def("__len__", &obs_vector::size, "The length of the vector")
        .def("clear", &obs_vector::clear)
        .def("__getitem__", &const_std_random_access<obs_vector>::get,
             py::return_value_policy<py::copy_const_reference>())
        .def("__iter__", py::iterator<obs_vector>())
    ;
    //py::to_python_converter<observation, observation_to_typle_converter>();

    py::class_<clutter_t::observation_container>("ObservationContainer")
        .def("__len__", &clutter_t::observation_container::size, "The length of the container")
        .def("__getitem__", &const_std_random_access<clutter_t::observation_container>::get,
             py::return_value_policy<py::copy_const_reference>())
        .def("__iter__", py::iterator<clutter_t::observation_container>())
    ;

    py::class_<kalman_filter::state_covariance_pair>(
            "_KalmanFilterStateAndCovPair",
            "A state vector and covariance pair"
        )
        .add_property("state", &get_kfs)
        .add_property("covariance", &get_kfc)
    ;

    py::register_exception_translator<key_error>(&translate_key_error);
    py::register_exception_translator<index_error>(&translate_index_error);

    // only implement 'const' members
    typedef std::deque<uint64_t> deque_unit64;
    py::class_<kalman_filter::states_and_cov_deque>("_StatesAndCovDeque")
        .def("__len__", &kalman_filter::states_and_cov_deque::size)
        .def("__getitem__", &const_std_random_access<kalman_filter::states_and_cov_deque>::get,
             py::return_value_policy<py::copy_const_reference>())
        .def("__iter__", py::iterator<kalman_filter::states_and_cov_deque, py::return_internal_reference<> >())
    ;

    py::class_<kalman_filter, shared_ptr_hack<kalman_filter> >(
            "KalmanFilter",
            "The Kalman filter used to interpolate missing states."
        )
        .add_property("predictions", py::make_function(&kalman_filter::predictions, py::return_internal_reference<>()),
                      "Return an iterable collection of predicted states and covariances.")
        .add_property("corrections", py::make_function(&kalman_filter::corrections, py::return_internal_reference<>()),
                      "Return an iterable collection of backwards-smoothed states and covariances.")
    ;
    py::register_ptr_to_python<boost::shared_ptr<const kalman_filter> >();

    py::enum_<mh_moves::move_type>(
        "MoveType",
        "Move types used as proposals when sampling track configurations via Matropolis-Hastings."
        )
        .value("NONE", mh_moves::NONE)
        .value("BIRTH", mh_moves::BIRTH)
        .value("DEATH", mh_moves::DEATH)
        .value("EXTEND", mh_moves::EXTEND)
        .value("REDUCE", mh_moves::REDUCE)
        .value("SPLIT", mh_moves::SPLIT)
        .value("MERGE", mh_moves::MERGE)
        .value("UPDATE", mh_moves::UPDATE)
        .value("TRANSFER", mh_moves::TRANSFER)
        .value("CROSS_OVER", mh_moves::CROSS_OVER)
        .value("FLIP", mh_moves::FLIP)
        .value("IDENTITY", mh_moves::IDENTITY)
        .value("MOVE_COUNT", mh_moves::MOVE_COUNT)
    ;

    py::class_<model::parameters, boost::shared_ptr<model::parameters> >(
        "ModelParameters",
        "The parameters of the dynamic model used by Biggles."
        )
        .add_property("mean_new_tracks_per_frame", &get_mntpf, &set_mntpf,
                      "Floating point value giving the mean number of tracks which are born in each frame.")
        .add_property("birth_rate", &get_mntpf, &set_mntpf,
                      "Floating point value giving the mean number of tracks which are born in each frame.")
        .add_property("mean_false_observations_per_frame", &get_mfopf, &set_mfopf,
                      "Floating point value giving the mean number of clutter observations in each frame.")
        .add_property("clutter_rate", &get_mfopf, &set_mfopf,
                      "Floating point value giving the mean number of clutter observations in each frame.")
        .add_property("frame_to_frame_survival_probability", &get_ftfsp, &set_ftfsp,
                      "Floating point value giving the probability of a feature surviving from frame to frame.")
        .add_property("survival_prob", &get_ftfsp, &set_ftfsp,
                      "Floating point value giving the probability of a feature surviving from frame to frame.")
        .add_property("generate_observation_probability", &get_gop, &set_gop,
                      "Floating point value giving the probability a feature will be observed in a frame.")
        .add_property("observation_prob", &get_gop, &set_gop,
                      "Floating point value giving the probability a feature will be observed in a frame.")
        .add_property("observation_error_covariance", &get_obs_err, &set_obs_err,
                      "2x2 covariance matrix of the observation error, \"R\"")
        .add_property("process_noise_covariance", &get_process_cov, &set_process_cov,
                      "4x4 covariance matrix of the process noise, \"Q\"")
        .def("initQ", &init_process_cov,
            "initialises the process noise from , \"Q\", as a diagonal matrix from [x, vx, y, vy],"
            " where each diagonal element will be (q*speed_of_light)^2." )
        .def("to_json", &model_parameters_to_json,
             "Serialise the model parameters to a string containing a JSON py::object.")
        .def("from_json", &model_parameters_from_json,
             "Construct a new :py:class:`ModelParameters` py::object from a string containing a serialised JSON "
             "representation.")
        .staticmethod("from_json")
    ;

    py::class_<mh_moves::move_t> ( "Move", "Biggles Move Function" )
        .def(py::init<mh_moves::move_type>())
        .def("__call__", &exec_move)
        .def("__repr__", &mh_moves::move_t::repr)
        .def("__str__", &mh_moves::move_t::str)
        .add_property("type", &mh_moves::move_t::type, &mh_moves::move_t::set_type, "Type of the move.")
    ;

    py::class_<track, shared_ptr_hack<track> >(
        "Track",
        "An ordered set of observations with a specified birth and death time."
        )
        .def(py::init<const track &>())
        .def(py::init<const track &, const track &, float>())
        .def("__len__", &track::size)
        .def("__repr__", &track_to_str)
        .def(py::self == py::self)
        .add_property("duration", &track::duration, "The number of time points.")
        .add_property("observations", py::make_function(&track::observations, py::return_internal_reference<1>()),
            "The observations of the track.")
        .add_property("first_time_stamp", &track::first_time_stamp, "The first time stamp of the track.")
        .add_property("last_time_stamp", &track::last_time_stamp,
            "The the first time stamp after the track has died.")
        .def("make_kalman_filter", &track::make_kalman_filter,
             "Construct a Kalman filter instance for this track and a given set of parameters.")
        .def("insert", &track::insert, "Inserts an observation into a track")
        .def("time_stamps_from_observations", &track::time_stamps_from_observations,
            "Copies track time stamps from its observations. Returns nothing.")
        .def("new", &track_from_time_stamps_and_obs,
            "Construct a new track from time stamps and observation collection.")
        .staticmethod("new")
    ;

    py::class_<track_collection>(
        "TrackCollection",
        "A collection of :py:class:`Track`s which is iterable. "
        "The collection is *immutable* and does not support indexing."
        )
        .def("__len__", &track_collection::size)
        .def("__iter__", py::iterator<track_collection>())
        .def("__repr__", &track_coll_to_string)
        .def("__getitem__", &get_track)
        .def("insert", &insert_track_into_collection,
            "Insert a track into the track collection")
        .def("remove", &remove_track_from_collection,
            "Remove a track from the track collection")
    ;
    py::register_ptr_to_python<boost::shared_ptr<const track> >();

    py::class_<observation>(
        "Observation",
        "A single observation within a frame. It has an x and y spatial co-ordinate and a single temporal "
        "co=ordinate which is usually the index of the frame in which it appears."
        )
        .def("__repr__", &obs_to_string)
        .add_property("x", &get_obs_x, &set_obs_x)
        .add_property("y", &get_obs_y, &set_obs_y)
        .add_property("t", &get_obs_t, &set_obs_t)
        .def("new", &obs_from_coordinates, "constructs a new observation from coordiantes")
        .staticmethod("new")
    ;

    py::class_<observation_collection>(
        "ObservationCollection",
        "A collection of :py:class:`Observation`s which is iterable. "
        "The collection does not support indexing but may be mutated by way of the :py:meth:`clear` method."
        )
        .def("__len__", &observation_collection::size)
        .def("__iter__", py::iterator<observation_collection>())
        .def("__repr__", &obs_coll_to_string)
        .def("__eq__", &obs_coll_are_equal)

        .def("clear", &observation_collection::clear)

        .add_property("first_time_stamp", &observation_collection::first_time_stamp,
                      "The first time stamp which is represented in the collection.")
        .add_property("last_time_stamp", &observation_collection::last_time_stamp,
                      "The first time stamp *after* the latest represented in the collection.")
        .def("remove", &observation_collection::erase, "remove observation")
        .def("insert", &insert_obs_into_collection, "insert observation")
        .def("count_at", &observation_collection::count_at_time_stamp, "number of observations at time stamp")
        .def("contains", &observation_collection::contains, "check if the collection contains an observation")
        .def("obs_at", &obs_coll_obs_at_ts, "the observations at a time stamp")
    ;

    py::class_<ObservationReservoir, boost::shared_ptr<ObservationReservoir> >("_ObservationPool",
        "A collection of :py:class:`Observation`s. Use 'observation_pool' to initialise"
        )
        .def("__len__", &ObservationReservoir::size)
        .def("__getitem__", &ObservationReservoir::operator[], py::arg( "index" ), py::return_internal_reference<>())
        .def("__iter__", py::iterator<ObservationReservoir>())
        .def("push", push_coords, "Add observation by coordinates")
    ;

    py::class_<indexed_partition_t>("IndexedPartition", "A indexed form of the partition")
        .def("add_track", &indexed_partition_t::add_track, "Adds a single indexed track to the partition")
        .def("push_track", &indexed_partition_t::push_track, "Adds a single track from a python 'dict'")
        .def("set_tracks", &indexed_partition_t::set_tracks, "Sets the tracks.")
        .def("clutter", &indexed_partition_t::clutter)
        .def("from_partition", &indexed_partition_t::from_partition)
        .def("set_clutter", &indexed_partition_t::set_clutter, "Sets the clutter.")
        .def("to_partition", &indexed_partition_t::to_partition)
        .def("tracks", &indexed_partition_t::tracks)
    ;

    py::class_<partition, boost::shared_ptr<partition> >(
        "Partition",
        "A set of spurious observations (*clutter*) and tracks which form a particular tracking solution."
        )
        .def(py::init<const partition &>())
        .def("__repr__", &partition_to_string)
        .add_property("first_time_stamp", &partition::first_time_stamp)
        .add_property("last_time_stamp", &partition::last_time_stamp)
        .add_property("duration", &partition::duration)
        .add_property("clutter", py::make_function(&partition::clutter, py::return_internal_reference<1>()))
        .add_property("tracks", py::make_function(&partition::tracks, py::return_internal_reference<1>()))
        .def("to_json", &partition_to_json,
             "Serialise the partition to a string containing a JSON py::object.")
        .def("from_json", &partition_from_json,
             "Construct a new :py:class:`Partition` py::object from a string containing a serialised JSON "
             "representation.")
        //.def("observation_pool", &partition::pool, "The observation pool of the partition")
        .def("minimise_clutter", &maximum_creation, "Creates a new partition that has put as much clutter as "
             "possible into tracks.")
        .def("new", &partition_from_collections,
             "Construct a new partition from a track collection and clutter")
        .def("set_expansion", &partition::set_expansion,
            "Set the partition expansion x-min, xmax, y-min, ymax")
        .def("volume", &partition::volume,
            "The spatial size of the partition in squared length unit")
        .def("observation_count", &partition::observation_count,
            "return the number of observations in the partition")
        .def("as_string", &partition::as_string, "returns a long and unique string representation of the partition")
        .def("ged", &graph_edit_distance, "Graph edit distance to another partition")
        .def("recorder_size", &recorder_size, "Size of the capability recorder")
        .staticmethod("new")
        .staticmethod("from_json")
    ;

    py::class_<server::tracking_state, boost::shared_ptr<server::tracking_state> >(
        "TrackingState",
        "An object representing the current state of a tracking job along with the current 'best' (or "
        "modal) samples for the partition and model parameters."
        )
        .def_readwrite("sample_count", &server::tracking_state::sample_count)
        .def_readwrite("current_log_pdf", &server::tracking_state::current_log_pdf)
        .def_readwrite("current_model_parameters", &server::tracking_state::current_model_parameters)
        //.def_readwrite("current_partition", &server::tracking_state::current_partition)
        .def_readwrite("current_move_type", &server::tracking_state::current_move_type)
        .def_readwrite("best_log_pdf", &server::tracking_state::best_log_pdf)
        .def_readwrite("best_model_parameters", &server::tracking_state::best_model_parameters)
        //.def_readwrite("best_partition", &server::tracking_state::best_partition)
        .def_readwrite("acceptance_rate", &server::tracking_state::acceptance_rate,
            "The acceptance rate of the *all* Metropolis-Hastings samples")
        .def_readwrite("recent_acceptance_rate", &server::tracking_state::recent_acceptance_rate,
            "The acceptance rate of the *latest* Metropolis-Hastings samples")
        .def_readwrite("accepted", &server::tracking_state::accepted,
            "Was the last partition proposal accepted")
        //.def_readwrite("last_proposed_partition", &server::tracking_state::last_proposed_partition)
        .def_readwrite("last_proposed_move", &server::tracking_state::last_proposed_move)
        .def_readwrite("last_proposal_density", &server::tracking_state::last_proposal_density)
        .def_readwrite("last_sample_density", &server::tracking_state::last_sample_density)
        .def_readwrite("last_pdr", &server::tracking_state::last_pdr)
        .def_readwrite("move_histogram", &server::tracking_state::move_histogram)
        .def_readwrite("moves_accepted", &server::tracking_state::moves_accepted)
        .def_readwrite("moves_rejected", &server::tracking_state::moves_rejected)
        .def_readwrite("moves_identity", &server::tracking_state::moves_identity)
        .def_readonly("internals", &server::tracking_state::internals)
        .add_property("current_partition", &server::tracking_state::current_partition_fun)
        .add_property("best_partition", &server::tracking_state::best_partition_fun)
        .add_property("last_proposed_partition", &server::tracking_state::last_proposed_partition_fun)
    ;

    py::class_<server::engine>(
        "Engine",
        "A threaded tracking engine. It is initialised via the :py:meth:`set_initial_conditions` method. "
        "Once started tracking is performed in a separate thread. The current state of the tracker may by "
        "retrieved via the :py:attr:`tracking_state` attribute.")
        .def("set_initial_conditions", &server::engine::set_initial_conditions_part)
        .def("start", &server::engine::start)
        .def("stop", &server::engine::stop)
        .add_property("tracking_state", &server::engine::tracking_state_ptr)
    ;

    py::class_<server::stepper>(
        "Stepper",
        "A non-threaded step-by-step sampler. It is initialised via the constructor. "
        "The current state of the tracker may by retrieved via the :py:attr:`tracking_state` attribute.")
        .def(py::init<const model::parameters &, partition &>())
        .def("step", &server::stepper::step,
            "draw a number of (partition, parameter) samples. Default number of steps = 1.")
        .add_property("tracking_state", &server::stepper::tracking_state_ptr)
        .def("fix_observation_error", &server::stepper::fix_observation_error,
            "set the observation error to a fixed value and stop sampling it")
        .def("sample_observation_error", &server::stepper::sample_observation_error,
            "Sample the observation error every *lag* samples")
    ;
    py::class_<clutter_t>("Clutter","The clutter type")
        .def("__len__", size0)
        .def("observations", all_observations)
        .def("insert", insert0)
    ;

    py::def("log_tracks_density", model::log_tracks_given_parameters_density);
    py::def("log_clutter_density", model::log_clutter_given_parameters_density);
    py::def("log_partition_density", model::log_partition_given_parameters_density);
    py::def("log_prior_density", model::log_parameters_prior_density);
    //py::def("log_posterior_density", model::log_partition_given_parameters_and_data_density);
    py::def("sample_model_parameters_given_partition", sample_model_parameters_given_partition);

    py::def("move_name", mh_moves::move_name, "The name of the move type");
    py::def("move_sign", mh_moves::move_sign, "A 3 character represenatation of the move name");

    py::def("observation_pool", create_observation_pool, "Create an empty observation pool.");
    py::def("log_pdf", model::log_partition_given_parameters_and_data_density, "Calculate the log PDF");
}
