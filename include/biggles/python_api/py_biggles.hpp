#ifndef PY_BIGGLES_HPP__
#define PY_BIGGLES_HPP__

#include <deque>
#include <vector>
#include <sstream>

#include <boost/python.hpp>
#include <boost/foreach.hpp>

#include "biggles/io.hpp"
#include "biggles/kalman_filter.hpp"
#include "biggles/tracker.hpp"
#include "biggles/server/engine.hpp"
#include "biggles/server/stepper.hpp"
#include "biggles/tools/max_creation.hpp"
#include "biggles/tools/partition_graph_edit_distance.hpp"
#include "biggles/tools/clutter_transform.hpp"
#include "biggles/tools/track_coll_transform.hpp"
//#include "biggles/tools/sundries.hpp"

using namespace biggles;
//using namespace boost::python;
namespace py=boost::python;

#define FOREACH BOOST_FOREACH

namespace biggles {
namespace python {

// from http://wiki.python.org/moin/boost.python/StlContainers
inline void IndexError() { PyErr_SetString(PyExc_IndexError, "Index out of range"); }
template<class T>
struct const_std_random_access
{
    typedef typename T::value_type V;
    static const V& get(T const& x, size_t i)
    {
        if( i<0 ) i+=x.size();
        return x[i];
    }
};

template<class T>
std::deque<T> deque_from_list(const py::list& ns) {
    std::deque<T> result;
    for (int i = 0; i < len(ns); ++i) {
        result.push_back(py::extract<T>(ns[i]));
    }
    return result;
}

class indexed_partition_t {
    clutter_index_t clutter_;
    std::deque<indexed_track> tracks_;
    py::dict convert_track(const indexed_track& idx_track) const {
        py::dict result;
        result["dynamic_drag"] = 1.0f;
        py::list time_span;
        time_span.append(idx_track.first_time_stamp);
        time_span.append(idx_track.last_time_stamp);
        result["time_span"] = time_span;
        py::list observations;
        FOREACH(size_t s, idx_track.observations) observations.append(s);
        result["observations"] = observations;
        return result;
    }

public:
    void from_partition(const partition& source_part) {
        clutter_observation_to_index(source_part.pool(), source_part.clutter(), clutter_);
        track_coll_observation_to_index(source_part.pool(), source_part.tracks(), tracks_);
    }
    void add_track(time_stamp first_ts, time_stamp last_ts, const py::list& obs_idx) {
        indexed_track new_track;
        new_track.first_time_stamp = first_ts;
        new_track.last_time_stamp = last_ts;
        new_track.observations = deque_from_list<size_t>(obs_idx);
        tracks_.push_back(new_track);
    }
    void set_clutter(const py::list& obs_idx) {
        clutter_ = deque_from_list<size_t>(obs_idx);
    }
    void push_track(const py::dict& idx_track) {
        py::list time_span = py::extract<py::list>(idx_track["time_span"]);
        time_stamp first_ts = py::extract<time_stamp>(time_span[0]);
        time_stamp last_ts = py::extract<time_stamp>(time_span[1]);
        add_track(first_ts, last_ts, py::extract<py::list>(idx_track["observations"]));
    }
    void set_tracks(const py::list& idx_tracks) {
        for (int i = 0; i < len(idx_tracks); ++i) {
            push_track(py::extract<py::dict>(idx_tracks[i]));
        }
    }
    partition to_partition(const boost::shared_ptr<ObservationReservoir>& observation_pool) {
        clutter_ptr clutter_p(new clutter_t());
        boost::shared_ptr<track_collection> tracks_p(new track_collection);
        clutter_index_to_observation(observation_pool, clutter_, *clutter_p);
        track_coll_index_to_observation(observation_pool, tracks_, *tracks_p);
        return partition(observation_pool, tracks_p, clutter_p);
    }
    py::list clutter() const {
        py::list ret;
        FOREACH(size_t s, clutter_) ret.append(s);
        return ret;
    }
    py::list tracks() const {
        py::list ret;
        FOREACH(const indexed_track& tr, tracks_) ret.append(convert_track(tr));
        return ret;
    }
};


class key_error : public std::runtime_error {
public:
    explicit key_error( const std::string& what_arg ) : std::runtime_error(what_arg) {}
};
class index_error : public std::runtime_error {
public:
    explicit index_error( const std::string& what_arg ) : std::runtime_error(what_arg) {}
};

/*
/// \brief a proxy for the model::parameters observation_error_covariance
class observation_covariance {
    float _00;
    float _01;
    float _11;
public:
    float a00() const { return _00; }
    float a01() const { return _01; }
    float a11() const { return _11; }
    void set_a00(const float value) { _00 = value; }
    void set_a01(const float value) { _01 = value; }
    void set_a11(const float value) { _11 = value; }
    py::tuple get(size_t index) const {
        if (index > 1)
            throw index_error("list index out of range");
        if (index == 0)
            return py::make_tuple(_00, _01);
        return py::make_tuple(_01, _11);
    }
};
*/

/*
/// \brief Convert observations from the ObservationReservoir into a python list [[x, y, t], ...]
py::list extract_observations_from_pool(const ObservationReservoir& pool) {
    py::list result;
    FOREACH (observation o, pool) result.append(py::make_tuple(x(o), y(o), t(o)));
    return result;
}
*/

struct observation_to_tuple_converter {
    static PyObject* convert(const observation& o) {
        return py::incref(py::make_tuple(x(o), y(o), t(o)).ptr());
    }
};

template<typename T>
struct Deque_to_python_list {

  static PyObject* convert(std::deque<T> const& v) {
    py::list l;
    typename std::deque<T>::const_iterator p;
    for(p=v.begin();p!=v.end();++p){
      l.append(object(*p));
    }
    return py::incref(l.ptr());
  }
};

py::tuple exec_move(const mh_moves::move_t& move, const partition& start_partition) {
    mh_moves::move_t::result_t proposal = move.exec(start_partition);
    if (proposal.success)
        return py::make_tuple(proposal.success, proposal.pmr, *proposal.proposed_partition);
    else
        return py::make_tuple(proposal.success, proposal.pmr, partition());
}

size_t recorder_size(const partition& start_partition) {
    if (not start_partition.get_capability_recorder()) {
        //std::cout << "nothing" << std::endl;
        return 0;
    }
    return start_partition.get_capability_recorder()->size();
}

/*
/// \brief biggles move
struct move {
    mh_moves::move_fun _move; ///< \brief the actual function to be called
    mh_moves::move_type _type;
    move() : _move(0), _type(mh_moves::NONE) {}
    move(mh_moves::move_type type) : _move(mh_moves::get_move(type)), _type(type) {}
    py::tuple exec(const partition& start_partition) {
        if (_type == mh_moves::NONE)
            throw std::runtime_error("move type not defined");
        partition_ptr_t start_part_ptr(new partition(start_partition));
        partition_ptr_t end_part_ptr;
        capability_recorder_ptr cap_rec_ptr = new_cap_recorder(*start_part_ptr);
        start_part_ptr->set_capability_recorder(cap_rec_ptr);
        float pdr = 0.f;
        bool result = _move(start_part_ptr, end_part_ptr, pdr);
        if (result)
            return py::make_tuple(result, pdr, *end_part_ptr);
        else
            return py::make_tuple(result, pdr, partition());
    }
    void set_type(mh_moves::move_type type) {
        _move = mh_moves::get_move(type);
        _type = type;
    }
    mh_moves::move_type type() const { return _type; }
    std::string repr() const { return std::string("biggles.") + mh_moves::move_name(_type) + "Move"; }
    std::string str() const { return mh_moves::move_name(_type); }
};
*/

shared_const_track_ptr get_track(const track_collection &tc, size_t i) {
    if (i>= tc.size())
        throw key_error("index out of range");
    track_collection::iterator it = tc.begin();
    std::advance(it, i);
    return *it;
}

shared_const_track_ptr remove_track_from_collection(track_collection &tc, size_t i) {
    if (i>= tc.size())
        throw key_error("index out of range");
    track_collection::iterator it = tc.begin();
    std::advance(it, i);
    shared_const_track_ptr it2(*it);
    tc.remove(*it);
    return it2;
}

void translate_key_error(const key_error &e) {  PyErr_SetString(PyExc_KeyError, e.what()); }
void translate_index_error(const index_error &e) {  PyErr_SetString(PyExc_IndexError, e.what()); }

void insert_track_into_collection(track_collection &tc, const track &t) { tc.insert(t); }
void insert_obs_into_collection(observation_collection &oc, const observation &o) { oc.insert(o); }

inline bool obs_coll_are_equal(const observation_collection& oc1, const observation_collection& oc2) {
    time_stamp first_ts = oc1.first_time_stamp();
    time_stamp last_ts = oc1.last_time_stamp();
    if (first_ts != oc2.first_time_stamp()) return false;
    if (last_ts != oc2.last_time_stamp()) return false;
    observation_collection::const_range obs1_at_t, obs2_at_t;
    typedef std::set<observation> obs_set;
    for (; first_ts < last_ts; ++first_ts) {
        obs1_at_t = oc1.at_time_stamp(first_ts);
        obs2_at_t = oc2.at_time_stamp(first_ts);
        if (obs_set(obs1_at_t.first, obs1_at_t.second) != obs_set(obs2_at_t.first, obs2_at_t.second)) return false;
    }
    return true;
}

// pseudo constructor for track
track track_from_time_stamps_and_obs(const time_stamp &first_ts, const time_stamp &last_ts,
    const observation_collection& oc)
{
    return track(first_ts, last_ts, oc.begin(), oc.end(), 1.0f);
}

std::vector<observation> obs_coll_obs_at_ts(const observation_collection& oc, const time_stamp ts) {
    observation_collection::const_range cr(oc.at_time_stamp(ts));
    return std::vector<observation>(cr.first, cr.second);
}

// accessors for model parameter properties
float get_mntpf(const model::parameters& p) { return model::mean_new_tracks_per_frame(p); }
void set_mntpf(model::parameters& p, float v) { model::mean_new_tracks_per_frame(p) = v; }
float get_mfopf(const model::parameters& p) { return model::mean_false_observations_per_frame(p); }
void set_mfopf(model::parameters& p, float v) { model::mean_false_observations_per_frame(p) = v; }
float get_ftfsp(const model::parameters& p) { return model::frame_to_frame_survival_probability(p); }
void set_ftfsp(model::parameters& p, float v) { model::frame_to_frame_survival_probability(p) = v; }
float get_gop(const model::parameters& p) { return model::generate_observation_probability(p); }
void set_gop(model::parameters& p, float v) { model::generate_observation_probability(p) = v; }

void init_process_cov(model::parameters& p, py::list& oc) {
    model::process_noise_covariance(p) = detail::initQ();
    for (int i = 0; i < 4; ++i) {
        float item = py::extract<float>(oc[i]);
        model::process_noise_covariance(p)(i, i) = item * detail::speed_of_light_ * item * detail::speed_of_light_;
    }
}

void set_process_cov(model::parameters& p, py::list& oc) {
    matrix4f Q;
    for (int row = 0; row < 4; row++) {
        for (int col = 0; col < 4; col++) {
            Q(row, col) = py::extract<float>(oc[row][col]);
        }
    }
    model::process_noise_covariance(p) = Q;
}

py::list get_process_cov(const model::parameters& p) {
    const matrix4f &q(model::process_noise_covariance(p));
    py::list pc;
    for (int row = 0; row < 4; row++) {
        py::list rlist;
        for (int col = 0; col < 4; col++) rlist.append(q(row, col));
        pc.append(rlist);
    }
    return pc;
}

void set_obs_err(model::parameters& p, py::list& oc) {
    model::observation_error_covariance(p) <<
        py::extract<float>(oc[0][0]),
        py::extract<float>(oc[0][1]),
        py::extract<float>(oc[1][0]),
        py::extract<float>(oc[1][1]);
}
py::list get_obs_err(const model::parameters& p) {
    py::list oc;
    py::list row0;
    py::list row1;
    row0.append(model::observation_error_covariance(p)(0,0));
    row0.append(model::observation_error_covariance(p)(0,1));
    row1.append(model::observation_error_covariance(p)(1,0));
    row1.append(model::observation_error_covariance(p)(1,1));
    oc.append(row0);
    oc.append(row1);
    return oc;
}

// conversion to/from JSON
std::string model_parameters_to_json(const model::parameters& p)
{
    std::stringstream ss;
    io::write_model_parameters_to_json_stream(ss, p);
    return ss.str();
}

boost::shared_ptr<model::parameters> model_parameters_from_json(const std::string& s)
{
    std::stringstream ss(s);
    return io::read_model_parameters_from_json_stream(ss);
}


//bool equal_partitions(const partition& p, const partition& q) { return p == q; }

std::string partition_to_json(const partition& p)
{
    std::stringstream ss;
    io::write_partition_to_json_stream(ss, p);
    return ss.str();
}

float graph_edit_distance(const partition& p, const partition& q) { return partition_ged(p, q); }

partition_ptr_t maximum_creation(const partition& p) {
    partition_ptr_t target_partition_ptr;
    partition_ptr_t origin_partition_ptr(new partition(p));
    max_creation(origin_partition_ptr, target_partition_ptr);
    return target_partition_ptr;
}

boost::shared_ptr<partition> partition_from_json(const std::string& s)
{
    std::stringstream ss(s);
    return io::read_partition_from_json_stream(ss);
}

track_collection copy_track_collection(const track_collection& tracks) {
    track_collection tc;
    BOOST_FOREACH(const boost::shared_ptr<const track>& t, tracks) {
        tc.insert(track(t->first_time_stamp(), t->last_time_stamp(),
            t->observations().begin(), t->observations().end(), 1.f));
    }
    return tc;
}

partition partition_from_collections(const track_collection& tracks, const py::list& clutter_list) {
    track_collection tc;
    BOOST_FOREACH(const boost::shared_ptr<const track>& t, tracks) {
        tc.insert(track(t->first_time_stamp(), t->last_time_stamp(),
            t->observations().begin(), t->observations().end(), 1.f));
    }
    partition::track_collection_ptr tracks_ptr(new track_collection(tc));
    std::deque<observation> clutter(deque_from_list<observation>(clutter_list));
    clutter_ptr clutter_ptr(new clutter_t(clutter.begin(), clutter.end()));
    return partition(tracks_ptr, clutter_ptr, expansion_from_observations(*tracks_ptr, *clutter_ptr));
}



// pseudo constructor
observation obs_from_coordinates(const float& x_, const float& y_, const time_stamp & t_) {
    return new_obs(x_, y_, t_);
}

std::string track_to_str(const track& t) {
    std::stringstream ss;
    ss << "biggles.Track( time [" << t.first_time_stamp() << ", " << t.last_time_stamp() << "), "
        << t.size()
        << " observations [" << t.observations().first_time_stamp() << ", "
        << t.observations().last_time_stamp() << ") )";
    return ss.str();
}

std::string obs_to_string(const observation& o) {
    std::stringstream ss;
    ss << "biggles.Observation(" << x(o) << ", " << y(o) << ", " << t(o) << ")";
    return ss.str();
}

std::string obs_coll_to_string(const observation_collection& oc) {
    std::stringstream ss;
    ss << "biggles.ObservationCollection( time [" << oc.first_time_stamp() << ", " << oc.last_time_stamp()
        << "), " << oc.size() << " observations )";
    return ss.str();
}

std::string partition_to_string(const partition& p) {
    std::stringstream ss;
    ss << "biggles.Partition( time [" << p.first_time_stamp() << ", " << p.last_time_stamp()
        << "), " << p.tracks().size() << " tracks, " << p.clutter().size() << " clutter observations )";
    return ss.str();
}

std::string track_coll_to_string(const track_collection& tc) {
    std::stringstream ss;
    ss << "biggles.TrackCollection( " << tc.size() << " tracks)";
    return ss.str();
}

// accessors for observation properties
float get_obs_x(const observation& o) { return x(o); }
float get_obs_y(const observation& o) { return y(o); }
time_stamp get_obs_t(const observation& o) { return t(o); }

// accessors for observation properties
void set_obs_x(observation& o, float val)      { o->set_x(val); }
void set_obs_y(observation& o, float val)      { o->set_y(val); }
void set_obs_t(observation& o, time_stamp val) { o->set_t(val); }

boost::shared_ptr<ObservationReservoir> create_observation_pool() {
    return boost::shared_ptr<ObservationReservoir>(new ObservationReservoir);
}

// HACK: a nasty shared_ptr sub-class which implicitly removes the 'const-ness' when initialised
// with a const pointer. Used to work around boost::python not supporting boost::shared_ptr<const T>.
template<typename T>
struct shared_ptr_hack : public boost::shared_ptr<T> {
    typedef boost::shared_ptr<T> base;

    // default and copy constructors
    shared_ptr_hack() : base() { }
    shared_ptr_hack(const shared_ptr_hack& h) : base(h) { }

    // dangerous direct pointer constructors
    shared_ptr_hack(const T* t) : base(const_cast<T*>(t)) { }
    shared_ptr_hack(T* t) : base(t) { }

    // dirty const-removing constructor
    shared_ptr_hack(const boost::shared_ptr<const T>* p) : base(boost::const_pointer_cast<T>(p)) { }
};

// This is an ungoldy horrible hack which explicitly re-implements some boost
// internals in order to support implicit const-casting of boost::shared_ptr.
// This appears to be the only approach on earlier boosts.
#if BOOST_VERSION < 15000
namespace boost{ namespace python { namespace objects {

template <>
void* pointer_holder<boost::shared_ptr<const track>, const track>::holds(type_info dst_t, bool null_ptr_only)
{
    typedef boost::shared_ptr<const track> Pointer;
    typedef track Value;

    if (dst_t == python::type_id<Pointer>()
        && !(null_ptr_only && get_pointer(this->m_p))
    )
        return &this->m_p;

    Value* p
#  if BOOST_WORKAROUND(__SUNPRO_CC, BOOST_TESTED_AT(0x590))
        = static_cast<Value*>( get_pointer(this->m_p) )
#  else
        = const_cast<Value*>( get_pointer(this->m_p) )
#  endif
        ;

    if (p == 0)
        return 0;

    if (void* wrapped = holds_wrapped(dst_t, p, p))
        return wrapped;

    type_info src_t = python::type_id<Value>();
    return src_t == dst_t ? p : find_dynamic_type(p, src_t, dst_t);
}

template <>
void* pointer_holder<boost::shared_ptr<const kalman_filter>, const kalman_filter>::holds(type_info dst_t, bool null_ptr_only)
{
    typedef boost::shared_ptr<const kalman_filter> Pointer;
    typedef kalman_filter Value;

    if (dst_t == python::type_id<Pointer>()
        && !(null_ptr_only && get_pointer(this->m_p))
    )
        return &this->m_p;

    Value* p
#  if BOOST_WORKAROUND(__SUNPRO_CC, BOOST_TESTED_AT(0x590))
        = static_cast<Value*>( get_pointer(this->m_p) )
#  else
        = const_cast<Value*>( get_pointer(this->m_p) )
#  endif
        ;

    if (p == 0)
        return 0;

    if (void* wrapped = holds_wrapped(dst_t, p, p))
        return wrapped;

    type_info src_t = python::type_id<Value>();
    return src_t == dst_t ? p : find_dynamic_type(p, src_t, dst_t);
}

} } }
#endif

// accessors for kalman filter state and covariance pair
py::object get_kfs(const kalman_filter::state_covariance_pair& sc) {
    const kalman_filter::state_vector& sv(sc.first);
    return py::make_tuple(sv[0], sv[1], sv[2], sv[3]);
}

py::object get_kfc(const kalman_filter::state_covariance_pair& sc) {
    const kalman_filter::covariance_matrix& cm(sc.second);
    return py::make_tuple(
        py::make_tuple(cm(0,0), cm(0,1), cm(0,2), cm(0,3)),
        py::make_tuple(cm(1,0), cm(1,1), cm(1,2), cm(1,3)),
        py::make_tuple(cm(2,0), cm(2,1), cm(2,2), cm(2,3)),
        py::make_tuple(cm(3,0), cm(3,1), cm(3,2), cm(3,3))
    );
}

/*
py::object get_proc_noise(const kalman_filter &kf) {
    const kalman_filter::covariance_matrix &q(kf.Q());
    return py::make_tuple(
        py::make_tuple(q(0,0), q(0,1), q(0,2), q(0,3)),
        py::make_tuple(q(1,0), q(1,1), q(1,2), q(1,3)),
        py::make_tuple(q(2,0), q(2,1), q(2,2), q(2,3)),
        py::make_tuple(q(3,0), q(3,1), q(3,2), q(3,3))
    );
}
*/

size_t (clutter_t::*size0)() const = &clutter_t::size;
void (clutter_t::*insert0)(const observation& o) = &clutter_t::insert;
void (ObservationReservoir::*push_coords)(const float, const float, const time_stamp) = &ObservationReservoir::push;
clutter_t::observation_container (clutter_t::*all_observations)() const = &clutter_t::observations;

} } // namespace biggles::python

#endif // PY_BIGGLES_HPP__
