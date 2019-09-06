/// @file utility.hpp Various utility functions common to the test programs

#ifndef BIGGLES_TEST_UTILITY_HPP__
#define BIGGLES_TEST_UTILITY_HPP__

#include <boost/program_options.hpp>
#include <boost/format.hpp>

#include "../partition.hpp"
#include "../observation_collection.hpp"
#include "../observation.hpp"
#include "../track_collection.hpp"
#include "../track.hpp"
#include "../detail/random.hpp"
#include "../mh_moves/mh_moves.hpp"
#include "../mh_moves/utility.hpp"
#include "../tools/partition_graph_edit_distance.hpp"
#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include "debug.hpp"

const std::string DEF("\033[0m");
const std::string BOLD("\033[1m");

const std::string BLACK("\033[30m");
const std::string RED("\033[31m");
const std::string GREEN("\033[32m");
const std::string BROWN("\033[33m");
const std::string BLUE("\033[34m");
const std::string MAGENTA("\033[35m");
const std::string CYAN("\033[36m");
const std::string WHITE("\033[37m");
const std::string DEF_FG("\033[22;39m");

const std::string BLACK_BG("\033[40m");
const std::string RED_BG("\033[41m");
const std::string GREEN_BG("\033[42m");
const std::string BROWN_BG("\033[43m");
const std::string BLUE_BG("\033[44m");
const std::string MAGENTA_BG("\033[45m");
const std::string CYAN_BG("\033[46m");
const std::string WHITE_BG("\033[47m");
const std::string DEF_BG("\033[49m");

typedef std::deque<float>::const_iterator data_iter;
inline std::pair<float, float> mean_stdev(data_iter b, data_iter e) {
    float n = 0.0f;
    float mean = 0.0f;
    float m2 = 0.0;
    for (data_iter it = b; it != e; ++it) {
        float delta = *it - mean;
        mean += delta / ++n;
        m2 += delta*(*it - mean);
    }
    return std::make_pair(mean, std::sqrt(m2/(n-1)));
}

/// \brief error on the probability of the binomial distribution, if estimated through sampling
float binom_prob_error_est(size_t successes, size_t n_samples) {
    float prob = float(successes)/float(n_samples);
    return std::sqrt( prob*(1.f- prob)/float(n_samples));
}

/// \brief the error of the fraction \e est_numer / \e est_denom if their errors are \e err_numer resp. \e est_denom
float error_of_fraction(float est_numer, float err_numer, float est_denom, float err_denom) {
    float est_frac = est_numer/est_denom;
    return std::abs(est_frac) * ( err_numer/std::abs(est_numer) +  err_denom/std::abs(est_denom) );
}

float auto_correlation(data_iter b, data_iter e, size_t lag = 1);
float correlation(data_iter b1, data_iter e1, data_iter b2);

typedef std::map<std::string, bool> test_t;

/// @brief Write a partition out to disk in GNU plotutils and MATLAB-friendly ways.
///
/// MATLAB data is stored as a table with each row being observation x- and y-co-ordinate, time stamp and track id.
/// Track id 0 corresponds to the clutter.
///
/// @param prefix
/// @param p
/// @param name
inline void write_partition(const char* prefix, const biggles::partition& p, const char* name)
{
    using namespace biggles;

    std::ofstream ds_xy((std::string("biggles_") + prefix + std::string("_") + name + std::string("_xy.txt")).c_str());
    std::ofstream ds_tx((std::string("biggles_") + prefix + std::string("_") + name + std::string("_tx.txt")).c_str());
    std::ofstream ds_ty((std::string("biggles_") + prefix + std::string("_") + name + std::string("_ty.txt")).c_str());
    std::ofstream ds_matlab((std::string("biggles_") + prefix + std::string("_") + name + std::string("_matlab.txt")).c_str());

    int m(0);
    BOOST_FOREACH(const boost::shared_ptr<const biggles::track>& p_t, p.tracks())
    {
        ds_xy << "#m=" << 1+(m%5) << ",S=1\n";
        ds_tx << "#m=" << 1+(m%5) << ",S=1\n";
        ds_ty << "#m=" << 1+(m%5) << ",S=1\n";

        ++m;

        BOOST_FOREACH(const biggles::observation& o, *p_t)
        {
            ds_xy << std::setw(16) << x(o) << std::setw(16) << y(o) << '\n';
            ds_tx << std::setw(16) << t(o) << std::setw(16) << x(o) << '\n';
            ds_ty << std::setw(16) << t(o) << std::setw(16) << y(o) << '\n';
            ds_matlab << std::setw(16) << x(o)
                      << std::setw(16) << y(o)
                      << std::setw(16) << t(o)
                      << std::setw(16) << m
                      << '\n';
        }

        ds_xy << '\n';
        ds_tx << '\n';
        ds_ty << '\n';
    }

    ds_xy << "#m=0,S=1\n";
    ds_tx << "#m=0,S=1\n";
    ds_ty << "#m=0,S=1\n";

    for (clutter_t::const_iterator it = p.clutter().begin(); it not_eq p.clutter().end(); ++it) {
        BOOST_FOREACH(const biggles::observation& o, it->second)
        {
            ds_xy << std::setw(16) << x(o) << std::setw(16) << y(o) << '\n';
            ds_tx << std::setw(16) << t(o) << std::setw(16) << x(o) << '\n';
            ds_ty << std::setw(16) << t(o) << std::setw(16) << y(o) << '\n';

            ds_matlab << std::setw(16) << x(o)
                      << std::setw(16) << y(o)
                      << std::setw(16) << t(o)
                      << std::setw(16) << 0
                      << '\n';
        }
    }
}


/// \brief coverts a track into a string (unpolished method)
std::string track_to_str(const biggles::track & track0, biggles::time_stamp end_ts = 0);

/// \brief Fills params with some values
inline bool generate_parameters(biggles::model::parameters& params) {
    using namespace biggles;
    model::mean_new_tracks_per_frame(params)           = 0.2f;
    model::mean_false_observations_per_frame(params)   = 0.2f;
    model::frame_to_frame_survival_probability(params) = 0.8f;
    model::generate_observation_probability(params)    = 0.8f;
    float sigma_R = 0.3f;
    model::observation_error_covariance(params)       << sigma_R * sigma_R, 0, 0, sigma_R * sigma_R;
    model::constraint_radius(params) = 1.0f;
    model::process_noise_covariance(params) = detail::initQ();
    return true;
}

/// \brief Finds the best parameters for partition from 'nloops' samples
inline float find_good_parameters(const biggles::partition& const_partition,
    biggles::model::parameters& good_parameters, const int nloops) {
    using namespace biggles;
    float best_pdf = -std::numeric_limits<int>::max();
    model::parameters target_parameters;
    generate_parameters(target_parameters);
    for (int i = 0; i < nloops; ++i) {
        sample_model_parameters_given_partition(const_partition, target_parameters);
        float log_pdf = model::log_partition_given_parameters_and_data_density(const_partition, target_parameters);
        if (log_pdf > best_pdf) {
            best_pdf = log_pdf;
            good_parameters = target_parameters;
        }
    }
    return best_pdf;
}

inline std::string partition_to_string(const biggles::partition& part) {
    biggles::track_collection::const_iterator it1;
    std::deque<std::string> track_str;
    biggles::time_stamp final = part.last_time_stamp();
    biggles::track_collection tc1 = part.tracks();
    for (it1 = tc1.begin(); it1 not_eq tc1.end(); ++it1)
        track_str.push_back(track_to_str(**it1, final));
    std::sort(track_str.begin(), track_str.end());
    std::stringstream res;
    res << "|";
    for (size_t i = 0; i < track_str.size(); ++i) {
        res << track_str[i] << "|";
    }
    return res.str();
}

inline std::deque<std::string> partition_to_string_list(const biggles::partition& part) {
    biggles::track_collection::const_iterator it1;
    std::deque<std::string> track_str;
    biggles::time_stamp final = part.last_time_stamp();
    biggles::track_collection tc1 = part.tracks();
    const std::string lim("|");
    for (it1 = tc1.begin(); it1 not_eq tc1.end(); ++it1)
        track_str.push_back(lim + track_to_str(**it1, final) + lim);
    std::sort(track_str.begin(), track_str.end());
    return track_str;
}



inline std::string parameters_to_string(const biggles::model::parameters& p) {
    return boost::str(
        boost::format("b = %.3f, c = %.3f, o = %.3f, s = %.3f")
            % biggles::model::mean_new_tracks_per_frame(p)
            % biggles::model::mean_false_observations_per_frame(p)
            % biggles::model::generate_observation_probability(p)
            % biggles::model::frame_to_frame_survival_probability(p)
    );
}

/*
inline std::string crc16(const std::string &data, bool uppercase) {
    std::stringstream str;
    boost::crc_16_type  result;
    result.process_bytes( data.c_str(), data.size() );
    if   (uppercase) str << boost::format("%04X") % result.checksum();
    else str << boost::format("%04x") % result.checksum();
    return str.str();
}
*/


bool are_equal(const biggles::clutter_t& oc1, const biggles::clutter_t& oc2);

inline bool operator==(const biggles::partition& part1, const biggles::partition& part2) {
    if (not are_equal(part1.clutter(), part2.clutter())) {
        return false;
    }
    biggles::track_collection tc1 = part1.tracks();
    biggles::track_collection tc2 = part2.tracks();
    if (tc1.size() not_eq tc2.size()) {
        return false;
    }
    biggles::track_collection::const_iterator it1, it2;
    for (it1 = tc1.begin(); it1 not_eq tc1.end(); ++it1) {
        bool match = false;
        const biggles::track& track_to_find = **it1;
        for (it2 = tc2.begin(); it2 not_eq tc2.end(); ++it2) {
            if (track_to_find == **it2) {
                match = true;
                break;
            }
        }
        if (not match) {
            return false;
        }
    }
    return true;
}

std::string clutter_to_str(const biggles::clutter_t& clutter, biggles::time_stamp end_ts = 0);

/// \brief proposal move proposal-density-ratio test
std::deque<std::string> move_pdr_test(const biggles::partition_ptr_t& orignal_part_ptr, const std::string& forward_move,
    const int total, test_t& test_results);
/// \brief proposal move proposal-density-ratio test
std::deque<std::string> move_pdr_test_adv(const biggles::partition_ptr_t& orignal_part_ptr,
    const std::string& forward_move, const int total, test_t& test_results);

struct random_gen_paras {
    int min_tracks;
    int max_tracks;
    float p_yes; ///< \brief for clutter
    float p_no; ///< \brief for clutter
    float p_tr; ///< \brief for track
    float lambda;
    std::string input;
    std::string output;
};

struct general_paras {
    int total;
    unsigned int seed;
    random_gen_paras rgp;
    std::string opts;
};

bool parse_args(int argc, char** argv, general_paras & gen_pars);

biggles::partition_ptr_t demote(biggles::partition_ptr_t& part_sample);

/// \brief creates a random partition for test and debugging purposes
biggles::partition_ptr_t random_debug_partition(const random_gen_paras& paras);

/// \brief creates a random partition for test and debugging purposes. This version creates a larger variety
biggles::partition_ptr_t random_debug_partition_adv(const random_gen_paras& paras);

class biggles_move {
    struct move {
        virtual bool exec(const biggles::partition_ptr_t& start_partition, biggles::partition_ptr_t& end_partition, float& pdr) = 0;
    };

    struct birth : public move {
        bool exec(const biggles::partition_ptr_t& start_partition, biggles::partition_ptr_t& end_partition, float& pdr) {
            return biggles::mh_moves::birth(start_partition, end_partition, pdr);
        }
    };
    struct death : public move {
        bool exec(const biggles::partition_ptr_t& start_partition, biggles::partition_ptr_t& end_partition, float& pdr) {
            return biggles::mh_moves::death(start_partition, end_partition, pdr);
        }
    };
    struct merge : public move {
        bool exec(const biggles::partition_ptr_t& start_partition, biggles::partition_ptr_t& end_partition, float& pdr) {
            return biggles::mh_moves::merge(start_partition, end_partition, pdr);
        }
    };
    struct split : public move {
        bool exec(const biggles::partition_ptr_t& start_partition, biggles::partition_ptr_t& end_partition, float& pdr) {
            return biggles::mh_moves::split(start_partition, end_partition, pdr);
        }
    };
    struct extend : public move {
        bool exec(const biggles::partition_ptr_t& start_partition, biggles::partition_ptr_t& end_partition, float& pdr) {
            return biggles::mh_moves::extend(start_partition, end_partition, pdr);
        }
    };
    struct reduce : public move {
        bool exec(const biggles::partition_ptr_t& start_partition, biggles::partition_ptr_t& end_partition, float& pdr) {
            return biggles::mh_moves::reduce(start_partition, end_partition, pdr);
        }
    };
    struct update : public move {
        bool exec(const biggles::partition_ptr_t& start_partition, biggles::partition_ptr_t& end_partition, float& pdr) {
            return biggles::mh_moves::update(start_partition, end_partition, pdr);
        }
    };

    struct transfer : public move {
        bool exec(const biggles::partition_ptr_t& start_partition, biggles::partition_ptr_t& end_partition, float& pdr) {
            return biggles::mh_moves::transfer(start_partition, end_partition, pdr);
        }
    };

    struct cross_over : public move {
        bool exec(const biggles::partition_ptr_t& start_partition, biggles::partition_ptr_t& end_partition, float& pdr) {
            return biggles::mh_moves::cross_over(start_partition, end_partition, pdr);
        }
    };

    boost::shared_ptr<move> move_ptr;
public:
    biggles_move(const biggles::mh_moves::move_type &move_type) {
        if (move_type == biggles::mh_moves::BIRTH)
            move_ptr = boost::shared_ptr<move>(new birth);
        else if (move_type == biggles::mh_moves::DEATH)
            move_ptr = boost::shared_ptr<move>(new death);
        else if (move_type == biggles::mh_moves::MERGE)
            move_ptr = boost::shared_ptr<move>(new merge);
        else if (move_type == biggles::mh_moves::SPLIT)
            move_ptr = boost::shared_ptr<move>(new split);
        else if (move_type == biggles::mh_moves::EXTEND)
            move_ptr = boost::shared_ptr<move>(new extend);
        else if (move_type == biggles::mh_moves::REDUCE)
            move_ptr = boost::shared_ptr<move>(new reduce);
        else if (move_type == biggles::mh_moves::UPDATE)
            move_ptr = boost::shared_ptr<move>(new update);
        else if (move_type == biggles::mh_moves::TRANSFER)
            move_ptr = boost::shared_ptr<move>(new transfer);
        else if (move_type == biggles::mh_moves::CROSS_OVER)
            move_ptr = boost::shared_ptr<move>(new cross_over);
        else
            throw std::logic_error("move type not recognised: " + move_type);
    }
    biggles_move(const std::string& move_type) {
        if (move_type == "Birth")
            move_ptr = boost::shared_ptr<move>(new birth);
        else if (move_type == "Death")
            move_ptr = boost::shared_ptr<move>(new death);
        else if (move_type == "Merge")
            move_ptr = boost::shared_ptr<move>(new merge);
        else if (move_type == "Split")
            move_ptr = boost::shared_ptr<move>(new split);
        else if (move_type == "Extend")
            move_ptr = boost::shared_ptr<move>(new extend);
        else if (move_type == "Reduce")
            move_ptr = boost::shared_ptr<move>(new reduce);
        else if (move_type == "Update")
            move_ptr = boost::shared_ptr<move>(new update);
        else if (move_type == "Transfer")
            move_ptr = boost::shared_ptr<move>(new transfer);
        else if (move_type == "Cross-over")
            move_ptr = boost::shared_ptr<move>(new cross_over);
        else
            throw std::logic_error("move type not recognised: " + move_type);
    }
    bool operator()(const biggles::partition_ptr_t& start_partition, biggles::partition_ptr_t& end_partition, float& pdr) {
        return move_ptr->exec(start_partition, end_partition, pdr);
    }
};

template<class NUM>
class TickTimer {
    NUM lag_;
    NUM last_;
public:
    TickTimer(NUM lag) : lag_(lag), last_(0) {}
    bool passed(const NUM time) {
        if (time < last_ + lag_)
            return false;
        last_ = time;
        return true;
    }
};

/// \brief the standard deviation on the samples proposal density ratio
///
/// \param n_forth number of occurance of the forward partition
/// \param total total number of forward samples
/// \param n_back number of occurance of the backward partition (given the forward partition)
float pdr_sample_stdev(const float& n_forth, const float& total, const float& n_back);

struct dump_not_push {
    template<class DATA>
    void push_back(const DATA& data) { std::cout << data << std::endl; }
};

namespace biggles {
    namespace sampling {
        float sample_beta(float alpha, float beta);
    }
}



#endif // BIGGLES_TEST_UTILITY_HPP__
