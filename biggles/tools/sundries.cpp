#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/timer.hpp>
#include "biggles/io.hpp"
#include "tools/sundries.hpp"
#include "detail/random.hpp"
#include "detail/physics.hpp"
#include "tools/debug.hpp"

using namespace biggles;

const float MAX_X = 8.f;
const float MAX_Y = 8.f;
const time_stamp END_TIME = 20;

int encode3_obs(const observation& o) {
    return (int(x(o)) % 3 + 3 * (int(y(o)) % 3));
}

float auto_correlation(data_iter b, data_iter e, size_t lag) {
    size_t len = std::distance(b, e);
    BOOST_ASSERT(lag>0);
    BOOST_ASSERT(len > lag+1);
    data_iter b2 = b;
    data_iter e1 = b;
    std::advance(b2, lag);
    std::advance(e1, len-lag);
    std::pair<float, float> ms1 = mean_stdev(b, e1);
    std::pair<float, float> ms2 = mean_stdev(b2, e);
    std::deque<float> products(len-1);
    std::transform(b, e1, b2, products.begin(), std::multiplies<float>());
    float sum = std::accumulate(products.begin(), products.end(), 0.f);
    return (sum - (len-lag)*ms1.first*ms2.first)/((len-lag-1)*ms1.second*ms2.second);
}

float correlation(data_iter b1, data_iter e1, data_iter b2) {
    size_t len = std::distance(b1, e1);
    data_iter e2 = b2;
    std::advance(e2, len);
    std::deque<float> products(len);
    std::pair<float, float> ms1 = mean_stdev(b1, e1);
    std::pair<float, float> ms2 = mean_stdev(b2, e2);
    std::transform(b1, e1, b2, products.begin(), std::multiplies<float>());
    float sum = std::accumulate(products.begin(), products.end(), 0.f);
    return (sum - len*ms1.first*ms2.first)/((len-1)*ms1.second*ms2.second);
}

std::string track_to_str(const track & track0, time_stamp end_ts) {

    time_stamp ts=track0.first_time_stamp();
    std::stringstream track_ss;
    track_ss << std::string(ts, ' ');
    for (; ts < track0.last_time_stamp(); ++ts) {
        observation_collection::const_range obs_at_t(track0.observations().at_time_stamp(ts));
        if (obs_at_t.first != obs_at_t.second)
            track_ss << encode3_obs(*obs_at_t.first);
        else
           track_ss << "_";
    }
    while (ts++ < end_ts) {
        track_ss << " ";
    }
    return track_ss.str();
}

float pdr_sample_stdev(const float& n_forth, const float& total, const float& n_back) {
    return sqrtf( n_forth*n_forth/(total*total*n_back*n_back) * (
        n_forth*(1.f - n_forth/total) + n_forth/n_back*n_forth * (1.f - n_back/n_forth)
        ) );
}

partition_ptr_t demote(partition_ptr_t& part_sample) {
    clutter_ptr clutter(new clutter_t(part_sample->clutter()));
    boost::shared_ptr<track_collection> tracks(new track_collection());
    BOOST_FOREACH(const boost::shared_ptr<const track>& t, part_sample->tracks()) {
        BOOST_FOREACH(const observation& o, t->observations()) {
            clutter->insert(o);
        }
    }
    return partition_ptr_t(new partition(tracks, clutter, part_sample->expansion()));
}

bool parse_args(int argc, char** argv, general_paras& gen_pars) {
    namespace po = boost::program_options;
    po::options_description desc("Options");
    random_gen_paras& paras(gen_pars.rgp);
    if (gen_pars.seed == 0)
        gen_pars.seed = std::time(0);

    desc.add_options()
        ("samples,n", po::value< int >(&gen_pars.total)->default_value(int(gen_pars.total)),
            "number of samples")
        ("min-tracks,m", po::value< int >(&paras.min_tracks)->default_value(int(paras.min_tracks)),
            "min. number of tracks")
        ("max-tracks,x", po::value< int >(&paras.max_tracks)->default_value(int(paras.max_tracks)),
            "max. number of tracks")
        ("p-yes", po::value< float >(&paras.p_yes)->default_value(float(paras.p_yes)), "p-yes")
        ("p-no", po::value< float >(&paras.p_no)->default_value(float(paras.p_no)), "p-no")
        ("p-tr", po::value< float >(&paras.p_tr)->default_value(float(paras.p_tr)), "p-tr")
        ("lambda,l", po::value< float >(&paras.lambda)->default_value(float(paras.lambda)), "lambda for track length")
        ("seed,s", po::value< std::string >()->default_value(boost::str(boost::format("%x") % gen_pars.seed)),
            "random seed in hex format")
        ("input,i", po::value< std::string >(&paras.input)->default_value(""), "partition input file")
        ("output,o", po::value< std::string >(&paras.output)->default_value(""), "partition output file")
        ("opts", po::value< std::string >(&gen_pars.opts), "additional options string")
        ("help,h", "print help")
    ;

    po::options_description argsdesc("Command line arguments");
    argsdesc.add(desc);

    po::variables_map args;
    po::store(po::command_line_parser(argc, argv).options(argsdesc).run(), args);
    po::notify(args);
    std::stringstream interpreter;
    interpreter << std::hex << args["seed"].as<std::string>();
    interpreter >> gen_pars.seed;
    if (args.count("help") > 0) {
        std::cout << argsdesc << std::endl;
        return false;
    }
    return true;
}


bool are_equal(const biggles::clutter_t& oc1, const biggles::clutter_t& oc2) {
    if ((oc1.size() == 0) xor (oc2.size() == 0)) {
        return false;
    }
    biggles::time_stamp first_ts = oc1.first_time_stamp();
    biggles::time_stamp last_ts = oc1.last_time_stamp();
    if (first_ts != oc2.first_time_stamp()) return false;
    if (last_ts != oc2.last_time_stamp()) return false;
    for (biggles::time_stamp ts = first_ts; ts < last_ts; ++ts) {
        if (oc1.size(ts) != oc2.size(ts)) return false;
        BOOST_FOREACH(const biggles::observation& o, oc1.observations(ts)) {
            if (not oc2.contains(o))
                return false;
        }
    }
    return true;
}



std::string clutter_to_str(const clutter_t& clutter, time_stamp end_ts) {
    time_stamp ts=clutter.first_time_stamp();
    std::stringstream track_ss;
    for (int i = 0; i < int(ts); ++i)
        track_ss << " :";
    for (; ts < clutter.last_time_stamp(); ++ts) {
        std::deque<observation> obs_at_t(clutter.observations(ts));
        if (not obs_at_t.empty()) {
            for (std::deque<observation>::iterator it = obs_at_t.begin(); it != obs_at_t.end(); ++it) {
                track_ss << encode3_obs(*it);
            }
            track_ss << ":";
        }
        else {
           track_ss << " :";
        }
    }
    while (ts++ < end_ts) {
        track_ss << " :";
    }
    return track_ss.str();
}

partition_ptr_t random_debug_partition(const random_gen_paras& paras) {
    if (paras.input.size() > 0) {
        std::ifstream fh(paras.input.c_str());
        partition_ptr_t part(io::read_partition_from_json_stream(fh));
        return part;
    }
    BOOST_ASSERT(paras.min_tracks <= paras.max_tracks);
    BOOST_ASSERT(7 > paras.max_tracks);
    BOOST_ASSERT(paras.min_tracks >= 0);
    boost::shared_ptr<track_collection> tracks(new track_collection());
    std::deque<observation> clutter_obs;
    const float p_yes = paras.p_yes;
    const float p_no = paras.p_no;
    float p = p_yes;
    for (int y = 0; y < 3; ++y) {
        for (int t = 0 ; t < 10; ++t) {
            if (sampling::uniform_real(0, 1) < p) {
                p = p_yes;
                clutter_obs.push_back(new_obs(1, y, t));
            } else {
                p = p_no;
            }
        }
    }
    int num_tracks = sampling::uniform_int(paras.min_tracks, paras.max_tracks + 1);
    int i_track = 0;
    float p_tr = paras.p_tr;
    for (int x = 0; x<3; x+=2) {
        for (int y = 0; y < 3; ++y) {
            if (i_track >= num_tracks) break;
            observation_collection tr_obs;
            time_stamp tlo, thi;
            bool done = false;
            while (not done) {
                tr_obs.clear();
                tlo = sampling::uniform_int(0, 8);
                thi = std::min(sampling::poisson(paras.lambda) + 1 + tlo, time_stamp(10));
                for (time_stamp t = tlo; t < thi; ++t)
                    if (sampling::uniform_real(0, 1) < p_tr) tr_obs.insert(
                        new_obs(x, y, t));
                done = tr_obs.size() > 1;
            }
            tracks->insert(track(tlo, thi, tr_obs.begin(), tr_obs.end(), 1.0f));
            ++i_track;
        }
    }
    clutter_ptr clutter(new clutter_t(clutter_obs.begin(), clutter_obs.end()));
    partition_ptr_t part(new partition(tracks, clutter, partition::expansion_2d(0.f, 2.f, 0.f, 2.f)));
    if (paras.output.size() > 0) {
        std::ofstream fh(paras.output.c_str());
        io::write_partition_to_json_stream(fh, *part);
    }
    return part;
}

observation sample_obs(const observation& last_obs, const time_stamp time, const float speed_limit) {
    float angle = sampling::uniform_real() * 2.f * 3.141592653589793f;
    float radius = sampling::uniform_real() * std::abs(t(last_obs)-time) * speed_limit;
    return new_obs(radius*std::cos(angle) + x(last_obs), radius*std::sin(angle) + y(last_obs), time);
}

observation sample_obs(const float x, const float y, const time_stamp time) {
    return new_obs(sampling::uniform_real() * x, sampling::uniform_real() * y, time);
}

partition_ptr_t random_debug_partition_adv(const random_gen_paras& paras) {
    if (paras.input.size() > 0) {
        std::ifstream fh(paras.input.c_str());
        partition_ptr_t part(io::read_partition_from_json_stream(fh));
        return part;
    }
    BOOST_ASSERT(paras.min_tracks <= paras.max_tracks);
    BOOST_ASSERT(7 > paras.max_tracks);
    BOOST_ASSERT(paras.min_tracks >= 0);
    boost::shared_ptr<track_collection> tracks(new track_collection());
    std::deque<observation> clutter_obs;
    const float p_yes = paras.p_yes;
    const float p_no = paras.p_no;
    float p = p_yes;
    for (int y = 0; y < 3; ++y) {
        for (int t = 0 ; t < 10; ++t) {
            if (sampling::uniform_real(0, 1) < p) {
                p = p_yes;
                clutter_obs.push_back(new_obs(1, y, t));
            } else {
                p = p_no;
            }
        }
    }
    int num_tracks = sampling::uniform_int(paras.min_tracks, paras.max_tracks + 1);
    int i_track = 0;
    float p_tr = paras.p_tr;
    while (i_track < num_tracks) {
        observation_collection tr_obs;
        time_stamp tlo, thi;
        bool done = false;
        while (not done) {
            bool started = false;
            tr_obs.clear();
            tlo = sampling::uniform_int(0, 8);
            thi = std::min(sampling::poisson(paras.lambda) + 1 + tlo, time_stamp(END_TIME));
            for (time_stamp t = tlo; t < thi; ++t) {
                if (sampling::uniform_real() < p_tr) {
                    if (not started) {
                        tr_obs.insert(sample_obs(MAX_X, MAX_Y, t));
                        started = true;
                    } else {
                        tr_obs.insert(sample_obs(tr_obs.back(), t, detail::speed_of_light_));
                    }
                }
            }
            done = tr_obs.size() > 1;
        }
        tracks->insert(track(tlo, thi, tr_obs.begin(), tr_obs.end(), 1.0f));
        ++i_track;
    }
    clutter_ptr clutter(new clutter_t(clutter_obs.begin(), clutter_obs.end()));
    partition_ptr_t part(new partition(tracks, clutter, partition::expansion_2d(0.f, 2.f, 0.f, 2.f)));
    if (paras.output.size() > 0) {
        std::ofstream fh(paras.output.c_str());
        io::write_partition_to_json_stream(fh, *part);
    }
    return part;
}

std::deque<std::string> move_pdr_test(const partition_ptr_t& orignal_part_ptr, const std::string& forward_move,
        const int total, test_t& test_results)
{
    std::string backward_move;
    if (forward_move == "Birth") backward_move = "Death";
    if (forward_move == "Death") backward_move = "Birth";
    if (forward_move == "Extend") backward_move = "Reduce";
    if (forward_move == "Reduce") backward_move = "Extend";
    if (forward_move == "Merge") backward_move = "Split";
    if (forward_move == "Split") backward_move = "Merge";
    if (forward_move == "Update") backward_move = "Update";
    if (forward_move == "Transfer") backward_move = "Transfer";
    if (forward_move == "Cross-over") backward_move = "Cross-over";

    time_stamp global_fst = orignal_part_ptr->first_time_stamp();
    time_stamp global_lst = orignal_part_ptr->last_time_stamp();


    dump_not_push messages;
    //std::deque<std::string> messages;
    messages.push_back("Note: the test has two elements:");
    messages.push_back("1) the sampled 'F/B' ratio needs to be the similar to the calculated.");
    messages.push_back("   if increasing the number of samples does not help, then there is something WRONG!");
    messages.push_back("2) F-PDR and B-PDR need to be indentical. The F-PDR + B-PDR is reported, if there is a problem.");

    messages.push_back(boost::str(
        boost::format("original '%s'") % partition_to_string(*orignal_part_ptr)
    ));

    messages.push_back(boost::str(
        boost::format("clutter '%s'") % clutter_to_str(orignal_part_ptr->clutter(), orignal_part_ptr->last_time_stamp())
    ));


    messages.push_back("Pos 0123456789");
    BOOST_FOREACH (const std::string& s, partition_to_string_list(*orignal_part_ptr)) {
        messages.push_back(boost::str( boost::format("   %s") % s ));
    }

    messages.push_back(boost::str(
        boost::format("Forward move (F): '%s', Backward move (B): '%s'") % forward_move % backward_move
    ));

    biggles_move forward(forward_move);
    biggles_move backward(backward_move);


    std::map<std::string, float> forward_freq;
    std::map<std::string, float> move_pdr;
    std::map<std::string, std::set<float> > backward_pdr;
    std::map<std::string, float> undo;
    int fail_count = 0;
    const int back_counts = total;
    bool done = false;
    float count_pdr_add_to_0 = 0.f;
    float reversed = 0;
    //OK(orignal_part_ptr->get_capability_recorder()->transfer_pairs().size());
    for (int i = 0; i < total ; i++) {
        float forth_prob=0.f;
        partition_ptr_t target_part_ptr = orignal_part_ptr;
        done = forward(orignal_part_ptr, target_part_ptr, forth_prob);
        if (not done) ++fail_count;
        else {
            std::string part_str(partition_to_string(*target_part_ptr));
            if (forward_freq.find(part_str) == forward_freq.end()) {
                target_part_ptr->set_first_time_stamp(global_fst);
                target_part_ptr->set_last_time_stamp(global_lst);
                float back_prob=0.f;
                capability_recorder_ptr cap_rec_ptr(new_cap_recorder(*target_part_ptr));
                target_part_ptr->set_capability_recorder(cap_rec_ptr);
                //OK(part_str);
                //OK(cap_rec_ptr->transfer_pairs().size());
                for (int i2 = 0; i2 < back_counts; ++i2) {
                    partition_ptr_t backward_part_ptr = target_part_ptr;
                    done = backward(target_part_ptr, backward_part_ptr, back_prob);
                    if (done) {
                        backward_part_ptr->set_first_time_stamp(global_fst);
                        backward_part_ptr->set_last_time_stamp(global_lst);
                    }
                    if (*orignal_part_ptr == *backward_part_ptr) {
                        undo[part_str]++;
                        count_pdr_add_to_0 += (forth_prob + back_prob == 0.f) ? 1.f : 0.f;
                        backward_pdr[part_str].insert(back_prob);
                        reversed ++;
                    }
                }
            }
            forward_freq[part_str]++;
            move_pdr[part_str] = forth_prob;
        }
    }
    float event_count = 0.f;
    //size_t in_limits = 0;
    std::map<int, int> error_limits;
    bool pdr_test_ok = true;
    for (std::map<std::string, float >::const_iterator it = forward_freq.begin(); it != forward_freq.end(); ++it) {
        std::string forward_partition(it->first);
        float forward_part_freq = it->second;
        const float pdr = move_pdr[forward_partition];
        std::set<float> back_pdrs = backward_pdr[forward_partition];
        const float ratio = expf(-pdr);
        float forward_ratio = forward_part_freq/float(total);
        float backward_part_freq = undo[forward_partition];
        float backward_ratio = backward_part_freq/back_counts;
        float first_back_pdr = 0.f;
        //float error_est = pdr_sample_stdev(forward_part_freq, total, backward_part_freq);
        float forw_prob = forward_part_freq/total;
        float forw_err = binom_prob_error_est(forward_part_freq, total);
        float back_prob = backward_part_freq/back_counts;
        float back_err = binom_prob_error_est(backward_part_freq, back_counts);
        float error_est = error_of_fraction(forw_prob, forw_err, back_prob, back_err);
        if (not back_pdrs.empty())
            first_back_pdr = *back_pdrs.begin();
        float est_diff = fabsf(ratio - (forward_ratio/backward_ratio));
        if (error_est > 0)
            error_limits[int(est_diff/ error_est)]++;
        else if (est_diff == 0)
            error_limits[0]++;
        else
            error_limits[9]++;
        if (fabsf(pdr + first_back_pdr) < 5.e-6) {
            messages.push_back(boost::str(
                boost::format("'%s', F = %6.2f%%, B|F = %6.2f%%, F/B = %.2f +/- %.2f (calc: %.2f), F-PDR = %.2f")
                % forward_partition % (100.f * forward_ratio) % (100.f * backward_ratio) % (forward_ratio/backward_ratio)
                % error_est % ratio % pdr
            ));
        } else {
            messages.push_back(boost::str(
                boost::format("'%s', F = %6.2f%%, B|F = %6.2f%%, F/B = %.2f +/- %.2f (calc: %.2f), F+B = %.2f + %.2f = %g")
                % forward_partition % (100.f * forward_ratio) % (100.f * backward_ratio) % (forward_ratio/backward_ratio)
                % error_est % ratio % pdr % first_back_pdr % (pdr+first_back_pdr)
            ));
            pdr_test_ok = false;
        }
        event_count += forward_part_freq;
        if (back_pdrs.size() > 1) {
            messages.push_back(boost::str(boost::format("  num B-PDR = %d") % back_pdrs.size()));
            BOOST_FOREACH( const float& pdr, back_pdrs) {
                messages.push_back(boost::str(boost::format("  B-PDR = %.2f") % pdr));
            }
        }
    }
    messages.push_back(boost::str(
        boost::format("failed, Freq = %.2f%%%%") % (float(fail_count)/float(total)*100.f)
    ));
    messages.push_back(boost::str(
        boost::format("reversed, Freq = %.2f%%%%") % (reversed/float(total)*100.f)
    ));
    messages.push_back(boost::str(
        boost::format("fail-count %d, event-count %d, sum = %d") % fail_count % int(event_count) % (fail_count + int(event_count))
    ));
    messages.push_back(boost::str(
        boost::format("count_pdr_add_to_0/reversed =  %.2f") % (count_pdr_add_to_0/reversed)
    ));
    for (std::map<int, int>::iterator it = error_limits.begin(); it != error_limits.end(); ++it) {
        messages.push_back(boost::str(
            boost::format("%d/%d (%6.2f%%) samples are in %d sigma range of calculated ratio")
                % (it->second) % forward_freq.size()
                % ( 100.f * float(it->second)/ float(forward_freq.size()))
                % (it->first + 1)
        ));
    }
    test_results["fail count + event count == total"] = float(fail_count) + event_count == float(total);
    test_results["All PDR tests passed"] = pdr_test_ok;
    //test_results["count of PDRs that add to 0 == undone moves"] = count_pdr_add_to_0 == reversed;
    return std::deque<std::string>();
    //return messages;
}

std::deque<std::string> move_pdr_test_adv(const partition_ptr_t& orignal_part_ptr, const std::string& forward_move,
        const int total, test_t& test_results)
{
    std::string backward_move;
    if (forward_move == "Birth") backward_move = "Death";
    if (forward_move == "Death") backward_move = "Birth";
    if (forward_move == "Extend") backward_move = "Reduce";
    if (forward_move == "Reduce") backward_move = "Extend";
    if (forward_move == "Merge") backward_move = "Split";
    if (forward_move == "Split") backward_move = "Merge";
    if (forward_move == "Update") backward_move = "Update";
    if (forward_move == "Transfer") backward_move = "Transfer";
    if (forward_move == "Cross-over") backward_move = "Cross-over";

    time_stamp global_fst = orignal_part_ptr->first_time_stamp();
    time_stamp global_lst = orignal_part_ptr->last_time_stamp();

    dump_not_push messages;
    //std::deque<std::string> messages;
    messages.push_back("Note: the test has two elements:");
    messages.push_back("1) the sampled 'F/B' ratio needs to be the similar to the calculated.");
    messages.push_back("   if increasing the number of samples does not help, then there is something WRONG!");
    messages.push_back("2) F-PDR and B-PDR need to be indentical. The F-PDR + B-PDR is reported, if there is a problem.");

    messages.push_back(boost::str(
        boost::format("original '%s'") % (orignal_part_ptr->as_string())
    ));

    messages.push_back(boost::str(
        boost::format("clutter '%s'") % clutter_to_str(orignal_part_ptr->clutter(), orignal_part_ptr->last_time_stamp())
    ));


    /*
    messages.push_back("Pos 0123456789");
    BOOST_FOREACH (const std::string& s, partition_to_string_list(*orignal_part_ptr)) {
        messages.push_back(boost::str( boost::format("   %s") % s ));
    }
    */

    messages.push_back(boost::str(
        boost::format("Forward move (F): '%s', Backward move (B): '%s'") % forward_move % backward_move
    ));

    biggles_move forward(forward_move);
    biggles_move backward(backward_move);


    std::map<std::string, float> forward_freq;
    std::map<std::string, float> move_pdr;
    std::map<std::string, std::set<float> > backward_pdr;
    std::map<std::string, float> undo;
    int fail_count = 0;
    const int back_counts = total;
    bool done = false;
    float count_pdr_add_to_0 = 0.f;
    float reversed = 0;
    //OK(orignal_part_ptr->get_capability_recorder()->transfer_pairs().size());
    for (int i = 0; i < total ; i++) {
        float forth_prob=0.f;
        partition_ptr_t target_part_ptr = orignal_part_ptr;
        done = forward(orignal_part_ptr, target_part_ptr, forth_prob);
        if (not done) ++fail_count;
        else {
            std::string part_str(target_part_ptr->as_string());
            if (forward_freq.find(part_str) == forward_freq.end()) {
                target_part_ptr->set_first_time_stamp(global_fst);
                target_part_ptr->set_last_time_stamp(global_lst);
                float back_prob=0.f;
                capability_recorder_ptr cap_rec_ptr(new_cap_recorder(*target_part_ptr));
                target_part_ptr->set_capability_recorder(cap_rec_ptr);
                //OK(part_str);
                //OK(cap_rec_ptr->transfer_pairs().size());
                boost::timer timer;
                for (int i2 = 0; i2 < back_counts; ++i2) {
                    partition_ptr_t backward_part_ptr = target_part_ptr;
                    timer.restart();
                    done = backward(target_part_ptr, backward_part_ptr, back_prob);
                    if (done) {
                        backward_part_ptr->set_first_time_stamp(global_fst);
                        backward_part_ptr->set_last_time_stamp(global_lst);
                    }
                    if (*orignal_part_ptr == *backward_part_ptr) {
                        undo[part_str]++;
                        count_pdr_add_to_0 += (forth_prob + back_prob == 0.f) ? 1.f : 0.f;
                        backward_pdr[part_str].insert(back_prob);
                        reversed ++;
                    }
                }
            }
            forward_freq[part_str]++;
            move_pdr[part_str] = forth_prob;
        }
    }
    float event_count = 0.f;
    //size_t in_limits = 0;
    std::map<int, int> error_limits;
    bool pdr_test_ok = true;
    for (std::map<std::string, float >::const_iterator it = forward_freq.begin(); it != forward_freq.end(); ++it) {
        std::string forward_partition(it->first);
        float forward_part_freq = it->second;
        const float pdr = move_pdr[forward_partition];
        std::set<float> back_pdrs = backward_pdr[forward_partition];
        const float ratio = expf(-pdr);
        float forward_ratio = forward_part_freq/float(total);
        float backward_part_freq = undo[forward_partition];
        float backward_ratio = backward_part_freq/back_counts;
        float first_back_pdr = 0.f;
        //float error_est = pdr_sample_stdev(forward_part_freq, total, backward_part_freq);
        float forw_prob = forward_part_freq/total;
        float forw_err = binom_prob_error_est(forward_part_freq, total);
        float back_prob = backward_part_freq/back_counts;
        float back_err = binom_prob_error_est(backward_part_freq, back_counts);
        float error_est = error_of_fraction(forw_prob, forw_err, back_prob, back_err);
        if (not back_pdrs.empty())
            first_back_pdr = *back_pdrs.begin();
        float est_diff = fabsf(ratio - (forward_ratio/backward_ratio));
        if (error_est > 0)
            error_limits[int(est_diff/ error_est)]++;
        else if (est_diff == 0)
            error_limits[0]++;
        else
            error_limits[9]++;
        if (fabsf(pdr + first_back_pdr) < 5.e-6) {
            messages.push_back(boost::str(
                boost::format("'%s', F = %6.2f%%, B|F = %6.2f%%, F/B = %.2f +/- %.2f (calc: %.2f), F-PDR = %.2f")
                % forward_partition % (100.f * forward_ratio) % (100.f * backward_ratio) % (forward_ratio/backward_ratio)
                % error_est % ratio % pdr
            ));
        } else {
            messages.push_back(boost::str(
                boost::format("'%s', F = %6.2f%%, B|F = %6.2f%%, F/B = %.2f +/- %.2f (calc: %.2f), F+B = %.2f + %.2f = %g")
                % forward_partition % (100.f * forward_ratio) % (100.f * backward_ratio) % (forward_ratio/backward_ratio)
                % error_est % ratio % pdr % first_back_pdr % (pdr+first_back_pdr)
            ));
            pdr_test_ok = false;
        }
        event_count += forward_part_freq;
        if (back_pdrs.size() > 1) {
            messages.push_back(boost::str(boost::format("  num B-PDR = %d") % back_pdrs.size()));
            BOOST_FOREACH( const float& pdr, back_pdrs) {
                messages.push_back(boost::str(boost::format("  B-PDR = %.2f") % pdr));
            }
        }
    }
    messages.push_back(boost::str(
        boost::format("failed, Freq = %.2f%%%%") % (float(fail_count)/float(total)*100.f)
    ));
    messages.push_back(boost::str(
        boost::format("reversed, Freq = %.2f%%%%") % (reversed/float(total)*100.f)
    ));
    messages.push_back(boost::str(
        boost::format("fail-count %d, event-count %d, sum = %d") % fail_count % int(event_count) % (fail_count + int(event_count))
    ));
    messages.push_back(boost::str(
        boost::format("count_pdr_add_to_0/reversed =  %.2f") % (count_pdr_add_to_0/reversed)
    ));
    int max_std = 0;
    for (std::map<int, int>::iterator it = error_limits.begin(); it != error_limits.end(); ++it) {
        max_std = std::max(max_std, it->first);
        messages.push_back(boost::str(
            boost::format("%d/%d (%6.2f%%) samples are in %d sigma range of calculated ratio")
                % (it->second) % forward_freq.size()
                % ( 100.f * float(it->second)/ float(forward_freq.size()))
                % (it->first + 1)
        ));
    }
    test_results["fail count + event count == total"] = float(fail_count) + event_count == float(total);
    test_results["All PDR tests passed"] = pdr_test_ok;
    //test_results["Sampled F/B is similar to calculated F/B"] = max_std < 8;
    //test_results["count of PDRs that add to 0 == undone moves"] = count_pdr_add_to_0 == reversed;
    return std::deque<std::string>();
    //return messages;
}

namespace biggles { namespace sampling {
    float sample_beta(float alpha, float beta) {
        typedef boost::gamma_distribution<float> dist_t;
        dist_t dist_x(alpha);
        dist_t dist_y(beta);
        boost::variate_generator<boost::mt19937&, dist_t > sample_x(detail::random_generator, dist_x);
        boost::variate_generator<boost::mt19937&, dist_t > sample_y(detail::random_generator, dist_y);
        float x(sample_x());
        float y(sample_y());
        return x/(x+y);
    }
}}
