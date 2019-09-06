#include "tools/partition_graph_edit_distance.hpp"

#include <set>
#include <deque>
#include <algorithm>
#include <boost/foreach.hpp>
#include <utility> // for std::pair
#include "tools/debug.hpp"

namespace biggles {

std::ostream& operator<<(std::ostream& os, const link_t& l) {
    os
        << "(" << x(l.first) << ", " << y(l.first) << ", " << t(l.first) << ")"
        << " -> "
        << "(" << x(l.second) << ", " << y(l.second) << ", " << t(l.second) << ")";
    return os;
}

bool operator<(const link_t& l1, const link_t& l2) {
    return l1.first < l2.first or (l1.first == l2.first and l1.second < l2.second);
}
bool operator==(const link_t& l1, const link_t& l2) {
    return l1.first == l2.first and l1.second == l2.second;
}

link_set_t get_links(const partition& part) {
    link_set_t result;
    BOOST_FOREACH (const shared_const_track_ptr& tp, part.tracks()) {
        const observation_collection& obs_coll = tp->observations();
        BOOST_ASSERT(obs_coll.size() > 1);
        observation_collection::const_iterator fir = obs_coll.begin();
        observation_collection::const_iterator sec = fir;
        ++sec;
        for (; sec != obs_coll.end(); ++fir, ++sec) {
            insert_result_t insert_result = result.insert(link_t(*fir, *sec));
            BOOST_ASSERT(insert_result.second);
        }
    }
    return result;
}

float partition_ged(const partition& part1, const partition& part2) {
    link_set_t links1 = get_links(part1);
    link_set_t links2 = get_links(part2);
    std::deque<link_t> sym_diff(links1.size() + links2.size());
    std::deque<link_t>::iterator end_iter = std::set_symmetric_difference(
        links1.begin(), links1.end(), links2.begin(), links2.end(), sym_diff.begin());
    return std::distance(sym_diff.begin(), end_iter);
}

float map_dot_product(numeric_partition_t::iterator b1, numeric_partition_t::iterator e1,
    numeric_partition_t::iterator b2)
{
    float sum = 0.f;
    while (b1 != e1) {
        sum += (b1++->second)*(b2++->second);
    }
    return sum;
}

void PartitionMeanVariance::add_partition(const partition& part1) {
    numeric_partition_t npart1;
    BOOST_FOREACH (const link_t link, get_links(part1)) {
        npart1.insert(npart1.end(), std::make_pair(link, 1.f));
    }
    add_partition(npart1);
}

void PartitionMeanVariance::rem_partition(const partition& part1) {
    numeric_partition_t npart1;
    BOOST_FOREACH (const link_t link, get_links(part1)) {
        npart1.insert(npart1.end(), std::make_pair(link, 1.f));
    }
    rem_partition(npart1);
}

void PartitionMeanVariance::add_partition(const numeric_partition_t& part1) {
    _n++;
    numeric_partition_t delta(part1);
    for (numpart_iter it = _mean.begin(); it != _mean.end(); ++it)
        delta[it->first] -= it->second; // missing keys will be created;
    for (numpart_iter it = delta.begin(); it != delta.end(); ++it)
        _mean[it->first] += (it->second)/_n;
    numeric_partition_t delta2(part1);
    for (numpart_iter it = _mean.begin(); it != _mean.end(); ++it)
        delta2[it->first] -= it->second; // missing keys will be created
    _m2 += map_dot_product(delta.begin(), delta.end(), delta2.begin());
}

void PartitionMeanVariance::rem_partition(const numeric_partition_t& part1) {
    numeric_partition_t delta2(part1);
    for (numpart_iter it = _mean.begin(); it != _mean.end(); ++it) {
        delta2[it->first] -= it->second; // missing keys will be created
        it->second *= _n; // n*mean
    }
    for (numpart_citer it = part1.begin(); it != part1.end(); ++it)
        _mean[it->first] -= it->second;
    _n--;
    numeric_partition_t delta(part1);
    for (numpart_iter it = _mean.begin(); it != _mean.end(); ++it) {
        it->second /= _n;
        delta[it->first] -= it->second;
    }
    _m2 -= map_dot_product(delta.begin(), delta.end(), delta2.begin());
}

float PartitionGelmanRubin::get() const {
    PartitionMeanVariance between_chain;
    float m_chains = _within_chain_mv.size();
    float n_samples = _within_chain_mv[0].n();
    for (size_t i = 1; i < m_chains; ++i)
        BOOST_ASSERT(n_samples == _within_chain_mv[i].n());
    float sum = 0.f;
    for (pmv_iter it = _within_chain_mv.begin(); it != _within_chain_mv.end(); ++it) {
        between_chain.add_partition(it->mean());
        sum += it->var();
    }
    float W = sum/m_chains;
    float B_by_n = between_chain.var();
    float var_estimate = (n_samples - 1.f)/n_samples * W + B_by_n;
    return sqrt(var_estimate/W);
}

} // namespace biggles
