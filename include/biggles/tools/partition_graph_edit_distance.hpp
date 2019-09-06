/// @file partition_graph_edit_distance.hpp A graph edit distance for track partitions

#ifndef PARTITION_GRAPH_EDIT_DISTANCE_HPP__
#define PARTITION_GRAPH_EDIT_DISTANCE_HPP__

#include "../partition.hpp"

namespace biggles {

typedef std::pair<observation, observation> link_t;
bool operator<(const link_t& l1, const link_t& l2);
bool operator==(const link_t& l1, const link_t& l2);

typedef std::set<link_t> link_set_t;
typedef link_set_t::iterator link_set_iter;
typedef std::pair<link_set_iter, bool> insert_result_t;

/** \brief calculates the graph edit distance between partitions
 *
 *  graph definition: \n
 *  \b vertex - observation \n
 *  \b edge - link \n
 *
 *  edit operations: \n
 *  <b> vertex insertion </b> - <em> should not appear </em>; weight 1 \n
 *  <b> vertex deletion </b> - <em> should not appear </em>; weight 1 \n
 *  <b> vertex substition </b> - <em> not applicable </em> \n
 *
 *  <b> edge insertion </b> - weight 1 \n
 *  <b> edge deletion </b> - weight 1 \n
 *  <b> edg substition </b> - not applicable \n
 *
 */
float partition_ged(const partition& p1, const partition&p2);

typedef std::map<link_t, float> numeric_partition_t; ///< \brief used to do statistics for partitions
typedef numeric_partition_t::iterator numpart_iter;
typedef numeric_partition_t::const_iterator numpart_citer;

/** Mean and sample variance for tracking partitions
 *
 */
class PartitionMeanVariance {
    float _n;
    numeric_partition_t _mean;
    float _m2;
public:
    PartitionMeanVariance() {}
    /// \brief adds a partition to the samples
    void add_partition(const partition& p1);
    /// \brief adds a "numeric" partition to the samples
    void add_partition(const numeric_partition_t& p1);
    /// \brief removes a partition to the samples
    void rem_partition(const partition& p1);
    /// \brief removes a "numeric" partition to the samples
    void rem_partition(const numeric_partition_t& p1);
    /// \brief the number of partitions added
    size_t n() const { return _n; }
    const numeric_partition_t& mean() const { return _mean; }
    /// \brief the sample variance
    float var() const {
        if (_n < 2.f)
            return 1.f/0.f;
        return _m2/(_n-1);
    }
};

class PartitionGelmanRubin {
    std::deque<PartitionMeanVariance> _within_chain_mv;
    typedef std::deque<PartitionMeanVariance>::const_iterator pmv_iter;
public:
    explicit PartitionGelmanRubin(int num_chains = 2) : _within_chain_mv(num_chains) {}
    void add_partition(int index, const partition& part) { _within_chain_mv.at(index).add_partition(part); }
    void rem_partition(int index, const partition& part) { _within_chain_mv.at(index).rem_partition(part); }
    float get() const;
};

} // namespace biggles;

#endif
