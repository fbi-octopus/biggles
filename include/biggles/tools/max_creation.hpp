#ifndef BIGGLES_MAX_CREATION_HPP__
#define BIGGLES_MAX_CREATION_HPP__

#include "../partition.hpp"

namespace biggles {
    /// \brief Puts as much clutter as possible into tracks using a randomised greedy algorithm
    bool max_creation(const partition_ptr_t& start_partition_ptr, partition_ptr_t& end_partition_ptr);
    /// \brief Puts as much clutter as possible into tracks using a randomised greedy algorithm. Start from the end.
    bool rev_creation(const partition_ptr_t& start_partition_ptr, partition_ptr_t& end_partition_ptr);
}

#endif // BIGGLES_MAX_CREATION_HPP__
