/// @file io.hpp Stream input/output support for Biggles data structures.

#ifndef BIGGLES_IO_HPP__
#define BIGGLES_IO_HPP__

#include "model.hpp"
#include "partition.hpp"
#include <iostream>

namespace biggles {

/// @brief Serialisation of Biggles data structures to disk
namespace io {

/// @brief Serialise a partition as a JSON document.
///
/// The JSON has no extraneous white space and so you may want to run it through a separate pretty-printer if you want
/// to make it suitable for human consumption.
///
/// @param os The stream to write JSON to.
/// @param partition The partition to serialised.
void write_partition_to_json_stream(std::ostream& os,
                const biggles::partition& partition);

/// @brief Serialise a partition from a JSON formatted stream
///
/// @param is The stream to serialise from
///
/// @return a shared pointer to the new partition
boost::shared_ptr<biggles::partition> read_partition_from_json_stream(std::istream& is);

/// @brief Serialise a set of model parameters as a JSON document.
///
/// The JSON has no extraneous white space and so you may want to run it through a separate pretty-printer if you want
/// to make it suitable for human consumption.
///
/// @param os The stream to write JSON to.
/// @param model_parameters The model_parameters to serialised.
void write_model_parameters_to_json_stream(std::ostream& os,
                const biggles::model::parameters& model_parameters);

/// @brief Serialise a set of model parameters from a JSON formatted stream
///
/// @param is The stream to serialise from
///
/// @return a shared pointer to the new model parameters
boost::shared_ptr<biggles::model::parameters>
    read_model_parameters_from_json_stream(std::istream& is);

}

}

#endif // BIGGLES_IO_HPP__
