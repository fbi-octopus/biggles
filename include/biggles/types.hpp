/// @file types.hpp Simple typedefs used by Biggles

#ifndef BIGGLES_TYPES_HPP__
#define BIGGLES_TYPES_HPP__

#include <Eigen/Dense>

namespace biggles {

typedef Eigen::Matrix2f matrix2f;
typedef Eigen::Matrix4f matrix4f;
typedef Eigen::Matrix<float, 2, 4> matrix2x4f;
typedef Eigen::Matrix<float, 4, 2> matrix4x2f;

typedef Eigen::Vector4f state_t; // the underlying state vector type: [x, x', y, y']'





namespace model {

/// @brief A tuple representing the model parameters.
///
/// The model parameters are represented by the tuple \f$ ( \lambda_b, \lambda_f, p_s, p_d, R ) \f$ where
///
/// - \f$ \lambda_b \f$: the mean number of new tracks appearing per frame.
/// - \f$ \lambda_f \f$: the mean number of false observations per frame.
/// - \f$ p_s \f$: the probability that a target will survive from frame \f$ t \f$ to \f$ t+1 \f$.
/// - \f$ p_d \f$: the probability that a target will generate an observation.
/// - \f$ r \f$: the constraint radius.
/*
typedef boost::tuples::tuple<float, float, float, float, matrix2f, float, matrix4f> parameters;
*/

struct parameters {
    float birth_rate;
    float clutter_rate;
    float survival_probability;
    float observation_probability;
    matrix2f observation_error_covariance;
    float constraint_radius;
    matrix4f process_noise_covariance;
    bool operator==(const parameters &other) const {
        return
        other.birth_rate == birth_rate and
        other.clutter_rate == clutter_rate and
        other.survival_probability == survival_probability and
        other.observation_probability == observation_probability and
        other.constraint_radius == constraint_radius and
        other.observation_error_covariance == observation_error_covariance and
        other.process_noise_covariance == process_noise_covariance;
    }
};

}

}

#endif // BIGGLES_MODEL_HPP__
