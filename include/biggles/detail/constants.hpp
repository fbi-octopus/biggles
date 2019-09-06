#ifndef BIGGLES_DETAIL_CONSTANTS_HPP___
#define BIGGLES_DETAIL_CONSTANTS_HPP___

#include <limits>

namespace biggles { namespace detail {
    /// \brief A very large floating point number
    const float myriad(std::numeric_limits<float>::max());
    /// \brief exp(-0.5)
    const float exp_one_sigma(0.60653065971263342);
    /// \brief exp(-2.0)
    const float exp_two_sigma(0.1353352832366127);
    /// \brief 1/sqrt(2*pi)
    const float gauss_factor(0.3989422804014327f);
}} // biggles::detail

#endif // BIGGLES_DETAIL_CONSTANTS_HPP___
