#ifndef FUN_HPP__
#define FUN_HPP__

#include <string>
#include <deque>
#include <functional>
#include <boost/assert.hpp>

//#include <stdlib.h>
namespace biggles {
    /// \brief x^2
    template<class DTYPE> DTYPE square(const DTYPE& x) { return x*x; }
    template<class DTYPE> DTYPE cubed(const DTYPE& x) { return x*x*x; }
    template<class DTYPE> DTYPE abs_diff(const DTYPE& x, const DTYPE& y) {
        return x>y ? x-y : y-x;
    }

    std::string local_time(const std::string& fmt = std::string("%y%m%d_%T"));

    std::string count2str(size_t cnt);
    std::string float2str(float num);
    std::string time2str(float seconds);

    /// \brief a class to hold a fixed number of values
    ///
    /// if a new value is inserted into the full container the last (and oldest) value is removed
    template <class DATATYPE>
    class fqueue {
        typedef typename std::deque<DATATYPE> container_t;
        container_t data_;
        size_t max_size_;
    public:
        explicit fqueue(size_t max_size) : max_size_(max_size) { BOOST_ASSERT(max_size_ > 0); }
        typedef typename container_t::const_iterator const_iterator;
        const_iterator begin() const { return data_.begin(); }
        const_iterator end() const { return data_.end(); }
        void push(const DATATYPE& value) {
            if (max_size_ == data_.size()) data_.pop_back();
            data_.push_front(value);
        }
        bool is_full() const { return max_size_ <= data_.size(); }
        size_t size() const { return data_.size(); }
        DATATYPE last() const { return data_.front(); }
        DATATYPE first() const { return data_.back(); }
    };

    /// \brief this is copied form cplusplus.com and it not required in C++11
    template <class ForwardIterator>
    bool is_sorted (ForwardIterator first, ForwardIterator last) {
      if (first==last) return true;
      ForwardIterator next = first;
      while (++next!=last) {
        if (*next < *first)
          return false;
        ++first;
      }
      return true;
    }

    template<class ITERATOR>
    bool any_weight(ITERATOR b, ITERATOR e) {
        typedef typename ITERATOR::value_type value_type;
        while (b != e) {
            if (*b++ > value_type(0))
                return true;
        }
        return false;
    }

    /// \brief (sum(\em fun(*\em it)) for \em it in [\em b, \em e)) + init
    template<class ITERATOR, class FUN>
    typename FUN::result_type fun_sum(ITERATOR b, ITERATOR e, const FUN& fun, typename FUN::result_type init) {
        while (b != e) init += fun(*b++);
        return init;
    }

    template<class INPUT_ITER, class OUTPUT_ITER, class FUN>
    void fun_partial_sum(INPUT_ITER b, INPUT_ITER e, OUTPUT_ITER ob, const FUN& fun, typename FUN::result_type init) {
        while (b != e) *ob++ = init += fun(*b++);
    }

    /// \brief wrapper of unary function \em FUN for std::partial_sum
    template <class FUN> class fun_sum_operator :
        public std::binary_function<typename FUN::result_type, typename FUN::argument_type, typename FUN::result_type>
    {
        typedef typename FUN::result_type result_type;
        typedef typename FUN::argument_type argument_type;
        const FUN& fun_;
    public:
        fun_sum_operator(const FUN& fun) : fun_(fun) {}
        result_type operator()(const result_type& prev, const argument_type& arg) const { return prev + fun_(arg); }
    }; // fun_sum_operator

template<class EIGEN_MATRIX> bool is_symmetric(const EIGEN_MATRIX &matrix) {
    return matrix == matrix.transpose();
}

template<class EIGEN_MATRIX> float measure_asymmetry(const EIGEN_MATRIX &matrix) {
    return (matrix - matrix.transpose()).determinant();
}

template<class EIGEN_MATRIX> EIGEN_MATRIX enforce_symmetry(const EIGEN_MATRIX &matrix) {
    return (matrix + matrix.transpose())/2.0;
}


} // namespace biggles

#endif
