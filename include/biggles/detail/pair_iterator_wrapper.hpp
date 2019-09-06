/// @file pair_iterator_wrapper.hpp A wrapper iterator which transforms pairs into their second value

#ifndef BIGGLES_DETAIL_PAIR_ITERATOR_WRAPPER_HPP__
#define BIGGLES_DETAIL_PAIR_ITERATOR_WRAPPER_HPP__

namespace biggles { namespace detail
{

/// @brief A BidirectionalIterator wrapper around an iterator yielding pairs.
///
/// This iterator hides the pair's first value and presents only the second.
template<typename BidirectionalIterator, typename Value>
class pair_iterator_wrapper : public std::iterator<std::bidirectional_iterator_tag, Value>
{
public:
    pair_iterator_wrapper() { }
    pair_iterator_wrapper(const pair_iterator_wrapper& i) : it_(i.it_) { }
    explicit pair_iterator_wrapper(const BidirectionalIterator& it) : it_(it) { }
    const pair_iterator_wrapper& operator = (const pair_iterator_wrapper& i) { it_ = i.it_; return *this; }
    bool operator == (const pair_iterator_wrapper& i) { return it_ == i.it_; }
    bool operator != (const pair_iterator_wrapper& i) { return it_ != i.it_; }
    Value& operator * () { return it_->second; }
    Value* operator -> () { return &(it_->second); }
    pair_iterator_wrapper<BidirectionalIterator, Value>& operator ++ () { ++it_; return *this; }
    pair_iterator_wrapper<BidirectionalIterator, Value>& operator -- () { --it_; return *this; }

    pair_iterator_wrapper<BidirectionalIterator, Value> operator ++ (int) {
        pair_iterator_wrapper<BidirectionalIterator, Value> rv(*this);
        ++it_; return rv;
    }

    pair_iterator_wrapper<BidirectionalIterator, Value> operator -- (int) {
        pair_iterator_wrapper<BidirectionalIterator, Value> rv(*this);
        --it_; return rv;
    }

protected:
    BidirectionalIterator it_;
};

} }

#endif // BIGGLES_DETAIL_PAIR_ITERATOR_WRAPPER_HPP__
