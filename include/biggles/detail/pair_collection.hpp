/// @file pair_collection.hpp defintions for the pair collection that are used for mergable tracks

#ifndef BIGGLES_DETAIL_PAIR_COLLECTION_HPP___
#define BIGGLES_DETAIL_PAIR_COLLECTION_HPP___

#include <boost/assert.hpp>
#include <boost/foreach.hpp>
#include <deque>
#include <iterator>
#include <map>
#include <set>

namespace biggles {
namespace detail {

template<typename Value>
class pair_collection_iterator : public std::iterator<std::forward_iterator_tag, const std::pair<Value, Value> >
{
public:
    typedef Value value;

    typedef std::pair<value, value> value_pair;

    pair_collection_iterator() { }

    pair_collection_iterator(const pair_collection_iterator& i)
        : map_iterator_(i.map_iterator_)
        , map_end_(i.map_end_)
        , set_iterator_(i.set_iterator_)
        , value_pair_(i.value_pair_)
    {
        bump_up();
    }

    pair_collection_iterator(typename std::map<value, std::set<value> >::const_iterator map_iterator,
                             typename std::map<value, std::set<value> >::const_iterator map_end,
                             typename std::set<value>::const_iterator set_iterator)
        : map_iterator_(map_iterator)
        , map_end_(map_end)
        , set_iterator_(set_iterator)
    {
        bump_up();
    }

    const pair_collection_iterator& operator = (const pair_collection_iterator& i)
    {
        map_iterator_ = i.map_iterator_;
        map_end_ = i.map_end_;
        set_iterator_ = i.set_iterator_;
        bump_up();
    }

    const pair_collection_iterator& operator ++ ()
    {
        if(map_iterator_ == map_end_)
            return *this;

        ++set_iterator_;
        bump_up();

        return *this;
    }

    bool operator == (const pair_collection_iterator& i)
    {
        return (i.map_iterator_ == map_iterator_) && (i.set_iterator_ == set_iterator_);
    }

    bool operator != (const pair_collection_iterator& i)
    {
        return (i.map_iterator_ != map_iterator_) || (i.set_iterator_ != set_iterator_);
    }

    const std::pair<value, value>& operator * ()
    {
        if(map_iterator_ != map_end_)
            value_pair_ = std::make_pair(map_iterator_->first, *set_iterator_);
        return value_pair_;
    }

    const std::pair<value, value>* operator -> ()
    {
        if(map_iterator_ != map_end_)
            value_pair_ = std::make_pair(map_iterator_->first, *set_iterator_);
        return &value_pair_;
    }

protected:
    typename std::map<value, std::set<value> >::const_iterator map_iterator_;
    typename std::map<value, std::set<value> >::const_iterator map_end_;
    typename std::set<value>::const_iterator set_iterator_;
    std::pair<value, value> value_pair_;

    void bump_up()
    {
        while(1)
        {
            if(map_iterator_ == map_end_)
                return;

            if(set_iterator_ == map_iterator_->second.end())
            {
                ++map_iterator_;
                if(map_iterator_ == map_end_)
                {
                    set_iterator_ = typename std::set<value>::const_iterator();
                }
                else
                {
                    set_iterator_ = map_iterator_->second.begin();
                }
            }
            else
            {
                return;
            }
        }
    }
};

template<typename Value>
class pair_collection
{
public:
    typedef Value value;

    typedef std::pair<const value, const value> value_pair;

private:
    struct pair_compare {
        bool operator() (const value_pair& a, const value_pair& b) const {
            return a.first < b.first or (a.first == b.first and a.second < b.second);
        }
    };

protected:
    typedef std::multiset<value_pair, pair_compare> value_pair_multiset;

public:
    typedef typename value_pair_multiset::const_iterator const_iterator;

    typedef typename value_pair_multiset::const_iterator iterator;

    pair_collection() { }

    pair_collection(const pair_collection& pc)
        : pairs_(pc.pairs_)
    { }

    const pair_collection& operator = (const pair_collection& pc)
    {
        pairs_ = pc.pairs_;
    }

    template<typename InputIterator>
    pair_collection(InputIterator first,
                    InputIterator last);

    template<typename InputIterator>
    void insert(InputIterator first, InputIterator last);

    void insert(const value_pair& pair);

    void erase_pairs_with_element(const value& v);

    bool contains_pair_with_element(const value& v) const;

    size_t size() const { return pairs_.size(); }

    bool empty() const { return pairs_.empty(); }

    const_iterator begin() const { return pairs_.begin(); }

    const_iterator end() const { return pairs_.end(); }

protected:
    value_pair_multiset pairs_;
};

template<typename Value>
template<typename InputIterator>
pair_collection<Value>::pair_collection(InputIterator first,
                                        InputIterator last)
{
    insert(first, last);
}

template<typename Value>
template<typename InputIterator>
void pair_collection<Value>::insert(InputIterator first,
                                    InputIterator last)
{
    BOOST_FOREACH(const value_pair& p, std::make_pair(first, last))
    {
        insert(p);
    }
}

template<typename Value>
void pair_collection<Value>::insert(const pair_collection<Value>::value_pair& pair)
{
    value_pair pair_to_insert(
        (pair.first < pair.second)
        ? std::make_pair(pair.first, pair.second)
        : std::make_pair(pair.second, pair.first));

    BOOST_FOREACH(const value_pair& p, pairs_.equal_range(pair_to_insert))
    {
        if(p == pair_to_insert)
            return;
    }
    pairs_.insert(pair_to_insert);
}

template<typename Value>
bool pair_collection<Value>::contains_pair_with_element(const pair_collection<Value>::value& v) const
{
    BOOST_FOREACH(const value_pair& p, pairs_)
    {
        if(p.first == v)
            return true;
        if(p.second == v)
            return true;
    }

    return false;
}

template<typename Value>
void pair_collection<Value>::erase_pairs_with_element(const value& v)
{
    for(;;) {
        // This erases a single pair with v as first or second
        const_iterator i(pairs_.begin());
        for(; (i != pairs_.end()) && (i->first != v) && (i->second != v); ++i) {
            /* nop */
        }
        if(i != pairs_.end()) {
            pairs_.erase(i);
        } else {
            break;
        }
    }
}

} }

#endif // BIGGLES_DETAIL_PAIR_COLLECTION_HPP___
