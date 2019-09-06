/// @file jansson_wrapper.hpp
/// @brief A light RAII wrapper for Jansson

#include <cstdarg>
#include <iostream>

#include <jansson.h>

namespace jansson
{

/// @brief A RAII wrapper around Jansson
///
/// The Jansson JSON library had a C interface with a reference-counting memory management solution. This is all well
/// and good for C but isn't very exception safe. This wrapper is a thin shim which re-casts the reference count
/// solution as a Resource Acquisition Is Initialisation pattern. Like <code>boost::shared_ptr<></code>, one has to be
/// careful about how one constructs objects of this class; the default constructor will 'steal' the reference. This is
/// to allow one to keep the same mental model with <code>json_t*</code> values as one has with naked pointers
/// constructed by <code>new</code>. For example, to construct a new JSON array, one can do:
///
/// \code{.cpp}
/// jansson::value empty_array(json_array());
/// jansson::value array(json_pack("[ii]", 1, 2));
/// \endcode
///
/// Jansson APIs which return 'borrowed' references, such as <code>json_array_get()</code> should be wrapped in
/// <code>json_incref</code>:
///
/// \code{.cpp}
/// jansson::value array(json_pack("[ii]", 1, 2));
/// jansson::value second_elem(json_incref(json_array_get(array.get(), 1)));
/// \endcode
///
/// As always, it is best to use a tool like Valgrind to check that you are not leaving dangling references.
class value
{
private:
    mutable json_t* wrapped_;
public:
    value() : wrapped_(NULL) {}
    explicit value(json_t* v) : wrapped_(v) { }
    value(const value& v) : wrapped_(v.wrapped_) { json_incref(wrapped_); }
    ~value() { clear(); }

    const value& operator = (const value& v)
    {
        clear();
        wrapped_ = v.wrapped_;
        json_incref(wrapped_);
        return *this;
    }

    void swap(value& other)
    {
        std::swap(other.wrapped_, wrapped_);
    }

    void clear() throw ()
    {
        if(NULL != wrapped_) { json_decref(wrapped_); }
        wrapped_ = NULL;
    }

    operator bool () const { return (NULL != wrapped_); }

    json_t* get() const { return wrapped_; }
};

value read(std::istream& is);

void write(std::ostream& os, const value& v);

void vunpack(value& v, const char* fmt, va_list ap);

void unpack(value& v, const char* fmt, ...);

}

