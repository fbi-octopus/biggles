#include "jansson_wrapper.hpp"

#include <cstdarg>
#include <cstdlib>
#include <stdexcept>
#include <sstream>
#include <utility>

#include <jansson.h>

namespace jansson
{

value read(std::istream& is)
{
    json_error_t err;
    std::ostringstream ss;
    ss << is.rdbuf();

    json_t* v = json_loads(ss.str().c_str(), 0, &err);
    if(NULL == v) {
        throw std::runtime_error(std::string("Error reading JSON: ") + err.text);
    }

    return value(v);
}

void write(std::ostream& os, const value& v)
{
    if(!v || !json_is_object(v.get())) {
        throw std::invalid_argument("expected a JSON object");
    }

    char* str = json_dumps(v.get(), 0);
    os << str;
    free(str);
}

void vunpack(value& v, const char* fmt, va_list ap)
{
    if(!v) {
        throw std::invalid_argument("Passed empty value to unpack");
    }

    json_error_t err;
    if(0 != json_vunpack_ex(v.get(), &err, JSON_STRICT, fmt, ap)) {
        throw std::runtime_error(std::string("Error unpacking JSON: ") + err.text);
    }
}

void unpack(value& v, const char* fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vunpack(v, fmt, ap);
}

}

