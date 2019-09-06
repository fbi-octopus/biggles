#include <boost/format.hpp>
#include <cmath>
#include "biggles/detail/fun.hpp"
#include "tools/debug.hpp"
#include <time.h>
namespace biggles {

    using boost::format;
    using boost::str;

    const std::string SI_PREFIX(" kMGTPEZY");

    int num_digits(const float number) {
        return int(std::floor(std::log10(std::abs(number))+1));
    }

    /// \brief converts a size_t to a string using SI prefixes
    std::string count2str(size_t cnt) {
        if (cnt < 1000) {
            return str(format("%6d") % cnt);
        }
        size_t ind = 0;
        float num(cnt);
        while (num >= 1000.f) {
            num /= 1000.f;
            ++ind;
        }
        return str(format("%5.1f%c") % num % SI_PREFIX[ind]);
    }

    /// \brief converts a float to a string using SI prefixes
    std::string float2str(float num) {
        if (std::abs(num) < 1000.f) {
            return str(format("%7.1f") % num);
        }
        size_t ind = 0;
        while (std::abs(num) >= 1000.f) {
            num /= 1000.f;
            ++ind;
        }
        std::string formstr(str(format("%%6.%df%%c") % (4 - num_digits(num))));
        return str(format(formstr) % num % SI_PREFIX[ind]);
    }

    std::string time2str(float seconds) {
        size_t sec = static_cast<size_t>(seconds);
        size_t minutes = sec/size_t(60);
        size_t hours = minutes/size_t(60);
        return str(format("%3d:%02d:%02d") % hours % (minutes % 60) % (sec % 60));
    }

    /// \brief returns the local time as a string formated by \e strftime
    std::string local_time(const std::string& fmt) {
        time_t rawtime;
        struct tm * timeinfo;
        char buffer [160];
        time (&rawtime);
        timeinfo = localtime (&rawtime);
        size_t n = strftime (buffer,160, fmt.c_str(),timeinfo);
        return std::string(buffer).substr(0, n);
    }
}
