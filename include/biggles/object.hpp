#ifndef BIGGLES_OBJECT_HPP__
#define BIGGLES_OBJECT_HPP__

#include <string>

namespace biggles {
    /// \brief the base class for all biggles classes
    class biggles_object {
    public:
        /// \brief a python-style descritpion string of the class. It is also used to sort objects
        virtual const std::string str() const { return "biggles"; };
        /// @name Comparison operators
        /// @{
        /// \brief operator less than based on biggles_object::str()
        bool operator< (const biggles_object& other) {  return str() <  other.str();  }
        /// \brief operator less than or equal based on biggles_object::str()
        bool operator<=(const biggles_object& other) {  return str() <= other.str();  }
        /// \brief operator equal based on biggles_object::str()
        //bool operator==(const biggles_object& other) {  return str() == other.str();  }
        /// \brief operator unequal based on biggles_object::str()
        bool operator!=(const biggles_object& other) {  return str() != other.str();  }
        /// @}
    };
}

inline std::ostream& operator<<(std::ostream& os, const biggles::biggles_object& obj) {
  os << obj.str();
  return os;
}

#endif
