/*  -*- c++ -*-  */
#ifndef PARAMETER_H
#define PARAMETER_H

#include <protomol/config/Value.h>
#include <protomol/type/SimpleTypes.h>

#include <ostream>

namespace ProtoMol {
  //________________________________________________________ Parameter
  struct Parameter {
    /**
     * Container struct for parameters providing wide range
     * of constructors.
     */

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Parameter();

    Parameter(const std::string & k, const Value & val);
    Parameter(const std::string & k, const Value & val, const Value & def);
    template<typename T>
    Parameter(const std::string & k, const Value & val,
      T def) : keyword(k), value(val), defaultValue(val) {
      defaultValue.set(def);
    }
    Parameter(const char *k, const Value & val);
    Parameter(const char *k, const Value & val, const Value & def);
    template<typename T>
    Parameter(const char *k, const Value & val,
      T def) : keyword(std::string(k)), value(val), defaultValue(val) {
      defaultValue.set(def);
    }

    Parameter(const std::string & k, const Value & val, const Text & t);
    Parameter(const std::string & k, const Value & val, const Value & def,
      const Text & t);
    template<typename T>
    Parameter(const std::string & k, const Value & val, T def, const Text &
      t) : keyword(k), value(val), defaultValue(val), text(t.text) {
      defaultValue.set(def);
    }
    Parameter(const char *k, const Value & val, const Text & t);
    Parameter(const char *k, const Value & val, const Value & def, const Text &
      t);
    template<typename T>
    Parameter(const char *k, const Value & val, T def, const Text &
      t) : keyword(std::string(k)), value(val), defaultValue(val),
      text(t.text) {
      defaultValue.set(def);
    }

    std::ostream &print(std::ostream &stream) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /// the keyword of the parameter
    std::string keyword;
    /// the value of the parameter
    Value value;
    /// optional default value of the parameter
    Value defaultValue;
    /// optional help text of the parameter
    std::string text;
  };
}
#endif /* PARAMETER_H */
