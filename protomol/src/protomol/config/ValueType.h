/* -*- c++ -*- */
#ifndef VALUETYPE_H
#define VALUETYPE_H

#include <protomol/type/AbstractEnumType.h>
#include <protomol/type/Vector3D.h>

#include <protomol/config/ConstraintValueType.h>

namespace ProtoMol {
  //_____________________________________________________ ValueEnum
  /**
   * Implements a map of enum and string for the types of Value.
   */
  class ValueEnum {
public:
    virtual ~ValueEnum() {}
    enum Enum {
      FIRST = 0,        // Used internally only
      UNDEFINED = 0,   // Value returned when no string matches
      STRING,
      INT,
      UINT,
      REAL,
      BOOL,
      VECTOR3D,
      VECTOR,
      INTEGRATOR,       // Extended types start here, which act like a STRING
      FORCE,
      LAST              // Used internally only
    };
public:
    /**
     * Defines types for extended types, which do not have
     * a corresponding and natural/fundamental type.
     * The mapping includes a string and the std::string conversion,
     * which enables such objects to act like a string.@n
     *
     * Modern C++ Design, p. 29.
     */
    template<Enum e>
    struct Enum2ValueType {
      explicit Enum2ValueType(std::string s) : str(s) {}
      explicit Enum2ValueType(char const *c) : str(std::string(c)) {}

      operator std::string() const {return str;}
      operator Enum() const {return e;}

      std::string str;
      enum {value = e};
    };

public:
    typedef Enum2ValueType<INTEGRATOR> Integrator;
    typedef Enum2ValueType<FORCE> Force;



protected:
    static const std::string str[];
  };

  //_____________________________________________________ ValueType
  typedef AbstractEnumType<ValueEnum> ValueType;

  //_____________________________________________ ValueTraits & Enum2ValueTraits
  /**
   * ValueTraits and Enum2ValueTraits define mapping and type implementation
   * behavior used by Value.
   */
  template<typename T>
  struct ValueTraits;

  /**
   * Trait class defining the actual type, enum value, the representation type,
   * possible constraints, and implements the conversions.@n
   *
   * Mapping enum value to type (ValueTraits<>)
   * Modern C++ Design, p. 29.
   */
  template<ValueType::Enum e>
  struct Enum2ValueTraits;

  /// STRING
  template<>
  struct ValueTraits<std::string> {
    typedef std::string Type;          // Type, template parameter
    typedef std::string RepType;       // Representation type of the value
    enum {value = ValueType::STRING};  // Enum value

    static RepType init() {return "";}
    static Type empty() {return "";}
    static unsigned int size(const RepType &v) {return v.size();}
    static bool read(std::istream &is, RepType &v) {is >> v; return true;}

    template<typename T> static bool convert(T v, RepType &r)
    {r = toString(v); return true;}                    // !NB everything is ok
    template<typename T> static RepType convert(T v)
    {RepType tmp; convert(v, tmp); return tmp;}
    template<typename T> static bool convertTest(T v)
    {RepType tmp; return convert(v, tmp);}

    // Constraints, defined in ConstraintValueType
    static bool check(const ConstraintValueType::NoConstraints &,
                      const RepType &) {return true;}
    static bool check(const ConstraintValueType::NotEmpty &,
                      const RepType &v) {return v != "";}
  };
  template<>
  struct Enum2ValueTraits<ValueType::STRING> {
    // Mapping enum value to type (ValueTraits<>)
    typedef ValueTraits<std::string> Type;
  };

  template<int n>
  struct ValueTraits < const char[n] > : public ValueTraits<std::string> {};
  template<int n>
  struct ValueTraits < char[n] > : public ValueTraits<std::string> {};


  /// INT
  template<>
  struct ValueTraits<int> {
    typedef int Type;
    typedef int RepType;
    enum {value = ValueType::INT};

    static RepType init() {return 0;}
    static Type empty() {return 0;}
    static unsigned int size(RepType) {return 1;}
    static bool read(std::istream &is, RepType &v)
    {std::string tmp; is >> tmp; return toInt(tmp, v);}

    template<typename T> static bool convert(T, RepType &) {return false;}
    template<typename T> static RepType convert(T v)
    {RepType tmp; convert(v, tmp); return tmp;}
    template<typename T> static bool convertTest(T v)
    {RepType tmp; return convert(v, tmp);}

    static bool convert(RepType v, RepType &r) {r = v; return true;}
    static bool convert(std::string v, RepType &r) {return toInt(v, r);}
    static bool convert(bool v, RepType &r) {r = (v ? 1 : 0); return true;}
    static bool convert(unsigned int v, RepType &r) {
      int c = static_cast<int>(v);
      if (static_cast<unsigned int>(c) == v) {
        r = c;
        return true;
      }
      return false;
    }
    static bool convert(Real v, RepType &r) {
      int c = static_cast<int>(v);
      if (static_cast<Real>(c) == v) {
        r = c;
        return true;
      }
      return false;
    }

    static bool check(const ConstraintValueType::NoConstraints &,
                      RepType) {return true;}
    static bool check(const ConstraintValueType::Zero &,
                      RepType v) {return v == 0;}
    static bool check(const ConstraintValueType::NotZero &,
                      RepType v) {return v != 0;}
    static bool check(const ConstraintValueType::Positive &,
                      RepType v) {return v > 0;}
    static bool check(const ConstraintValueType::Negative &,
                      RepType v) {return v < 0;}
    static bool check(const ConstraintValueType::NotPositive &,
                      RepType v) {return v <= 0;}
    static bool check(const ConstraintValueType::NotNegative &,
                      RepType v) {return v >= 0;}
  };

  template<>
  struct Enum2ValueTraits<ValueType::INT> {
    typedef ValueTraits<int> Type;
  };


  /// UINT
  template<>
  struct ValueTraits<unsigned int> {
    typedef unsigned int Type;
    typedef unsigned int RepType;
    enum {value = ValueType::UINT};

    static RepType init() {return 0;}
    static Type empty() {return 0;}
    static unsigned int size(RepType) {return 1;}
    static bool read(std::istream &is, RepType &v)
    {std::string tmp; is >> tmp; return toUInt(tmp, v);}

    template<typename T> static bool convert(T, RepType &) {return false;}
    template<typename T> static RepType convert(T v) 
    {RepType tmp; convert(v, tmp); return tmp;}
    template<typename T> static bool convertTest(T v)
    {RepType tmp; return convert(v, tmp);}

    static bool convert(RepType v, RepType &r) {r = v; return true;}
    static bool convert(std::string v, RepType &r) {return toUInt(v, r);}
    static bool convert(bool v, RepType &r) {r = (v ? 1 : 0); return true;}
    static bool convert(int v, RepType &r) {
      if (v >= 0) {
        r = static_cast<unsigned int>(v);
        return true;
      }
      return false;
    }
    static bool convert(Real v, RepType &r) {
      unsigned int c = static_cast<unsigned int>(v);
      if (static_cast<Real>(c) == v) {
        r = c;
        return true;
      }
      return false;
    }

    static bool check(const ConstraintValueType::NoConstraints &,
                      RepType) {return true;}
    static bool check(const ConstraintValueType::NotZero &,
                      RepType v) {return 0U != v;}
    static bool check(const ConstraintValueType::Zero &,
                      RepType v) {return 0U == v;}
    static bool check(const ConstraintValueType::Positive &,
                      RepType v) {return v > 0U;}
  };

  template<>
  struct Enum2ValueTraits<ValueType::UINT> {
    typedef ValueTraits<unsigned int> Type;
  };


  /// REAL
  template<>
  struct ValueTraits<Real> {
    typedef Real Type;
    typedef Real RepType;
    enum {value = ValueType::REAL};

    static RepType init() {return 0.0;}
    static Type empty() {return 0.0;}
    static unsigned int size(RepType) {return 1;}
    static bool read(std::istream &is, RepType &v)
    {std::string tmp; is >> tmp; return toReal(tmp, v);}

    template<typename T> static bool convert(T, RepType &) {return false;}
    template<typename T> static RepType convert(T v)
    {RepType tmp; convert(v, tmp); return tmp;}
    template<typename T> static bool convertTest(T v)
    {RepType tmp; return convert(v, tmp);}

    static bool convert(RepType v, RepType &r) {r = v; return true;}
    static bool convert(std::string v, RepType &r) {return toReal(v, r);}
    static bool convert(int v, RepType &r)
    {r = static_cast<Real>(v); return true;}
    static bool convert(unsigned int v, RepType &r)
    {r = static_cast<Real>(v); return true;}
    static bool convert(bool v, RepType &r) 
    {r = (v ? 1.0 : 0.0); return true;}

    static bool check(const ConstraintValueType::NoConstraints &,
                      RepType) {return true;}
    static bool check(const ConstraintValueType::Zero &,
                      RepType v) {return v == 0.0;}
    static bool check(const ConstraintValueType::NotZero &,
                      RepType v) {return v != 0.0;}
    static bool check(const ConstraintValueType::Positive &,
                      RepType v) {return v > 0.0;}
    static bool check(const ConstraintValueType::Negative &,
                      RepType v) {return v < 0.0;}
    static bool check(const ConstraintValueType::NotPositive &,
                      RepType v) {return v <= 0.0;}
    static bool check(const ConstraintValueType::NotNegative &,
                      RepType v) {return v >= 0.0;}
  };

  template<>
  struct Enum2ValueTraits<ValueType::REAL> {
    typedef ValueTraits<Real> Type;
  };


  /// BOOL
  template<>
  struct ValueTraits<bool> {
    typedef bool Type;
    typedef bool RepType;
    enum {value = ValueType::BOOL};

    static RepType init() {return true;}
    static Type empty() {return true;}
    static unsigned int size(RepType) {return 1;}
    static bool read(std::istream &is, RepType &v)
    {std::string tmp; is >> tmp; return toBool(tmp, v);}

    template<typename T> static bool convert(T, RepType &) {return false;}
    template<typename T> static RepType convert(T v)
    {RepType tmp; convert(v, tmp); return tmp;}
    template<typename T> static bool convertTest(T v)
    {RepType tmp; return convert(v, tmp);}

    static bool convert(RepType v, RepType &r) {r = v; return true;}
    static bool convert(int v, RepType &r) {r = (v != 0); return true;}
    static bool convert(unsigned int v, RepType &r) {r = (v != 0); return true;}
    static bool convert(Real v, RepType &r) {r = (v != 0.0); return true;}
    static bool convert(std::string v, RepType &r) {return toBool(v, r);}

    static bool check(const ConstraintValueType::NoConstraints &, RepType)
    {return true;}
  };

  template<>
  struct Enum2ValueTraits<ValueType::BOOL> {
    typedef ValueTraits<bool> Type;
  };


  /// VECTOR3D
  template<>
  struct ValueTraits<Vector3D> {
    typedef Vector3D Type;
    typedef Vector3D RepType;
    enum {value = ValueType::VECTOR3D};

    static RepType init() {return Vector3D(0.0, 0.0, 0.0);}
    static Type empty() {return Vector3D(0.0, 0.0, 0.0);}
    static unsigned int size(const RepType &) {return 3;}
    static bool read(std::istream &is, RepType &v) {
      std::string x, y, z;
      is >> x >> y >> z;
      return toVector3D(x + " " + y + " " + z, v);
    }

    template<typename T> static bool convert(T, RepType &) {return false;}
    template<typename T> static RepType convert(T v)
    {RepType tmp; convert(v, tmp); return tmp;}
    template<typename T> static bool convertTest(T v)
    {RepType tmp; return convert(v, tmp);}

    static bool convert(RepType v, RepType &r) {r = v; return true;}
    static bool convert(std::string v, RepType &r) {return toVector3D(v, r);}

    static bool check(const ConstraintValueType::NoConstraints &,
                      const RepType &) {return true;}
    static bool check(const ConstraintValueType::NotZero &,
                      const RepType &v) {return v.normSquared() != 0;}
    static bool check(const ConstraintValueType::Zero &,
                      const RepType &v) {return v.normSquared() == 0;}
  };

  template<>
  struct Enum2ValueTraits<ValueType::VECTOR3D> {
    typedef ValueTraits<Vector3D> Type;
  };


  /// VECTOR
  template<>
  struct ValueTraits<std::vector<Real> > {
    typedef std::vector<Real> Type;
    typedef std::vector<Real> RepType;
    enum {value = ValueType::VECTOR};

    static RepType init() {return std::vector<Real>();}
    static Type empty() {return std::vector<Real>();}
    static unsigned int size(const RepType &v) {return v.size();}
    static bool read(std::istream &is, RepType &v) {
      v.clear();
      std::string str;
      is >> str;
      if (isReal(str)) {
        v.push_back(toReal(str));
        while (is >> str) {
          if (!isReal(str))
            break;
          v.push_back(toReal(str));
        }
      } else if (str.size() > 2 && str[0] == '-' && str[1] == '-' &&
                 isUInt(str.substr(2))) {
        unsigned int n = toUInt(str.substr(2));
        for (unsigned int i = 0; i < n; i++) {
          if (!(is >> str))
            return false;
          if (!isReal(str)) {
            is.seekg((-1) * static_cast<int>(str.size()), std::ios::cur);
            is.clear();
            return false;
          }
          v.push_back(toReal(str));
        }

        return true;
      }
      is.seekg((-1) * static_cast<int>(str.size()), std::ios::cur);
      is.clear();
      return true;
    }

    template<typename T> static bool convert(T, RepType &) {return false;}
    template<typename T> static RepType convert(T v)
    {RepType tmp; convert(v, tmp); return tmp;}
    template<typename T> static bool convertTest(T v)
    {RepType tmp; return convert(v, tmp);}

    static bool convert(RepType v, RepType &r) {r = v; return true;}
    static bool convert(std::string v, RepType &r) {return toVector(v, r);}

    static bool check(const ConstraintValueType::NoConstraints &,
                      const RepType &) {return true;}
    static bool check(const ConstraintValueType::NotEmpty &,
                      const RepType &v) {return v.size() > 0;}
    static bool check(const ConstraintValueType::Empty &,
                      const RepType &v) {return v.empty();}
    static bool check(const ConstraintValueType::Zero &, const RepType &v) {
      for (unsigned int i = 0; i < v.size(); i++)
        if (v[i] != 0.0) return false;
      return true;
    }
    static bool check(const ConstraintValueType::NotZero &, const RepType &v) {
      for (unsigned int i = 0; i < v.size(); i++)
        if (v[i] == 0.0) return false;
      return true;
    }
    static bool check(const ConstraintValueType::Positive &,
                      const RepType &v) {
      for (unsigned int i = 0; i < v.size(); i++)
        if (v[i] <= 0.0) return false;
      return true;
    }
    static bool check(const ConstraintValueType::Negative &,
                      const RepType &v) {
      for (unsigned int i = 0; i < v.size(); i++)
        if (v[i] >= 0.0) return false;
      return true;
    }
    static bool check(const ConstraintValueType::NotPositive &,
                      const RepType &v) {
      for (unsigned int i = 0; i < v.size(); i++)
        if (v[i] > 0.0) return false;
      return true;
    }
    static bool check(const ConstraintValueType::NotNegative &,
                      const RepType &v) {
      for (unsigned int i = 0; i < v.size(); i++)
        if (v[i] < 0.0) return false;
      return true;
    }
  };

  template<>
  struct Enum2ValueTraits<ValueType::VECTOR> {
    typedef ValueTraits<std::vector<Real> > Type;
  };


  /// INTEGRATOR
  template<>
  struct ValueTraits<ValueType::Integrator> {
    typedef ValueType::Integrator Type;
    typedef std::string RepType;
    enum {value = ValueType::INTEGRATOR};

    static RepType init() {return "";}
    static Type empty() {return Type("");}
    static unsigned int size(const RepType &v) {return v.size();}
    static bool read(std::istream &is, RepType &v) {
      v = "";
      std::string str;
      is >> str;
      if (str != "{") return false;

      int count = 1;
      while (is >> str) {
        if (str == "{") count++;
        else if (str == "}") count--;

        if (count == 0) return true;

        v += (v.empty() ? "" : " ") + str;
      }

      return false;
    }

    template<typename T> static bool convert(T, RepType &) {return false;}
    template<typename T> static RepType convert(T v)
    {RepType tmp; convert(v, tmp); return tmp;}
    template<typename T> static bool convertTest(T v)
    {RepType tmp; return convert(v, tmp);}

    static bool convert(RepType v, RepType &r) {r = v; return true;}

    static bool check(const ConstraintValueType::NoConstraints &,
                      const RepType &) {return true;}
    static bool check(const ConstraintValueType::NotEmpty &,
                      const RepType &v) {return v != "";}
  };

  template<>
  struct Enum2ValueTraits<ValueType::INTEGRATOR> {
    typedef ValueTraits<ValueType::Integrator> Type;
  };

  /// FORCE
  template<>
  struct ValueTraits<ValueType::Force> {
    typedef ValueType::Force Type;
    typedef std::string RepType;
    enum {value = ValueType::FORCE};

    static RepType init() {return "";}
    static Type empty() {return Type("");}
    static unsigned int size(const RepType &v) {return v.size();}
    static bool read(std::istream &is, RepType &v) {
      v = "";
      std::string str;
      bool first = true;
      while (is >> str) {
        if (equalNocase(str, "Force"))  // TOBEFIXED !!!
          if (first) {
            first = false;
            continue;
          } else {
            is.seekg((-1) * static_cast<int>(str.size()), std::ios::cur);
            return true;
          }
        else if (!str.empty() && str[str.size() - 1] == ',')
          if (first) {
            first = false;
            continue;
          } else {
            is.seekg(-1, std::ios::cur);
            v += (v.empty() ? "" : " ") + str;
            v.resize(v.size() - 1);
            return true;
          }
        else if (str == "}") {
          is.seekg(-1, std::ios::cur);
          return true;
        }
        if (first) {
          is.clear(std::ios::failbit);
          return false;
        }
        v += (v.empty() ? "" : " ") + str;
      }

      if (v.empty())
        is.clear(std::ios::failbit);
      else
        is.clear(std::ios::goodbit);
      return !is.fail();
    }

    template<typename T> static bool convert(T, RepType &) {return false;}
    template<typename T> static RepType convert(T v)
    {RepType tmp; convert(v, tmp); return tmp;}
    template<typename T> static bool convertTest(T v)
    {RepType tmp; return convert(v, tmp);}

    static bool convert(RepType v, RepType &r) {r = v; return true;}

    static bool check(const ConstraintValueType::NoConstraints &,
                      const RepType &) {return true;}
    static bool check(const ConstraintValueType::NotEmpty &,
                      const RepType &v) {return v != "";}
  };

  template<>
  struct Enum2ValueTraits<ValueType::FORCE> {
    typedef ValueTraits<ValueType::Force> Type;
  };
}
#endif /* VALUETYPE_H */
