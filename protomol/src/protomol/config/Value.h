/*  -*- c++ -*-  */
#ifndef VALUE_H
#define VALUE_H

#include <protomol/config/ValueType.h>
#include <protomol/config/ConstraintValueType.h>
#include <typeinfo>

namespace ProtoMol {
  //________________________________________________________ Value
  /** Represents a value with a type and associated with a constraint.
   *  A Value object can behave as 8 different types,
   *   or be undefined:
   *
   * - string
   * - int
   * - unsigned int
   * - Real
   * - bool
   * - Vector3D
   * - vector<Real>
   * - Special, extended types; acts like string with a type (ValueType),
   *   i.e., integrator, force etc.
   */
  class Value {
    //________________________________________________________ PlaceHolder
    /**
     * Abstract class of the place holder
     */
    class PlaceHolder {
    public:
      virtual ~PlaceHolder() {}

      /// Returns the std::type_info of the value
      virtual const std::type_info &type() const = 0;
      /// Clones the value
      virtual PlaceHolder *clone() const = 0;
      /// True if the value is valid and has a type
      virtual bool valid() const = 0;
      /// Makes the value non-valid/cleared and returns previous state of valid.
      virtual bool clear() = 0;
      /// Initialized the value with the default value of its representation
      /// type and applying its constraints
      virtual bool init() = 0;
      /// Returns the number of elements, non-1 for string, Vector3D, vector
      virtual unsigned int size() const = 0;
      /// Reads the value from std::istream according its type
      virtual void read(std::istream &is) = 0;

      // Type methods
      /// Returns the type
      virtual ValueType::Enum getType() const = 0;
      /// Returns the type of the constraint
      virtual ConstraintValueType::Enum getConstraintType() const = 0;
      /// Returns the type as a string
      virtual const std::string &getConstraintTypeString() const = 0;
      /// Return the type of the constraint as a string
      virtual const std::string &getTypeString() const = 0;

      // Getter
      virtual bool get(std::string &v) const = 0;
      virtual bool get(int &v) const = 0;
      virtual bool get(unsigned int &v) const = 0;
      virtual bool get(Real &v) const = 0;
      virtual bool get(bool &v) const = 0;
      virtual bool get(Vector3D &v) const = 0;
      virtual bool get(std::vector<Real> &v) const = 0;

      // Setter
      bool set(const char *v) {return set(std::string(v));}
      virtual bool set(const std::string &v) = 0;
      virtual bool set(int v) = 0;
      virtual bool set(unsigned int v) = 0;
      virtual bool set(Real v) = 0;
      virtual bool set(bool v) = 0;
      virtual bool set(const Vector3D &v) = 0;
      virtual bool set(const std::vector<Real> &v) = 0;
      virtual bool set(const Value &val) = 0;
    };

    //________________________________________________________ PlaceHolder
    /**
     * Implementation of the placerholder taking a type defined
     * in ValueTrait and a constraint (Constraint)
     */
    template<typename ValueTrait, typename Constraint =
             ConstraintValueType::NoConstraints>
    class Holder : public PlaceHolder {
    public:
      Holder(const typename ValueTrait::RepType &v, bool def = true) :
        holder(def ? v : ValueTrait::init()),
        ok(ValueTrait::check(Constraint(), v) && def) {}
      Holder(const Holder &rhs) : holder(rhs.holder), ok(rhs.ok) {}
      virtual const std::type_info &type() const {return typeid(ValueTrait);}
      virtual PlaceHolder *clone() const {return new Holder(*this);}
      virtual bool valid() const {return ok;}
      virtual bool clear() {bool tmp = ok; ok = false; return tmp;}
      virtual bool init() {
        bool tmp = ok; holder = ValueTrait::init();
        ok = ValueTrait::check(Constraint(), holder);  return tmp;
      }
      virtual unsigned int size() const {return ValueTrait::size(holder);}
      virtual void read(std::istream &is) {
        ok = ValueTrait::read(is, holder) &&
          ValueTrait::check(Constraint(), holder);
      }

      // Type methods
      virtual ValueType::Enum getType() const {
        return (ValueType::Enum)ValueTrait::value;
      }
      virtual ConstraintValueType::Enum getConstraintType() const {
        return (ConstraintValueType::Enum)Constraint::value;
      }
      virtual const std::string &getTypeString() const {
        return ValueType::getString((ValueType::Enum)ValueTrait::value);
      }
      virtual const std::string &getConstraintTypeString() const {
        return ConstraintValueType::
          getString((ConstraintValueType::Enum)Constraint::value);
      }

      // Getter
      virtual bool get(std::string &v) const {
        return ValueTraits<std::string>::convert(holder, v) && ok;
      }
      virtual bool get(int &v) const {
        return ValueTraits<int>::convert(holder, v) && ok;
      }
      virtual bool get(unsigned int &v) const {
        return ValueTraits<unsigned int>::convert(holder, v) && ok;
      }
      virtual bool get(Real &v) const {
        return ValueTraits<Real>::convert(holder, v) && ok;
      }
      virtual bool get(bool &v) const {
        return ValueTraits<bool>::convert(holder, v) && ok;
      }
      virtual bool get(Vector3D &v) const {
        return ValueTraits<Vector3D>::convert(holder, v) && ok;
      }
      virtual bool get(std::vector<Real> &v) const {
        return ValueTraits<std::vector<Real> >::convert(holder, v) && ok;
      }

      // Setter
      virtual bool set(const std::string &v) {
        ok = ValueTrait::convert(v, holder) &&
          ValueTrait::check(Constraint(), holder);
        return ok;
      }
      virtual bool set(int v) {
        ok = ValueTrait::convert(v, holder) &&
          ValueTrait::check(Constraint(), holder);
        return ok;
      }
      virtual bool set(unsigned int v) {
        ok = ValueTrait::convert(v, holder) &&
          ValueTrait::check(Constraint(), holder);
        return ok;
      }
      virtual bool set(Real v) {
        ok = ValueTrait::convert(v, holder) &&
          ValueTrait::check(Constraint(), holder);
        return ok;
      }
      virtual bool set(bool v) {
        ok = ValueTrait::convert(v, holder) &&
          ValueTrait::check(Constraint(), holder);
        return ok;
      }
      virtual bool set(const Vector3D &v) {
        ok = ValueTrait::convert(v, holder) &&
          ValueTrait::check(Constraint(), holder);
        return ok;
      }
      virtual bool set(const std::vector<Real> &v) {
        ok = ValueTrait::convert(v, holder) &&
          ValueTrait::check(Constraint(), holder);
        return ok;
      }
      virtual bool set(const Value &val) {
        ok = val.valid();
        if (ok) {
          typename ValueTrait::RepType v;
          ok = val.get(v);
          if (ok)
            ok = ValueTrait::convert(v, holder) &&
              ValueTrait::check(Constraint(), holder);
        }
        return ok;
      }
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Data members
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    private:
      typename ValueTrait::RepType holder;
      bool ok;
    };

    //________________________________________________________ Value

  public:
    struct Undefined {};
    static const ProtoMol::Value::Undefined *undefined; ///< Undefined value
    //
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Value();
    ~Value();

    // Not all compilers like this ...
    //
    //template<typename T>
    //explicit Value(const T& value):
    //  myValue(new Holder<ValueTraits<T> >(value)){}
    //
    //template<typename T>
    //Value(const T& value, const Undefined*):
    //  myValue(new Holder<ValueTraits<T> >(value,false)){}
    //
    // ... so we expand ourself
    explicit Value(const std::string &value);
    explicit Value(const char *value);
    explicit Value(int value);
    explicit Value(unsigned int value);
    explicit Value(Real value);
    explicit Value(bool value);
    explicit Value(const Vector3D &value);
    explicit Value(const std::vector<Real> &value);
    explicit Value(const ValueType::Integrator &value);
    explicit Value(const ValueType::Force &value);

    Value(const Value &other);

    Value(const std::string &value, const Undefined *);
    Value(const char *value, const Undefined *);
    Value(int value, const Undefined *);
    Value(unsigned int value, const Undefined *);
    Value(Real value, const Undefined *);
    Value(bool value, const Undefined *);
    Value(const Vector3D &value, const Undefined *);
    Value(const std::vector<Real> &value, const Undefined *);
    Value(const ValueType::Integrator &value, const Undefined *);
    Value(const ValueType::Force &value, const Undefined *);

    Value(const Value &value, const Undefined *);

    template<typename T, typename C>
    Value(const T &value, const C &) :
      myValue(new Holder<ValueTraits<T> , C>(value)) {}

    template<typename T, typename C>
    Value(const T &value, const C &, const Undefined *) :
      myValue(new Holder<ValueTraits<T> , C>(value, false)) {}

    template<typename T, typename C>
    Value(const T &, const C &, const Value &value) :
      myValue(new Holder<ValueTraits<T>, C>(ValueTraits<T>::init())) {
      typename ValueTraits<T>::RepType v;
      if (value.get(v)) set(v);
    }


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    template<typename T>
    Value &operator=(const T &value) {
      set(value);
      return *this;
    }

    Value &operator=(const Value &rhs);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Conversion
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    //     template<typename T>
    //     operator T() const {
    //       T val = T();
    //       if (myValue != NULL) myValue->get(val);
    //       return val;
    //     }
    operator std::string() const;
    operator int() const;
    operator unsigned int() const;
    operator Real() const;
    operator bool() const;
    operator Vector3D() const;
    operator std::vector < Real>() const;


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Value
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool valid() const;
    bool isDefined() const;
    bool clear();
    bool init();
    unsigned int size() const;
    const std::type_info &type() const;
    std::string getString() const;
    /// Returns a string containing the type, value (only if valid) and
    /// constraint (if defined)
    std::string getDefinitionTypeString() const;
    /// Returns a complete debug string of the value, validity, type and
    /// constraint
    std::string debug() const;

    // Type
    ValueType::Enum getType() const;
    ConstraintValueType::Enum getConstraintType() const;
    const std::string &getTypeString() const;
    const std::string &getConstraintTypeString() const;
    /// Test on type
    bool equalType(const Value &v) const;
    /// Test on constraint type
    bool equalConstraint(const Value &v) const;
    /// Test on type and constraint type
    bool equalTypeAndConstraint(const Value &v) const;

    /// Sets the value, type and constraints are unchanged
    bool set(Value v);

    /// Sets the value, type and constraints are unchanged
    template<typename T>
    bool set(const T &value) {
      if (isDefined()) myValue->set(value);
      else {
        Value v(value);
        Value(v).swap(*this);
      }

      return valid();
    }


    /// Safe retrieve of the value, (trying) converting to the passed type
    template<typename T>
    bool get(T &value) const {
      return myValue != NULL ? myValue->get(value) : false;
    }

    /// Retrieve of the value by conversions
    template<typename T>
    T get() const {
      T tmp = T();
      get(tmp);
      return tmp;
    }

  private:
    Value &swap(Value &rhs);

    template<typename T>
    static bool equal(const Value &v1, const Value &v2);

    template<typename T>
    static bool less(const Value &v1, const Value &v2);

    template<typename T>
    static bool lessEqual(const Value &v1, const Value &v2);

    template<typename T>
    static bool greater(const Value &v1, const Value &v2);

    template<typename T>
    static bool greaterEqual(const Value &v1, const Value &v2);


    template<typename T>
    static bool equalType(const Value &v1, const Value &v2);

    void read(std::istream &is);

  public:
    friend Report::MyStreamer &operator<<(Report::MyStreamer &os,
                                          const Value &v);
    friend std::istream &operator>>(std::istream &is, Value &v) {
      v.read(is);
      return is;
    }

    // Boolean ==
    friend bool operator==(const Value &v1, const Value &v2);

    template<typename T>
    friend bool operator==(const Value &v1, const T &v2) {
      return v1 == Value(v2);
    }

    template<typename T>
    friend bool operator==(const T &v1, const Value &v2) {
      return Value(v1) ==
             v2;
    }

    // Boolean !=
    friend bool operator!=(const Value &v1, const Value &v2) {
      return !(v1 == v2);
    }

    template<typename T>
    friend bool operator!=(const Value &v1, const T &v2) {
      return !(v1 == Value(v2));
    }

    template<typename T>
    friend bool operator!=(const T &v1, const Value &v2) {
      return !(Value(v1) == v2);
    }

    //Boolean <
    friend bool operator<(const Value &v1, const Value &v2);

    template<typename T>
    friend bool operator<(const Value &v1, const T &v2) {
      return v1 < Value(v2);
    }

    template<typename T>
    friend bool operator<(const T &v1, const Value &v2) {
      return Value(v1) < v2;
    }

    //Boolean <=
    friend bool operator<=(const Value &v1, const Value &v2);

    template<typename T>
    friend bool operator<=(const Value &v1, const T &v2) {
      return v1 <= Value(v2);
    }

    template<typename T>
    friend bool operator<=(const T &v1, const Value &v2) {
      return Value(v1) <=
             v2;
    }

    //Boolean >
    friend bool operator>(const Value &v1, const Value &v2);

    template<typename T>
    friend bool operator>(const Value &v1, const T &v2) {
      return v1 > Value(v2);
    }

    template<typename T>
    friend bool operator>(const T &v1, const Value &v2) {
      return Value(v1) > v2;
    }

    //Boolean >=
    friend bool operator>=(const Value &v1, const Value &v2);

    template<typename T>
    friend bool operator>=(const Value &v1, const T &v2) {
      return v1 >= Value(v2);
    }

    template<typename T>
    friend bool operator>=(const T &v1, const Value &v2) {
      return Value(v1) >=
             v2;
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // private data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    PlaceHolder *myValue;
  };

  //_________________________________________________________________INLINES
  inline bool Value::valid() const {
    return myValue ? myValue->valid() : false;
  }

  inline ValueType::Enum Value::getType() const {
    return myValue ? myValue->getType() : ValueType::UNDEFINED;
  }

  inline bool Value::isDefined() const {
    return getType() != ValueType::UNDEFINED;
  }

  inline bool Value::clear() {
    return myValue ? myValue->clear() : false;
  }

  inline bool Value::init() {
    return myValue ? myValue->init() : false;
  }

  inline unsigned int Value::size() const {
    return valid() ? myValue->size() : 0;
  }

  inline const std::type_info &Value::type() const {
    return myValue ? myValue->type() : typeid(void);
  }

  inline ConstraintValueType::Enum Value::getConstraintType() const {
    return myValue ? myValue->getConstraintType() :
      ConstraintValueType::UNDEFINED;
  }

  inline const std::string &Value::getTypeString() const {
    return myValue ? myValue->getTypeString() :
      ValueType::getString(ValueType::UNDEFINED);
  }

  inline const std::string &Value::getConstraintTypeString() const {
    return myValue ? myValue->getConstraintTypeString() :
      ConstraintValueType::getString(ConstraintValueType::UNDEFINED);
  }

  inline bool Value::equalType(const Value &v) const {
    return v.getType() == getType();
  }

  inline bool Value::equalConstraint(const Value &v) const {
    return v.getConstraintType() == getConstraintType();
  }

  inline bool Value::equalTypeAndConstraint(const Value &v) const {
    return v.getType() == getType() &&
      v.getConstraintType() == getConstraintType();
  }
}
#endif
