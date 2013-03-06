#include <protomol/config/Value.h>

#include <protomol/base/StringUtilities.h>
#include <protomol/base/MathUtilities.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;
//____ Value

const Value::Undefined *Value::undefined(NULL);

Value::Value() :
  myValue(NULL) {}

Value::~Value() {
  if (myValue != NULL) {
    delete myValue;
    myValue = NULL;
  }
}

Value::Value(const string &value) :
  myValue(new Holder<ValueTraits<string> >(value)) {}

Value::Value(const char *value) :
  myValue(new Holder<ValueTraits<string> > (value)) {}

Value::Value(int value) :
  myValue(new Holder<ValueTraits<int> >(value)) {}

Value::Value(unsigned int value) :
  myValue(new Holder<ValueTraits<unsigned int> >(value)) {}

Value::Value(Real value) :
  myValue(new Holder<ValueTraits<Real> >(value)) {}

Value::Value(bool value) :
  myValue(new Holder<ValueTraits<bool> >(value)) {}

Value::Value(const Vector3D &value) :
  myValue(new Holder<ValueTraits<Vector3D> > (value)) {}

Value::Value(const vector<Real> &value) :
  myValue(new Holder<ValueTraits<vector<Real> > >(value)) {}

Value::Value(const ValueType::Integrator &value) :
  myValue(new Holder<ValueTraits<ValueType::Integrator> >(value)) {}

Value::Value(const ValueType::Force &value) :
  myValue(new Holder<ValueTraits<ValueType::Force> >(value)) {}

Value::Value(const Value &other) :
  myValue(other.myValue ? other.myValue->clone() :
            NULL) {}

Value::Value(const string &value, const Undefined *) :
  myValue(new Holder<ValueTraits<string> >(value, false)) {}

Value::Value(const char *value, const Undefined *) :
  myValue(new Holder<ValueTraits<string> >(value, false)) {}

Value::Value(int value, const Undefined *) :
  myValue(new Holder<ValueTraits<int> >(value, false)) {}

Value::Value(unsigned int value, const Undefined *) :
  myValue(new Holder<ValueTraits<unsigned int> >(value, false)) {}

Value::Value(Real value, const Undefined *) :
  myValue(new Holder<ValueTraits<Real> >(value, false)) {}

Value::Value(bool value, const Undefined *) :
  myValue(new Holder<ValueTraits<bool> >(value, false)) {}

Value::Value(const Vector3D &value, const Undefined *) :
  myValue(new Holder<ValueTraits<Vector3D> >(value, false)) {}

Value::Value(const vector<Real> &value, const Undefined *) :
  myValue(new Holder<ValueTraits<vector<Real> > >(value, false)) {}

Value::Value(const ValueType::Integrator &value, const Undefined *) :
  myValue(new Holder<ValueTraits<ValueType::Integrator> >(value, false)) {}

Value::Value(const ValueType::Force &value, const Undefined *) :
  myValue(new Holder<ValueTraits<ValueType::Force> >(value, false)) {}

Value::Value(const Value &value, const Undefined *) :
  myValue(value.myValue ? value.myValue->clone() :
            NULL) {
  clear();
}

Value &Value::operator=(const Value &rhs) {
  Value(rhs).swap(*this);
  return *this;
}

Value::operator string() const {
  string val;
  if (myValue != NULL) myValue->get(val);
  return val;
}

Value::operator int() const {
  int val = 0;
  if (myValue != NULL) myValue->get(val);
  return val;
}

Value::operator unsigned int() const {
  unsigned int val = 0;
  if (myValue != NULL) myValue->get(val);
  return val;
}

Value::operator Real() const {
  Real val = 0.0;
  if (myValue != NULL) myValue->get(val);
  return val;
}

Value::operator bool() const {
  bool val = false;
  if (myValue != NULL) myValue->get(val);
  return val;
}

Value::operator Vector3D() const {
  Vector3D val(0.0, 0.0, 0.0);
  if (myValue != NULL) myValue->get(val);
  return val;
}

Value::operator vector < Real> () const {
  vector<Real> val;
  if (myValue != NULL) myValue->get(val);
  return val;
}

string Value::getString() const {
  string tmp("");
  if (myValue) myValue->get(tmp);
  return tmp;
}

string Value::getDefinitionTypeString() const {
  string res("<");
  res += getTypeString();
  if (valid()) res += "=" + getString();

  if (getConstraintType() != ConstraintValueType::NOCONSTRAINTS &&
      getConstraintType() != ConstraintValueType::UNDEFINED)
    res += "," + getConstraintTypeString();
  res += ">";

  return res;
}

bool Value::set(Value v) {
  if (isDefined()) myValue->set(v);
  else Value(v).swap(*this);

  return valid();
}

Value &Value::swap(Value &rhs) {
  std::swap(myValue, rhs.myValue);
  return *this;
}

string Value::debug() const {
  return getString() + "," + (valid() ? "true" : "false") + "," +
         getTypeString() + "," + getConstraintTypeString();
}

void Value::read(istream &is) {
  if (isDefined()) myValue->read(is);
  else {
    string str;
    is >> str;
    set(str);
  }
}

template<typename T>
inline bool Value::equal(const Value &v1, const Value &v2) {
  T u1, u2;
  if (v1.get(u1) && v2.get(u2)) return u1 == u2;
  return false;
}

template<typename T>
inline bool Value::less(const Value &v1, const Value &v2) {
  T u1, u2;
  if (v1.get(u1) && v2.get(u2)) return u1 < u2;
  return false;
}

template<typename T>
inline bool Value::lessEqual(const Value &v1, const Value &v2) {
  T u1, u2;
  if (v1.get(u1) && v2.get(u2)) return u1 <= u2;
  return false;
}

template<typename T>
inline bool Value::greater(const Value &v1, const Value &v2) {
  T u1, u2;
  if (v1.get(u1) && v2.get(u2)) return u1 < u2;
  return false;
}

template<typename T>
inline bool Value::greaterEqual(const Value &v1, const Value &v2) {
  T u1, u2;
  if (v1.get(u1) && v2.get(u2)) return u1 <= u2;
  return false;
}

template<typename T>
inline bool Value::equalType(const Value &v1, const Value &v2) {
  T u1, u2;
  return v1.get(u1) && v2.get(u2) &&
         (v1.getType() == (ValueType::Enum)ValueTraits<T>::value ||
          v2.getType() ==
          (ValueType::Enum)ValueTraits<T>::value);
}

namespace ProtoMol {
  template<>
  inline bool Value::greaterEqual<Vector3D>(const Value &, const Value &) {
    return false;
  }

  template<>
  inline bool Value::greater<Vector3D>(const Value &, const Value &) {
    return false;
  }

  template<>
  inline bool Value::less<Vector3D>(const Value &, const Value &) {
    return false;
  }

  template<>
  inline bool Value::lessEqual<Vector3D>(const Value &, const Value &) {
    return false;
  }

  template<>
  inline bool Value::greaterEqual<vector<Real> >(const Value &,
                                                 const Value &) {
    return false;
  }

  template<>
  inline bool Value::greater<vector<Real> >(const Value &, const Value &) {
    return false;
  }

  template<>
  inline bool Value::less<vector<Real> >(const Value &, const Value &) {
    return false;
  }

  template<>
  inline bool Value::lessEqual<vector<Real> >(const Value &, const Value &) {
    return false;
  }

  bool operator==(const Value &v1, const Value &v2) {
    if (!v1.isDefined() && !v2.isDefined()) return true;
    if (v1.valid() != v2.valid()) return false;
    else if (!v1.valid()) return true;

    if (Value::equalType<Vector3D>(v1, v2))
      return Value::equal<Vector3D>(v1, v2);
    else if (Value::equalType<Real>(v1, v2))
      return Value::equal<Real>(v1, v2);
    else if (Value::equalType<int>(v1, v2))
      return Value::equal<int>(v1, v2);
    else if (Value::equalType<unsigned int>(v1, v2))
      return Value::equal<unsigned int>(v1, v2);
    else if (Value::equalType<bool>(v1, v2))
      return Value::equal<bool>(v1, v2);

    return v1.getString() == v2.getString();
  }

  bool operator<(const Value &v1, const Value &v2) {
    if (!v1.isDefined() && !v2.isDefined()) return true;
    if (v1.valid() != v2.valid()) return false;
    else if (!v1.valid()) return true;

    if (Value::equalType<Vector3D>(v1, v2))
      return Value::less<Vector3D>(v1, v2);
    else if (Value::equalType<Real>(v1, v2))
      return Value::less<Real>(v1, v2);
    else if (Value::equalType<int>(v1, v2))
      return Value::less<int>(v1, v2);
    else if (Value::equalType<unsigned int>(v1, v2))
      return Value::less<unsigned int>(v1, v2);
    else if (Value::equalType<bool>(v1, v2))
      return Value::less<bool>(v1, v2);

    return v1.getString() < v2.getString();
  }

  bool operator<=(const Value &v1, const Value &v2) {
    if (!v1.isDefined() && !v2.isDefined()) return true;
    if (v1.valid() != v2.valid()) return false;
    else if (!v1.valid()) return true;

    if (Value::equalType<Vector3D>(v1, v2))
      return Value::lessEqual<Vector3D>(v1, v2);
    else if (Value::equalType<Real>(v1, v2))
      return Value::lessEqual<Real>(v1, v2);
    else if (Value::equalType<int>(v1, v2))
      return Value::lessEqual<int>(v1, v2);
    else if (Value::equalType<unsigned int>(v1, v2))
      return Value::lessEqual<unsigned int>(v1, v2);
    else if (Value::equalType<bool>(v1, v2))
      return Value::lessEqual<bool>(v1, v2);

    return v1.getString() <= v2.getString();
  }

  bool operator>(const Value &v1, const Value &v2) {
    if (!v1.isDefined() && !v2.isDefined()) return true;
    if (v1.valid() != v2.valid()) return false;
    else if (!v1.valid()) return true;

    if (Value::equalType<Vector3D>(v1, v2))
      return Value::greater<Vector3D>(v1, v2);
    else if (Value::equalType<Real>(v1, v2))
      return Value::greater<Real>(v1, v2);
    else if (Value::equalType<int>(v1, v2))
      return Value::greater<int>(v1, v2);
    else if (Value::equalType<unsigned int>(v1, v2))
      return Value::greater<unsigned int>(v1, v2);
    else if (Value::equalType<bool>(v1, v2))
      return Value::greater<bool>(v1, v2);

    return v1.getString() > v2.getString();
  }

  bool operator>=(const Value &v1, const Value &v2) {
    if (!v1.isDefined() && !v2.isDefined()) return true;
    if (v1.valid() != v2.valid()) return false;
    else if (!v1.valid()) return true;

    if (Value::equalType<Vector3D>(v1, v2))
      return Value::greaterEqual<Vector3D>(v1, v2);
    else if (Value::equalType<Real>(v1, v2))
      return Value::greaterEqual<Real>(v1, v2);
    else if (Value::equalType<int>(v1, v2))
      return Value::greaterEqual<int>(v1, v2);
    else if (Value::equalType<unsigned int>(v1, v2))
      return Value::greaterEqual<unsigned int>(v1, v2);
    else if (Value::equalType<bool>(v1, v2))
      return Value::greaterEqual<bool>(v1, v2);

    return v1.getString() >= v2.getString();
  }

  MyStreamer &operator<<(MyStreamer &OS, const Value &v) {
    OS << v.getString();
    return OS;
  }
}
