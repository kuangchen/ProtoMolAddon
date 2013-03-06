#include <protomol/config/Parameter.h>
#include <protomol/base/StringUtilities.h>

#include <iomanip>

using namespace std;
using namespace ProtoMol;

//____ Parameter
Parameter::Parameter() :
  keyword(""), value(Value()), defaultValue(Value()) {}

Parameter::Parameter(const string &k, const Value &val) :
  keyword(k), value(val), defaultValue(val, Value::undefined) {}

Parameter::Parameter(const string &k, const Value &val, const Value &def) :
  keyword(k), value(val), defaultValue(def) {}

Parameter::Parameter(const char *k, const Value &val) :
  keyword(string(k)), value(val), defaultValue(val, Value::undefined) {}

Parameter::Parameter(const char *k, const Value &val, const Value &def) :
  keyword(string(k)), value(val), defaultValue(def) {}

Parameter::Parameter(const string &k, const Value &val, const Text &t) :
  keyword(k), value(val), defaultValue(val, Value::undefined), text(t.text) {}

Parameter::Parameter(const string &k, const Value &val, const Value &def,
                     const Text &t) :
  keyword(k), value(val), defaultValue(def), text(t.text) {}

Parameter::Parameter(const char *k, const Value &val, const Text &t) :
  keyword(string(k)), value(val),
  defaultValue(val, Value::undefined), text(t.text) {}

Parameter::Parameter(const char *k, const Value &val, const Value &def,
                     const Text &t) :
  keyword(string(k)), value(val), defaultValue(def), text(t.text) {}

ostream &Parameter::print(ostream &stream) const {
  stream.setf(ios::left);

  if (!keyword.empty())
    stream << Constant::PRINTINDENT << Constant::PRINTINDENT
           << setw(25) << keyword << " ";

  stream << setw(30);
  if (defaultValue.valid())
    stream << defaultValue.getDefinitionTypeString();
  else stream << value.getDefinitionTypeString();

  if (!text.empty()) {
    stream << " # ";
    int col = 2 * string(Constant::PRINTINDENT).length() + 59;
    fillFormat(stream, text, col, col, 120);

  } else stream << endl;

  return stream;
}
