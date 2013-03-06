#include <protomol/config/ValueType.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ValueType
const string ValueEnum::str[static_cast<int>(LAST) -
                            static_cast<int>(FIRST)] = {
  // Order is essential, must be in relation to Enum ordering
  string("undefined"),  // Returned when no enum matches
  string("string"),
  string("int"),
  string("uint"),
  string("real"),
  string("boolean"),
  string("coordinates"),
  string("vector"),
  string("integrator"),
  string("force")
};
