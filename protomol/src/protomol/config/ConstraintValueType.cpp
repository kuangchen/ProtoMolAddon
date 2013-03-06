#include <protomol/config/ConstraintValueType.h>

using namespace std;
using namespace ProtoMol;

//____ConstraintValueType
const string ConstraintValueEnum::str[static_cast<int>(LAST) -
                                      static_cast<int>(FIRST)] = {
  // Order is essential, must be in relation to Enum ordering
  string("undefined"),  // Returned when no enum matches
  string("no-constraints"),
  string("empty"),
  string("non-empty"),
  string("zero"),
  string("non-zero"),
  string("positive"),
  string("negative"),
  string("non-positive"),
  string("non-negative")
};
