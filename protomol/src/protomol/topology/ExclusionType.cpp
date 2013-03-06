#include <protomol/topology/ExclusionType.h>

using namespace std;
using namespace ProtoMol;
//____ExclusionType

const string ExclusionEnum::str[static_cast<int>(LAST) -
                                static_cast<int>(FIRST)] = {
  // Order is essential, must be in relation to Enum ordering
  string("undefined"),  // Returned when no enum matches
  string("none"),
  string("1-2"),
  string("1-3"),
  string("1-4"),
  string("scaled1-4")
};
