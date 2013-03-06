#include <protomol/force/born/BornRadii.h>

using namespace std;
using namespace ProtoMol;
//____ Born Radii calculation only

const string BornRadii::keyword("BornRadii");

void BornRadii::getParameters(vector<Parameter> &parameters) const {
  parameters.push_back
    (Parameter("-bornswitch", Value(myBornSwitch, ConstraintValueType::NotNegative()), 3,
               Text("Born Switch type")));
}

BornRadii BornRadii::make(const vector<Value> &values) {
  return BornRadii(values[0]);
}
