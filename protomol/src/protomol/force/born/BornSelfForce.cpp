#include <protomol/force/born/BornSelfForce.h>
using namespace std;
using namespace ProtoMol;
//____ BornSelfForce

const string BornSelfForce::keyword("BornSelf");

void BornSelfForce::getParameters(vector<Parameter> &parameters) const {
  parameters.push_back
    (Parameter("-bornswitch", Value(myBornSwitch, ConstraintValueType::NotNegative()), 3,
               Text("Born Switch type")));
  parameters.push_back
    (Parameter("-D", Value(myDielecConst, ConstraintValueType::NotNegative()), 80,
               Text("Bulk solvent dielectric constant.")));
}

BornSelfForce BornSelfForce::make(const vector<Value> &values) {
  return BornSelfForce(values[0],values[1]);
}
