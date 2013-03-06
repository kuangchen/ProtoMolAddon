#include <protomol/force/GB/GBForce.h>

using namespace std;
using namespace ProtoMol;

const string GBForce::keyword("GBForce");

void GBForce::getParameters(vector<Parameter> &parameters) const {

   parameters.push_back
     (Parameter("-soluteDielec", Value(soluteDielec, ConstraintValueType::NoConstraints()), 1.0, Text("Solute Dielectric")));
   parameters.push_back
     (Parameter("-solventDielec", Value(solventDielec, ConstraintValueType::NoConstraints()), 80.0, Text("Solvent Dielectric")));

}

GBForce GBForce::make(const vector<Value> &values) {
   return (GBForce(values[0],values[1]));
}
