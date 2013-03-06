#include <protomol/force/GB/GBACEForce.h>

using namespace std;
using namespace ProtoMol;

const string GBACEForce::keyword("GBACEForce");

void GBACEForce::getParameters(vector<Parameter> &parameters) const {

     parameters.push_back(
        Parameter("-solvationparam",Value(sigma, ConstraintValueType::NoConstraints()),2.26 / 418.4, Text("solvation parameter")));
        //kJ nm^{-2} -> KCal \AA^{-2} 1/418.4

     parameters.push_back(
        Parameter("-watersphereradius",Value(rho_s, ConstraintValueType::NoConstraints()), 1.4, Text("solvation parameter")));

}

GBACEForce GBACEForce::make(const vector<Value> &values)
{
   return (GBACEForce(values[0], values[1]));
}
