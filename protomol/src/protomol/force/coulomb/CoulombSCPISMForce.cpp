#include <protomol/force/coulomb/CoulombSCPISMForce.h>

using namespace std;
using namespace ProtoMol;

//____ CoulombForce with DiElectric Term for Implicit Solvation
const string CoulombSCPISMForce::keyword("CoulombSCPISM");
const string CoulombSCPISMForce::C1::keyword("C1");
const string CoulombSCPISMForce::C2::keyword("C2");
const string CoulombSCPISMForce::C3::keyword("C3");
const string CoulombSCPISMForce::C4::keyword("C4");

// Default constructor
CoulombSCPISMForce::CoulombSCPISMForce() : D(80) {}

// Constructor with parameters
CoulombSCPISMForce::CoulombSCPISMForce(Real Dval) : D(Dval) {}

void CoulombSCPISMForce::getParameters(vector<Parameter> &parameters) const {
  parameters.push_back
    (Parameter("-D", Value(D, ConstraintValueType::NotNegative()),
               80, Text("Bulk solvent dielectric")));
}

CoulombSCPISMForce CoulombSCPISMForce::make(const vector<Value> &values) {
  return CoulombSCPISMForce(values[0]);
}
