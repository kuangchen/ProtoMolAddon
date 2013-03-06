#include <protomol/force/coulomb/CoulombForceDiElec.h>

using namespace std;
using namespace ProtoMol;

//____ CoulombForce with DiElectric Term for Implicit Solvation

const string CoulombForceDiElec::keyword("CoulombDiElec");
const string CoulombForceDiElec::C1::keyword("C1");
const string CoulombForceDiElec::C2::keyword("C2");
const string CoulombForceDiElec::C3::keyword("C3");
const string CoulombForceDiElec::C4::keyword("C4");

// Default constructor
CoulombForceDiElec::CoulombForceDiElec() : EPS(1), D(53), S(0.25) {}

// Constructor with parameters
CoulombForceDiElec::CoulombForceDiElec(Real EPSval, Real Dval, Real Sval) :
  EPS(EPSval), D(Dval), S(Sval) {}

void CoulombForceDiElec::getParameters(vector<Parameter> &parameters)
const {
  parameters.push_back
    (Parameter("-EPS", Value(EPS, ConstraintValueType::NotNegative()), 1,
               Text("DiElect EPSILON")));
  parameters.push_back
    (Parameter("-D", Value(D, ConstraintValueType::NotNegative()),
               78, Text("DiElect Dval")));
  parameters.push_back
    (Parameter("-S", Value(S, ConstraintValueType::NotNegative()),
               0.3, Text("DiElect Sval")));
}

CoulombForceDiElec CoulombForceDiElec::make(const vector<Value> &values) {
  return CoulombForceDiElec(values[0], values[1], values[2]);
}
