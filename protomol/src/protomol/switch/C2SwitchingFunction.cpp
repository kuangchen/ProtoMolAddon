#include <protomol/switch/C2SwitchingFunction.h>
#include <protomol/base/MathUtilities.h>
#include <protomol/base/Exception.h>

using namespace std;
using namespace ProtoMol;
//____ C2SwitchingFunction

C2SwitchingFunction::C2SwitchingFunction() :
  mySwitchon(-1.0), myCutoff(0.0) {}

C2SwitchingFunction::C2SwitchingFunction(Real switchon,
                                         Real cutoff) :
  mySwitchon(switchon), mySwitchon2(switchon * switchon), myCutoff(cutoff),
  myCutoff2(cutoff * cutoff),
  mySwitch1(1.0 / power<3>(cutoff * cutoff - switchon * switchon)),
  mySwitch2(cutoff * cutoff - 3.0 * switchon * switchon),
  mySwitch3(4.0 / power<3>(cutoff * cutoff - switchon * switchon)) {}

void C2SwitchingFunction::getParameters(vector<Parameter> &parameters) const {
  parameters.push_back
    (Parameter("-switchon",
               Value(mySwitchon, ConstraintValueType::NotNegative()),
               Text("C2 swf switch on")));
  parameters.push_back
    (Parameter("-cutoff", Value(myCutoff, ConstraintValueType::Positive()),
               Text("C2 swf cutoff")));
}

C2SwitchingFunction C2SwitchingFunction::make(vector<Value> values) {
  Real switchon, cutoff;
  values[0].get(switchon);
  values[1].get(cutoff);
  if (!values[0].valid() || !values[0].valid() || switchon < 0.0 || cutoff <=
      0.0 || switchon >= cutoff)
    THROW(getId() + " switching function: 0 <= switchon (=" +
          values[0].getString() + ") < cutoff (=" +
          values[1].getString() + ").");

  return C2SwitchingFunction(switchon, cutoff);
}

Matrix3By3 C2SwitchingFunction::hessian(const Vector3D &rij, Real a) const {
  Real tm1 = myCutoff2 - mySwitchon2;
  Real tm2 = tm1 * tm1 * tm1;
  Real tm3 = 12.0 * (a - myCutoff2) * (a - mySwitchon2) / tm2;
  Real tm4 = 24.0 * (2.0 * a - myCutoff2 - mySwitchon2) / tm2;

  return
    Matrix3By3(tm3, 0, 0, 0, tm3, 0, 0, 0, tm3) +
    Matrix3By3(rij, rij) * tm4;
}
