#include <protomol/switch/C1SwitchingFunction.h>
#include <protomol/base/StringUtilities.h>

using namespace std;
using namespace ProtoMol;
//____ C1SwitchingFunction

C1SwitchingFunction::C1SwitchingFunction() :
  myCutoff(0.0), myCutoff2(0.0), myCutoff_3(0.0), my15Cutoff_1(0.0),
  my05Cutoff_3(0.0), my15Cutoff_3(0.0) {}

C1SwitchingFunction::C1SwitchingFunction(Real cutoff) :
  myCutoff(cutoff), myCutoff2(cutoff * cutoff),
  myCutoff_3(1.0 / (cutoff * cutoff * cutoff)), my15Cutoff_1(1.5 / cutoff),
  my05Cutoff_3(0.5 / (cutoff * cutoff * cutoff)),
  my15Cutoff_3(1.5 / (cutoff * cutoff * cutoff)) {}

void C1SwitchingFunction::getParameters(vector<Parameter> &parameters) const {
  parameters.push_back
    (Parameter("-cutoff", Value(myCutoff, ConstraintValueType::Positive()),
               Text("C1 swf cutoff")));
}

C1SwitchingFunction C1SwitchingFunction::make(vector<Value> values) {
  return C1SwitchingFunction(values[0]);
}

Matrix3By3 C1SwitchingFunction::hessian(const Vector3D &rij, Real a) const {
  Real rijNorm = sqrt(a);

  Real tm1 = 1.5 / (myCutoff * rijNorm);
  Real tm2 = 1.5 * myCutoff_3 * rijNorm;
  Real tm3 = -tm1 + tm2;
  Real tm4 = tm1 + tm2;

  return
    Matrix3By3(tm3, 0, 0, 0, tm3, 0, 0, 0, tm3) +
    Matrix3By3(rij, rij / a) * tm4;
}
