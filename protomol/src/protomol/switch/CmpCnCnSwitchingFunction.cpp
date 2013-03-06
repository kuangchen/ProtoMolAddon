#include <protomol/switch/CmpCnCnSwitchingFunction.h>
#include <protomol/base/MathUtilities.h>
#include <protomol/base/Exception.h>

using namespace std;
using namespace ProtoMol;
//____ CmpCnCnSwitchingFunction

CmpCnCnSwitchingFunction::CmpCnCnSwitchingFunction() :
  mySwitchon(-1.0), myCutoff(0.0), myCompoff(0.0), myOrder(0.0) {}

CmpCnCnSwitchingFunction::CmpCnCnSwitchingFunction(Real switchon, Real cutoff,
                                                   Real switchoff,
                                                   Real order) :
  mySwitchon(switchon), mySwitchon2(switchon * switchon), myCutoff(cutoff),
  myCutoff2(cutoff * cutoff), myCompoff(switchoff),
  myCompoff2(switchoff * switchoff), myOrder(order) {
  if (myOrder == 5.0) myOrder = 2.0;   // fix as 5 dosn't exist
  ordInt = (int)myOrder;
  ordIdx = ordInt - 2;   // /2-1;
  myIRange[0] = 1.0 / pow(switchon - switchoff, ordInt + 1);
  myItRange[0] = 1.0 / pow(switchoff - cutoff, ordInt + 1);

  for (int i = 1; i < ordInt + 1; i++) {
    myIRange[i] = myIRange[i - 1] / (switchon - switchoff);
    myItRange[i] = myItRange[i - 1] / (switchoff - cutoff);
  }
}

void CmpCnCnSwitchingFunction::getParameters(vector<Parameter> &parameters)
const {
  parameters.push_back
    (Parameter("-switchon",
               Value(mySwitchon, ConstraintValueType::NotNegative()),
               Text("CmpCnCn swf switch on")));
  parameters.push_back
    (Parameter("-switchoff", Value(myCompoff, ConstraintValueType::Positive()),
               Text("CmpCnCn swf switchoff")));
  parameters.push_back
    (Parameter("-cutoff", Value(myCutoff, ConstraintValueType::Positive()),
               Text("CmpCnCn swf cutoff")));
  parameters.push_back
    (Parameter("-n", Value(myOrder, ConstraintValueType::Positive()), 2.0,
               Text("CmpCnCn swf smoothness")));
}

CmpCnCnSwitchingFunction CmpCnCnSwitchingFunction::make(vector<Value> values) {
  Real switchon, cutoff, switchoff, order;
  values[0].get(switchon);
  values[2].get(cutoff);
  values[1].get(switchoff);
  values[3].get(order);
  if (!values[0].valid() || !values[1].valid() || !values[2].valid() ||
      !values[3].valid() || switchon < 0.0 || cutoff <= 0.0 ||
      switchoff <= 0.0 || switchon >= switchoff || switchoff > cutoff ||
      order < 2.0 || order > 6.0)
    THROW(getId() + " switching function: 0 <= switchon (=" +
          values[0].getString() + ") < cutoff (=" + values[1].getString() +
          ") < switchoff (=" + values[2].getString() + ")." +
          ", order 2,3,4,6 (=" + values[3].getString() + ").");

  return CmpCnCnSwitchingFunction(switchon, cutoff, switchoff, order);
}

Matrix3By3 CmpCnCnSwitchingFunction::hessian(const Vector3D &rij,
                                             Real a) const {
  Real dswcoef[NUMSWCM][MAXEQNCM] = {
    {30.0, -60.0, 30.0, 0, 0, 0, 0},
    {140.0, -420.0, 420.0, -140.0, 0, 0, 0},
    {630.0, -2520.0, 3780.0, -2520.0, 630.0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0},
    {12012.0, -72072.0, 180180.0, -240240.0, 180180.0, -72072.0, 12012.0}
  };

  Real d2swcoef[NUMSWCM][MAXEQNCM] = {
    {60.0, -180.0, 120.0, 0, 0, 0, 0},
    {420.0, -1680.0, 2100.0, -840.0, 0, 0, 0},
    {2520.0, -12600.0, 22680.0, -17640.0, 5040.0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0},
    {72072.0, -504504.0, 1441440.0, -2162160.0, 1801800.0, -792792.0, 144144.0}
  };

  Real sqrta = sqrt(a);
  Real c[MAXEQNCM + 1];
  Real swDiff;
  if (a > myCompoff2) {
    swDiff = sqrta - myCutoff;
    c[0] = pow(swDiff, ordInt - 1);
    for (int i = 1; i < ordInt + 2; i++) c[i] = c[i - 1] * swDiff;

    Real tm3, tm4;
    tm3 = tm4 = 0.0;
    for (int i = 0; i < ordInt + 1; i++) {
      tm3 += dswcoef[ordIdx][i] * c[i + 1] * myItRange[i];
      tm4 += d2swcoef[ordIdx][i] * c[i] * myItRange[i];
    }

    tm3 /= sqrta;
    tm4 /= a;
    tm4 -= tm3 / a;
    return Matrix3By3(tm3, 0, 0, 0, tm3, 0, 0, 0,
                      tm3) + Matrix3By3(rij, rij) * tm4;
  } else {
    swDiff = sqrta - myCompoff;
    c[0] = pow(swDiff, ordInt - 1);
    for (int i = 1; i < ordInt + 2; i++) c[i] = c[i - 1] * swDiff;

    Real tm3, tm4;
    tm3 = tm4 = 0.0;
    for (int i = 0; i < ordInt + 1; i++) {
      tm3 += dswcoef[ordIdx][i] * c[i + 1] * myIRange[i];
      tm4 += d2swcoef[ordIdx][i] * c[i] * myIRange[i];
    }

    tm3 /= sqrta;
    tm4 /= a;
    tm4 -= tm3 / a;
    return 
      -((Matrix3By3(tm3, 0, 0, 0, tm3, 0, 0, 0, tm3) +
         Matrix3By3(rij, rij) * tm4));
  }
}
