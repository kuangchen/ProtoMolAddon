/* -*- c++ -*- */
#ifndef CMPCNCNSWITCHINGFUNCTION_H
#define CMPCNCNSWITCHINGFUNCTION_H

#include <vector>
#include <string>

#include <protomol/config/Parameter.h>
#include <protomol/type/Matrix3By3.h>

namespace ProtoMol {
  //____ CmpCnCnSwitchingFunction

  /**
   * The switching function provide C2,3,4 or 6 with compliment continuous
   */
  class CmpCnCnSwitchingFunction {

  public:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types and Enums
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    enum {USE = 1};
    enum {MODIFY = 1};
    enum {CUTOFF = 1};
    enum {MAXEQNCM = 7};
    enum {NUMSWCM = 5};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    CmpCnCnSwitchingFunction();
    CmpCnCnSwitchingFunction(Real switchon, Real cutoff, Real switchoff,
                             Real order);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CmpCnCnSwitchingFunction
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool roughTest(Real distSquared) const {return distSquared <= myCutoff2;}

    Real cutoffSquared() const {return myCutoff2;}

    Real cutoff() const {return myCutoff;}

    void operator()(Real &value, Real &deriv, Real distSquared) const {
      Real swcoef[NUMSWCM][MAXEQNCM] = {
        {10.0, -15.0, 6.0, 0, 0, 0, 0},
        {35.0, -84.0, 70.0, -20.0, 0, 0, 0},
        {126.0, -420.0, 540.0, -315.0, 70.0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0},
        {1716.0, -9009.0, 20020.0, -24024.0, 16380.0, -6006.0, 924.0}
      };

      Real dswcoef[NUMSWCM][MAXEQNCM] = {
        {30.0, -60.0, 30.0, 0, 0, 0, 0},
        {140.0, -420.0, 420.0, -140.0, 0, 0, 0},
        {630.0, -2520.0, 3780.0, -2520.0, 630.0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0},
        {12012.0, -72072.0, 180180.0, -240240.0, 180180.0, -72072.0, 12012.0}
      };
      Real c[MAXEQNCM + 1];
      Real swDiff;
      deriv = 0.0;
      if (distSquared > myCutoff2) value = 0.0;
      else if (distSquared > myCompoff2) {
        Real sqrtd = sqrt(distSquared);
        swDiff = sqrtd - myCutoff;
        c[0] = pow(swDiff, ordInt);
        for (int i = 1; i < ordInt + 2; i++) c[i] = c[i - 1] * swDiff;

        value = deriv = 0.0;
        for (int i = 0; i < ordInt + 1; i++) {
          value += swcoef[ordIdx][i] * c[i + 1] * myItRange[i];
          deriv += dswcoef[ordIdx][i] * c[i] * myItRange[i];
        }

        deriv /= sqrtd;
      } else if (distSquared >= mySwitchon2) {
        Real sqrtd = sqrt(distSquared);
        swDiff = sqrtd - myCompoff;
        c[0] = pow(swDiff, ordInt);
        for (int i = 1; i < ordInt + 2; i++) c[i] = c[i - 1] * swDiff;

        value = deriv = 0.0;
        for (int i = 0; i < ordInt + 1; i++) {
          value += swcoef[ordIdx][i] * c[i + 1] * myIRange[i];
          deriv -= dswcoef[ordIdx][i] * c[i] * myIRange[i];
        }

        value = 1.0 - value;
        deriv /= sqrtd;

      } else value = 0.0;
    }

    Matrix3By3 hessian(const Vector3D &rij, Real distSquared) const;

    static const std::string getId() {return "ComplimentCnCn";}
    void getParameters(std::vector<Parameter> &parameters) const;
    static unsigned int getParameterSize() {return 4;}
    static CmpCnCnSwitchingFunction make(std::vector<Value> values);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real mySwitchon, mySwitchon2, myCutoff, myCutoff2, myCompoff, myCompoff2,
         myOrder;
    Real myIRange[MAXEQNCM];
    Real myItRange[MAXEQNCM];
    int ordInt, ordIdx;
  };
  //____ INLINES
}
#endif /* CMPCNCNSWITCHINGFUNCTION_H */
