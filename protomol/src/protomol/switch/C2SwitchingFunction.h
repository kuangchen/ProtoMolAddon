/* -*- c++ -*- */
#ifndef C2SWITCHINGFUNCTION_H
#define C2SWITCHINGFUNCTION_H

#include <vector>
#include <string>

#include <protomol/config/Parameter.h>
#include <protomol/type/Matrix3By3.h>

namespace ProtoMol {
  //____ C2SwitchingFunction

  /**
   * The switching function provide C2 continues
   */
  class C2SwitchingFunction {
  public:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types and Enums
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    enum {USE = 1};
    enum {MODIFY = 1};
    enum {CUTOFF = 1};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    C2SwitchingFunction();
    C2SwitchingFunction(Real switchon, Real cutoff);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class C2SwitchingFunction
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool roughTest(Real distSquared) const {return distSquared <= myCutoff2;}

    Real cutoffSquared() const {return myCutoff2;}

    Real cutoff() const {return myCutoff;}

    void operator()(Real &value, Real &deriv, Real distSquared) const {
      deriv = 0.0;
      if (distSquared > myCutoff2)
        value = 0.0;
      else if (distSquared >= mySwitchon2) {
        Real c2 = myCutoff2 - distSquared;
        Real c4 = c2 * (mySwitch2 + 2.0 * distSquared);
        value = mySwitch1 * (c2 * c4);
        deriv = mySwitch3 * (c2 * c2 - c4);
      } else
        value = 1.0;
    }

    Matrix3By3 hessian(const Vector3D &rij, Real distSquared) const;

    static const std::string getId() {return "C2";}
    void getParameters(std::vector<Parameter> &parameters) const;
    static unsigned int getParameterSize() {return 2;}
    static C2SwitchingFunction make(std::vector<Value> values);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real mySwitchon, mySwitchon2, myCutoff, myCutoff2, mySwitch1, mySwitch2,
         mySwitch3;
  };
  //____ INLINES
}
#endif /* C2SWITCHINGFUNCTION_H */
