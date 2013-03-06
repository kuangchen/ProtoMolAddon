/* -*- c++ -*- */
#ifndef CnSWITCHINGFUNCTION_H
#define CnSWITCHINGFUNCTION_H

#include <vector>
#include <string>

#include <protomol/config/Parameter.h>
#include <protomol/type/Matrix3By3.h>

namespace ProtoMol {
  //____ CnSwitchingFunction

  /**
   * The switching function provide C2, 3, 4 0r 6 continuous
   */
  class CnSwitchingFunction {

  public:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types and Enums
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    enum {USE = 1};
    enum {MODIFY = 1};
    enum {CUTOFF = 1};
    enum {MAXEQNN = 7};
    enum {NUMSW = 5};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    CnSwitchingFunction();
    CnSwitchingFunction(Real order, Real switchon, Real switchoff, Real cutoff );
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CnSwitchingFunction
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool roughTest(Real distSquared) const {return distSquared <= myCutoff2;}

    Real cutoffSquared() const {return myCutoff2;}

    Real cutoff() const {return myCutoff;}

    void operator()(Real &value, Real &deriv, Real distSquared) const;
    Matrix3By3 hessian(const Vector3D &rij, Real distSquared) const;

    static const std::string getId() {return "Cn";}
    void getParameters(std::vector<Parameter> &parameters) const;
    static unsigned int getParameterSize() {return 4;}
    static CnSwitchingFunction make(std::vector<Value> values);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real mySwitchon, mySwitchon2, myCutoff, myCutoff2, myOrder, mySwitchoff;
    Real myIRange[MAXEQNN];
    int ordInt, ordIdx;

    static Real swcoef[][MAXEQNN], dswcoef[][MAXEQNN], d2swcoef[][MAXEQNN];
  };
  //____ INLINES
}
#endif /* CnSWITCHINGFUNCTION_H */
