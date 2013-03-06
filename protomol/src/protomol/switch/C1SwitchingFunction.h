/* -*- c++ -*- */
#ifndef C1SWITCHINGFUNCTION_H
#define C1SWITCHINGFUNCTION_H

#include <vector>

#include <protomol/config/Parameter.h>
#include <protomol/type/Matrix3By3.h>

namespace ProtoMol {
  //____ C1SwitchingFunction

  /**
   * The switching function provide C1 continues
   */
  class C1SwitchingFunction {
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
    C1SwitchingFunction();
    C1SwitchingFunction(Real cutoff);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class C1SwitchingFunction
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool roughTest(Real distSquared) const {return distSquared <= myCutoff2;}

    Real cutoffSquared() const {return myCutoff2;}

    Real cutoff() const {return myCutoff;}

    void operator()(Real &value, Real &derivOverD, Real distSquared) const {
      if (distSquared > myCutoff2) {
        value = 0.0;
        derivOverD = 0.0;
        return;
      }
      Real dist = sqrt(distSquared);
      value = 1.0 - dist * (my15Cutoff_1 - distSquared * my05Cutoff_3);
      derivOverD = -my15Cutoff_1 / dist + dist * my15Cutoff_3;
    }

    Matrix3By3 hessian(const Vector3D &rij, Real distSquared) const;

    static const std::string getId() {return "C1";}
    void getParameters(std::vector<Parameter> &parameters) const;
    static unsigned int getParameterSize() {return 1;}
    static C1SwitchingFunction make(std::vector<Value> values);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real myCutoff, myCutoff2, myCutoff_3, my15Cutoff_1, my05Cutoff_3,
         my15Cutoff_3;
  };

  //____ INLINES
}
#endif /* C1SWITCHINGFUNCTION_H */
