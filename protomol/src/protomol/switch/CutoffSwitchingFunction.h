/* -*- c++ -*- */
#ifndef CUTOFFSWITCHINGFUNCTION_H
#define CUTOFFSWITCHINGFUNCTION_H

#include <vector>

#include <protomol/config/Parameter.h>

namespace ProtoMol {
  //____ CutoffSwitchingFunction

  /**
   * Cutoff switching function, implements a simple truncation of the potential.
   */
  class CutoffSwitchingFunction {
  public:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types and Enums
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    enum {USE = 1};
    enum {MODIFY = 0};
    enum {CUTOFF = 1};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    CutoffSwitchingFunction() : myCutoff(0.0), myCutoff2(0.0) {}
    CutoffSwitchingFunction(Real cutoff) :
      myCutoff(cutoff), myCutoff2(cutoff * cutoff) {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CutoffSwitchingFunction
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /// simple and fast test if we should apply the switching function
    bool roughTest(Real distSquared) const {return distSquared <= myCutoff2;}
    Real cutoffSquared() const {return myCutoff2;}
    Real cutoff() const {return myCutoff;}

    void operator()(Real &value, Real &derivOverD, Real distSquared) const {
      derivOverD = 0;
      value = (distSquared > myCutoff2 ? 0.0 : 1.0);
    }

    static const std::string getId() {return "Cutoff";}
    void getParameters(std::vector<Parameter> &parameters) const {
      parameters.push_back
        (Parameter("-cutoff", Value(myCutoff, ConstraintValueType::Positive()),
                   Text("cutoff swf cutoff")));
    }
    static unsigned int getParameterSize() {return 1;}

    static CutoffSwitchingFunction make(std::vector<Value> values) {
      return CutoffSwitchingFunction(values[0]);
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real myCutoff, myCutoff2;
  };
}
#endif /* CUTOFFSWITCHINGFUNCTION_H */
