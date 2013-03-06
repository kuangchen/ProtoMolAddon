/* -*- c++ -*- */
#ifndef BORNSWITCH_H
#define BORNSWITCH_H

#include <vector>

#include <protomol/config/Parameter.h>
#include <protomol/type/Matrix3By3.h>

namespace ProtoMol {
  //____ BornSwitch

  /**
   * The switching function for born calculations
   * Equation numbers from "Notes on SCPISM" by C.R.Sweet, based on Hassan's papers.
   */
  class BornSwitch {
  public:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types and Enums
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    enum {CUTOFF = 5};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    BornSwitch(){sw = 3;}; //default switch
    BornSwitch(int s) {sw = s;};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class BornSwitch
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Real switchValue(Real r) const {  //Eq. (3), (23)
      Real f;
      if (r >= CUTOFF) return 0;      
      if (sw == 3) {
        f = (1 - .0016 * r * r * r * r);
        return f * f * f * f;
      } else if (sw == 1) {
        f = (1 - .04 * r * r);
        return f * f;
      } else if (sw == 2) {
        f = (1 - 0.008 * r * r * r);
        return f * f * f;
      } else {
        return 1;
      }
    }

    Real switchDerivative(Real r) const {  //Eq. (6), (24)
      Real f, fr;
      if (r >= CUTOFF) return 0;      
      if (sw == 3) {
        f = (1 - .0016 * r * r * r * r);
        fr = -16 * .0016 * r * r * r * f * f * f;
        return fr;
      } else if (sw == 1) {
        fr = -0.16 * (1 - 0.04 * r * r) * r;
        return fr;
      } else if (sw == 2) {
        fr = -9 * .008 * r * r * (1 - .008 * r * r * r) *
          (1 - .008 * r * r * r);
        return fr;
      } else {
        return 0;
      }
    }

    Real switchSecondDerivative(Real r) const {  //Eq. (20), (25)
      Real f, frr;
      if (r >= CUTOFF) return 0;      
      if (sw == 3) {
        f = (1 - .0016 * r * r * r * r);
        frr = 0.00049152 * f * f * r * r * r * r * r * r - 0.0768 * f * f * f *
              r * r;
        return frr;
      } else if (sw == 1) {
        frr = 0.0192 * r * r - 0.16;
        return frr;
      } else if (sw == 2) {
        f = (1 - 0.008 * r * r * r);
        frr = 0.003456 * f * r * r * r * r - 0.144 * f * f * r;
        return frr;
      } else {
        return 0;
      } 
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    int sw;
  };

  //____ INLINES
}
#endif /* BornSwitch_H */
