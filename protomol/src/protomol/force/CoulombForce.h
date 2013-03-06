/* -*- c++ -*- */
#ifndef COULOMBFORCE_H
#define COULOMBFORCE_H

#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/config/Parameter.h>
#include <string>

namespace ProtoMol {
  //____ CoulombForce
  class CoulombForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // This uses the weighted charges on each atom, so the Coulomb
    // constant here is one.
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  public:
    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CoulombForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, Real &force, Real /*distSquared*/,
                    Real rDistSquared, const Vector3D &,
                    const GenericTopology *topo, int atom1, int atom2,
                    ExclusionClass excl) const {
      energy = topo->atoms[atom1].scaledCharge *
               topo->atoms[atom2].scaledCharge *
               sqrt(rDistSquared);
      if (excl == EXCLUSION_MODIFIED)
        energy *= topo->coulombScalingFactor;

      force = energy * rDistSquared;
    }

    static void accumulateEnergy(ScalarStructure *energies, Real energy) {
      (*energies)[ScalarStructure::COULOMB] += energy;
    }

    static Real getEnergy(const ScalarStructure *energies) {
      return (*energies)[ScalarStructure::COULOMB];
    }

    // Parsing
    static std::string getId() {return keyword;}
    static unsigned int getParameterSize() {return 0;}
    void getParameters(std::vector<Parameter> &) const {}

    static CoulombForce make(const std::vector<Value> &) {
      return CoulombForce();
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Sub Classes
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    class C1 {
    public:
      static Real kernel(Real r) {return 1.0 / r;}

      static Real dKernel(Real r) {r = 1.0 / r; return -r *
                                                              r;}

      static Real kernelR(Real rr) {return rr;}

      static Real dKernelR(Real rr) {return -rr * rr;}

      static Real smooth(Real r, Real /*c*/, Real cr) {
        return cr * (1.5 - 0.5 * r * r * cr * cr);
      }

      static Real smooth0(Real /*c*/, Real cr) {return 1.5 * cr;}

      static Real dSmooth(Real r, Real /*c*/, Real cr) {
        return -r * cr * cr * cr;
      }

      static Real smoothKernel(Real r, Real c, Real cr) {
        return r < c ?  smooth(r, c, cr) : kernel(r);
      }

      static Real dSmoothKernel(Real r, Real c, Real cr) {
        return r < c ? dSmooth(r, c, cr) : dKernel(r);
      }

      static std::string getKeyword() {return keyword;}

      static std::string getForceKeyword() {return CoulombForce::keyword;}

    public:
      static const std::string keyword;
    };
  public:
    class C2 {
    public:
      static Real kernel(Real r) {return 1.0 / r;}

      static Real dKernel(Real r) {r = 1.0 / r; return -r * r;}

      static Real kernelR(Real rr) {return rr;}

      static Real dKernelR(Real rr) {return -rr * rr;}

      static Real smooth(Real r, Real /*c*/, Real cr) {
        r = r * r * cr * cr;
        return cr * (1.875 - r * (1.25 - 0.375 * r));
      }

      static Real smooth0(Real /*c*/, Real cr) {return 1.875 * cr;}

      static Real dSmooth(Real r, Real c, Real cr) {
        c = r * cr * cr;
        return cr * c * (1.5 * r * c - 2.5);
      }

      static Real smoothKernel(Real r, Real c, Real cr) {
        return r < c ?  smooth(r, c, cr) : kernel(r);
      }

      static Real dSmoothKernel(Real r, Real c, Real cr) {
        return r < c ? dSmooth(r, c, cr) : dKernel(r);
      }

      static std::string getKeyword() {return keyword;}

      static std::string getForceKeyword() {return CoulombForce::keyword;}

    public:
      static const std::string keyword;
    };
  public:
    class C3 {
    public:
      static Real kernel(Real r) {return 1.0 / r;}

      static Real dKernel(Real r) {r = 1.0 / r; return -r *
                                                              r;}

      static Real kernelR(Real rr) {return rr;}

      static Real dKernelR(Real rr) {return -rr * rr;}

      static Real smooth(Real r, Real /*c*/, Real cr) {
        r = r * r * cr * cr;
        return 0.0625 * cr * (35.0 - r * (35.0 - r * (21.0 - 5.0 * r)));
      }

      static Real smooth0(Real /*c*/, Real cr) {return 2.1875 * cr;}

      static Real dSmooth(Real r, Real c, Real cr) {
        c = r * r * cr * cr;
        return r * cr * cr * cr * (-4.375 + c * (5.25 - 1.875 * c));
      }

      static Real smoothKernel(Real r, Real c, Real cr) {
        return r < c ?  smooth(r, c, cr) :  kernel(r);
      }

      static Real dSmoothKernel(Real r, Real c, Real cr) {
        return r < c ? dSmooth(r, c, cr) : dKernel(r);
      }

      static std::string getKeyword() {return keyword;}

      static std::string getForceKeyword() {return CoulombForce::keyword;}

    public:
      static const std::string keyword;
    };
  public:
    class C4 {
    public:
      static Real kernel(Real r) {return 1.0 / r;}

      static Real dKernel(Real r) {r = 1.0 / r; return -r * r;}

      static Real kernelR(Real rr) {return rr;}

      static Real dKernelR(Real rr) {return -rr * rr;}

      static Real smooth(Real r, Real /*c*/, Real cr) {
        r = r * r * cr * cr;
        return 0.0078125 * cr *
          (315.0 - r * (420.0 - r * (378.0 - r * (180.0 - r * 35.0))));
      }

      static Real smooth0(Real /*c*/, Real cr) {return 2.4609375 * cr;}

      static Real dSmooth(Real r, Real c, Real cr) {
        c = r * r * cr * cr;
        return -r * cr * cr * cr *
          (6.5625 - c * (11.8125 - c * (8.4375 - c * 2.1875)));
      }

      static Real smoothKernel(Real r, Real c, Real cr) {
        return r < c ?  smooth(r, c, cr) :  kernel(r);
      }

      static Real dSmoothKernel(Real r, Real c, Real cr) {
        return r < c ? dSmooth(r, c, cr) : dKernel(r);
      }

      static std::string getKeyword() {return keyword;}

      static std::string getForceKeyword() {return CoulombForce::keyword;}

    public:
      static const std::string keyword;
    };

  public:
    static const std::string keyword;
  private:
  };
}
#endif /* COULOMBFORCE_H */
