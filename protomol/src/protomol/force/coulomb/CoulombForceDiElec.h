/* -*- c++ -*- */
#ifndef COULOMBFORCEDIELEC_H
#define COULOMBFORCEDIELEC_H

#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/config/Parameter.h>
#include <string>
#include <math.h>

namespace ProtoMol {
  //____ CoulombForce
  class CoulombForceDiElec {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // This uses the weighted charges on each atom.
    // A DiElectric value is added to the electrostatic computation according
    // to Shen and Freed 2002 "All-atom fast protein folding simulations:The
    // Villin Headpiece" PROTEINS:Structure, Function, and Genetics 49:439-445
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  public:
    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    // Default constructor
    CoulombForceDiElec();
    // Constructor with parameters
    CoulombForceDiElec(Real EPSval, Real Dval, Real Sval);
    // Dielectric function constant terms D = 53.0 and S = 0.25
    // To change the value of eps(0)= 1 given with DiElecCom
    // Set X in DiElecCom = (D-X)*0.5 to your chosen value 2,3,4,etc.. 
    // which yields eps(0)=2,3,4,etc...
    // eps(0)=1 is insufficient to model solvent polarization effects between
    // near molecules
    // Try 2,3,or 4

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CoulombForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy,
                    Real &force,
                    Real /*distSquared*/, Real rDistSquared, const Vector3D &,
                    const GenericTopology *topo,
                    int atom1, int atom2, ExclusionClass excl) const {
      Real rDist = sqrt(rDistSquared);
      Real Dist = 1.0 / rDist;

      Real DiElecCom = (D - EPS) * 0.5;

      // Dielectric term for energy calculation
      Real rDiElecEXP = exp(-S * Dist);
      Real rDiElecP1 = DiElecCom * (S * S * Dist * Dist + 2 * S * Dist + 2);
      Real rDiElec = 1.0 / (D - (rDiElecP1 * rDiElecEXP));

      energy = topo->atoms[atom1].scaledCharge *
               topo->atoms[atom2].scaledCharge *
               rDist * rDiElec;

      if (excl == EXCLUSION_MODIFIED)
        energy *= topo->coulombScalingFactor;

      // reciprocal force term calculation includes dielectric derivative
      Real drDistrDiElecPoly = S * S * Dist * Dist * Dist + 2 * S * Dist *
                               Dist + 2 * Dist;
      Real drDistrDiElecP1 = DiElecCom * rDiElecEXP *
                             (3 * S * S * Dist * Dist + 4 * S * Dist + 2);
      Real drDistrDiElecP2 = DiElecCom * rDiElecEXP * S * drDistrDiElecPoly;
      Real drDistrDiElecP3 = D * Dist - DiElecCom * rDiElecEXP *
                             drDistrDiElecPoly;

      Real drDistrDiElec =
        (D - drDistrDiElecP1 + drDistrDiElecP2) / 
        (drDistrDiElecP3 * drDistrDiElecP3);
      force = energy / rDiElec * drDistrDiElec;

      // no negative term - magnitude only
      // additional rDist term normalized force for later mult by directional
      // vec components
    }

    static void accumulateEnergy(ScalarStructure *energies, Real energy) {
      (*energies)[ScalarStructure::COULOMB] += energy;
    }

    static Real getEnergy(const ScalarStructure *energies) {
      return (*energies)[ScalarStructure::COULOMB];
    }

    // Parsing
    static std::string getId() {return keyword;}
    void getParameters(std::vector<Parameter> &) const;
    static unsigned int getParameterSize() {return 3;}
    static CoulombForceDiElec make(const std::vector<Value> &);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Sub Classes
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    class C1 {
    public:
      static Real kernel(Real r)                      {return 1.0 / r;}
      static Real dKernel(Real r)                {r = 1.0 / r; return -r * r;}
      static Real kernelR(Real rr)                    {return rr;}
      static Real dKernelR(Real rr)                   {return -rr * rr;}
      static Real smooth(Real r, Real /*c*/, Real cr) {
        return cr * (1.5 - 0.5 * r * r * cr * cr);
      }

      static Real smooth0(Real /*c*/, Real cr)         {return 1.5 * cr;}

      static Real dSmooth(Real r, Real /*c*/, Real cr) {
        return -r * cr * cr * cr;
      }

      static Real smoothKernel(Real r, Real c, Real cr) {
        return r < c ?  smooth(r, c, cr) : kernel(r);
      }

      static Real dSmoothKernel(Real r, Real c, Real cr) {
        return r < c ? dSmooth(r, c, cr) : dKernel(r);
      }

      static std::string getKeyword()                 {return keyword;}

      static std::string getForceKeyword() {return CoulombForceDiElec::keyword;}

    public:
      static const std::string keyword;
    };
  public:
    class C2 {
    public:
      static Real kernel(Real r)                      {return 1.0 / r;}
      static Real dKernel(Real r)                  {r = 1.0 / r; return -r * r;}
      static Real kernelR(Real rr)                    {return rr;}
      static Real dKernelR(Real rr)                   {return -rr * rr;}
      static Real smooth(Real r, Real /*c*/, Real cr) {
        r = r * r * cr * cr;
        return cr * (1.875 - r * (1.25 - 0.375 * r));
      }

      static Real smooth0(Real /*c*/, Real cr)         {return 1.875 * cr;}

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

      static std::string getKeyword()                 {return keyword;}
      static std::string getForceKeyword() {return CoulombForceDiElec::keyword;}

    public:
      static const std::string keyword;
    };
  public:
    class C3 {
    public:
      static Real kernel(Real r)                      {return 1.0 / r;}
      static Real dKernel(Real r)                  {r = 1.0 / r; return -r * r;}

      static Real kernelR(Real rr)                    {return rr;}

      static Real dKernelR(Real rr)                   {return -rr * rr;}

      static Real smooth(Real r, Real /*c*/, Real cr) {
        r = r * r * cr * cr;
        return 0.0625 * cr * (35.0 - r * (35.0 - r * (21.0 - 5.0 * r)));
      }

      static Real smooth0(Real /*c*/, Real cr)         {return 2.1875 * cr;}

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

      static std::string getKeyword()                 {return keyword;}

      static std::string getForceKeyword() {return CoulombForceDiElec::keyword;}

    public:
      static const std::string keyword;
    };
  public:
    class C4 {
    public:
      static Real kernel(Real r)                      {return 1.0 / r;}
      static Real dKernel(Real r)                  {r = 1.0 / r; return -r * r;}
      static Real kernelR(Real rr)                    {return rr;}
      static Real dKernelR(Real rr)                   {return -rr * rr;}
      static Real smooth(Real r, Real /*c*/, Real cr) {
        r = r * r * cr * cr;
        return 0.0078125 * cr *
          (315.0 - r * (420.0 - r * (378.0 - r * (180.0 - r * 35.0))));}

      static Real smooth0(Real /*c*/, Real cr)         {return 2.4609375 * cr;}

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

      static std::string getKeyword()                 {return keyword;}
      static std::string getForceKeyword() {return CoulombForceDiElec::keyword;}

    public:
      static const std::string keyword;
    };
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real EPS;
    Real D;
    Real S;
  };
}
#endif /* COULOMBFORCEDIELEC_H */
