/* -*- c++ -*- */
#ifndef COULOMBSCPISMFORCE_H
#define COULOMBSCPISMFORCE_H

#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/config/Parameter.h>
#include <protomol/base/PMConstants.h>
#include <protomol/base/Report.h>
#include <string>
#include <math.h>

#include <iostream>
using namespace std;

using namespace ProtoMol::Report;

namespace ProtoMol {
  //____ CoulombForce
  class CoulombSCPISMForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Electrostatic interaction term for SCP ISM for MD of
    // Hassan et al. PROTEINS (51):109 (2003)
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    // Default constructor
    CoulombSCPISMForce();
    // Constructor with parameters
    CoulombSCPISMForce(Real Dval);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CoulombSCPISForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, Real &force,
                    Real /*distSquared*/, Real rDistSquared, const Vector3D &,
                    const GenericTopology *topo,
                    int atom1, int atom2, ExclusionClass excl) const {

      //SCPISM valid?
      if(!topo->doSCPISM)
        report << error << "CoulombSCPISMForce requires SCPISM parameters." << endr;
      // If either molecule belongs to a water, do nothing.
      // Won't happen in most simulations, but could in the 
      // case of comparing forces.
      if (topo->molecules[topo->atoms[atom1].molecule].water ||
	      topo->molecules[topo->atoms[atom2].molecule].water)
	        return;

      //std::cout << "EPS: " << EPS << std::endl;
      //std::cout << "D: " << D << std::endl;
      //std::cout << "S: " << S << std::endl;

      Real rDist = sqrt(rDistSquared);
      Real Dist = 1.0 / rDist;

      // Screening term for interaction energy calculation
      Real alpha_ij = topo->atoms[atom1].mySCPISM_A->sqrtalphaSCPISM *
                      topo->atoms[atom2].mySCPISM_A->sqrtalphaSCPISM;
      Real Dp1 = D + 1.0;
      Real rDp1 = 1.0 / Dp1;
      Real Dm1 = D - 1.0;
      Real k = Dm1 * 0.5;
      Real rDiElecEXP = k * exp(-alpha_ij * Dist);
      Real DiElec = Dp1 / (1.0 + rDiElecEXP) - 1.0;
      Real rDiElec = 1.0 / DiElec;

      energy = topo->atoms[atom1].scaledCharge *
               topo->atoms[atom2].scaledCharge *
               rDist * rDiElec *
               (excl == EXCLUSION_MODIFIED ? topo->coulombScalingFactor : 1);

      // reciprocal force term calculation
      Real drDiElec = alpha_ij * rDp1 * (1.0 + DiElec) * (D - DiElec);
      force = (rDiElec * drDiElec + rDist) * energy * rDist;
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
    static unsigned int getParameterSize() {return 1;}
    static CoulombSCPISMForce make(const std::vector<Value> &);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Sub Classes
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    class C1 {
    public:
      static Real kernel(Real r)                   {return 1.0 / r;}
      static Real dKernel(Real r)                  {r = 1.0 / r; return -r * r;}
      static Real kernelR(Real rr)                 {return rr;}
      static Real dKernelR(Real rr)                {return -rr * rr;}
      static Real smooth(Real r, Real /*c*/, Real cr) {
        return cr * (1.5 - 0.5 * r * r * cr * cr);
      }

      static Real smooth0(Real /*c*/, Real cr)         {return 1.5 * cr;}
      static Real dSmooth(Real r, Real /*c*/, Real cr) {
        return -r * cr * cr * cr;
      }

      static Real smoothKernel(Real r, Real c, Real cr) {
        return r < c ?  smooth(r, c, cr) : kernel( r);
      }

      static Real dSmoothKernel(Real r, Real c, Real cr) {
        return r < c ? dSmooth(r, c, cr) : dKernel(r);
      }

      static std::string getKeyword()                 {return keyword;}

      static std::string getForceKeyword() {
        return CoulombSCPISMForce::keyword;
      }

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
        return r < c ?  smooth(r, c, cr) :  kernel(r);
      }

      static Real dSmoothKernel(Real r, Real c, Real cr) {
        return r < c ? dSmooth(r, c, cr) : dKernel(r);
      }

      static std::string getKeyword()                 {return keyword;}

      static std::string getForceKeyword() {return CoulombSCPISMForce::keyword;}

    public:
      static const std::string keyword;
    };
  public:
    class C3 {
    public:
      static Real kernel(Real r)                      {return 1.0 / r;}
      static Real dKernel(Real r)                 {r = 1.0 / r; return -r * r;}
      static Real kernelR(Real rr)                    {return rr;}
      static Real dKernelR(Real rr)                   {return -rr * rr;}

      static Real smooth(Real r, Real /*c*/, Real cr)   {
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

      static std::string getForceKeyword() {return CoulombSCPISMForce::keyword;}

    public:
      static const std::string keyword;
    };
  public:
    class C4 {
    public:
      static Real kernel(Real r)                      {return 1.0 / r;}
      static Real dKernel(Real r)                {r = 1.0 / r; return -r * r;}
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

      static std::string getForceKeyword() {return CoulombSCPISMForce::keyword;}

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
    Real D;
  };

  //____ INLINES
}
#endif /* COULOMBSCPISMFORCE_H */
